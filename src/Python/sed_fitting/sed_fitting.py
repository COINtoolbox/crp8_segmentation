import os
import glob
import pandas as pd
import numpy as np
from multiprocessing import Pool
from sed_utils import process_sed, sedpy_filter, build_model, build_sps, sfr_exp

HOME = os.environ["HOME"]

# galaxies csv files with fluxes and errors per region
data_dir = f"{HOME}/flux_galaxies"
csv_files = sorted(glob.glob(os.path.join(data_dir, "*.csv")))

# csv with redshifts
redshift_csv = f"{HOME}/galaxy_redshifts.csv"
df_z = pd.read_csv(redshift_csv)  # columns: 'gal_name', 'z'
gal_redshifts = dict(zip(df_z.gal_name, df_z.z))

# parameters 
run_params = {
    "object_redshift": None,  # will be updated according to each galaxy
    "fixed_metallicity": True,
    "add_duste": True,
    "add_agn": False,
    "verbose": True,
    "zcontinuous": 1,
    "dynesty": False,
    "emcee": False,
    "optimize": True,
    "min_method": 'lm',
    "nmin": 16,
    "nwalkers": 100,
    "random_seed": 42
}

if __name__ == "__main__":
    np.random.seed(42)

    for csv_file in csv_files:
        # galaxy name/id
        gal_name = os.path.splitext(os.path.basename(csv_file))[0]

        # outputs 
        output_dir = f"{HOME}/teste_sed_coin/{gal_name}"
        os.makedirs(output_dir, exist_ok=True)

        # redshift
        zgal = gal_redshifts.get(gal_name, None)
        if zgal is None:
            print(f"Warning: redshift not found for galaxy {gal_name}. Using object_redshift=None")

        run_params_gal = run_params.copy()
        run_params_gal["object_redshift"] = zgal

        # read fluxes and errors 
        df = pd.read_csv(csv_file)
        flux_cols = sorted([c for c in df.columns if c.startswith("F") and not c.endswith("_err") and not c.endswith("_n_eff")])
        err_cols  = sorted([c for c in df.columns if c.endswith("_err")])

        df["fluxes"] = df[flux_cols].values.tolist()
        df["errors"] = df[err_cols].values.tolist()

        # filters in sedpy
        filters = [sedpy_filter(c) for c in flux_cols]

        # input for each region
        region_inputs = [(i, df.fluxes[i], df.errors[i]) for i in range(len(df))]

        # build model
        print(f"Building model for {gal_name} with z = {zgal}...")
        model = build_model(**run_params_gal)
        sps = build_sps(**run_params_gal)
        print("done")

        # process
        final_vals = []
        with Pool(processes=4) as pool:  # you can change
            results = pool.starmap(
                process_sed,
                [(i, flux, err, model, sps, filters, output_dir, run_params_gal) for i, flux, err in region_inputs]
            )
        final_vals.extend(results)

        # save results
        df_results = pd.DataFrame(final_vals, columns=model.free_params)

        # add region column (starting at 1, not 0)
        df_results.insert(0, "region", np.arange(1, len(df_results) + 1))
        df_results["sfr"] = sfr_exp(tage=df_results["tage"].values, tau=df_results["tau"].values, mass=df_results["mass"].values)

        out_csv = os.path.join(output_dir, f"{gal_name}_sed_results.csv")
        df_results.to_csv(out_csv, index=False)



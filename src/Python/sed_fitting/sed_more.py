#more properties

import os
import glob
import pandas as pd
import numpy as np
import astropy.units as u
from astropy.cosmology import LambdaCDM

DEFAULT_COSMO = LambdaCDM(H0=70, Om0=0.3, Ode0=0.7)

# pixel scale (arcsec/pixel)
pixel_scale_arcsec = 0.02


def kpc_per_pixel(z, pix_scale_arcsec, cosmo=DEFAULT_COSMO):
    """
    Convert pixel scale from arcsec/pixel to kpc/pixel at a given redshift.
    """
    DA = cosmo.angular_diameter_distance(z)
    kpc_per_arcsec = (DA * u.arcsec.to(u.rad)).to(u.kpc).value
    return pix_scale_arcsec * kpc_per_arcsec


# paths results from sed_fitting.py
data_dir = f"{HOME}/sed_results"
gal_dirs = sorted(glob.glob(os.path.join(data_dir, "*")))

# redshifts 
redshift_csv = f"{HOME}/galaxy_redshifts.csv"
df_z = pd.read_csv(redshift_csv)
gal_redshifts = dict(zip(df_z.gal_name, df_z.z))


for gal_dir in gal_dirs:

    gal_name = os.path.basename(gal_dir)
    csv_file = os.path.join(gal_dir, f"{gal_name}_sed_results.csv")

    if not os.path.exists(csv_file):
        print(f"Skipping {gal_name}: no SED results found")
        continue

    print(f"\nProcessing {gal_name}...")

    # read results
    df_results = pd.read_csv(csv_file)

    # region
    if "region" not in df_results.columns:
        print(f"Warning: 'region' column not found in {gal_name}, creating default index")
        df_results.insert(0, "region", np.arange(1, len(df_results) + 1))

    # redshift
    zgal = gal_redshifts.get(gal_name, None)

    if "zred" in df_results.columns:
        z_used = df_results["zred"].median()
        print(f"Using median fitted redshift: z = {z_used:.4f}")
    else:
        if zgal is None:
            print(f"Warning: no redshift available for {gal_name}, skipping")
            continue
        z_used = zgal
        print(f"Using input redshift: z = {z_used:.4f}")

    # physical scale
    pixel_scale_kpc = kpc_per_pixel(z_used, pixel_scale_arcsec)

    # more properties
    df_results["z_used"] = z_used
    df_results["kpc_per_pixel"] = pixel_scale_kpc

    df_results["area_kpc2"] = df_results["n_pix"] * (pixel_scale_kpc ** 2)

    df_results["sigma_sfr_10myr"] = df_results["sfr_10myr"] / df_results["area_kpc2"]
    df_results["sigma_sfr_100myr"] = df_results["sfr_100myr"] / df_results["area_kpc2"]

    # save
    out_csv = os.path.join(gal_dir, f"{gal_name}_sed_results_area.csv")
    df_results.to_csv(out_csv, index=False)

    print(f"Saved: {out_csv}")

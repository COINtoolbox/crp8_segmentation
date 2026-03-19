# sed_utils.py

import os
HOME = os.environ["HOME"]
os.environ['SPS_HOME'] = f"{HOME}/fsps"
import fsps
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from prospect.models import priors
from prospect.models.templates import TemplateLibrary
from prospect.fitting import fit_model, lnprobfn
from prospect.likelihood import chi_spec, chi_phot
from prospect.utils.obsutils import fix_obs
from prospect.sources import CSPSpecBasis
from prospect.models.sedmodel import SedModel
import sedpy


# --- Helper functions ---
def sedpy_filter(colname):
    """
    Map a flux column name (e.g. 'F090W') to a sedpy filter name (e.g. 'jwst_f090w')
    """
    base = colname.strip().replace('"', '').replace("'", '').lower()
    
    # filters HST
    hst_filters = ['f225w','f275w','f336w','f435w','f606w','f814w','f105w','f125w','f140w','f160w']
    
    if base in hst_filters:
        return 'hst_' + base
    else:
        # assumes JWST if is not in HST filters
        return 'jwst_' + base


def build_obs(flux, flux_errs=None, filters=None, snr=10, **extras):
    """
    Build a dictionary of observational data.
    
    Parameters
    ----------
    flux : array_like
        Observed fluxes in Jy, one per filter.
    flux_errs : array_like, optional
        Flux uncertainties in Jy. If None, uncertainties are estimated
        assuming a constant signal-to-noise ratio given by `snr`.
    filters : list of str
        List of sedpy filter names corresponding to the flux measurements.
    snr : float, optional
        Signal-to-noise ratio used to estimate flux uncertainties when
        `flux_errs` is not provided. Default is 10.
    **extras : dict
        Additional keyword arguments.

    Returns
    -------
    obs : dict
        Observational dictionary.
    """
    
    if filters is None:
        raise ValueError("filters must be provided as a list of sedpy filter names")
    
    flux = np.array(flux)
    
    if flux_errs is None or len(flux_errs) != len(flux):
        flux_errs = (1./snr) * flux
    flux_errs = np.array(flux_errs)
    
    obs = {}
    obs["filters"] = sedpy.observate.load_filters(filters)
    obs["maggies"] = flux
    obs["maggies_unc"] = flux_errs
    obs["phot_wave"] = np.array([f.wave_effective for f in obs["filters"]])
    obs["wavelength"] = None
    obs["spectrum"] = None
    obs['unc'] = None
    obs['mask'] = None

    return fix_obs(obs)


def age_of_universe_gyr_at_z(z, H0_kms_Mpc=70.0, Om=0.3, Ol=0.7):
    """
    Numerical integral to compute age of universe at redshift z (Gyr)
    using a flat LambdaCDM with H0 (km/s/Mpc), Om, Ol.
    """
    # Convert constants
    pc = 3.085677581e16       # meters
    Mpc = pc * 1e6
    H0_SI = H0_kms_Mpc * 1000.0 / Mpc   # s^-1
    sec_per_Gyr = 1e9 * 365.25 * 24 * 3600.0
    H0_Gyr = H0_SI * sec_per_Gyr        # in Gyr^-1

    # integrate from z -> zmax
    zmax = 1e5
    zs = np.linspace(z, zmax, 200000)
    integrand = 1.0 / ((1.0 + zs) * np.sqrt(Om * (1.0 + zs)**3 + Ol))
    integral = np.trapz(integrand, zs)   # in units of 1/H0
    age_gyr = integral / H0_Gyr
    return float(age_gyr)


def build_model(object_redshift=None, fixed_metallicity=None, add_duste=False, add_agn=False, **extras):
    """Build a Prospector SedModel object with optional free redshift"""
    model_params = TemplateLibrary["parametric_sfh"]

    # Initial guesses
    model_params["dust2"]["init"] = 1
    model_params["logzsol"]["init"] = 1.0
    model_params["mass"]["init"] = 1e8

    # Priors
    model_params["dust2"]["prior"] = priors.TopHat(mini=0.0, maxi=5.0)
    model_params["tau"]["prior"] = priors.LogUniform(mini=0.5, maxi=10.0)
    model_params["mass"]["prior"] = priors.TopHat(mini=10**5, maxi=10**13.5)
    model_params["logzsol"]["prior"] = priors.TopHat(mini=-2, maxi=0.3)

    # IMF: Chabrier (2003)
    model_params["imf_type"]["init"] = 1
    model_params["imf_type"]["isfree"] = False

    # Dust attenuation law: Calzetti (2000)
    model_params["dust_type"]["init"] = 2
    model_params["dust_type"]["isfree"] = False

    # Metallicity
    if fixed_metallicity is not None:
        model_params["logzsol"]["isfree"] = True
        model_params["logzsol"]["init"] = fixed_metallicity

    # Redshift
    if object_redshift is not None:
        model_params["zred"]["isfree"] = False
        model_params["zred"]["init"] = object_redshift
        age_at_z = age_of_universe_gyr_at_z(object_redshift)
        model_params["tage"]["prior"] = priors.TopHat(mini=0.01, maxi=max(0.01, age_at_z))
        model_params["tage"]["init"] = min(6.0, age_at_z*0.8)
    else:
        model_params["zred"]["isfree"] = True
        model_params["zred"]["init"] = 0.7
        model_params["zred"]["prior"] = priors.Uniform(mini=0.0, maxi=2.0)
        model_params["tage"]["prior"] = priors.TopHat(mini=0.01, maxi=13.8)
        model_params["tage"]["init"] = 4.0

    # Dust emission
    if add_duste:
        model_params.update(TemplateLibrary["dust_emission"])
        
    # if add_duste:
    #     model_params['add_dust_emission'] = TemplateLibrary['alpha']['add_dust_emission']
    #     model_params['duste_umin'] = TemplateLibrary['alpha']['duste_umin']
    #     model_params['duste_qpah'] = TemplateLibrary['alpha']['duste_qpah']
    #     model_params['duste_gamma'] = TemplateLibrary['alpha']['duste_gamma']
    #     model_params['dust_type']['init'] = 4
    #     model_params["duste_qpah"]["prior"] = priors.TopHat(mini=0.1, maxi=10)

    # AGN
    if add_agn:
        model_params['fagn'] = TemplateLibrary['alpha']['fagn']
        model_params['agn_tau'] = TemplateLibrary['alpha']['agn_tau']
        model_params['add_agn_dust'] = TemplateLibrary['alpha']['add_agn_dust']
        model_params['fagn']['isfree'] = True
        model_params['agn_tau']['isfree'] = True
        model_params['add_agn_dust']['init'] = True

    return SedModel(model_params)


def build_sps(zcontinuous=1, **extras):
    """Build SPS object for Prospector"""
    return CSPSpecBasis(zcontinuous=zcontinuous)


def plot_sed(wspec, pspec, wphot, obs, output_file):
    """Plot observed and model SED"""
    ymin, ymax = obs["maggies"].min()*0.8, obs["maggies"].max()/0.4
    cmap = plt.get_cmap('gist_rainbow_r')
    norm = mcolors.Normalize(vmin=0, vmax=len(obs["filters"])-1)

    plt.figure(figsize=(16,8))
    plt.loglog(wspec, pspec, label='Fitted spectrum', lw=1, color='mediumpurple', alpha=1)
    plt.errorbar(
        wphot, obs['maggies'], yerr=obs['maggies_unc'],
        label='Observed photometry', marker='o', markersize=10, alpha=0.8,
        ls='', lw=2, ecolor='black', markerfacecolor='none',
        markeredgecolor='black', markeredgewidth=3
    )

    for i, f in enumerate(obs['filters']):
        w, t = f.wavelength.copy(), f.transmission.copy()
        t = t / t.max()
        t = 10**(0.2*(np.log10(ymax/np.min(obs["maggies"])))) * t * np.min(obs["maggies"])
        plt.loglog(w, t, lw=3, color=cmap(norm(i)), alpha=0.7)

    plt.xlabel(r'Wavelength [$\AA$]', fontsize=15)
    plt.ylabel('Flux Density [maggies]', fontsize=15)
    #xmin = 3e3
    #xmax = 1e7
    #plt.xlim([xmin, xmax])
    plt.xlim([np.min(wphot)*0.8, np.max(wphot)/0.8])
    plt.ylim([ymin, ymax])
    plt.tick_params(axis='both', which='major', labelsize=15)
    plt.legend(loc='best', fontsize=15)
    plt.tight_layout()
    plt.savefig(output_file)
    plt.close()


def process_sed(region_index, flux, flux_err=None, model=None, sps=None, filters=None, output_dir=None, run_params={}):
    """Fit SED for a single region and save plot"""
    if output_dir is None:
        output_dir = os.getcwd()
    if filters is None:
        raise ValueError("filters must be provided as a list of sedpy filter names")

    flux = np.array(flux)/3631
    flux_err = np.array(flux_err)/3631 if flux_err is not None else None

    obs = build_obs(flux, flux_errs=flux_err, filters=filters)
    
    print(f"Fitting SED for region {region_index}...")
    output = fit_model(obs, model, sps, lnprobfn=lnprobfn, **run_params)
    print(f"Done optimization for region {region_index}")

    results, _ = output.get("optimization", ([], 0))
    if len(results) == 0:
        raise RuntimeError(f"Optimization failed for region {region_index}!")

    ind_best = np.argmin([r.cost for r in results])
    theta_best = results[ind_best].x.copy()

    pspec, pphot, _ = model.mean_model(theta_best, obs=obs, sps=sps)

    # --- reconstruir wspec após ajuste ---
    a = 1.0 + model.params.get('zred', 0.0)
    wspec = sps.wavelengths * a
    wphot = obs["phot_wave"]

    plot_file = os.path.join(output_dir, f"sedfit_region_{region_index+1}.png") # to start in 1, not 0
    plot_sed(wspec, pspec, wphot, obs, plot_file)

    return theta_best

def sfr_exp(tage, tau, mass):
    """
    Instantaneous SFR for an exponential SFH.

    Parameters
    ----------
    tage : float
        Time since the onset of star formation (Gyr), i.e. the duration over which the star formation history is defined
    tau : float
        SFH timescale in Gyr
    mass : float
        Total formed stellar mass in Msun

    Returns
    -------
    sfr_yr : float
        Instantaneous SFR in Msun/yr
    """

    A = mass / (tau * (1.0 - np.exp(-tage / tau)))  # Msun / Gyr
    sfr_gyr = A * np.exp(-tage / tau)               # Msun / Gyr
    sfr_yr = sfr_gyr / 1e9                          # Msun / yr
    
    return sfr_yr
    
    
def sfr_avg_exp(tage, tau, mass, dt):
    """
    Average SFR over the last dt Gyr for an exponential SFH.

    Parameters
    ----------
    tage : float or array
        Time since the onset of star formation (Gyr), i.e. the duration over which the star formation history is defined
    tau : float or array
        SFH timescale in Gyr
    mass : float or array
        Total formed stellar mass in Msun
    dt : float
        Time window in Gyr (e.g. 0.01 = 10 Myr, 0.1 = 100 Myr)

    Returns
    -------
    sfr : float or array
        Average SFR in Msun/yr
    """
    A = mass / (tau * (1.0 - np.exp(-tage / tau)))
    t1 = np.maximum(0.0, tage - dt)
    sfr_avg = (A * tau / dt) * (np.exp(-t1 / tau) - np.exp(-tage / tau)) # Msun / Gyr
    sfr_avg_yr = sfr_avg / 1e9                                           # Msun / yr
    
    return sfr_avg_yr

def mass_weighted_age_exp(tage, tau):
    """
    Mass-weighted stellar age for an exponential SFH.

    Parameters
    ----------
    tage : float or array
        Time since the onset of star formation (Gyr), i.e. the duration over which the star formation history is defined
    tau : float or array
        Exponential SFH timescale in Gyr

    Returns
    -------
    t_mw : float or array
        Mass-weighted stellar age in Gyr
    """
    exp_term = np.exp(-tage / tau)
    t_mw = (tau - (tage + tau) * exp_term) / (1.0 - exp_term)
    
    return t_mw
    

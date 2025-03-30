#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct  6 15:18:27 2022
Optimized Radial Velocity Correction Script

@author: maria.tsantaki
"""
from __future__ import print_function, division
import numpy as np
from astropy.io import fits as fits
from PyAstronomy import pyasl
from observations import read_observations, local_norm, read_raw_observations, snr
import matplotlib.pyplot as plt
import glob
from astroquery.gaia import Gaia
import astropy.units as u
from astropy.coordinates import SkyCoord
from astroquery.simbad import Simbad

# Set Gaia Data Release to DR3
Gaia.MAIN_GAIA_TABLE = "gaiadr3.gaia_source"
# Constants
RV_RANGE = (-200., 200.)  # Radial velocity range for cross-correlation
RV_STEP = 0.2  # Step size in km/s
BROADENING_RES = 100000  # Instrumental broadening resolution
GAIA_SEARCH_WIDTH = u.Quantity(0.03, u.deg)
GAIA_SEARCH_HEIGHT = u.Quantity(0.03, u.deg)
SPEED_OF_LIGHT = 299792.458  # km/s

def read_template(fname):
    """Reads and processes the template spectrum."""
    intervals = [5350., 6800.]
    tw, tf, _ = read_observations(fname, intervals[0], intervals[1])
    tf = pyasl.instrBroadGaussFast(tw, tf, BROADENING_RES, fullout=False, maxsig=None)
    tf = 1.0 - tf
    return tw, tf

def read_spec(fname):
    """Reads and normalizes observed spectra."""
    SNR = snr(fname)
    if SNR == None:
        print('SNR is None')
        return [np.nan], [np.nan]
    intervals = [5500., 6700.]
    dw, df, d = local_norm(fname, intervals, SNR)
    df = pyasl.instrBroadGaussFast(dw, df, BROADENING_RES, fullout=False, maxsig=None)
    df = 1.0 - df
    return dw, df

def velocity_shift(dw, df, tw, tf, name='lol', plot=False):
    """Computes the velocity shift using cross-correlation."""
    
    rv, cc = pyasl.crosscorrRV(dw, df, tw, tf, *RV_RANGE, RV_STEP, skipedge=10)
    # Find the index of maximum cross-correlation function
    maxind = np.argmax(cc)

    if abs(rv[maxind])>99.:
        # search Gaia DR3 for RV limits
        result_table = Simbad.query_object(name)

        if result_table:
            print(result_table)        
            ra = result_table['RA'][0]
            dec = result_table['DEC'][0]

            coord = SkyCoord(ra=ra, dec=dec, unit=(u.hourangle, u.deg), frame='icrs')       
            r = Gaia.query_object_async(coordinate=coord, width=GAIA_SEARCH_WIDTH, height=GAIA_SEARCH_HEIGHT,)
            if r and not np.isnan(r['radial_velocity'][0]):
                rv_range = (r['radial_velocity'][0] - 5, r['radial_velocity'][0] + 5)
                rv, cc = pyasl.crosscorrRV(dw, df, tw, tf, *rv_range, RV_STEP, skipedge=10)
                maxind = np.argmax(cc)
                print(f"Cross-correlation function is maximized at dRV = {rv[maxind]:.2f} km/s")
                v = rv[maxind]
            else:
                print("No RV shift found in Gaia.")
                v = 0.0
    print("Cross-correlation function is maximized at dRV = ", rv[maxind], " km/s")
    v = rv[maxind]

    if plot:
        plt.plot(rv, cc, 'bp-')
        plt.plot(rv[maxind], cc[maxind], 'ro')
        plt.title('CCF')
        plt.show()
    return v

def rv_correction(filename, velocity):
    """Applies radial velocity correction to FITS file header."""

    hdulist = fits.open(filename)
    header = hdulist[0].header #read header
    
    relat_wave = header['CRVAL1'] * np.sqrt((1. - velocity/SPEED_OF_LIGHT) / (1. + velocity/SPEED_OF_LIGHT))
    relat_step = header['CDELT1'] * np.sqrt((1. - velocity/SPEED_OF_LIGHT) / (1. + velocity/SPEED_OF_LIGHT))
    
    #Change the header with the new values
    header['CRVAL1'], header['CDELT1'], header['CD1_1'] = relat_wave, relat_step, relat_step

    #Save results
    hdulist.writeto(filename[:-5]+"_rv.fits", overwrite=True)
    hdulist.close()
    return

def plot_rv(observed, rv):

    wave, flux, _, _ = read_raw_observations(observed)
    plt.plot(wave, flux, label='original')
    wave_rv, flux, _, _ = read_raw_observations(observed[:-5]+'_rv.fits')
    plt.plot(wave_rv, flux, label=f'RV corrected ({rv:.2f} km/s)')
    plt.title(str(observed))
    plt.legend()
    plt.show()
    return

if __name__ == '__main__':

    template = 'HARPS.Archive_Sun-4_norm.fits'
    tw, tf = read_template(template)
    ffiles = glob.glob("data/*.fits")
    
    for observed in ffiles:
        name = observed.split('/')[-1].split('_')[0]
        print(f"Processing: {name}")
        dw, df = read_spec(observed)
        if dw[0] == np.nan:
            print('Something wrong with the spectra!')
        # # Plot template and data
        # plt.title("Template (blue) and data (red)")
        # plt.plot(dw, df)
        # plt.plot(tw, tf, 'b.-')
        # plt.show()
    else:
            rv = velocity_shift(dw, df, tw, tf, name=name, plot=False)
            # Apply RV correction and plot if significant
            rv_correction(observed, round(rv, 2))
            if abs(rv)>10:
                # plot this file onced it is saved!
                plot_rv(observed, rv)

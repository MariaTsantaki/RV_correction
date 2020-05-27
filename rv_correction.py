#!/usr/bin/python
from __future__ import print_function, division
import numpy as np
from astropy.io import fits as fits
from PyAstronomy import pyasl
from observations import read_obs_intervals, read_observations, snr
import matplotlib.pyplot as plt

def rv_correction(filename, velocity, output = "auto"):
    c = 299792.458 #km s-1
    hdulist = fits.open(filename)
    header = hdulist[0].header #read header
    start_wave = header['CRVAL1'] #initial wavelenght
    #print start_wave
    step = header['CDELT1'] #step in wavelenght
    #print step
    # New starting wavelenght and step
    relat_wave = start_wave * np.sqrt((1. - velocity/c)/(1. + velocity/c))
    relat_step = step * np.sqrt((1. - velocity/c)/(1. + velocity/c))
    #Change the header with the new values
    header['CRVAL1'] = relat_wave
    header['CDELT1'] = relat_step
    header['CD1_1'] = relat_step
    #print relat_wave
    #print relat_step
    #Save results
    if output == "auto":
        outputfile = filename[:-5]+"_rv.fits"
    else: 
        outputfile = output 
    hdulist.writeto(outputfile, overwrite=True)
    hdulist.close()
    return

def read_template(fname):
    '''Choose the template.'''
    SNR = snr(fname)
    intervals = [5339.8, 5500.40]
    intervals = [6500., 6800.]
    tw, tf, d = read_observations(fname, intervals[0]-1.0, intervals[1]+1.0)
    tf = 1.0 - tf
    return tw, tf

def read_spec(fname):
    '''Read spectra.'''
    SNR = snr(fname)
    intervals = [6500., 6800.]
    dw, df, d = read_obs_intervals(fname, [intervals], SNR)
    df = 1.0 - df
    return dw, df

def velocity_shift(dw, df, tw, tf):
    '''Carry out the cross-correlation.
    The RV-range is -30 - +30 km/s in steps of 0.6 km/s.
    The first and last 20 points of the data are skipped.'''
    rv, cc = pyasl.crosscorrRV(dw, df, tw, tf, -40., 40., 0.01, skipedge=100)
    # Find the index of maximum cross-correlation function
    maxind = np.argmax(cc)
    print("Cross-correlation function is maximized at dRV = ", rv[maxind], " km/s")
    #plt.plot(rv, cc, 'bp-')
    #plt.plot(rv[maxind], cc[maxind], 'ro')
    #plt.show()
    return rv[maxind]

if __name__ == '__main__':
    import argparse

    args = argparse.ArgumentParser(description='Calculate for velocity shifts and correct them.')
    args.add_argument('observed', type=str, help='Spectrum')
    args.add_argument('-p', '--plot', help='plot flag', action='store_true')
    args.add_argument('-t', '--template', type=str, help='Template', choices=['initial.spec', 'HARPS.Archive_Sun-4_norm.fits', 'HARPS.Archive_Arcturus_norm.fits', 'HARPS.GBOG_Procyon_norm.fits'], default='HARPS.Archive_Sun-4_norm.fits')
    args.add_argument('-o', '--output', type=str, help="Name of output file", default="auto")
    args.add_argument('-r', '--rv', type=float, help="radial velocity to correct", default=-99999.99)
    args = args.parse_args()

    template = args.template
    observed = args.observed
    plot_flag = args.plot
    output = args.output
    rv = args.rv

    print(template, observed, plot_flag, output)

    dw, df, d = read_observations(observed, 4000, 6900)
    if plot_flag:
        plt.plot(dw, df)
        plt.show()

    tw, tf = read_template(template)
    dw, df = read_spec(observed)
    if plot_flag:

        # Plot template and data
        plt.title("Template (blue) and data (red)")
        plt.plot(tw, tf, 'b.-')
        plt.plot(dw, df, 'r.-')
        plt.show()

    if rv == -99999.99:
        rv = velocity_shift(dw, df, tw, tf)
    print("Correction for rv:", rv)
    rv_correction(observed, round(rv,3), output=output)

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct  6 15:18:27 2022

@author: maria.tsantaki
"""
from __future__ import print_function, division
import numpy as np
from astropy.io import fits as fits
from PyAstronomy import pyasl
from observations import read_observations, local_norm, read_raw_observations, eso_fits, snr
import matplotlib.pyplot as plt
import argparse

def read_template(fname):
    '''Choose the template.'''
    intervals = [5350., 6800.]
    tw, tf, d = read_observations(fname, intervals[0], intervals[1])
    tf = 1.0 - tf
    return tw, tf

def read_spec(fname):
    '''Read spectra.'''
    SNR = snr(fname)
    intervals = [5500., 6700.]
    dw, df, d = local_norm(fname, intervals, SNR)
    df = 1.0 - df
    return dw, df

def velocity_shift(dw, df, tw, tf, plot=False):
    '''Carry out the cross-correlation.
    The RV-range is -200 - +200 km/s in steps of 0.2 km/s.
    The first and last 10 points of the data are skipped.'''
    
    rv, cc = pyasl.crosscorrRV(dw, df, tw, tf, -200., 200., 0.2, skipedge=10)
    # Find the index of maximum cross-correlation function
    maxind = np.argmax(cc)
    print("Cross-correlation function is maximized at dRV = ", rv[maxind], " km/s")
    v = rv[maxind]

    if plot:
        plt.plot(rv, cc, 'bp-')
        plt.plot(rv[maxind], cc[maxind], 'ro')
        plt.title('CCF')
        plt.show()
    return v

def rv_correction(filename, velocity):
    c = 299792.458 #km s-1
    hdulist = fits.open(filename)
    header = hdulist[0].header #read header
    if 'PRODCATG' in header:
        # ESO Product
        if eso_fits(hdulist) is None:
            hdulist.writeto(filename[:-5]+"_rv.fits", overwrite=True)
            hdulist.close()
        else:
             wave, flux = eso_fits(hdulist)
             start_wave = wave[0]
             step = wave[1] - wave[0]
             # New starting wavelenght and step
             relat_wave = start_wave * np.sqrt((1. - velocity/c)/(1. + velocity/c))
             relat_step = step * np.sqrt((1. - velocity/c)/(1. + velocity/c))
             #Change the header with the new values
             header['CRVAL1'] = relat_wave
             header['CDELT1'] = relat_step
             header['CD1_1'] = relat_step
             hdulist.writeto(filename[:-5]+"_rv.fits", overwrite=True)
             hdulist.close()

    elif 'CRVAL1' in header:
        start_wave = header['CRVAL1'] #initial wavelenght
        step = header['CDELT1'] #step in wavelenght

        # New starting wavelenght and step
        relat_wave = start_wave * np.sqrt((1. - velocity/c)/(1. + velocity/c))
        relat_step = step * np.sqrt((1. - velocity/c)/(1. + velocity/c))
        #Change the header with the new values
        header['CRVAL1'] = relat_wave
        header['CDELT1'] = relat_step
        header['CD1_1'] = relat_step
        hdulist.writeto(filename[:-5]+"_rv.fits", overwrite=True)
        hdulist.close()

    else:
        #Save results
        hdulist.writeto(filename[:-5]+"_rv.fits", overwrite=True)
        hdulist.close()
    return

def plot_rv(observed, rv):

    wave, flux, delta_l, start_wave = read_raw_observations(observed)
    plt.plot(wave, flux, label='original')
    wave_rv, flux, delta_l, start_wave = read_raw_observations(observed[:-5]+'_rv.fits')
    plt.plot(wave_rv, flux, label='RV cor=' + str(rv))
    plt.title(str(observed))
    plt.legend()
    plt.show()
    return

if __name__ == '__main__':

    args = argparse.ArgumentParser(description='Calculate for velocity shifts and correct them.')
    args.add_argument('observed', type=str, help='Spectrum')
    args = args.parse_args()

    print(args.observed)
    dw, df = read_spec(args.observed)

    template = 'data/HARPS.Archive_Sun-4_norm.fits'
    tw, tf = read_template(template)

    rv = velocity_shift(dw, df, tw, tf, plot=False)
    # Save fits file
    rv_correction(args.observed, round(rv, 2))
    # plot this file onced it is saved!
    # plot_rv(args.observed, rv)

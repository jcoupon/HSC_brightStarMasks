#!/usr/bin/env python

"""
Jean coupon - 2017
script to plot results for project Bright star masks
"""

import os, sys
import numpy as np
import re

import collections

import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

import halomodel

from scipy import interpolate
from scipy.special import erf
#from scipy import linalg

from astropy.io import ascii,fits
import astropy.wcs       as wcs
#import astropy.stats     as astats
#import photutils         as pu


import plot_utils

# from   astropy.cosmology import FlatLambdaCDM
# cosmo = FlatLambdaCDM(H0=72.0, Om0=0.258)
# h = cosmo.H0.value/100.0

# ----------------------------------------------------- #
# import buildMasks routines
# ----------------------------------------------------- #

sys.path.append(os.path.join(os.path.dirname(__file__)))
import buildMasks as BM

# ----------------------------------------------------- #
# global variables
# ----------------------------------------------------- #

EPS = 1.0e-8
COLOR = [ 'blue', 'green', 'red', 'orange', 'magenta', 'purple', 'lightblue', \
    'pink', 'Cyan', 'Brown', 'DarkRed', 'Indigo', 'DarkGreen', 'Salmon']

MARKER = ["o", "^", "s", "p"]

FIGX = 6.0
FIGY = 4.0

MARKER_SIZE = 2
FONT_SIZE = 16

# ----------------------------------------------------- #
# main
# ----------------------------------------------------- #

def main(args):
    function = getattr(sys.modules[__name__], args.option)(args, show=False)
    return

# ----------------------------------------------------- #
# Main functions
# ----------------------------------------------------- #


def skyMap(args, show=False):

    """ options
    """
    # bins_mag = np.linspace(0.0, 24.0, 50)
    title = args.title
    fileInName = args.input.split(",")
#    keys=['ra', 'dec', 'G_Gaia']
    keys=['ra', 'dec']

    """ read input file
    """
    cols, _ = getCols(fileInName[0], keys, select=args.select, dictionary=True)
    #cols = {'ra': [0.0], 'dec': [0.0] }

    """ initialise figure
    """
    # initialise figure
    fig = plt.figure(figsize=(FIGX*1.5, FIGY*1.7)); size = MARKER_SIZE

    fields={}
    fields['Spring'] = {'raMin': 123.0, 'raMax': 229.0, 'decMin': -8.0, 'decMax':8.0, 'position': [0.1, 0.20, 0.7, 0.2] }
    fields['Fall'] = {'raMin': 324.0, 'raMax': 47.0, 'decMin': -10.0, 'decMax':10.0, 'position': [0.1, 0.48, 0.5, 0.2] }
    fields['Northern'] = {'raMin': 195.0, 'raMax': 255.0, 'decMin': 40.0, 'decMax':50.0, 'position': [0.1, 0.7, 0.5, 0.3] }
    fields['Elais-N1'] = {'raMin': 238.0, 'raMax': 247.0, 'decMin': 52.0, 'decMax':58.0, 'position': [0.72, 0.42, 0.2, 0.25] }
    fields['Aegis'] = {'raMin': 211.0, 'raMax': 218.0, 'decMin': 50.0, 'decMax':55.0, 'position': [0.72, 0.74, 0.2, 0.22] }

    vmin = 1

    if "density" in args.labels:
        cmap = plt.cm.jet
        vmax = 461
        size = 0.5
        density =  1.0
        label=r'Star density ($G_\mathrm{Gaia}<18$) [deg$^{-2}$]'

    if "masked" in args.labels:
        cmap = plt.cm.viridis
        vmax = 14450.0
        size = 1.0
        density = 80000.0/100.0
        label=r'Masked fraction in %'

    ax = {}
    hb = {}
    for f in fields:
        ax[f], hb[f] = doMapPlot(fig, fields[f], cols['ra'], cols['dec'], size=size, vmin=vmin, vmax=vmax, cmap=cmap)
        ax[f].set_title(f)
        print f, hb[f].norm.vmin, hb[f].norm.vmax

    cmap_min = vmin / ( size * size ) / density
    cmap_max = vmax / ( size * size ) / density

    gradient = np.linspace(0.0, 1.0, 256)
    gradient = np.vstack((gradient, gradient))

    #cmap = plt.cm.jet

    ax2 = fig.add_axes([0.1, 0.08, 0.8, 0.04]) # xmin, ymin, dx, and dy
    ax2.imshow(gradient, aspect='auto', cmap=cmap, extent=[cmap_min,cmap_max,0,1])
    ax2.set_xlabel(label)

    ax2.get_yaxis().set_visible(False)

    fig.savefig(args.output)
    if show:
        plt.show()

    return

    #for (ra, dec, depth) in zip(data[key["ra"]], data[key["dec"]], data[key["depth"]]):
    #    plotPatch(pixscale, ra, dec, ax, cmap(norm(depth)))

    return


def doMapPlot(fig, field, ra, dec, size=0.5, vmin=1, vmax=500, cmap = plt.cm.viridis):

    """ Routine to draw the density of object and return
    a pointer to the axes

    INPUT
    fig: figure objects
    field: dictionnary that contains field limits
    and position on the image
    ra: object r.a.
    dec: object declination
    size: size of one cell on the side in degree

    OUTPUT
    axes object ax

    """

    if field['raMax'] < field['raMin']:
        twoPi = 360.0
    else:
        twoPi = 0.0

    width = field['raMax'] + twoPi - field['raMin']
    height = field['decMax'] - field['decMin']

    center=[(field['raMin']+field['raMax']+twoPi)/2.0, (field['decMin']+field['decMax'])/2.0]

    pixelScale = 0.1   # degrees per pixel
    w = createWCS(center, pixelScale)

    ax = fig.add_axes(field['position'], projection=w, title=args.title)
    setAxes(ax, field, pixelScale)

    Npix = [ax.get_xlim()[0] - ax.get_xlim()[1], ax.get_ylim()[1] - ax.get_ylim()[0]]

    if field['raMax'] < field['raMin']:
        select = ((ra < field['raMin']) | (ra > field['raMax'])) & (field['decMin'] < dec) & (dec < field['decMax'])
    else:
        select = (field['raMin'] < ra) & (ra < field['raMax']+twoPi) & (field['decMin'] < dec) & (dec < field['decMax'])

    [x, y] = zip(*w.wcs_world2pix(zip(ra[select], dec[select]), 1))

    extent=[ax.get_xlim()[0], ax.get_xlim()[1], ax.get_ylim()[0], ax.get_ylim()[1]]
    hb = ax.hexbin(x, y, gridsize=(int(Npix[0]/(size/pixelScale)), int(Npix[1]/(size/pixelScale))), mincnt=1, cmap=cmap, extent=extent, vmin=vmin, vmax=vmax)

    return ax, hb


def createWCS(center, pixscale):

    # Create a new WCS object.  The number of axes must be set
    # from the start
    w = wcs.WCS(naxis=2)

    # Set up an tangential projection
    # Vector properties may be set with Python lists, or Numpy arrays
    w.wcs.crpix = [1.0, 1.0]
    w.wcs.cdelt = [pixscale, pixscale]
    w.wcs.crval = center
    #w.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    w.wcs.ctype = ["RA---AIT", "DEC--AIT"]
    #w.wcs.ctype = ["RA---SIN", "DEC--SIN"]
    #w.wcs.ctype = ["RA---MOL", "DEC--MOL"]

    return w

def setAxes(ax, field, pixelScale):


    if field['raMax'] < field['raMin']:
        twoPi = 360.0
    else:
        twoPi = 0.0

    width = field['raMax'] + twoPi - field['raMin']
    height = field['decMax'] - field['decMin']

    ax.set_xlim(-width/pixelScale/2.0, width/pixelScale/2.0)
    ax.set_ylim(-height/pixelScale/2.0, height/pixelScale/2.0)

    ax.set_autoscale_on(False)
    ax.set_aspect('equal')

    lon=ax.coords[0]
    lat=ax.coords[1]

    lon.set_axislabel('R.A.')
    lat.set_axislabel('Dec. [deg]')

    lon.grid(color='black', alpha=0.5, linestyle='solid')
    lat.grid(color='black', alpha=0.5, linestyle='solid')

    lon.set_major_formatter('hh:mm')
    lat.set_major_formatter('dd')

    #lon.set_separator(('d', "'", '"'))
    #lat.set_separator("dms")
    lat.set_separator("")

    ax.invert_xaxis()

    return

def CAMIRA(args, show=False):

    """ Options
    """
    fileInName = args.input.split(",")
    labels = args.labels.split(",")

    #title = os.path.basename(args.input)
    title = args.title


    #fig = plt.figure(figsize=(FIGX, FIGY)); ax = plt.gca()

    fig, ax = plt.subplots(2, 1, sharex=True, sharey=True, figsize=(FIGX*0.8, FIGX*0.8)); size = MARKER_SIZE


    for i,f in enumerate(fileInName):

        if i == 0:
            plot_utils.axes(ax[i], '', r'$\frac{\^N_\mathrm{m}-\^N_\mathrm{m, w/masks}}{\^N_\mathrm{m, w/masks}}$', [0.0, 1.2], [-0.5, 1.0], title=title, fontsize = FONT_SIZE)
        else:
            plot_utils.axes(ax[i], r'Cluster redshift' , r'$\frac{\^N_\mathrm{m}-\^N_\mathrm{m, w/masks}}{\^N_\mathrm{m, w/masks}}$', [0.0, 1.2], [-0.5, 1.0], title=title, fontsize = FONT_SIZE)


        ax[i].plot([0.0, 1.20], [0.0, 0.0], lw=size, color='black', ls=':')

        data = ascii.read("{0:s}".format(f))

        z = data['z_nomask']
        eta = (data['N_nomask']-data['N_mask'])/data['N_mask']

        if "sirius" in f:
            alpha = 1.0
            color='red'
        else:
            alpha = 1.0
            color='blue'

        plot_utils.markers(ax[i], z, eta, None, size/3.0, color, labels[i], alpha=alpha)


        ax[i].locator_params(axis='y', nbins=4)


        if not args.nolegend:
            ax[i].legend( frameon=False, numpoints=1, loc='upper right')

    fig.set_tight_layout(True)
    fig.savefig(args.output)
    if show:
        plt.show()

    return



    return

def plotWoftheta(args, show=False):


    """ Options
    """
    fileInName = args.input.split(",")
    labels = args.labels.split(",")

    #title = os.path.basename(args.input)
    title = args.title

    fig = plt.figure(figsize=(FIGX, FIGY)); ax = plt.gca(); size = MARKER_SIZE

    plot_utils.axes(ax, r'$\theta$ [deg]', r'$w(\theta)$', [5.e-5, 3.e0], [4e-5, 8.e1], title=title, xlog=True, ylog=True, xexp=True, yexp=True)

    IC=0.0

    for i,f in enumerate(fileInName):
        data = ascii.read("{0:s}".format(f), format="no_header")
        select = (data['col2'] > 0.0) & (data['col3'] > 0.0)

        if "CFHTLenS" in f:
            ax.plot(data['col1'][select], (data['col2'][select]+IC), lw=size, color=COLOR[i], label=labels[i])
            #ax.fill_between(data['col1'][select], data['col2'][select]-data['col3'][select], data['col2'][select]+data['col3'][select], color=COLOR[i], alpha=0.2, label=labels[i])
        else:

            alpha = 1.0
            if "oldMask" in f or "insideNewMask" in f:
                alpha = 0.5

            plot_utils.markers(ax, data['col1'][select], (data['col2'][select]+IC), data['col3'][select], size, COLOR[i-1], labels[i], alpha=alpha)

    if not args.nolegend:
        plt.legend( frameon=False, numpoints=1, loc='lower left')

    fig.set_tight_layout(True)
    fig.savefig(args.output)
    if show:
        plt.show()

    return


def plotGalFrac(args, show=False):

    """ options
    """
    bins_mag = np.linspace(5.0, 24.0, 100.0)
    title = args.title

    fileInName = args.input.split(",")

    if args.select is not None:
        sel = args.select.split(",")
    else:
        sel = ["" for i in len(keys)]


    """ read input files
    """

    # HSC galaxies
    cols, _ = getCols(fileInName[0], ['G_Gaia', 'extended_HSC', 'extended_SDSS'], select=sel[0], dictionary=True)
    extended = cols['extended_HSC'] == 1
    all_hist, bin_edges = np.histogram(cols['G_Gaia'], bins=bins_mag)
    extended_hist, bin_edges = np.histogram(cols['G_Gaia'][extended] , bins=bins_mag)
    nonZero = (all_hist > 0) & (extended_hist > 0)

    galFrac_from_HSC = np.array(extended_hist[nonZero], dtype='f')/np.array(all_hist[nonZero], dtype='f')
    galFrac_from_HSC_err = np.sqrt(np.array(extended_hist[nonZero], dtype='f'))/np.array(all_hist[nonZero], dtype='f')
    bins_HSC = ((bin_edges[1:]+bin_edges[:-1])/2.0)[nonZero]

    # SDSS galaxies in HSC sample
    extended = cols['extended_SDSS'] == 1
    all_hist, bin_edges = np.histogram(cols['G_Gaia'], bins=bins_mag)
    extended_hist, bin_edges = np.histogram(cols['G_Gaia'][extended] , bins=bins_mag)
    nonZero = (all_hist > 0) & (extended_hist > 0)

    galFrac_from_HSC_SDSS = np.array(extended_hist[nonZero], dtype='f')/np.array(all_hist[nonZero], dtype='f')
    galFrac_from_HSC_SDSS_err = np.sqrt(np.array(extended_hist[nonZero], dtype='f'))/np.array(all_hist[nonZero], dtype='f')
    bins_HSC_SDSS = ((bin_edges[1:]+bin_edges[:-1])/2.0)[nonZero]

    # SDSS galaxies
    cols, _ = getCols(fileInName[1], ['G_Gaia', 'extended_SDSS'], select=sel[1], dictionary=True)
    extended = cols['extended_SDSS'] == 1
    all_hist, bin_edges = np.histogram(cols['G_Gaia'], bins=bins_mag)
    extended_hist, bin_edges = np.histogram(cols['G_Gaia'][extended] , bins=bins_mag)
    nonZero = (all_hist > 0) & (extended_hist > 0)

    galFrac_from_SDSS = np.array(extended_hist[nonZero], dtype='f')/np.array(all_hist[nonZero], dtype='f')
    galFrac_from_SDSS_err = np.sqrt(np.array(extended_hist[nonZero], dtype='f'))/np.array(all_hist[nonZero], dtype='f')
    bins_SDSS = ((bin_edges[1:]+bin_edges[:-1])/2.0)[nonZero]



    """ initialise figure
    """
    fig = plt.figure(figsize=(FIGX, FIGY)); ax = plt.gca(); size = MARKER_SIZE
    plot_utils.axes(ax, r'$G_\mathrm{Gaia}$', r'Galaxy fraction', [10.0, 22.0], [0.00, 0.30], ylog=False, title=title)

    ax.fill_between(bins_HSC, galFrac_from_HSC-galFrac_from_HSC_err, galFrac_from_HSC+galFrac_from_HSC_err, color=COLOR[0], alpha=0.2, label='HSC-SSP extended')
    ax.plot(bins_HSC, galFrac_from_HSC, color=COLOR[0], lw=size)

    ax.fill_between(bins_SDSS, galFrac_from_SDSS-galFrac_from_SDSS_err, galFrac_from_SDSS+galFrac_from_SDSS_err, color=COLOR[1], alpha=0.2, label='SDSS extended')
    ax.plot(bins_SDSS, galFrac_from_SDSS, color=COLOR[1], lw=size)

    ax.plot([18.0, 18.0], [0.0, 0.30], lw=size, ls='--', color='black')


    # ax.fill_between(bins_HSC_SDSS, galFrac_from_HSC_SDSS-galFrac_from_HSC_SDSS_err, galFrac_from_HSC_SDSS-galFrac_from_HSC_SDSS_err, color=COLOR[2], alpha=0.2, label='SDSS extended (in HSC-SSP)')
    # ax.plot(bins_HSC_SDSS, galFrac_from_HSC_SDSS, color=COLOR[2], lw=size)

    if not args.nolegend:
        plt.legend(frameon=False, numpoints=1, ncol=1, loc='upper left')


    fig.set_tight_layout(True)

    fig.savefig(args.output)
    if show:
        plt.show()


    return


def plotMagDiff(args, show=False):

    """ options
    """
    # bins_mag = np.linspace(0.0, 24.0, 50)
    title = args.title
    fileInName = args.input.split(",")
    keys_mags=['G_Gaia', 'phot_g_mean_mag']

    """ read input file
    """
    cols, _ = getCols(fileInName[0], keys_mags, select=args.select, dictionary=True)

    """ initialise figure
    """
    fig = plt.figure(figsize=(FIGX, FIGX)); ax = plt.gca(); size = MARKER_SIZE
    plot_utils.axes(ax, r'Gaia', r'Tycho-2 (emulated)', [3.0, 14.0], [3.0, 14.0], title=title)

    hb = ax.hexbin(cols[keys_mags[0]], cols[keys_mags[1]], gridsize=80, mincnt=1, cmap='viridis')

    x = np.linspace(0.0, 20.0, 50)
    ax.plot(x, x, color='black')


    fig.set_tight_layout(True)

    fig.savefig(args.output)
    if show:
        plt.show()

    return





def plotMagDist(args, show=False):

    """ options
    """
    #bins_mag = np.linspace(0.0, 22.5, 200)
    bins_mag = np.linspace(0.0, 22.5, 100)
    title = args.title
    fileInName = args.input.split(",")
    keys = args.keys.split(",")
    if args.select is not None:
        sel = args.select.split(",")
    else:
        sel = ["" for i in range(len(keys))]

    if args.labels is None:
        labels = fileInName
    else:
        labels = args.labels.split(",")

    """ read data
    """
    data = {}
    for l,f,k,s in zip(labels, fileInName, keys, sel):
        if s == "":
            s = None
        data[l], _ = getCols(f, [k], select=s)

    """ initialise figure
    """
    fig = plt.figure(figsize=(FIGX, FIGY)); ax = plt.gca(); size = MARKER_SIZE
    plot_utils.axes(ax, r'$G_\mathrm{Gaia}$', r'N', [8.0, 22.0], [0.5, 1e4], ylog=True, yexp=True, title=title)
    # plot_utils.axes(ax, r'$i$', r'N', [14.0, 23.0], [0.5, 1e4], ylog=True, yexp=True, title=title)

    """ loop over input files
    """
    ymax = 0.0
    for i,l in enumerate(labels):
        #hist, bin_edges = np.histogram(data[f], bins=50, range=[0.0, 20.0], density=False)
        #bins = (bin_edges[1:]+bin_edges[:-1])/2.0
        #ax.fill_between(bins, 0.0*hist, hist, color=COLOR[i], alpha=0.2, label=f)

        if i == 2:
            hist, bin_edges, patches = ax.hist(data[l], bins_mag, histtype='step', color='black', fill=False, label=labels[i], lw=size)
        else:
            hist, bin_edges, patches = ax.hist(data[l], bins_mag, histtype='step', color=COLOR[i], alpha=0.3, fill=True, label=labels[i])

        ymax=max(ymax, np.max(hist))

    ax.set_ylim([0.0, ymax*1.2])

    ax.locator_params(axis='x', nbins=5)

    #ax.set_ylim([0.0, ymax*1.2])

    if not args.nolegend:
        plt.legend(frameon=False, numpoints=1, ncol=1, loc='upper left')

    fig.set_tight_layout(True)

    fig.savefig(args.output)
    if show:
        plt.show()

    return

def plotSaturVsSeeing(args, show=False):

    """ options
    """
    bins_seeing = np.linspace(0.3, 1.3, 15)
    title = args.title
    fileInName = args.input.split(",")

    filters = ['g', 'r', 'i', 'z', 'Y']
    keys_PSF=['gpsf_pix', 'rpsf_pix', 'ipsf_pix', 'zpsf_pix', 'ypsf_pix']
    keys_mags=['gmag_psf', 'rmag_psf', 'imag_psf', 'zmag_psf', 'ymag_psf']

    """ read input file
    """
    cols, _ = getCols(fileInName[0], keys_PSF+keys_mags, select=args.select, dictionary=True)

    """ convert moments to FWHM
    """
    for k in keys_PSF:
        cols[k] *= 2.0 * np.sqrt(2) # object db is in arcsec


    """ initialise figure
    """
    fig = plt.figure(figsize=(FIGX, FIGY)); ax = plt.gca(); size = MARKER_SIZE
    plot_utils.axes(ax, r'PSF size [arcsec]', r'PSF mag', [0.3, 1.40], [16.0, 22.0], title=title)

    """ loop over filters
    """
    for i,f in enumerate(filters):

        """ loop over seeing bins
        """
        x=[]
        lower=[]
        upper=[]
        median=[]
        for l,r in zip(bins_seeing[:-1], bins_seeing[1:]):
            select = (l < cols[keys_PSF[i]]) & (cols[keys_PSF[i]] < r) & (cols[keys_mags[i]] > 0.0)
            if select.any():
                x.append(np.mean(cols[keys_PSF[i]][select]))
                lower.append(np.percentile(cols[keys_mags[i]][select], 25.0))
                upper.append(np.percentile(cols[keys_mags[i]][select], 75.0))
                median.append(np.percentile(cols[keys_mags[i]][select], 50.0))

        print f, median[0]

        ax.fill_between(x, lower, upper, color=COLOR[i], alpha=0.2, label=f)
        ax.plot(x, median, color=COLOR[i], lw=size)


    """ legend
    """
    if not args.nolegend:
        plt.legend( frameon=False, numpoints=1, ncol=1)


    """ output file
    """
    fig.set_tight_layout(True)

    fig.savefig(args.output)
    if show:
        plt.show()

    return


def plotMaskRadius(args, show=False):

    from scipy.optimize import curve_fit

    title = args.title
    brightest = BM.MAG_BRIGHTEST
    limit = BM.MAG_LIMIT

    # Initialise figure
    fig = plt.figure(figsize=(FIGX, FIGY)); ax = plt.gca(); size = MARKER_SIZE

    plot_utils.axes(ax, r'$G_{\rm Gaia}$ [mag]', r'r [arcsec]', [2.0, 18.5], [3.0, 1.e4], yexp=True, ylog=True, title=title)
    ax.locator_params(axis='x', nticks=4)

    ax.fill_between([brightest, limit], [3.0, 3.0], [1.e4, 1.e4], color=COLOR[0], alpha=0.1)

    ax.text((brightest+limit)/2.0, 10.0, 'Extended \n halo',
            horizontalalignment='center',
            verticalalignment='center',
            fontsize=FONT_SIZE, color=COLOR[0])

    G_Gaia = np.linspace(brightest, 18.0, num=100)

    fileInName = args.input.split(",")

    r_Czakon = 200.0*pow(10.0, 0.25*(7.0 - G_Gaia)) + 12.0*pow(10.0, 0.05*(16.0 - G_Gaia))
    ax.plot(G_Gaia, r_Czakon, color=COLOR[1], lw=size, label='Sirius (old)', ls='--')


    maskRadiusPara = ascii.read(fileInName[1], header_start=-1)

    bright = maskRadiusPara["mask_type"] == "bright"
    faint = maskRadiusPara["mask_type"] == "faint"

    allStars = ascii.read(fileInName[0], header_start=-1)
    x = allStars['mag']
    y = allStars['radiusIsotropic']
    #ax.plot(x, y, color=COLOR[i], lw=size, label='Mask radius')
    plot_utils.markers(ax, x, y, y*0.0, size, COLOR[0], '')

    ax.plot(G_Gaia[G_Gaia < limit], BM.r_vs_mag(G_Gaia[G_Gaia < limit], maskRadiusPara["a"][bright], maskRadiusPara["b"][bright]), color=COLOR[0], lw=size, label='Arcturus (this work)')
    ax.plot(G_Gaia[G_Gaia > limit], BM.r_vs_mag(G_Gaia[G_Gaia > limit], maskRadiusPara["a"][faint], maskRadiusPara["b"][faint]), color=COLOR[0], lw=size, label='')

    bestSeeing = ascii.read(fileInName[2], header_start=-1)
    select = np.isfinite(bestSeeing['seeing'])
    x = bestSeeing['mag'][select]
    y = bestSeeing['radiusIsotropic'][select]
    #ax.plot(x, y, color=COLOR[i], lw=size, label='Mask radius')
    plot_utils.markers(ax, x, y, None, size, COLOR[2], '10% best PSF', alpha=1.0, marker='^', fillstyle='none')

    #ax.scatter(x, y, label='10% best seeing', color=COLOR[2], marker=MARKER[1],  fillstyle='none')

    print np.mean(bestSeeing['radiusIsotropic']/allStars['radiusIsotropic'])

    worstSeeing = ascii.read(fileInName[3], header_start=-1)
    select = np.isfinite(worstSeeing['seeing'])
    x = worstSeeing['mag'][select]
    y = worstSeeing['radiusIsotropic'][select]
    #ax.plot(x, y, color=COLOR[i], lw=size, label='Mask radius')
    plot_utils.markers(ax, x, y, None, size, COLOR[3], '10% worst PSF', alpha=1.0, marker='s', fillstyle='none')
    # ax.scatter(x, y, label='10% worst seeing', color=COLOR[3], marker=MARKER[0], s=5*size)

    print np.mean(worstSeeing['radiusIsotropic']/allStars['radiusIsotropic'])



    # plot_utils.markers(ax, x, data['radiusSpikes'], data['radiusIsotropic']*0.0, size, COLOR[2], 'Spikes')
    # plot_utils.markers(ax, x, data['radiusBleedTrail'], data['radiusBleedTrail']*0.0, size, COLOR[3], 'Bleedtrails')


    # ax.plot([limit, limit], [3.0, 1.e4], color='black', lw=size, ls='--', label='')

    ax.plot([brightest, brightest], [3.0, 1.e4], color='black', lw=size, ls=':', label='')

    if not args.nolegend:
        plt.legend( frameon=False, numpoints=1, ncol=1)

    fig.set_tight_layout(True)

    fig.savefig(args.output)
    if show:
        plt.show()

    return


#def r_vs_mag(m, a, b):
#
#    return a*np.exp(-m/b)


def plotSeeingDist(args, show=False):

    title = args.title

    # input data points
    fileInName = args.input.split(",")

    keys=['gpsf_pix', 'rpsf_pix', 'ipsf_pix', 'zpsf_pix', 'ypsf_pix']
    filtName = dict(zip(keys, ['g', 'r', 'i', 'z', 'Y']))

    PSF, _ = getCols(fileInName[0], keys, dictionary=True)

    # convert to FWHM
    for k,v in PSF.items():
        v *= 0.17 * 2.0 * np.sqrt(2) # randoms db is in pixel

    # Initialise figure
    fig = plt.figure(figsize=(FIGX, FIGY)); ax = plt.gca(); size = MARKER_SIZE

    plot_utils.axes(ax, r'PSF size [arcsec]', r'n', [0.3, 1.5], [0.0, 5.0], title=title)

    ymax=0.0
    for i,k in enumerate(keys):
        hist, bin_edges = np.histogram(PSF[k], bins=50, range=[0.3, 2.0], density=True)
        bins = (bin_edges[1:]+bin_edges[:-1])/2.0
        ax.fill_between(bins, 0.0*hist, hist, color=COLOR[i], alpha=0.2, label='')
        ax.plot(bins,  hist, color=COLOR[i],  lw=size, label=filtName[k])

        ymax=max(ymax, np.max(hist))

    ax.set_ylim([0.0, ymax*1.2])

    if not args.nolegend:
        plt.legend( frameon=False, numpoints=1, ncol=1)

    fig.set_tight_layout(True)

    fig.savefig(args.output)
    if show:
        plt.show()

    return



def plotStackedStar(args, show=False):
    """
    Plot fits image
    """
    import astropy.wcs as wcs


    """ Options
    """

    brightest = BM.MAG_BRIGHTEST
    limit = BM.MAG_LIMIT
    fileInName = args.input.split(",")

    """ read image and its properties
    """
    IMAGE_EXT = 1
    fileIn = fits.open(fileInName[0])
    image = fileIn[IMAGE_EXT].data
    w = wcs.WCS(fileIn[IMAGE_EXT].header)
    pixelScale = wcs.utils.proj_plane_pixel_scales(w)[0]
    fileIn.close()

    rmax = [float(r) for r in args.rmax.split(",")]
    mag = rmax[2]

    maskRadiusPara = ascii.read(fileInName[1], header_start=-1)
    bright = maskRadiusPara["mask_type"] == "bright"
    faint = maskRadiusPara["mask_type"] == "faint"

    radiusBright = BM.r_vs_mag(mag, maskRadiusPara["a"][bright][0], maskRadiusPara["b"][bright][0])
    radiusFaint = BM.r_vs_mag(mag, maskRadiusPara["a"][faint][0], maskRadiusPara["b"][faint][0])

    """ draw figure
    """
    fig = plt.figure(figsize=(FIGY, FIGY)); ax = plt.gca(); size = MARKER_SIZE

    vmin = 0.0
    vmax = 0.4

    # side = 2.0 * rmax[0]
    if mag < limit:
        side = radiusBright*1.2
    else:
        side = radiusFaint*2.0

    nx = ny = int(2.0*side/pixelScale/3600.0)
    image = crop(image, nx, ny)
    N = image.shape

    im = ax.imshow(image, cmap=plt.cm.viridis, interpolation='none', origin='lower',  vmin=vmin, vmax=vmax, aspect='auto')


    #radius_iso = rmax[0]/pixelScale/3600.0
    if mag < limit:
        radiusIsotropic = radiusBright/pixelScale/3600.0
    else:
        radiusIsotropic = radiusFaint/pixelScale/3600.0
    ax.add_artist(plt.Circle((N[0]/2.0, N[1]/2.0), radiusIsotropic, color='white', fill=False, linewidth=size))

    if mag > limit :
        radiusSpike = 1.50*radiusIsotropic
        h = 2.0*radiusSpike
        w = h * 0.1
        ax.add_artist(plt.Rectangle((N[0]/2.0 - w/2, N[1]/2.0-h/2), w, h, color='white', fill=False, linewidth=size))

        radiusBleedTrail = 1.50*radiusIsotropic
        w = 2.0*radiusBleedTrail
        h = w * 0.1
        ax.add_artist(plt.Rectangle((N[0]/2.0 - w/2, N[1]/2.0-h/2), w, h, color='white', fill=False, linewidth=size))

    ax.text(0.05*N[0], (1-0.1)*N[1], args.labels, fontsize=FONT_SIZE, color='white')
    ax.text(0.75*N[0], (1-0.1)*N[1], args.title, fontsize=FONT_SIZE, color='white')


    ax.add_artist(plt.Rectangle((0.05*N[0], 0.04*N[1]), radiusIsotropic, 0.0, color='white', fill=False, linewidth=size))


    # ax.plot([0.05*N[0], 0.05*N[0] + radiusIsotropic], [0.04*N[1], 0.04*N[1]], color='white', linestyle='-', linewidth=2)
    ax.text(0.05*N[0], 0.07*N[1], r"{0:.0f} arcsec".format(side), fontsize=FONT_SIZE, color='white')

    ax.set_xticklabels([])
    ax.set_yticklabels([])
    ax.axis('off')

    ax.set_rasterized(True)
    #fig.savefig(args.output,  dpi=60, bbox_inches='tight', pad_inches=0.05) #)
    fig.savefig(args.output,  dpi=200, bbox_inches='tight', pad_inches=0.05) #)


    return


def plotSourceDensity(args, show=False):

    from scipy.optimize import curve_fit

    title = args.title
    fileInName = args.input.split(",")
    rmax = float(args.rmax)

    if args.labels is None:
        label = fileInName
    else:
        label = args.labels.split(",")

    """ initialise figure
    """
    xmin = 3.0
    xmax = rmax

    #sys.stderr.write("{0:f} {1:f}".format(xmin, xmax))
    # return

    fig = plt.figure(figsize=(FIGX*0.7, FIGY*0.7)); ax = plt.gca(); size = MARKER_SIZE
    plot_utils.axes(ax, r'$r\,\,\mathrm{[arcsec]}$', r'Source density', [xmin, xmax], [0.0, 2.5], ylog=False, title=title)

    """ loop over input files
    """
    ymax = 0.0
    for i,(l,f) in enumerate(zip(label,fileInName)):

        data = ascii.read(f, header_start=-1)
        bins = data['r']
        histIsoptropic = data['nIsotropic']
        histBleedTrail = data['nBleedTrail']
        histSpikes = data['nSpikes']

        hist = histIsoptropic

        ax.fill_between(bins, 0.0*hist, hist, color=COLOR[i], alpha=0.2, label=l) #, drawstyle="steps")

        # first file = parent incompleteness
        if i == 0:
            histIsoptropicParent = histIsoptropic
            histBleedTrailParent = histBleedTrail
            histSpikesParent = histSpikes

            popt, pcov = curve_fit(BM.completness_func, bins, histIsoptropic, p0=[rmax/2.0, 5.0])
            ax.plot(bins, BM.completness_func(bins, popt[0], popt[1]), color=COLOR[i], lw=size, label='', ls='--')


        # second file = children overdensity
        if i == 1:
            histIsoptropicPrimary = histIsoptropic
            histBleedTrailPrimary = histBleedTrail
            histSpikesPrimary = histSpikes

    radiusIsotropic = BM.getRadius(bins, histIsoptropicParent, histIsoptropicPrimary, rmax, args.mag)
    radiusBleedTrail = BM.getRadius(bins, histBleedTrailParent, histBleedTrailPrimary, rmax, args.mag)
    radiusSpikeRadius = BM.getRadius(bins, histSpikesParent, histSpikesPrimary, rmax, args.mag)

    ax.arrow(radiusIsotropic, 1.70, 0.0, -0.20 , head_width=(xmax-xmin)*0.02, head_length=0.1, fc='black', ec='black')
    #ax.arrow(YSpikeRadius, 1.60, 0.0, -0.30 , head_width=(xmax-xmin)*0.02, head_length=0.1, fc='magenta', ec='magenta')

    if not args.nolegend:
        plt.legend(frameon=False, numpoints=1, ncol=1, loc='upper left')

    #ax.locator_params(axis='x', nticks=3)
    ax.set_xlim([xmin, xmax])

    ax.locator_params(axis='x', nbins=5)
    fig.subplots_adjust(bottom=0.18)
    fig.subplots_adjust(left=0.18)


    #fig.set_tight_layout(True)

    fig.savefig(args.output)
    if show:
        plt.show()

    return

    if not args.nolegend:
        plt.legend(frameon=False, numpoints=1, ncol=1, loc='upper left')

    #ax.locator_params(axis='x', nticks=3)
    ax.set_xlim([xmin, xmax])

    ax.locator_params(axis='x', nbins=5)
    fig.subplots_adjust(bottom=0.15)


    #fig.set_tight_layout(True)

    fig.savefig(args.output)
    if show:
        plt.show()

def PSFwGhost(args, show=False):


    rmin = 0.05
    rmax = 500
    cmax = 1.0
    cmin = 1.0E-12
    xlabel = 'r (arcsec)'
    ylabel = 'Contrast'

    r = np.arange(0, rmax, 0.05)
    k1 = kolmogorov(r, 0.25, 7, 0.86)
    k2 = kolmogorov(r, 0.25, 2, 0.14)
    k = []
    for i in range(len(k1)):
    	k.append(k1[i]+k2[i])

    ip = instpsf(r)

    cg1 = circ_ghost(r, 2.25*11.2, 3.53547E-08)
    cg2 = circ_ghost(r, 4.68*11.2, 4.1028E-09)
    cg3 = circ_ghost(r, 27.66*11.2, 1.1742E-09)
    cg4 = circ_ghost(r, 13.83*11.2, 9.40048E-10)
    cg5 = circ_ghost(r, 13.85*11.2, 4.68405E-10)
    cg6 = circ_ghost(r, 23.00*11.2, 3.39806E-10)
    cg7 = circ_ghost(r, 9.18*11.2, 2.13051E-10)
    cg8 = circ_ghost(r, 25.42*11.2, 1.39032E-10)
    cg9 = circ_ghost(r, 11.59*11.2, 1.33885E-10)
    cg10 = circ_ghost(r, 20.76*11.2, 4.17137E-11)

    t = []
    for i in range(len(k)):
    	t.append(k[i]+cg1[i]+cg2[i]+cg3[i]+cg4[i]+cg5[i]+cg6[i]+cg7[i]+cg8[i]+cg9[i]+cg10[i])


    fig = plt.figure(figsize=(FIGX, FIGX)); ax = plt.gca(); size = MARKER_SIZE

    plot_utils.axes(ax, r'$r\,\,\mathrm{[arcsec]}$', r'Contrast', [rmin, rmax], [cmin, cmax], ylog=True, xlog=True, yexp=True, title=args.title)

    #ax.yscale('log')

    ax.plot(r, ip, color='green', marker=".", markersize=0, linewidth=size, linestyle="--", label="Instrumental")
    ax.plot(r, k, color='blue', marker=".", markersize=0, linewidth=size, linestyle="-", label="Atmospheric")
    ax.plot(r, k1, color='blue', marker=".", markersize=0, linewidth=size, linestyle=":", label=r"Moffat ($\beta=7$)")
    ax.plot(r, k2, color='blue', marker=".", markersize=0, linewidth=size, linestyle="-.", label=r"Moffat ($\beta=2$)")
    ax.plot(r, cg1, color='red', marker=".", markersize=0, linewidth=size/2.0, linestyle=":", label="")
    ax.plot(r, cg2, color='red', marker=".", markersize=0, linewidth=size/2.0, linestyle=":", label="")
    ax.plot(r, cg3, color='red', marker=".", markersize=0, linewidth=size/2.0, linestyle=":", label="")
    ax.plot(r, cg4, color='red', marker=".", markersize=0, linewidth=size/2.0, linestyle=":", label="")
    ax.plot(r, cg5, color='red', marker=".", markersize=0, linewidth=size/2.0, linestyle=":", label="")
    ax.plot(r, cg6, color='red', marker=".", markersize=0, linewidth=size/2.0, linestyle=":", label="")
    ax.plot(r, cg7, color='red', marker=".", markersize=0, linewidth=size/2.0, linestyle=":", label="")
    ax.plot(r, cg8, color='red', marker=".", markersize=0, linewidth=size/2.0, linestyle=":", label="")
    ax.plot(r, cg9, color='red', marker=".", markersize=0, linewidth=size/2.0, linestyle=":", label="")
    ax.plot(r, cg10, color='red', marker=".", markersize=0, linewidth=size/2.0, linestyle=":", label="")
    ax.plot(r, t, color='black', marker=".", markersize=0, linewidth=size, linestyle="-", label="Total")
    ax.text(0.1, 3.53547E-08, 'CCD-WinS2', fontsize=FONT_SIZE/1.5)
    ax.text(0.1, 4.1028E-09, 'FilS2-FilS1', fontsize=FONT_SIZE/1.5)
    ax.text(0.1, 1.1742E-09, 'CCD-FilS1', fontsize=FONT_SIZE/1.5)
    #ax.text(0.1, 9.40048E-10, 'CCD-WinS1')
    ax.text(0.1, 5.40048E-10, 'CCD-WinS1', fontsize=FONT_SIZE/1.5)
    #ax.legend((sptype,), loc=0)

    #ax.xticks()
    #ax.yticks()
    #ax.xlim(rmin, rmax)
    #ax.ylim(cmin, cmax)
    #ax.xscale('log')
    #ax.yscale('log')

    #ax.locator_params(axis='y', nbins=5)


    if not args.nolegend:
        plt.legend(frameon=False, numpoints=1, ncol=1, loc='upper right')

    ax.set_rasterized(True)

    fig.set_tight_layout(True)
    fig.savefig(args.output, dpi=300)
    if show:
        plt.show()


    return

def kolmogorov0(beta):
	return (2.0**(1.0/beta)-1.0)*(beta-1.0)


def kolmogorov(r, hw, beta, f):
	return (2.0**(1.0/beta)-1.0)*(beta-1.0)*f/(1.0+(2.0**(1.0/beta)-1.0)*(r/hw)**2)**beta/kolmogorov0(beta)


def instpsf(r):
	sigma = 0.357/2.35
	return np.exp(-r**2/2.0/sigma/sigma)


def circ_ghost(r, rmax, cont):
	p = []
	for i in range(len(r)):
		if r[i] <= rmax:
			p.append(cont)
		else:
			v = cont*(kolmogorov(r[i]-rmax, 0.25, 7, 0.8)+kolmogorov(r[i]-rmax, 0.25, 2, 0.2))
			p.append(v)
	return p


# ----------------------------------------------------- #
# Utils
# ----------------------------------------------------- #


def getCols(fileInName, keys, select=None, array=False, dictionary=False):
    '''
    Return column arrays from fileInName (fits file)

    INPUT
    fileInName: fits file name
    keys: column names list or expression (e.g. "col1 - col2")
    selection: selection (string)
    array: if set returns 2D array (N,K)

    OUTPUT
    list of column arrays
    (optional) selection array

    WARNING
    Numeric column names cannot be used
    '''

    # first open the fits file and read the data
    fileIn = fits.open(fileInName)
    data = fileIn[1].data
    cols = fileIn[1].columns
    fileIn.close()

    if select is not None:
        cmd = ""
        for s in re.split('(\W+)', select):
            if s in cols.names:
                cmd += "data['"+s+"']"
            else:
                cmd += s
        sel = eval(cmd)
        str_select = "[sel]"
    else:
        str_select = ""

    result=[]
    # replace key value with command
    for k in keys:
        cmd = ""
        for s in re.split('(\W+)', k):
            if s in cols.names:
                cmd += "data"+str_select+"['"+s+"']"
            else:
                cmd += s

        # print cmd
        result.append(eval(cmd))

    if array:
        K = len(result)
        N = len(result[0])
        result_tmp = np.zeros((N,K))
        for d in range(K):
            result_tmp[:,d] = result[d]

        result = result_tmp

    if dictionary:
        result = dict(zip(keys, result))

    # return result
    if select is not None:
        return result, sel
    else:
        return result, None





def crop(image, Nx_cropped, Ny_cropped):
    """
    crop image according to scales
    (between 0 and 1)
    """

    (Nx, Ny) = image.shape

    dx = max(0, Nx-Nx_cropped)
    dy = max(0, Ny-Ny_cropped)

    return image[0+dx/2:Nx-dx/2, 0+dy/2:Ny-dy/2]


# ----------------------------------------------------- #
# Main
# ----------------------------------------------------- #



if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument('option', help="Which quantity to plot")
    parser.add_argument('-i', '--input', default=None, help='input file')
    parser.add_argument('-o', '--output', default=None, help='output file')
    parser.add_argument('-l', '--labels', default=None, help='parameters labels')
    parser.add_argument('-t', '--title', type=str,   default="",   help='Title')
    parser.add_argument('--nolegend', action='store_true', help='Do not display legend')
    parser.add_argument('-s', '--select', type=str,   default=None, help='Selection')
    parser.add_argument('-k', '--keys', default=None, help='Keys to plot')
    parser.add_argument('-rmax', type=str, default="5.0", help='Rmax for source density plot [arcsec]')
    parser.add_argument('-mag', type=float, default="9.0", help='Star magnitude')

    args = parser.parse_args()
    main(args)

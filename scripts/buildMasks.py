#!/usr/bin/env python

"""
Jean coupon - 2017
script to build masks from a star catalogue

usage: buildMasks.py [-h] [-i INPUT] [-o OUTPUT] [-k KEYS] [-rmax RMAX]
                     [-mag MAG] [-s SELECT] [-worstSeeing] [-bestSeeing]
                     [-basename BASENAME]
                     option

positional arguments:
  option                Which quantity to plot

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        input file
  -o OUTPUT, --output OUTPUT
                        output file
  -k KEYS, --keys KEYS  Keys
  -rmax RMAX            Rmax for source density plot [arcsec]
  -mag MAG              Star magnitude
  -s SELECT, --select SELECT
                        Selection
  -worstSeeing          Compute radius around worst seeing stars
  -bestSeeing           Compute radius around best seeing stars
  -basename BASENAME    Base name for region files


"""

from __future__ import print_function, division


import os
import sys
import errno
import datetime
import re
import collections
from io import open

import pickle
import numpy as np

from astropy.io import ascii,fits
from astropy.coordinates import SkyCoord
import astropy.wcs as wcs

from scipy.special import erf
from scipy.optimize import curve_fit
from astropy.table import Table, Column

# ----------------------------------------------------- #
# global variables
# ----------------------------------------------------- #

EPS = 1.0e-8

""" separate bright stars, dominated by the halo,
from faint stars
"""
MAG_LIMIT = 9.0

""" brightest star in HSC-SSP footprint
"""
MAG_BRIGHTEST = 2.57

NOW = datetime.datetime.now()

# ----------------------------------------------------- #
# main
# ----------------------------------------------------- #

def main(args):
    function = getattr(sys.modules[__name__], args.option)(args)
    return

# ----------------------------------------------------- #
# Main functions
# ----------------------------------------------------- #

def maskRadius(args, show=False):

    pixelSize = 0.17 # in arcsec
    # for the horizontal bleed trails in arsec
    if args.mag is not None:
        ymax = 0.1*r_vs_mag(args.mag, 700, 4.04)
    else:
        ymax = 5.0

    # for the vertical spikes arsec
    if args.mag is not None:
        xmax = 0.1*r_vs_mag(args.mag, 700, 4.04)
    else:
        xmax = 2.0

    """ options
    """
    bins_r = np.linspace(3.0, args.rmax, 50)
    fileOutName = args.output.split(",")
    fileInName = args.input.split(",")
    keys = args.keys.split(",")

    if args.select is not None:
        sel = args.select.split(",")
    else:
        sel = [None] * len(fileInName)

    """ loop over input files
    """
    for fi,fo,k,s in zip(fileInName, fileOutName, keys, sel):

        [iPSF], _ = getCols(fi, ['user_ipsf_FWHM'], select=s)
        if args.bestSeeing:
            select = (iPSF > 0.0) & (np.isfinite(iPSF))
            bestSeeing10 = (iPSF+EPS < np.percentile(iPSF[select], 10.0))
            meanSeeing = np.mean(iPSF[(bestSeeing10) & (iPSF > 0.0)])
        elif args.worstSeeing:
            select = (iPSF > 0.0) & (np.isfinite(iPSF))
            worstSeeing10 = (iPSF+EPS > np.percentile(iPSF[select], 90.0))
            meanSeeing = np.mean(iPSF[(worstSeeing10) & (iPSF > 0.0)])
        else:
            meanSeeing = np.mean(iPSF[iPSF > 0.0])

        if False:

            [r, x, y, imag_psf, ymag_psf, iPSF], _ = getCols(fi, [k, 'x', 'y', 'imag_psf', 'ymag_psf', 'user_ipsf_FWHM'], select=s)

            #print np.mean(imag_psf[np.isfinite(imag_psf)])

            if args.bestSeeing:
                r = r[bestSeeing10]
                x = x[bestSeeing10]
                y = y[bestSeeing10]
                imag_psf = imag_psf[bestSeeing10]
                ymag_psf = ymag_psf[bestSeeing10]
                iPSF = iPSF[bestSeeing10]


            if args.worstSeeing:
                r = r[worstSeeing10]
                x = x[worstSeeing10]
                y = y[worstSeeing10]
                imag_psf = imag_psf[worstSeeing10]
                ymag_psf = ymag_psf[worstSeeing10]
                iPSF = iPSF[worstSeeing10]


            x = x * 3600.0 * pixelSize
            y = y * 3600.0 * pixelSize

            """ mean density in outer radius
            """
            # all objects
            rmax = min(max(r), args.rmax)
            select = (r > 0.8*rmax) & (r < rmax)
            area = np.pi*(pow(rmax, 2.0) - pow(0.8*rmax, 2.0))
            nsources = float(len(r[select]))/area

            # bright objects in the y band
            nsourcesYBright = nsources * float(len(r[imag_psf-ymag_psf > 1.0]))/float(len(r))

            """ histogram of isotropic density
            """
            histIsoptropic, bin_edges = np.histogram(r, bins=bins_r, density=False)
            bins = (bin_edges[:-1]+bin_edges[1:])/2.0
            area = np.pi*(pow(bin_edges[1:], 2.0)-pow(bin_edges[:-1], 2.0))

            histIsoptropic = histIsoptropic/(area*nsources)

            """ histogram of horizontal density
            """
            select = (abs(x) < rmax) & (abs(y) < ymax)

            histBleedTrail, bin_edges = np.histogram(abs(x[select]), bins=bins_r, density=False)
            areaBleedTrail = 4.0*ymax*(bin_edges[1:]-bin_edges[:-1])

            histBleedTrail = histBleedTrail/(areaBleedTrail*nsources)

            """ histogram of vertical density (for y-band spikes)
            """
            select = (abs(x) < xmax) & (abs(y) < rmax) & (imag_psf-ymag_psf > 1.0)

            histSpikes, bin_edges = np.histogram(abs(y[select]), bins=bins_r, density=False)
            areaSpikes = 4.0*xmax*(bin_edges[1:]-bin_edges[:-1])

            histSpikes = histSpikes/(areaSpikes*nsourcesYBright)

            """ write output file
            """

            out = Table([bins, histIsoptropic, histBleedTrail, histSpikes], names=['r', 'nIsotropic', 'nBleedTrail', 'nSpikes'])
            ascii.write(out, fo, format="commented_header", overwrite=True)

        else:

            data = ascii.read(fo, header_start=-1)
            bins = data['r']
            histIsoptropic = data['nIsotropic']
            histBleedTrail = data['nBleedTrail']
            histSpikes = data['nSpikes']

        # first file = parent incompleteness
        if fo == fileOutName[0]:
            histIsoptropicParent = histIsoptropic
            histBleedTrailParent = histBleedTrail
            histSpikesParent = histSpikes

        # second file = children overdensity
        if fo == fileOutName[1]:
            histIsoptropicPrimary = histIsoptropic
            histBleedTrailPrimary = histBleedTrail
            histSpikesPrimary = histSpikes

    radiusIsotropic = getRadius(bins, histIsoptropicParent, histIsoptropicPrimary, args.rmax, args.mag)
    radiusBleedTrail = getRadius(bins, histBleedTrailParent, histBleedTrailPrimary, args.rmax, args.mag)
    radiusSpikeRadius = getRadius(bins, histSpikesParent, histSpikesPrimary, args.rmax, args.mag)

    print(meanSeeing, radiusIsotropic, radiusBleedTrail, radiusSpikeRadius)

    return

def getRadius(bins, nParent, nPrimary, rmax, mag):

    incompleteRadius = 0.0
    try:
        popt, pcov = curve_fit(completness_func, bins, nParent, p0=[rmax/2.0, 5.0])
        completeness = completness_func(bins, popt[0], popt[1])
        incompleteRadius = bins[completeness > 0.95][0]
    except:
        pass

    overdensityRadius = 0.0
    try:
        overdensityRadius = bins[nPrimary > 1.20][-1]
    except:
        pass

    if mag < MAG_LIMIT:
        r = max(incompleteRadius, overdensityRadius)
    else:
        r = incompleteRadius

    return r


def get_tract_coord(args, show=False):

    # read input tract list
    with open(args.input, 'r') as f:
        tract_indices = [int(t) for t in f.readlines()]

    # remove faulty tracts
    for index in  [6248, 6478, 17865, 17975]:
        if index in tract_indices:
            tract_indices.remove(index)

    sys.stdout.write("Importing LSST pipeline libraries...")
    sys.stdout.flush()
    import lsst.daf.persistence as dafPersist
    import lsst.afw.cameraGeom as camGeom
    import lsst.afw.coord as afwCoord
    import lsst.afw.geom as afwGeom
    import lsst.afw.image as afwImage
    sys.stdout.write("Done\n")

    skyMap = pickle.load( open(args.skyMap, "rb" ) )

    coord = {'index':[], 'ra':[], 'dec':[]}
    for index in tract_indices:

        tract = skyMap.generateTract(index)

        for patch in tract:
            header = build_header(tract, patch)
            print(repr(header))
            return

        tract_ra, tract_dec = \
            tract.getCtrCoord().getRa().asDegrees(), \
            tract.getCtrCoord().getDec().asDegrees()

        coord['index'].append(index)
        coord['ra'].append(tract_ra)
        coord['dec'].append(tract_dec)

    Table(coord).write(args.output, overwrite=True)

    return


def build_header(tract, patch):

    from astropy import wcs

    tract_wcs = tract.getWcs()
    # print(patch.getOuterBBox())

    # Create a new WCS object.  The number of axes must be set
    # from the start
    w = wcs.WCS(naxis=2)

    # TODO
    # Set up an "Airy's zenithal" projection
    # Vector properties may be set with Python lists, or Numpy arrays
    # w.wcs.crpix = [-234.75, 8.3393]
    # w.wcs.cdelt = np.array([-0.066667, 0.066667])
    # w.wcs.crval = [0, -90]
    # w.wcs.ctype = ["RA---AIR", "DEC--AIR"]
    # w.wcs.set_pv([(2, 1, 45.0)])

    # Now, write out the WCS object as a FITS header
    return w.to_header()


def MakeMasterFile(args, show=False):

    """ make a master region file
    """

    fileInName = args.input.split(",")

    """ star catalogue
    """
    stars, _ = getCols(fileInName[0], ['source_id', 'ra', 'dec', 'G_Gaia'], dictionary=True)
    Nstars = len(stars['source_id'])

    bright = stars['G_Gaia'] < MAG_LIMIT
    faint = stars['G_Gaia'] >= MAG_LIMIT

    rReg = np.zeros(Nstars)

    """ parameters to set the regions radius
    """
    maskRadiusPara = ascii.read(fileInName[1], header_start=-1)
    a_bright =  maskRadiusPara["a"][maskRadiusPara["mask_type"] == "bright"][0]
    b_bright =  maskRadiusPara["b"][maskRadiusPara["mask_type"] == "bright"][0]
    a_faint =  maskRadiusPara["a"][maskRadiusPara["mask_type"] == "faint"][0]
    b_faint =  maskRadiusPara["b"][maskRadiusPara["mask_type"] == "faint"][0]

    rReg[bright] = r_vs_mag(stars['G_Gaia'][bright], a_bright, b_bright)
    rReg[faint] = r_vs_mag(stars['G_Gaia'][faint], a_faint, b_faint)

    rRegDeg = rReg/3600.0

    ra = stars['ra']
    dec = stars['dec']
    ID = stars['source_id']
    mag = stars['G_Gaia']

    fileOutName = args.output
    tract = -1
    patch =[-1,-1]

    writeMaskFile(ra, dec, rRegDeg, ID, fileOutName, "HSC-G", tract, patch, mag)

    return



def fitMaskRadius(args, show=False):

    from scipy.optimize import curve_fit

    G_Gaia = np.linspace(0.0, 20.0, num=100)

    fileInName = args.input

    mask_type = [ "bright", "faint" ]
    a = []
    b = []
    data = ascii.read(fileInName, header_start=-1)
    x = data['mag']
    y = data['radiusIsotropic']

    # bright part
    popt, pcov = curve_fit(r_vs_mag, x[x<MAG_LIMIT], y[x<MAG_LIMIT], p0=[8.e2, 4.0])
    a.append(popt[0])
    b.append(popt[1])

    # faint part
    popt, pcov = curve_fit(r_vs_mag, x[x>MAG_LIMIT], y[x>MAG_LIMIT], p0=[8.e2, 4.0])
    a.append(popt[0])
    b.append(popt[1])

    out = Table([mask_type, a, b], names=['mask_type', 'a', 'b'])
    ascii.write(out, args.output, format="commented_header", comment="# Parameters to draw regions for masking.\n# For bright and faint, the circle radius is defined as \n# r [arcsec] = a*np.exp(-mag/b) \n# ")

    return

def makePatchFiles(args, verbose=0):
    """ Makes patch-level region files given a
    skyMap, star catalogue and a recipe to draw
    each mask size according to the star
    brightness

    INPUT
    -i SKYMAP,STARCAT,RADIUSFILE

    OUTPUT
    -o directory where to put the region files

    """

    """ LSST pipeline libraries
    """
    sys.stderr.write("Importing LSST pipeline libraries...")
    import lsst.daf.persistence as dafPersist
    import lsst.afw.cameraGeom as camGeom
    import lsst.afw.coord as afwCoord
    import lsst.afw.geom as afwGeom
    import lsst.afw.image as afwImage
    sys.stderr.write("Done\n")

    fileInName = args.input.split(",")
    dirOutName = args.output

    filters=["HSC-G", "HSC-R", "HSC-I", "HSC-Z", "HSC-Y"]

    """ skyMap
    """
    sys.stderr.write("Loading sky map: {0:s}...".format(fileInName[0]))
    skyMap = pickle.load( open(fileInName[0], "rb" ) )
    sys.stderr.write("Done\n")

    NTracts = len(skyMap)

    """ star catalogue
    """
    sys.stderr.write("Reading star catalogue: {0:s}...".format(fileInName[1]))
    stars, _ = getCols(fileInName[1], ['source_id', 'ra', 'dec', 'G_Gaia'], dictionary=True)
    Nstars = len(stars['source_id'])

    starsTree = dataTree(stars['ra'], stars['dec'])

    bright = stars['G_Gaia'] < MAG_LIMIT
    faint = stars['G_Gaia'] >= MAG_LIMIT
    sys.stderr.write("Done\n")

    """ mask radius
    """
    sys.stderr.write("Computing masks radius from: {0:s}...".format(fileInName[2]))
    maskRadiusPara = ascii.read(fileInName[2], header_start=-1)
    a_bright =  maskRadiusPara["a"][maskRadiusPara["mask_type"] == "bright"][0]
    b_bright =  maskRadiusPara["b"][maskRadiusPara["mask_type"] == "bright"][0]
    a_faint =  maskRadiusPara["a"][maskRadiusPara["mask_type"] == "faint"][0]
    b_faint =  maskRadiusPara["b"][maskRadiusPara["mask_type"] == "faint"][0]
    sys.stderr.write("Done\n")

    rMask = np.zeros(Nstars)
    rMask[bright] = r_vs_mag(stars['G_Gaia'][bright], a_bright, b_bright) # in arcsec
    rMask[faint] = r_vs_mag(stars['G_Gaia'][faint], a_faint, b_faint) # in arsec

    rMask /= 3600.0

    """ largest possible radius in star catalogue
    """
    starMaxRadius = max(rMask)

    sys.stderr.write("Largest mask radius: {0:f} deg\n".format(starMaxRadius))

    sys.stderr.write("\r" + "Writing region files: {0:d}/{1:d} tracts done".format(0, NTracts))
    sys.stderr.flush()

    #sys.stderr.write("Writing region files\n".format(starMaxRadius))
    for i, tract in enumerate(skyMap):

        # if tract.getId() != 9376:
        #    continue

        wcs = tract.getWcs()
        pixelScale = wcs.pixelScale().asDegrees()

        """ maximum radii along the diagonal of the tract and patches
        """
        tractMaxRadius =  max(tract.getBBox().getMaxX(),tract.getBBox().getMaxY())*pixelScale*np.sqrt(2.) / 2.
        patchMaxRadius = (max(tract.getPatchInnerDimensions()) + tract.getPatchBorder())*pixelScale*np.sqrt(2.) / 2.

        """ find stars around tract
        """
        raTract, decTract = tract.getCtrCoord().getRa().asDegrees(), tract.getCtrCoord().getDec().asDegrees()
        starsTract = findNeighbors(starsTree, raTract, decTract, tractMaxRadius+starMaxRadius)

        Nstars = len(stars['ra'][starsTract])

        # print "\n", Nstars
        if Nstars == 0:
            if (i+1)%100 == 0:
                sys.stderr.write("\r" + "Writing region files: {0:d}/{1:d} tracts done".format(i+1, NTracts))
                sys.stderr.flush()
            continue

        mkdir_p(dirOutName+"/patches/{0:d}".format(tract.getId()))
        mkdir_p(dirOutName+"/tracts")

        starsRa = stars['ra'][starsTract]
        starsDec = stars['dec'][starsTract]
        starsRMasks = rMask[starsTract]
        starsID = stars['source_id'][starsTract]
        starsMag = stars['G_Gaia'][starsTract]

        """ write tract file
        """
        for f in filters:
            fileOutName=dirOutName+"/tracts/"+args.basename+"-{0:d}-{1:s}.reg".format(tract.getId(), f)
            if f == "HSC-I" :
                writeMaskFile(starsRa, starsDec, starsRMasks, starsID, fileOutName, f, tract.getId(), [-1,-1], starsMag)
            else:
                ln_sf(args.basename+"-{0:d}-{1:s}.reg".format(tract.getId(), "HSC-I"), fileOutName)


        """ loop over patches
        """
        for patch in tract:
            bbox = patch.getOuterBBox()
            p = afwGeom.Point2D((bbox.getMinX() + bbox.getMaxX())/2.0, (bbox.getMinY() + bbox.getMaxY())/2.0)
            center = wcs.pixelToSky(p).toIcrs()
            raPatch, decPatch = center.getRa().asDegrees(), center.getDec().asDegrees()

            D = distAngSpherDeg(raPatch, decPatch, starsRa, starsDec)

            inPatch = ( D < patchMaxRadius + 1.5*starsRMasks )



            if len(starsDec[inPatch]) == 0:
                continue

            for f in filters:
                fileOutName=dirOutName+"/patches/{0:d}/{1:s}-{0:d}-{2:d},{3:d}-{4:s}.reg".format(tract.getId(), args.basename, patch.getIndex()[0], patch.getIndex()[1], f)
                if f == "HSC-I" :
                    writeMaskFile(starsRa[inPatch], starsDec[inPatch], starsRMasks[inPatch], starsID[inPatch], fileOutName, f, tract.getId(), patch.getIndex(),starsMag[inPatch])
                else:
                    ln_sf(args.basename+"-{0:d}-{1:d},{2:d}-{3:s}.reg".format(tract.getId(), patch.getIndex()[0], patch.getIndex()[1], "HSC-I"), fileOutName)

        if (i+1)%100 == 0:
            sys.stderr.write("\r" + "Writing region files: {0:d}/{1:d} tracts done".format(i+1, NTracts))
            sys.stderr.flush()

    sys.stderr.write("\r" + "Writing region files: {0:d}/{1:d} tracts done\n".format(NTracts, NTracts))
    sys.stderr.flush()

    return

# ----------------------------------------------------- #
# Utils
# ----------------------------------------------------- #

def r_vs_mag(m, a, b):

    return a*np.exp(-m/b)

def completness_func(x, x0, dx):
    return 1.0/2.0 * (1.0+erf((x-x0)/dx))


def ln_sf(source_file, target_file):
    try:
        os.symlink(source_file, target_file)
    except OSError as e:
        if e.errno == errno.EEXIST:
            os.remove(target_file)
            os.symlink(source_file, target_file)

def mkdir_p(path):
    try:
        os.makedirs(path)
    except OSError as exc: # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else: raise


def dataTree(ra, dec):
    from scipy import spatial

    # Spherical to euclidean geometry -> 3D
    x = np.cos(dec*np.pi/180.0)*np.cos(ra*np.pi/180.0)
    y = np.cos(dec*np.pi/180.0)*np.sin(ra*np.pi/180.0)
    z = np.sin(dec*np.pi/180.0)

    return spatial.KDTree(list(zip(x, y, z)))

def findNeighbors(tree, ra, dec, radius_in_deg):

    # Spherical to euclidean geometry -> 3D
    x = np.cos(dec*np.pi/180.0)*np.cos(ra*np.pi/180.0)
    y = np.cos(dec*np.pi/180.0)*np.sin(ra*np.pi/180.0)
    z = np.sin(dec*np.pi/180.0)

    return tree.query_ball_point((x,y,z), radius_in_deg*np.pi / 180.0)

def distAngSpherDeg(ra1, dec1, ra2, dec2):

    sin2_ra  = 0.5*(1.0 - np.cos(ra1*np.pi/180.0)*np.cos(ra2*np.pi/180.0) - np.sin(ra1*np.pi/180.0)*np.sin(ra2*np.pi/180.0))
    sin2_dec = 0.5*(1.0 - np.cos(dec1*np.pi/180.0)*np.cos(dec2*np.pi/180.0) - np.sin(dec1*np.pi/180.0)*np.sin(dec2*np.pi/180.0))

    arg = sin2_dec + np.cos(dec1*np.pi/180.0)*np.cos(dec2*np.pi/180.0)*sin2_ra
    if hasattr(arg, '__iter__'):
        arg[arg < EPS] = EPS
    else:
        arg = max(EPS, arg)

    return 2.0*np.arcsin(np.sqrt(arg))*180.0/np.pi

def writeMaskFile(ra, dec, r, ID, fileOutName, filter, tract, patch, mag):

    f = open(fileOutName, "w+")

    f.write("# BRIGHT STAR CATALOG: Jean Coupon\n")
    f.write("# GENERATED ON: {0:s}\n".format(NOW.strftime('%Y-%m-%d')))
    f.write("# TRACT: {0:d}\n".format(tract))
    f.write("# PATCH: {0:d},{1:d}\n".format(patch[0], patch[1]))
    f.write("# FILTER {0:s}\n".format(filter))
    f.write("wcs; fk5\n")

    for i in range(len(ra)):
        f.write("circle({0:f},{1:f},{2:f}d) # ID: {3:d}, mag: {4:f}\n".format(ra[i], dec[i], r[i], ID[i], mag[i])) #ID[i].translate(None, '-')))
        if mag[i] > MAG_LIMIT:
            f.write("box({0:f},{1:f},{2:f}d,{3:f}d,0) # ID: {4:d}, mag: {5:f}\n".format(ra[i], dec[i], 0.1*3.0*r[i], 3.0*r[i], ID[i], mag[i])) #ID[i].translate(None, '-')))
            f.write("box({0:f},{1:f},{2:f}d,{3:f}d,0) # ID: {4:d}, mag: {5:f}\n".format(ra[i], dec[i], 3.0*r[i], 0.1*3.0*r[i], ID[i], mag[i])) #ID[i].translate(None, '-')))

    f.close()
    return

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





# ----------------------------------------------------- #
# Main
# ----------------------------------------------------- #



if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument('option', help="Which quantity to plot")
    parser.add_argument('-i', '--input', default=None, help='input file')
    parser.add_argument('-o', '--output', default=None, help='output file')

    parser.add_argument('-skyMap', default=None, help='input sky map')


    parser.add_argument('-k', '--keys', default=None, help='Keys')
    parser.add_argument('-rmax', type=float, default=5.0, help='Rmax for source density plot [arcsec]')
    parser.add_argument('-mag', type=float, default=None, help='Star magnitude')
    parser.add_argument('-s', '--select', type=str,   default=None, help='Selection')
    parser.add_argument('-worstSeeing', action='store_true', help='Compute radius around worst seeing stars')
    parser.add_argument('-bestSeeing', action='store_true', help='Compute radius around best seeing stars')
    parser.add_argument('-basename', type=str,   default="BrightStarMask", help='Base name for region files')

    args = parser.parse_args()

    main(args)

#!/usr/bin/env python

"""
Jean coupon - 2017
script to convert coordinates into pixels
"""

import os, sys
import numpy as np
import re

import collections


from astropy.io import ascii,fits
from astropy.coordinates import SkyCoord
import astropy.wcs       as wcs
#import astropy.stats     as astats
#import photutils         as pu


# ----------------------------------------------------- #
# global variables
# ----------------------------------------------------- #

EPS = 1.0e-8

# ----------------------------------------------------- #
# main
# ----------------------------------------------------- #

def main(args):
    function = getattr(sys.modules[__name__], args.option)(args)
    return

# ----------------------------------------------------- #
# Main functions
# ----------------------------------------------------- #


def toPixels(args):
    """
    Convert coordinates from pixels.
    Assumes HSC pixel size
    """

    verbose = True

    """ options
    """
    pixelScale = 0.17
    fileInName = args.input.split(",")

    """ read input file
    """
    cols, _ = getCols(fileInName[0], ['user_ra', 'user_dec', 'ra', 'dec'], dictionary=True)
    # coord = SkyCoord(cols['ra'], cols['dec'], unit="deg")
    N = len(cols['ra'])

    """ create wcs object
    """
    w = createZeroWcs(pixelScale, 0.0, 0.0)

    x = np.zeros(N)
    y = np.zeros(N)

    w.wcs.crval = [cols['user_ra'][0], cols['user_dec'][0]]
    count = 0

    """ loop over all objects
    but do the conversion in blocks sharing the
    same nearby star
    """

    # count = 0
    for i in range(N):
        w.wcs.crval = [cols['user_ra'][i], cols['user_dec'][i]]
        x[i], y[i] = w.wcs_world2pix(cols['ra'][i], cols['dec'][i], 1)

        # try:
        #    np.testing.assert_array_almost_equal(w.wcs.crval, [cols['user_ra'][i], cols['user_dec'][i]], decimal=4)
        # except:
        #    x[i-count:i+1], y[i-count:i+1] = w.wcs_world2pix(cols['ra'][i-count:i+1], cols['dec'][i-count:i+1], 1)
        #    w.wcs.crval = [cols['user_ra'][i], cols['user_dec'][i]]
        #    count = 0
        # count += 1

        #try:
        #    np.testing.assert_array_almost_equal(w.wcs.crval, [cols['user_ra'][i], cols['user_dec'][i]], decimal=4)
        #except:
        #    x[i-count:i+1], y[i-count:i+1] = wcs.utils.skycoord_to_pixel(coord[i-count:i+1], w, origin=1, mode='wcs')
        #    w.wcs.crval = [cols['user_ra'][i], cols['user_dec'][i]]
        #    count = 0

        #count += 1

        if verbose:
            if (i+1)%10000 == 0:
                sys.stderr.write("\r" + "Processed {0:d}/{1:d} objects ({2:.2f}%)".format(i+1, N, 100.0*float(i+1)/float(N)))
                sys.stderr.flush()

    # count -= 1
    # x[i-count:i+1], y[i-count:i+1] = w.wcs_world2pix(cols['ra'][i-count:N], cols['dec'][i-count:N], 1)

    #count -= 1
    #x[i-count:N], y[i-count:N] = wcs.utils.skycoord_to_pixel(coord[i-count:N], w, origin=1)

    if verbose: sys.stderr.write("\r" + "Processed {0:d}/{1:d} objects ({2:.2f}%)\n".format(i+1, N, 100.0))

    outCols = [fits.Column(name='x', format='E', array=x), fits.Column(name='y', format='E', array=y)]

    fileIn = fits.open(fileInName[0])
    tbhdu  = fits.HDUList([fileIn[0], fits.BinTableHDU.from_columns(fileIn[1].columns + fits.ColDefs(outCols))])
    tbhdu.writeto(args.output, clobber=True)

    fileIn.close()

    return

def createZeroWcs(pixelScale, ra, dec):
    """
    Create a new WCS object, centering ra,dec  on 0,0 pixel
    """

    w = wcs.WCS(naxis=2)
    w.wcs.crpix = [0.0, 0.0]
    w.wcs.cdelt = np.array([-pixelScale, pixelScale])
    w.wcs.crval = [ra, dec]
    w.wcs.ctype = ["RA---TAN", "DEC--TAN"]

    return w


def stackStar(args, show=False):

    """
    Stack fits images
    """

    import astropy.wcs as wcs

    verbose = True

    IMAGE_EXT = 1
    WEIGHT_EXT = 2


    N = len(os.listdir(args.input))

    """ first record smaller coordinates and pixelScale
    """
    nx, ny = np.inf, np.inf
    for i,f in enumerate(os.listdir(args.input)):

        # print "Reading", args.input+"/"+f
        fileIn = fits.open(args.input+"/"+f)

        nx = min(nx, fileIn[IMAGE_EXT].header['NAXIS2'])
        ny = min(ny, fileIn[IMAGE_EXT].header['NAXIS1'])

        if i == 0:
            w = wcs.WCS(fileIn[IMAGE_EXT].header)
            pixelScale = wcs.utils.proj_plane_pixel_scales(w)[0]

        fileIn.close()


    if verbose: sys.stderr.write("Final image size: [{0:d},{1:d}]\n".format(nx, ny))


    """ Main arrays
    """
    mean = np.zeros((nx, ny))
    stddev = np.zeros((nx, ny))

    """ compute mean and standard deviation
    """
    if verbose:
        sys.stderr.write("\r" + "Mean: {0:d}/{1:d} images".format(0, N))
        sys.stderr.flush()

    for i,f in enumerate(os.listdir(args.input)):
        fileIn = fits.open(args.input+"/"+f, memmap=True)
        data = fileIn[IMAGE_EXT].data
        data[data != data] = 0.0
        mean += crop(data, nx, ny)
        del data
        fileIn.close()
        if verbose:
            if (i+1)%1 == 0:
                sys.stderr.write("\r" + "Mean {0:d}/{1:d} images".format(i+1, N))
                sys.stderr.flush()

    if verbose: sys.stderr.write("\r" + "Mean {0:d}/{1:d} images\n".format(i+1, N))

    mean /= float(N)

    if verbose:
        sys.stderr.write("\r" + "Standard dev {0:d}/{1:d} images".format(0, N))
        sys.stderr.flush()

    for i,f in enumerate(os.listdir(args.input)):
        fileIn = fits.open(args.input+"/"+f, memmap=True)
        data = fileIn[IMAGE_EXT].data
        data[data != data] = 0.0
        stddev += pow(crop(data, nx, ny)-mean, 2.0)
        del data
        fileIn.close()

        if verbose:
            if (i+1)%1 == 0:
                sys.stderr.write("\r" + "Standard dev {0:d}/{1:d} images".format(i+1, N))
                sys.stderr.flush()

    if verbose: sys.stderr.write("\r" + "Standard dev {0:d}/{1:d} images\n".format(i+1, N))

    stddev = np.sqrt(1.0/float(N) * stddev)

    """ create output image
    """
    wStack, thduStack = createFitsImage([0.0, 0.0], pixelScale, nx, ny)

    """ last loop over files
    """
    weight = np.ones((nx, ny))

    if verbose:
        sys.stderr.write("\r" + "Stacking {0:d}/{1:d} images".format(0, N))
        sys.stderr.flush()

    for i,f in enumerate(os.listdir(args.input)):
        fileIn = fits.open(args.input+"/"+f, memmap=True)

        data = crop(fileIn[IMAGE_EXT].data, nx, ny)
        weight[:] = 1.0

        if N > 3:
            clip = (data != data) | (data > 3.0*stddev)
            data[clip] = 0.0
            weight[clip] = 0.0

        thduStack[IMAGE_EXT].data += data
        thduStack[WEIGHT_EXT].data += weight

        del fileIn[IMAGE_EXT].data
        fileIn.close()

        if verbose:
            if (i+1)%1 == 0:
                sys.stderr.write("\r" + "Stacking {0:d}/{1:d} images".format(i+1, N))
                sys.stderr.flush()

    if verbose: sys.stderr.write("\r" + "Stacking {0:d}/{1:d} images\n".format(i+1, N))

    thduStack[IMAGE_EXT].data /= thduStack[WEIGHT_EXT].data

    """ write output files
    """
    if verbose: sys.stderr.write("Writing stacked image...")
    thduStack.writeto(args.output, clobber=True)
    if verbose: sys.stderr.write("done\n")


    return

def stackStarOLD(args, show=False):

    """
    Stack fits image
    """
    import astropy.wcs as wcs
    import pickle

    doNoClip = False

    IMAGE_EXT = 1
    WEIGHT_EXT = 2

    N = len(os.listdir(args.input))

    """ first record smaller coordinates and pixelScale
    """
    nx, ny = np.inf, np.inf
    for i,f in enumerate(os.listdir(args.input)):

        print "Reading", args.input+"/"+f

        fileIn = fits.open(args.input+"/"+f)
        image = fileIn[IMAGE_EXT].data

        nx = min(nx, image.shape[0])
        ny = min(ny, image.shape[1])

        if i == 0:

            w = wcs.WCS(fileIn[IMAGE_EXT].header)
            pixelScale = wcs.utils.proj_plane_pixel_scales(w)[0]

        fileIn.close()

    """ main arrays
    """
    print N, nx, ny
    data = np.zeros((N, nx, ny))
    weight = np.zeros((N, nx, ny))

    """ stacking loop
    """
    if True:
        for i,f in enumerate(os.listdir(args.input)):

            fileIn = fits.open(args.input+"/"+f)
            image = fileIn[IMAGE_EXT].data

            image = crop(image, nx, ny)

            data[i] = image
            weight[i] = np.ones((nx,ny))

            fileIn.close()

        # pickle.dump([data, weight], open( "save.pickle", "wb" ) )

    # data, weight = pickle.load( open( "save.pickle", "rb" ) )


    """ create output images
    (sigmal clipped and non-clipped versions
    """
    wStack, thduStack = createFitsImage([0.0, 0.0], pixelScale, nx, ny)
    if doNoClip:
        wStackNoClip, thduStackNoClip = createFitsImage([0.0, 0.0], pixelScale, nx, ny)

    mean = np.mean(data, axis=0)
    std = np.std(data, axis=0)

    """ loop over files
    """
    for i in range(N):

        if doNoClip:
            thduStackNoClip[IMAGE_EXT].data += data[i]
            thduStackNoClip[WEIGHT_EXT].data += weight[i]

        if N > 3:
            clip = data[i] > 5.0*std
            data[i][clip] = 0.0
            weight[i][clip] = 0.0

        thduStack[IMAGE_EXT].data += data[i]
        thduStack[WEIGHT_EXT].data += weight[i]

    thduStack[IMAGE_EXT].data /= thduStack[WEIGHT_EXT].data
    if doNoClip:
        thduStackNoClip[IMAGE_EXT].data /= thduStackNoClip[WEIGHT_EXT].data

    """ write output files
    """
    thduStack.writeto(args.output, clobber=True)
    if doNoClip:
        thduStackNoClip.writeto(".".join(args.output.split(".")[:-1])+"_noClip."+args.output.split(".")[-1], clobber=True)

    return


def crop(image, Nx_cropped, Ny_cropped):
    """
    crop image according to scales
    (between 0 and 1)
    """

    (Nx, Ny) = image.shape

    edge_x = 0
    edge_y = 0
    if (Nx-Nx_cropped) % 2 !=0: edge_x = 1
    if (Ny-Ny_cropped) % 2 !=0: edge_y = 1

    return image[0+(Nx-Nx_cropped)/2:Nx-(Nx-Nx_cropped)/2-edge_x, 0+(Ny-Ny_cropped)/2:Ny-(Ny-Ny_cropped)/2-edge_y]


def createFitsImage(center, pixelScale, nx, ny):
    """ creates fits image object

    Input:
    - center: center of pointing,
    celestial coordinates  in decimal degree
    - pixelScale: size of pixel in deg^-1
    - nx, ny: size of the array

    Output:
    - wcs: wcs object
    - fits image object (hdu list)
    """

    # Create a new WCS object.
    w = wcs.WCS(naxis=2)

    # Set up tangential projection
    # reference pixel at the image center
    # (fits convention -> first pixel coord = 1,1)
    w.wcs.crpix = [(nx+1.0)/2.0, (ny+1.)/2.0]
    w.wcs.cdelt = [- pixelScale, pixelScale]
    w.wcs.crval = center
    w.wcs.ctype = ["RA---TAN", "DEC--TAN"]

    # header is an astropy.io.fits.Header object.
    # create PrimaryHDU object
    header = w.to_header()
    HDUs = [fits.PrimaryHDU(header=header)]

    # create data array
    HDUs.append(fits.ImageHDU(data=np.zeros((nx, ny)), name="DATA", header=header))

    # create exposure map (= weight)
    HDUs.append(fits.ImageHDU(data=np.zeros((nx, ny)), name="WEIGHT", header=header))

    # create background map
    # HDUs.append(fits.ImageHDU(data=np.zeros((nx, ny)), name="BGK", header=header))

    return w, fits.HDUList(HDUs)



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

    args = parser.parse_args()

    main(args)

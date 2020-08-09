# des_transient_masks

This repo contains Python code and data that allow one to create masks marking regions of the Dark Energy Camera
that contain statistically significant excesses of transient objects, which means that they probably contain
defects that are contaminating the DES images. True astrophysical transients should be randomly located on the array,
but spurious signals due to CCD defects will cluster in pixel coordinats.

The transient catalog and techniques to use them in this way were 
developed by Pedro Bernardinelli.  Gary Bernstein assisted with putting the code into its present form. Questions
may be directed to either.

## Contents
* `ops_epochs.fits` is a download of the `desoper:OPS_EPOCH` table on DESDM.  It gives names (like `Y2E1`) to
intervals of time during the DES survey.  These are the standard epochs used in choosing various calibration filesets
in DESDM processing of DES images.  The transient masks are similarly divided, since CCD defects can come and go over time.
* `y6a1_transient_maps.fits.gz` (19 MB) contains a series of FITS-format 3d images, one for each observing epoch.  Each epoch's 
image has shape (when read into numpy) of (62,256,256).  A pixel at location (_CCDNUM-1, Y, X_) in this image contains the
number of transients in Pedro's Y6A1 transient catalog which fall on the given _CCDNUM_, at pixel locations of

8*_X_ +1 < x < 8*_X_ + 9, 16*_Y_ +1 < y < 16*_Y_ + 17.

In other words it is a 3d histogram of the transients on DECam during that epoch, with bins of width 8x16 pixels in the
_(x,y)_ span of each 2048x4096 CCD image.

* `transient_masks.py` contains classes to manipulate the transient count maps file above to locate regions of 
excess transients on the array and create masks.  The tasks that this code can execute are described below.

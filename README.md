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

## Operations
The `transient_masks.py` file contains code to perform various masking tasks.  Some of the tasks are accessible by
running the file as a script from the command line.  For others, you need to import the file as a module and use the
classes / functions within Python.

* The `gather_data` function is the one that was used to generate the `y6a1_transient_maps.fits.gz` maps.  It uses
the `ops_epoch.fits` table to define the epoch, and the (very large, not included here) list of transients.  Users
should not have need for this.
* The `plot_epoch` function will generate a multipage PDF file plotting the 2d histogram of transients for each CCD in
a given epoch, and marking the bins of the histogram that are above threshold and considered masked.
* When run from the command line, if one provides the `--epoch` argument, then the output file of the program will be
a multiextension FITS file containing one 2048 x 4096 image for each CCD.  The datatype of the images are uint8, and
each pixel contains either 0 (valid) or 1 (excess transients). \[This task is done by the `EpochMask.write_masks` method.\]
* When run from the command line and given an "event table" at the `--input` argument, the output will be a new version of the
table that excludes any event that occurs in a CCD region that is masked during its epoch.  The input and output tables
can be in any format that the `astropy.table.Table.read/write` methods recognize.  They must contain columns corresponding
to the `EXPNUM, CCDNUM, X, Y` values for the time and pixel location of the event.

## Usage, options, and thresholds
The possible operations make use of several possible arguments:
* `x_shrink,y_shrink` or `x_bin,y_bin` can be set to integers 1,2,4,... if one wants to bin the transient counts 
this many times more coarsely
than the 8x16 (in x,y) that they are already stored in the `transient_maps` file.  A coarser binning could be helpful in reducing
shot noise in the count maps, making the masking more sensitive to defects, which might be good.  Another reason to potenially
increase `y_bin` is that the bad columns may be only partly masked, i.e. gaps in the masks, and larger `y_bin` might help fill
these gaps.  But the downside of increased binning is that if the defect is localized (single pixel or column), then larger bins
are just diluting the transient excess and it could get lost; and also the binned regions are larger, so good pixels could be
getting thrown away.  The default (keeping the 8x16 inputs) is a good compromise.
* `fp_rate` determines which count bins are masked.  For a given CCD, the modal transient count is used to define an expected
Poisson distribution of bin counts.  The threshold for masking is set such that a true Poisson distribution will trigger
masking for `fp_rate` of the bins.  The threshold is also forced to be at least 2x the median bin count.  
The default `fp_rate=1e-4` keeps spurious masking to O(1) bin per CCD and agrees pretty well with "eyeball" judgements about
which bins are real excesses.
* `default` or `keep_no_epoch` determine what will be done with events that occur on exposures that are not within one of
the tabulated epochs.
* `expnum_col,ccdnum_col,x_col,y_col` tell the code the column names of an input table that contain the required quantities
for filtering.

There are docstrings in the code, and `./transient_masks.py -h` will print some help.

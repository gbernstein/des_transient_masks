''' Routines and classes for use in building and using image masks built by
looking for excess transients in physical regions of CCD.  Assume that
transient counts will be stored in binned format at some nominal resolution.
Provide a way to build this stored data.  

Calling this script interactively will cause it to filter a specified
list of events.  Or can be used to output full-res mask files for a
chosen epoch.
'''

import os
import sys
import pickle
import astropy.table as tb 
import numpy as np 
import astropy.io.fits as pf
import matplotlib.pyplot as pl 
from matplotlib.backends.backend_pdf import PdfPages
import scipy.stats as st

import argparse

# DECam info
x_size = 2048
y_size = 4096
n_ccds = 62
half_ccdnum = 31  # The CCDNUM that is only half working

def drop_threshold(counts, fp_rate, ccdnum=0):
    ''' Given an array of counts per cell, return a boolean array
    marking the array elements that satisfy both
       * The counts are more than twice the mean, and
       * The counts would occur <fp_rate of the time in a Poisson
         distribution defined by the modal count.
    If the ccdnum is the half-chip, modes and medians
    are calculated accordingly.  If there are no counts,
    all cells are flagged.
    
    Returns: the boolean array; 
             the count threshold;
             the fraction of cells above threshold;
             the fraction of counts in high cells'''

    out = np.ones_like(counts, dtype=bool)
    if np.all(counts==0):
        return out, 0., 1., 1.
    if ccdnum==half_ccdnum:
        tmp = counts[:,:counts.shape[1]//2]
        mode = st.mode(tmp.flatten())[0] 
        mean = np.mean(tmp.flatten())
    else:
        mode = st.mode(counts.flatten())[0] 
        mean = np.mean(counts.flatten())

    poisson = st.poisson(mode + 1)
    threshold = max(np.ceil(2 * mean), poisson.isf(fp_rate))

    if ccdnum==half_ccdnum:
        out[:,:counts.shape[1]//2] = tmp > threshold
        frac_cells = np.count_nonzero(out[:,:counts.shape[1]//2]) / out[:,:counts.shape[1]//2].size
        frac_counts = np.sum(counts[out]) / np.sum(counts)
    else:
        out = counts > threshold
        frac_cells = np.count_nonzero(out) / out.size
        frac_counts = np.sum(counts[out]) / np.sum(counts)
    return out, threshold, frac_cells, frac_counts

def bin_number(x,width=8.,xmin=0.,xmax=2047.9):
    ''' Place x values into equal-width bins of given width and spanning
    (xmin,xmax).  Values are clipped to xmin/xmax.'''
    ix = np.floor((np.clip(x,xmin,xmax) - xmin) / width)
    return ix.astype(int,casting='unsafe')
        

class EpochMask:
    '''Holds the mask arrays for all CCDs for one DECam observing epoch'''
    def __init__(self, name, min_expnum, max_expnum, counts, fp_rate=1e-4,
                     x_shrink=None, y_shrink=None):
        ''' Build a mask set for the specified epoch using:
         * the `count` array of transient counts per input cell
         * the `fp_rate` to use for setting threshold
         * optionally, bin the `counts` by factors `x_shrink,y_shrink`
           to make a coarser mask
        '''
        self.name = name
        self.min_expnum = min_expnum
        self.max_expnum = max_expnum
        # Determine the bin width in use for the input counts
        xw_in = x_size // counts.shape[2]  
        yw_in = y_size // counts.shape[1]

        if x_shrink is None:
            # use the counts as is.
            tmp = counts
            self.xw = xw_in
        else:
            # Bin up the count data in x
            if counts.shape[2] % x_shrink != 0:
                raise ValueError('x_shrink must divide count dimensions equally')
            tmp = np.sum( counts.reshape(counts.shape[0],
                                         counts.shape[1],
                                         counts.shape[2]//x_shrink, x_shrink),
                              axis=3)
            self.xw = xw_in * x_shrink
        if y_shrink is None:
            # use the counts as is.
            self.yw = yw_in
        else:
            # Bin up the count data in y
            if counts.shape[1] % y_shrink != 0:
                raise ValueError('y_shrink must divide count dimension equally')
            tmp = np.sum( tmp.reshape(tmp.shape[0],
                                      tmp.shape[1]//y_shrink, y_shrink,
                                      tmp.shape[2]),
                              axis=2)
            self.yw = yw_in * y_shrink

        # Now build mask by thresholding, doing each CCD individually:
        frac_cells = []
        masked_counts = 0.
        self.mask = np.ones_like(tmp,dtype=bool)
        for c in range(n_ccds):
            m,thresh,f_cells, f_counts = drop_threshold(tmp[c], fp_rate,ccdnum=c+1)
            self.mask[c] = m
            if f_cells<1.:
                frac_cells.append(f_cells)
                masked_counts = masked_counts + f_counts * np.sum(counts[c])

        # Save away the overall cell /count fractions
        self.frac_cells = np.mean(frac_cells)
        self.frac_counts = masked_counts / np.sum(counts)

    def apply(self,ccdnum,x,y):
        ''' Return a boolean array giving the mask value at points
        defined by the arrays `ccdnum,x,y` (which must have matching shapes).
        x and y assume 1-indexed positions as per FITS conventions.
        '''
        ix = bin_number(x, width=self.xw, xmin=1., xmax=x_size+0.99)
        iy = bin_number(y, width=self.yw, xmin=1., xmax=y_size+0.99)
        return self.mask[ccdnum-1,iy,ix]

    def write_masks(self, fitsname):
        ''' Write full-res mask arrays to a FITS file.'''
        hdus = [pf.PrimaryHDU()]
        for i in range(n_ccds):
            m = self.mask[i]
            hdr = pf.Header()
            hdr['CCDNUM'] = i+1
            out = np.ones( (m.shape[0],int(self.yw),
                            m.shape[1],int(self.xw)), dtype=np.uint8 ) * \
              m[:,np.newaxis,:,np.newaxis]
            out = out.reshape(y_size,x_size)
            hdus.append(pf.ImageHDU(out, header = hdr))
        pf.HDUList(hdus).writeto(fitsname,overwrite=True)
        return
        
    def in_epoch(self,expnum):
        ''' Return boolean array to match input expnum array denoting
        which exposures are in this epoch.'''
        return np.logical_and(expnum >= self.min_expnum, expnum<= self.max_expnum)
    
class AllMasks:
    ''' Class containing mask set for all epochs. '''

    def __init__(self, fitsfile,lazy=True,
                fp_rate=1e-4, x_shrink=None, y_shrink=None):
        '''Initialize the class from a FITS file.  
        If `lazy=True`, then keep the file open and read the data for an
        epoch only when it's needed.
        '''

        self.epochs = {}  # Dictionary of loaded EpochMasks
        self.loaded_all_epochs = False
        self.ff = pf.open(fitsfile)
        if lazy:
            # Save the epoch names, expnum ranges, and extension numbers for
            # each epoch
            self.names = []
            mins = []
            maxes = []
            extns = []
            for i,hdu in enumerate(self.ff):
                if 'EPOCH' in hdu.header:
                    # This HDU holds an epoch.  Record its name and range
                    self.names.append(hdu.header['EPOCH'])
                    mins.append(hdu.header['MIN_EXP'])
                    maxes.append(hdu.header['MAX_EXP'])
                    extns.append(i)
            self.mins = np.array(mins)
            self.maxes = np.array(maxes)
            self.loaded = np.zeros_like(mins,dtype=bool)
            self.extns = np.array(extns)
            self.fp_rate = fp_rate
            self.x_shrink = x_shrink
            self.y_shrink = y_shrink
        else:
            # Not lazy, load all extensions now
            for hdu in self.ff:
                if 'EPOCH' in hdu.header:
                    # This HDU holds an epoch.  Save its data
                    name = hdu.header['EPOCH']
                    minexp = hdu.header['MIN_EXP']
                    maxexp = hdu.header['MAX_EXP']
                    self.epochs[name] =  EpochMask(name, minexp, maxexp, hdu.data,
                                                fp_rate=fp_rate, x_shrink=x_shrink, y_shrink=y_shrink)
                    # Report stats
                    print('Epoch {:s} masks {:.2f}% of area, {:.2f}% of transients'.format(name,
                                self.epochs[name].frac_cells*100,
                                self.epochs[name].frac_counts*100))
            # Done reading
            self.loaded_all_epochs = True
            self.ff.close()
                                     
    def load(self,epoch):
        'Insure mask is loaded for specified epoch'
        if self.loaded_all_epochs or epoch in self.epochs:
            # Already have it.
            return

        # Look through names for this epoch:
        found = False
        for i,n in enumerate(self.names):
            if n==epoch:
                found = True
                iext = self.extns[i]
                self.epochs[n] =  EpochMask(n, self.mins[i], self.maxes[i],
                                            self.ff[iext].data,
                                            fp_rate=self.fp_rate,
                                            x_shrink=self.x_shrink,
                                            y_shrink=self.y_shrink)
                self.loaded[i] = True
                # Report stats
                print('Epoch {:s} masks {:.2f}% of area, {:.2f}% of transients'.format(n,
                                self.epochs[n].frac_cells*100,
                                self.epochs[n].frac_counts*100))
                break
        if not found:
            raise ValueError('Non-existent epoch ' + epoch)
        if np.all(self.loaded):
            # If we've loaded last epoch, close the file
            self.loaded_all_epochs = True
            self.ff.close()
        return
        
    def apply(self, expnum, ccdnum, x, y, default=True):
        '''Return a boolean array of shape matching all arguments' shape,
        which states whether the given exposure/location is masked.
        Value of `default` is assigned when no epoch contains the expnum.
        '''
        # Initialize results
        if default:
            out = np.ones_like(expnum, dtype=bool)
        else:
            out = np.zeros_like(expnum, dtype=bool)

        # Keep track of whether each input has been done
        done = np.zeros_like(out)

        # Loop through each loaded epoch
        for name,m in self.epochs.items():
            use = e.in_epoch(expnum)
            if np.any(use):
                out[use] = e.apply(ccdnum[use],x[use],y[use])
                done = np.logical_or(done, use)

        if not (self.loaded_all_epochs or np.all(done)):
            # There might be inputs that need additional loads.
            # Step through unloaded epochs
            for i,l in enumerate(self.loaded):
                # Move along if this epoch is already loaded
                if l:
                    continue
                # See if this epoch is needed
                use = np.logical_and(expnum>=self.mins[i], expnum<=self.maxes[i])
                if np.any(use):
                    # Yes! Load it and use it as lookup
                    self.load(self.names[i])
                    out[use] = self.epochs[name].apply(ccdnum[use],x[use],y[use])
                    done = np.logical_or(done, use)
                    if np.all(done):
                        # No need to try other epochs
                        break

        return out


# This function will read the transient file and sum up counts, and save it all
# to a single FITS file.

def gather_data(transient_file='y6a1_transients.fits', epoch_file='ops_epoch.fits',
                out_file = 'y6a1_transient_maps.fits',
                    x_width=8, y_width=16):
    if x_size%x_width !=0:
        raise ValueError('x_width does not divide CCD size')
    if y_size%y_width !=0:
        raise ValueError('y_width does not divide CCD size')

    transients = tb.Table.read(transient_file)
    epochs = tb.Table.read(epoch_file)

    # Build an HDU for each epoch that has counts in the transients file
    hdus = [pf.PrimaryHDU()]
    for i in range(len(epochs)):
        print('Gathering epoch',epochs['NAME'][i])
        minexp = epochs['MINEXPNUM'][i]
        maxexp = epochs['MAXEXPNUM'][i]
        data = transients[(transients['EXPNUM'] <= maxexp) & \
                              (transients['EXPNUM'] >= minexp)]
        if len(data)==0:
            # No transients in this epoch, so no count array
            continue

        # Collect info on the transients and bin them
        cyx = np.vstack( [np.array(data['CCDNUM'],dtype=float),
                          data['YWIN_IMAGE'],
                          data['XWIN_IMAGE']])
        counts = np.histogramdd(cyx.T, bins=(n_ccds,y_size//y_width,x_size//x_width),
                                    range= ( (0.5,62.5), (1.,y_size+1), (1.,x_size+1.)))[0]

        # Now make an image HDU out of it
        hdr = pf.Header()
        hdr['EPOCH'] = epochs['NAME'][i]
        hdr['EXTNAME'] = epochs['NAME'][i]
        hdr['MIN_EXP'] = epochs['MINEXPNUM'][i]
        hdr['MAX_EXP'] = epochs['MAXEXPNUM'][i]
        hdus.append(pf.ImageHDU(counts.astype(np.uint16,casting='unsafe'),
                                              header=hdr))

    # Save the file
    pf.HDUList(hdus).writeto(out_file,overwrite=True)
    return

def plot_epoch(epoch, pdf_file=None, countfile='y6a1_transient_map.fits',
               fp_rate=1e-4, x_shrink=None, y_shrink=None):
    '''Make a PDF file showing the (binned) transient count maps,
    overlaid with marks of over-threshold pixels.  Multipage
    PDF output, one CCD per page.
    '''

    pp = PdfPages(pdf_file)

    # Find the epoch in the counts file
    ff = pf.open(countfile)
    # Not lazy, load all extensions now
    found = False
    for hdu in ff:
        if 'EPOCH' in hdu.header:
            # This HDU holds an epoch.  Is it ours?
            name = hdu.header['EPOCH']
            if name==epoch:
                counts = hdu.data
                found = True
                break
    if not found:
        raise ValueError("Epoch " + epoch + " counts not found")
    
    # Rebin counts if desired
    if x_shrink is None:
        # use the counts as is.
        tmp = counts
    else:
        # Bin up the count data in x
        if counts.shape[2] % x_shrink != 0:
            raise ValueError('x_shrink must divide count dimensions equally')
        tmp = np.sum( counts.reshape(counts.shape[0],
                                     counts.shape[1],
                                    counts.shape[2]//x_shrink, x_shrink),
                          axis=3)
    if y_shrink is not None:
        # Bin up the count data in y
        if tmp.shape[1] % y_shrink != 0:
                raise ValueError('y_shrink must divide count dimension equally')
        tmp = np.sum( tmp.reshape(tmp.shape[0],
                                  tmp.shape[1]//y_shrink, y_shrink,
                                  tmp.shape[2]),
                              axis=2)
    for i in range(n_ccds):
        pl.figure(figsize=(7.5,10))
        mask,thresh,frac_cells,frac_counts = drop_threshold(tmp[i],fp_rate)
        pl.imshow(tmp[i],origin='lower',interpolation='nearest',vmin=0,vmax=3*thresh,
                      cmap = 'cividis')
        yx = np.where(mask)
        pl.plot(yx[1],yx[0], 'wx', alpha=0.3)
        pl.colorbar()
        pl.gca().set_aspect('equal')
        pl.title('CCDNUM {:02d} Transient fraction: {:.3f}\% Area {:.3f}\%'.format(i+1,
                                frac_counts*100, frac_cells*100))
        pp.savefig()
        pl.close()

    pp.close()

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description=
        'If `epoch` argument is given, output is a FITS file holding\n' + \
        'mask images derived from counts file for selected epoch.\n' + \
        '  Otherwise, apply mask built from transient map to `input` table.  Output\n' + \
        'file is a table that retains only unmasked objects',
        formatter_class = argparse.RawTextHelpFormatter)
    parser.add_argument('-i','--input', help='Input table', type=str)
    parser.add_argument('-e','--epoch', help='Name of epoch to make masks for', type=str)
    parser.add_argument('-o','--output', help='Output mask file or filtered table', type=str)
    parser.add_argument('--counts', help='FITS file with transient count maps', type=str,
                        default='y6a1_transient_maps.fits.gz')
    parser.add_argument('--fp_rate', help='Poisson false-positive rate on masks', type=float,
                        default=1e-4)
    parser.add_argument('--x_bin', help='Additional x binning factor on counts', type=int)
    parser.add_argument('--y_bin', help='Additional y binning factor on counts', type=int)
    parser.add_argument('--keep_no_epoch', help='Retain events in undefined epoch',
                            action='store_true')
    parser.add_argument('--expnum_col', help='Column holding EXPNUM', type=str,
                        default='EXPNUM')
    parser.add_argument('--ccdnum_col', help='Column holding CCDNUM', type=str,
                        default='CCDNUM')
    parser.add_argument('--x_col', help='Column holding x coord', type=str,
                        default='XWIN_IMAGE')
    parser.add_argument('--y_col', help='Column holding y coord', type=str,
                        default='YWIN_IMAGE')

    args = parser.parse_args()

    # Load masks
    all = AllMasks(args.counts,lazy=True,
                fp_rate=args.fp_rate, x_shrink=args.x_bin, y_shrink=args.y_bin)

    if args.epoch:
        # Write mask files for this epoch
        all.load(args.epoch)
        all.epochs[args.epoch].write_masks(args.output)
        sys.exit(0)

    # Otherwise read and filter an input table
    events = tb.Table.read(args.input)

    # Look up masks
    mask = all.apply(events[args.expnum_col],
                         events[args.ccdnum_col],
                         events[args.x_col],
                         events[args.y_col],
                         default = not args.keep_no_epoch)

    # Filter table and save - need to invert mask
    events[~mask].write(args.output,overwrite=True)

    sys.exit(0)

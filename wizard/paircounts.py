"""
Do paircounts

.. module:: paircounts

"""

from __future__ import print_function, division

import numpy as np
import healpy as hp
from os import path

from .dataset import interpret_bins

# TODO: incorporate the healpixel commands into the commandline to filter HPIX there
def paircounts(reference_catalogs, reference_bins,
               unknown_catalogs, unknown_bins,
               healpix_nside=32,
               paircounts_path='',
               load_paircounts=True,
               load_paircounts_path_only=True,
               min_sep=0.01, max_sep=180, nbins=100,
               healpix_neighbors=-1,
               keep_empty=False, overwrite=False,
               time0=0,
               **kwargs):
    """Calculate paircounts

    Parameters
    ----------
    reference_catalogs : a list of dataframes. If more than one entry is
                         present, computes correlations within reference bins
    reference_bins : dictionary of lists, where key is the column, and then the
                     list is [min, max, step] or entries. Calculates paircounts
                     in each combination of all bins here with unknown_bins
    unknown_catalogs : a list of dataframes. If more than one entry is
                       present, computes correlations within unknown bins
    unknown_bins : dictionary of lists, where key is the column, and then the
                   list is [min, max, step] or entries. Calculates paircounts
                   in each combination of all bins here with reference_bins
    healpix_nside : int [Default 32] Size of healpixel regions we correlate
    paircounts_path : string, default ''. If not empty, location to save
    load_paircounts : bool [Default True] If True, just load paircounts.
    load_paircounts_path_only : bool [Default True] If True, when loading
                                paircounts, append only the string for the
                                loading. When you run a really big job, this
                                can make things easier in the combine stage.
    min_sep : float [Default 0.01] minimum separation we calcualte in arcmin
    max_sep : float [Default 180] maximum separation we calcualte in arcmin
    nbins : int [Default 100] Number of bins in log separation for correlation
    healpix_neighbors : int [Default 1] If > -1, looks this many pixels away
                        for doing pairs. So if 1, when correlating a reference
                        pixel, we only look at unknown of said pixel and the
                        immediate neighbors. If 2, then the neighbors of the
                        neighbors are also correlated, &c. This should speed
                        up calculations immensely.
    keep_empty: bool [Default False] If rcat_healpixel != ucat_healpixel,
                and sum(npairs) == 0, then don't append.

    Returns
    -------
    a dictionary of paircounts keyed by healpixel and reference / unknown bins:

    e.g.
    results[DuDr] = [[hpix_Du, hpix_Dr, arr]] where arr has all the different
    slicings from reference_bins and unknown_bins and paircounts. These are
    always ordered equivalently, such that the DuDr if you had just calculated
    the FULL survey with these bins would just be the sum of
    results[DuDr][:][2]

    Notes
    -----
    All catalogs are assumed to have RA and DEC and a healpixel column 'HPIX'
    that is created during the dataset phase.
    """

    if time0: from time import time


    # turn string into array of values
    reference_bins_integerized = interpret_bins_integerized(reference_bins, 'reference__{0}_region')
    unknown_bins_integerized = interpret_bins_integerized(unknown_bins, 'unknown__{0}_region')

    # generate combinations of all the possible values
    reference_bin_combinations = generate_combinations(reference_bins_integerized)
    unknown_bin_combinations = generate_combinations(unknown_bins_integerized)

    if load_paircounts:
        # load results from disk if load_paircounts is True
        results = {}
        for rcat_i, rcat in zip(['Dr', 'Rr'], reference_catalogs):
            rcat_healpixels = np.unique(rcat['HPIX'])

            for ucat_j, ucat in zip(['Du', 'Ru'], unknown_catalogs):
                results_key = '{0}{1}'.format(ucat_j, rcat_i)
                results[results_key] = []

                for rcat_healpixel_i, rcat_healpixel in enumerate(rcat_healpixels):
                    if time0:
                        print(rcat_i, ucat_j, rcat_healpixel, rcat_healpixel_i / len(rcat_healpixels), len(rcat_healpixels), time() - time0)

                    # determine the unknown healpixels
                    if healpix_neighbors >= 0:
                        ucat_healpixels = find_healpix_neighbors(
                            rcat_healpixel, healpix_nside, healpix_neighbors)
                    elif healpix_neighbors == 0:
                        ucat_healpixels = np.array([rcat_healpixel])
                    else:
                        ucat_healpixels = np.unique(ucat['HPIX'])
                    for ucat_healpixel in ucat_healpixels:
                        # if time0:
                        #     print(rcat_i, ucat_j, rcat_healpixel, ucat_healpixel, time() - time0)
                        pair_path = '{0}{1}{2}_refpix_{3}_unkpix_{4}.npy'.format(
                                paircounts_path, ucat_j, rcat_i, rcat_healpixel, ucat_healpixel)

                        if keep_empty:
                            # always try to load if we are keeping empty
                            do_load = True
                        else:
                            if path.exists(pair_path):
                                do_load = True
                            else:
                                do_load = False
                        if do_load:
                            if load_paircounts_path_only:
                                # just append the path
                                results[results_key].append([rcat_healpixel, ucat_healpixel, pair_path])
                            else:
                                pairs = np.load(pair_path)
                                # no point in adding empty entries. Need to check for
                                # r_hpix = u_hpix because n_reference and n_unknown get
                                # updated if you are in the same healpixel
                                if ((rcat_healpixel != ucat_healpixel) &
                                    (pairs['npairs'].sum() == 0)):
                                    continue
                                # append the array
                                results[results_key].append([rcat_healpixel, ucat_healpixel, pairs])

        # TODO: also load up the auto correlations
    else:
        results = {}
        # cross
        for rcat_i, rcat in zip(['Dr', 'Rr'], reference_catalogs):
            rcat_healpixels = np.unique(rcat['HPIX'])
            for ucat_j, ucat in zip(['Du', 'Ru'], unknown_catalogs):
                results_key = '{0}{1}'.format(ucat_j, rcat_i)
                results[results_key] = []
                for rcat_healpixel_i, rcat_healpixel in enumerate(rcat_healpixels):
                    rcat_filter_pixel = rcat[rcat['HPIX'] == rcat_healpixel]
                    # determine the unknown healpixels
                    if healpix_neighbors >= 0:
                        ucat_healpixels = find_healpix_neighbors(
                            rcat_healpixel, healpix_nside, healpix_neighbors)
                    elif healpix_neighbors == 0:
                        ucat_healpixels = np.array([rcat_healpixel])
                    else:
                        ucat_healpixels = np.unique(ucat['HPIX'])
                    for ucat_healpixel_i, ucat_healpixel in enumerate(ucat_healpixels):
                        ucat_filter_pixel = ucat[ucat['HPIX'] == ucat_healpixel]
                        # create array based on length of reference_bin_combinations and unknown_bin_combinations
                        dtype = [('npairs', '{0}i4'.format(nbins)),
                                 ('n_reference', 'i8'),
                                 ('n_unknown', 'i8')]
                        # add dtypes based on combos
                        for key in sorted(reference_bins_integerized.keys()):
                            dtype += [(key, 'i4')]
                        for key in sorted(unknown_bins_integerized.keys()):
                            dtype += [(key, 'i4')]

                        # make the actual array
                        ref_len = len(reference_bin_combinations)
                        unk_len = len(unknown_bin_combinations)
                        pairs = np.zeros((ref_len * unk_len), dtype=dtype)

                        if time0:
                            print(rcat_i, ucat_j, rcat_healpixel, ucat_healpixel, rcat_healpixel_i / len(rcat_healpixels), ucat_healpixel_i / len(ucat_healpixels), len(rcat_filter_pixel), len(ucat_filter_pixel), time() - time0)

                        save_path = '{0}{1}{2}_refpix_{3}_unkpix_{4}.npy'.format(
                                    paircounts_path, ucat_j, rcat_i, rcat_healpixel, ucat_healpixel)
                        if path.exists(save_path) and not overwrite:
                            print('already computed {0}'.format(save_path))
                            results[results_key].append([rcat_healpixel, ucat_healpixel, save_path])
                            continue

                        # calculate paircounts in healpixel
                        for entry_i, unknown_combos in enumerate(unknown_bin_combinations):
                            # filter unknown
                            ucat_filter_pixel_item = ucat_filter_pixel[ucat_filter_pixel[unknown_combos.keys()].isin(unknown_combos).all(axis=1)]
                            for entry_j, reference_combos in enumerate(reference_bin_combinations):
                                entry_ij = entry_i * ref_len + entry_j
                                # filter reference
                                rcat_filter_pixel_item = rcat_filter_pixel[rcat_filter_pixel[reference_combos.keys()].isin(reference_combos).all(axis=1)]
                                # no point in doing the calculation if no entries
                                if ((len(rcat_filter_pixel_item) > 0) &
                                    (len(ucat_filter_pixel_item) > 0)):
                                    ucat_ra = ucat_filter_pixel_item['RA']
                                    ucat_dec = ucat_filter_pixel_item['DEC']
                                    if 'W' in ucat_filter_pixel_item:
                                        ucat_w = ucat_filter_pixel_item['W']
                                    else:
                                        ucat_w = None
                                    rcat_ra = rcat_filter_pixel_item['RA']
                                    rcat_dec = rcat_filter_pixel_item['DEC']
                                    if 'W' in rcat_filter_pixel_item:
                                        rcat_w = rcat_filter_pixel_item['W']
                                    else:
                                        rcat_w = None
                                    # calculate treecorr
                                    npairs, n1, n2, logr = treecorr_correlation(
                                        rcat_ra, rcat_dec,
                                        ucat_ra, ucat_dec,
                                        rcat_w, ucat_w,
                                        max_sep, min_sep, nbins)
                                    # update pairs
                                    pairs[entry_ij]['npairs'] = npairs

                                # update pairs array
                                for key in sorted(unknown_combos):
                                    # [0] because each entry is a list
                                    pairs[entry_ij][key] = unknown_combos[key][0]
                                for key in sorted(reference_combos):
                                    # [0] because each entry is a list
                                    pairs[entry_ij][key] = reference_combos[key][0]
                                pairs[entry_ij]['n_reference'] = len(rcat_filter_pixel_item)
                                pairs[entry_ij]['n_unknown'] = len(ucat_filter_pixel_item)

                        if len(paircounts_path) > 0:
                            if keep_empty:
                                do_save = True
                            else:
                                if rcat_healpixel == ucat_healpixel:
                                    # we use other parameters later like
                                    # n_reference and n_unknown, so save
                                    do_save = True
                                else:
                                    if pairs['npairs'].sum() == 0:
                                        # seriously, just wasting space!
                                        do_save = False
                                    else:
                                        do_save = True

                            if do_save:
                                np.save(save_path, pairs)

                        # append to array to return to world
                        results[results_key].append([rcat_healpixel, ucat_healpixel, pairs])

        # auto references
        # TODO: CODE THIS UP
        for rcat_i, rcat in enumerate(reference_catalogs):
            # do auto


            for rcat_j, rcat2 in enumerate(reference_catalogs[rcat_i + 1:], rcat_i + 1):
                pass

        # auto unknowns
        # TODO: CODE THIS UP
        for ucat_i, ucat in enumerate(unknown_catalogs):
            # do auto


            for ucat_j, ucat2 in enumerate(unknown_catalogs[ucat_i + 1:], ucat_i + 1):
                pass

    return results

# treecorr cross
def treecorr_correlation(cat1_ra, cat1_dec,
                         cat2_ra, cat2_dec,
                         cat1_w=None, cat2_w=None,
                         max_sep=180, min_sep=0.01, nbins=100):
    """

    Parameters
    ----------
    catN_ra/dec : array of ra and dec values
    catN_w : array of weights (default is none)
    min_sep : in arcmin, minimum separation we consider [default: 0.01]
    max_sep : in arcmax, maximum separation we consider [default: 180]
    nbins   : number of bins used [default: 100]

    Returns
    -------
    npairs, n1, n2, logr

    """
    import treecorr
    # setup treecorr
    config = {'nbins': nbins,
              'min_sep': min_sep,
              'max_sep': max_sep,
              'sep_units': 'arcmin'}
    cat1 = treecorr.Catalog(ra=cat1_ra, dec=cat1_dec,
                            w=cat1_w,
                            ra_units='deg', dec_units='deg')
    cat2 = treecorr.Catalog(ra=cat2_ra, dec=cat2_dec,
                            w=cat2_w,
                            ra_units='deg', dec_units='deg')
    n2 = treecorr.NNCorrelation(**config)

    # run treecorr to get paircounts
    n2.process(cat1, cat2)

    return n2.npairs, len(cat1_ra), len(cat2_ra), n2.logr

# treecorr auto
def treecorr_autocorrelation(cat_ra, cat_dec,
                             cat_w=None,
                             max_sep=180, min_sep=0.01, nbins=100):
    """

    Parameters
    ----------
    cat_ra/dec : array of ra and dec values
    cat_w : array of weights (default is none)
    min_sep : in arcmin, minimum separation we consider [default: 0.01]
    max_sep : in arcmax, maximum separation we consider [default: 180]
    nbins   : number of bins used [default: 100]

    Returns
    -------
    npairs, n1, logr

    """
    import treecorr
    # setup treecorr
    config = {'nbins': nbins,
              'min_sep': min_sep,
              'max_sep': max_sep,
              'sep_units': 'arcmin'}
    cat = treecorr.Catalog(ra=cat_ra, dec=cat_dec,
                            w=cat_w,
                            ra_units='deg', dec_units='deg')
    n2 = treecorr.NNCorrelation(**config)

    # run treecorr to get paircounts
    n2.process(cat)

    return n2.npairs, len(cat_ra), n2.logr

# create combinations from dict of lists
def generate_combinations(values):
    from itertools import product
    keys = sorted(values.keys())
    lists = [values[key] for key in keys]
    biglist = product(*lists)
    items = []
    for item in biglist:
        items.append({keys[i]: [item[i]] for i in xrange(len(item))})
    return items

def find_healpix_neighbors(seed_pixel, nside, neighbors):
    out_pixels = [seed_pixel]
    # TODO: This is silly, this conversion between arrays and lists
    for ith in range(neighbors):
        out_pixels += list(hp.pixelfunc.get_all_neighbours(nside, out_pixels).flatten())
        out_pixels = np.unique(out_pixels)
        out_pixels = out_pixels[out_pixels != -1]
        out_pixels = list(out_pixels)
    out_pixels = np.array(out_pixels)
    return out_pixels

def interpret_bins_integerized(bins, string='{0}'):
    bins_dictionary = interpret_bins(bins)

    bins_integerized = {}
    for key in bins_dictionary:
        if bins_dictionary[key][0] == 'between':
            # nbins -> nbins - 1 values, +2 for above and below
            bins_integerized[string.format(key)] = range(len(bins_dictionary[key][1]) + 1)
        elif bins_dictionary[key][0] == 'equal':
            # still +1 because 0th entry is none of the above
            bins_integerized[string.format(key)] = range(len(bins_dictionary[key][1]) + 1)

    return bins_integerized

"""
Create datasets given fits files

.. module:: dataset
"""

from __future__ import print_function
import numpy as np
import pandas as pd
import healpy as hp
import fitsio

def dataset(reference_data_path, reference_random_path,
            unknown_data_path, unknown_random_path,
            load_dataset=False, dataset_path='',
            reference_columns=['RA', 'DEC', 'Z'],
            reference_ra_column='RA',
            reference_dec_column='DEC',
            reference_w_column='',
            unknown_columns=['RA', 'DEC', 'MEAN_Z', 'Z_MC'],
            unknown_ra_column='RA',
            unknown_dec_column='DEC',
            unknown_w_column='',
            max_objects=30000000,
            label_every=10000,
            healpix_nside=32,
            healpix_filter_size=20,
            healpix_filter_start=-1,
            unknown_bins={'MEAN_Z': ['between', [0, 1.5, 0.07]]},
            reference_bins={'Z': ['between', [0, 1.2, 0.1]]},
            time0=0):
    """Load up datasets into understandable form

    Parameters
    ----------
    [reference,unknown]_[data,random]_path : string, locations to fits files of
                                             catalogs. They are assumed to have
                                             the entries in reference_columns
                                             or unknown_columns in them.
                                             IF you pass a list of 2 items,
                                             then interpret as hdf5 file, with
                                             entry 1 being file path, entry 2
                                             the key of the dataset.
    [reference,unknown]_[ra,dec,w]_column : string, location of said columns.
                                            w means 'weight' and does NOT have
                                            to be present. These will be
                                            renamed in my dataset

    max_objects : int [default: 30,000,000] Maximum number of objects from
                  catalog to load. Useful if memory issues ariase.
    label_every : int [default: 10,000] Maximum number of labels calculated and
                  assigned when augmenting catalog with `region' values which
                  will be later used in the paircounts. Again, useful if memory
                  issues arise.
    healpix_nside : int [default: 32] nside we divide the paircounts into for
                    later paircounting. This therefore also controls the number
                    of files saved, and thus the harddisk space of the output


    Returns
    -------
    reference_catalogs, unknown_catalogs : list of pandas dataframes of
                                           reference and unknown catalogs

    Notes
    -----
    Cannot handle loading up any rows that are vectors!
    """
    if time0: from time import time

    if load_dataset:
        # load up our already augmented hdf5 files
        hdf = pd.HDFStore(dataset_path)
        reference_catalogs = []
        unknown_catalogs = []
        for ith, catalog_append, first_key in zip([0, 1], [reference_catalogs, unknown_catalogs], ['reference', 'unknown']):
            for jth, second_key in enumerate(['data', 'random']):
                if time0:
                    print(first_key, second_key, time() - time0)
                key = '/{0}/{1}'.format(first_key, second_key)
                if len([[reference_data_path, reference_random_path],
                        [unknown_data_path, unknown_random_path]][ith][jth]) == 0:
                    continue
                if key in hdf:
                    catalog = hdf[key]
                    # possibly filter on healpixel region for references
                    # TODO: Does this mess up the autocorrelation signal?
                    if ((healpix_filter_start >= 0) &
                        (first_key == 'reference')):
                        hpix = catalog['HPIX']
                        conds = hpix_filter(hpix, healpix_filter_start, healpix_filter_size)
                        catalog = catalog[conds]
                    catalog_append.append(catalog)
        hdf.close()
    else:
        # load up fits files, augment
        # TODO: add ability to NOT do randoms for one or the other...
        reference_catalogs = []
        unknown_catalogs = []

        # references
        # convert reference_bins to dictionary of arrays
        reference_bins_dictionary = interpret_bins(reference_bins)
        unknown_bins_dictionary = interpret_bins(unknown_bins)

        # TODO: a bit messy and unclear. Basically I didn't want to rewrite this for loop over and over...
        catalogs_paths = [[reference_data_path, reference_random_path], [unknown_data_path, unknown_random_path]]
        catalog_appends = [reference_catalogs, unknown_catalogs]
        columns = [reference_columns, unknown_columns]
        ra_columns = [reference_ra_column, unknown_ra_column]
        dec_columns = [reference_dec_column, unknown_dec_column]
        w_columns = [reference_w_column, unknown_w_column]

        augment_labels = ['reference', 'unknown']
        bins_dictionaries = [reference_bins_dictionary, unknown_bins_dictionary]

        # reference and unknown
        for catalog_append, catalogs_path, column, ra_column, dec_column, w_column, augment_label, bins_dictionary in zip(catalog_appends, catalogs_paths, columns, ra_columns, dec_columns, w_columns, augment_labels, bins_dictionaries):
            # data and randoms
            for catalog_path in catalogs_path:
                if time0:
                    print(augment_label, catalog_path, time() - time0)

                # load up catalog
                if len(catalog_path) == 2:
                    # h5 file
                    catalog = load_h5_file(catalog_path[0], catalog_path[1], column, max_objects)
                else:
                    catalog = load_fits_file(catalog_path, column, max_objects)

                # rename ra, dec, w columns
                if ra_column != 'RA':
                    catalog.rename(columns={ra_column: 'RA'}, inplace=True)
                if dec_column != 'DEC':
                    catalog.rename(columns={dec_column: 'DEC'}, inplace=True)
                if len(w_column) > 0:
                    if w_column != 'W':
                        catalog.rename(columns={w_column: 'W'}, inplace=True)

                # add healpixel region
                catalog['HPIX'] = radec_to_index(catalog['DEC'], catalog['RA'], healpix_nside)

                # possibly filter on healpixel region for reference samples
                # TODO: Does this mess up the autocorrelation signal?
                if ((healpix_filter_start >= 0) &
                    (augment_label == 'reference')):
                    hpix = catalog['HPIX']
                    conds = hpix_filter(hpix, healpix_filter_start, healpix_filter_size)
                    catalog = catalog[conds]

                # augment reference catalogs based on bins
                for key in bins_dictionary:
                    if bins_dictionary[key][0] == 'between':
                        labels = digitize(catalog[key], bins_dictionary[key][1])
                    elif bins_dictionary[key][0] == 'equal':
                        # assign unique label to each value. This is kludgey
                        labels = np.zeros(len(catalog[key]), dtype=int)
                        values = bins_dictionary[key][1]
                        for val_i, val in enumerate(values):
                            conds = np.where(catalog[key] == val)
                            labels[conds] = val_i + 1
                            # 0 == not in list
                        # for ith in range(len(catalog[key])):
                        #     if catalog[key][ith] in values:
                        #         labels[ith] = values.index(catalog[key][ith]) + 1
                    catalog['{0}__{1}_region'.format(augment_label, key)] = labels

                # add to list
                catalog_append.append(catalog)


        # save
        if len(dataset_path) > 0:
            reference_catalogs[0].to_hdf(dataset_path, '/reference/data')
            reference_catalogs[1].to_hdf(dataset_path, '/reference/random')
            unknown_catalogs[0].to_hdf(dataset_path, '/unknown/data')
            unknown_catalogs[1].to_hdf(dataset_path, '/unknown/random')

    return reference_catalogs, unknown_catalogs

def hpix_filter(hpix, healpix_filter_start, healpix_filter_size):
    hpix_unique = np.unique(hpix)  # also sorts
    hpix_min = healpix_filter_start * healpix_filter_size
    hpix_max = min([(healpix_filter_start + 1) * healpix_filter_size, len(hpix_unique) - 1])
    conds = (hpix >= hpix_unique[hpix_min]) & (hpix < hpix_unique[hpix_max])
    return conds

def digitize(cat_z, zbins, label_every=0):
    """Digitize catalog assignment

    Parameters
    ----------
    cat_z: z values
    zbins: range of z values

    Returns
    -------
    labels: zbin labels

    Notes
    -----
    label == 0 means OUTSIDE zbins to the left
    label == len(zbins) + 1 means OUTSIDE zbins to the right!
    This means that zbins[label] gives you the RIGHT side of the bin
    so if you want the CENTER of the bin, you actually must find:
        0.5 * (zbins[label - 1] + zbins[label])

    """
    if label_every > 0:
        labels = np.zeros(len(cat_z), dtype=int) - 1
        n_iter = int(np.ceil(len(cat_z) / label_every))
        for ith in xrange(n_iter):
            indx_lower = ith * label_every
            indx_upper = indx_lower + label_every
            cat_z_i = cat_z[indx_lower: indx_upper]
            labels_i = np.digitize(cat_z_i, zbins)
            labels[indx_lower: indx_upper] = labels_i
    else:
        labels = np.digitize(cat_z, zbins)
    return labels

def load_fits_file(filename, columns, max_objects=0, **kwargs):
    """Load up fits file with fitsio, returns pandas dataframe

    Parameters
    ----------
    filename : string that fitsio reads
    columns : columns we extract from file
    kwargs : whatever to pass to fitsio to read the filename

    Returns
    -------
    df: pandas dataframe with specified columns

    Notes
    -----

    """
    data = fitsio.read(filename, columns=columns, **kwargs)
    if max_objects > 0 and max_objects < len(data):
        # only take max_objects
        indx = np.random.choice(len(data), max_objects, replace=False)
        data = data[indx]
    df = pd.DataFrame({key: data[key].byteswap().newbyteorder()
                       for key in columns})
    return df

def load_h5_file(filename, tablename, columns, sample_num=0, **kwargs):
    """Load up h5 file, returns pandas dataframe

    Parameters
    ----------
    filename: string that fitsio reads
    columns : columns we extract from file
    kwargs: whatever to pass to pandas to read the filename

    Returns
    -------
    df: pandas dataframe with specified columns

    Notes
    -----

    """
    try:
        df = pd.read_hdf(filename, tablename, columns=columns)
    except TypeError:
        # problem with fixed pandas tables, so just load the full thing, drop
        df = pd.read_hdf(filename, tablename)
        drop_columns = [column for column in df.columns if column not in columns]
        if len(drop_columns) > 0:
            df.drop(drop_columns, axis=1, inplace=True)
    if sample_num > 0 and sample_num < len(df):
        # only take sample_num objects
        indx = np.random.choice(len(df), sample_num, replace=False)
        df = df.iloc[indx]
        # redo index
        df.index = np.arange(len(df))
    return df

# thanks someone on the internet
def index_to_radec(index, nside):
    theta, phi = hp.pixelfunc.pix2ang(nside, index)

    # dec, ra
    return -np.degrees(theta - np.pi / 2.), np.degrees(phi)

def radec_to_index(dec, ra, nside):
    return hp.pixelfunc.ang2pix(nside, np.radians(-dec + 90.), np.radians(ra))

def interpret_bins(bins):
    results = {}
    for key in bins:
        if bins[key][0] == 'between':
            # if 3 entries, then turn into arange
            results[key] = ['between', np.arange(*bins[key][1])]
        elif bins[key][0] == 'equal':
            # else leave in
            results[key] = ['equal', bins[key][1]]

    return results

def catalog_to_histogram(z, bins, weights=None, z_min=None, z_max=None):
    conds = (z > z_min) & (z < z_max)
    if sum(conds) == 0:
        bins_out = bins
        pdf = np.zeros(len(bins) - 1)
    else:
        # bin
        pdf, bins_out = np.histogram(z, bins, weights=weights)

    centers = 0.5 * (bins_out[1:] + bins_out[:-1])
    edges=bins_out
    conds_edges = (edges > z_min) & (edges < z_max)
    edges=edges[conds_edges]
    conds = (centers > z_min) & (centers < z_max)
    centers = centers[conds]
    pdf = pdf[conds]

    # normalize with trapezoidal
    norm = np.sum(0.5 * (centers[1:] - centers[:-1]) * (pdf[1:] + pdf[:-1]))
    if norm == 0:
        norm = 1

    return pdf / norm, centers ,edges

def _modify_indices_redshift(self_z, other_z, other_index_in=None):
    if not other_index_in:
        other_index_in = np.array([True] * len(other_z))
    # adjust indices based on redshifts.
    self_index = np.array([np.any(np.abs(zi - other_z[other_index_in]) < 0.001)
                             for zi in self_z])
    other_index = np.array([np.any(np.abs(zi - self_z[self_index]) < 0.001)
                             for zi in other_z])
    return self_index, other_index

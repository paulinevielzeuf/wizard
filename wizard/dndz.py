"""
Go from paircounts to phi(z)

.. module:: dndz
"""
from __future__ import print_function, division
import numpy as np
from astropy.cosmology import Planck15, WMAP9
import pandas as pd

def dndz(pairs, redshifts,
         dndz_path='', integrate='paircounts', w_estimator='LS',
         units='physical', redshift_column='Z', r_min=0.2, r_max=10, power=-1,
         z_min=0, z_max=1.0, cov_type='full',
         min_sep=0.01, max_sep=180, nbins=100,
         load_dndz=False, combine_keys=[],
         chunk_size=9999, time0=False):
    """Go from pairs to an estimate of the p(z) from clustering plus covariances

    Parameters
    ----------
    pairs : dictionary of dictionaries. pairs[DuDr][0] gets the 0th jackknife
            sample of the DuDr pairs in all the different bins.
    redshifts : redshifts of the reference samples, keyed so that
        redshifts[pairs['reference__{0}_region'.format(redshift_column)]]
        is the redshift that will go into phi(z) and any angular conversions
    dndz_path : string, location where we store the resulting dndz
    integrate : ['paircounts'=default, 'w'] what we integrate to marginalize
                out the theta term
    w_estimator : ['LS'=default, 'Natural'] estimator we use for going from
                  paircounts to w
    units : ['physical'=default, 'angular', 'comoving'] Units we do integration
    redshift_column : string [default: 'Z'], when using non-angular units, where
            do we look to find the central redshift
    r_min : float [default: 0.2] minimum range we integrate over
    r_max : float [default: 10] maximum range we integrate over
    power : float [default: -1] power by which we weight when doing integration
    z_min : float [default: 0] minimum redshift used
    z_max : float [default: 1] maximum redshift used
    cov_type : ['full'=default, 'tomo', 'diagonal']. If full, do jackknife for
               _all_ entries. If tomo, do jackknife covariance only for each
               bin that is not the reference redshift bin. If diagonal, only
               compute the diagonal component.
    min_sep : in arcmin, minimum separation we consider [default: 0.01]
    max_sep : in arcmax, maximum separation we consider [default: 180]
    nbins   : number of bins used [default: 100]
    load_dndz : bool [Default True] If True, just load dndz.
    combine_keys : list of str [default: empty] In addition to jackknife and
                   redshift region, we group by all other available reasonable
                   keys. Any entries here are summed over.
    chunk_size : int [Default: 9999] number of pairs we process at a time.
                 Useful if there are memory issues


    Returns
    -------
    df : pandas dataframe of the resultant parameters
    cov : covariance matrix
    stats : summary statistics
    tags : the names by which the paircounts were potentially grouped

    Notes
    -----
    This is pretty messy and needs to be cleaned up!

    TODO: double check the calculation of covariance, mean, width, errors on mean and width!
    """
    if time0: from time import time

    if load_dndz:
        df = pd.read_hdf(dndz_path, 'results')
        cov = pd.read_hdf(dndz_path, 'cov').values
        stats = pd.read_hdf(dndz_path, 'stats')
        tags = pd.read_hdf(dndz_path, 'tags').values
        return df, cov, stats, tags

    # TODO: Accomodate autocorrelation values

    redshift_str = 'reference__{0}_region'.format(redshift_column)
    groupby_ignore = [redshift_str]

    # TODO: This is pretty gross and could be cleaned up
    if integrate == 'paircounts':
        # integrate, then make w
        integrated_pairs = {}
        for kind in pairs:
            if time0:
                print(integrate, kind, time() - time0)
            # need to deal with the jackknife regions
            integrated_pairs[kind] = {}
            for jackknife in pairs[kind]:
                # if time0:
                #     print(integrate, kind, jackknife, time() - time0)

                arr = integrate_pairs(
                    pairs[kind][jackknife], redshifts[pairs[kind][jackknife][redshift_str]],
                    integrate_key='npairs',
                    r_min=r_min, r_max=r_max, power=power, units=units,
                    max_sep=max_sep, min_sep=min_sep, nbins=nbins,
                    chunk_size=chunk_size)
                integrated_pairs[kind][jackknife] = arr
                integrated_pairs[kind][jackknife]['jackknife'] = jackknife
            arr = pd.concat(integrated_pairs[kind], ignore_index=True)
            # if time0:
            #     print(integrate, kind, 'concated', time() - time0)

            # do combine_keys
            # these are always summed over
            sum_keys = ['npairs', 'n_reference', 'n_unknown']
            # these are always last
            save_keys = ['jackknife', redshift_str]
            if len(combine_keys) > 0:
                combine_groupby_keys = [key for key in arr.columns if key not in sum_keys + save_keys + combine_keys] + save_keys
                arr2 = arr.groupby(combine_groupby_keys).apply(combine_function, sum_keys=sum_keys)
            else:
                # in effect, combine_keys becomes every key that isn't a sum_key
                combine_groupby_keys = [key for key in arr.columns if key not in sum_keys + save_keys] + save_keys
                # but in that case we really ought to just reset the index
                arr2 = arr.set_index(combine_groupby_keys)

            # note that by having the indices, we make it easier to handle the
            # different paircounts (in case RuRr had a different number of
            # jackknives than DuDr for example)
            integrated_pairs[kind] = arr2
        # make w
        w = angular_correlation(integrated_pairs, w_estimator=w_estimator)

        # make df. we'll use DuDr since you almost always want that anyways
        df = integrated_pairs['DuDr'].copy().drop(
            ['npairs', 'n_reference', 'n_unknown'], axis=1)
        # add w
        df['w'] = w
        # add npairs and numbers (n_reference -> NRu or such)
        df['DuDr'] = integrated_pairs['DuDr']['npairs']
        df['NDu'] = integrated_pairs['DuDr']['n_unknown']
        df['NDr'] = integrated_pairs['DuDr']['n_reference']
        df['RuRr'] = integrated_pairs['RuRr']['npairs']
        df['NRu'] = integrated_pairs['RuRr']['n_unknown']
        df['NRr'] = integrated_pairs['RuRr']['n_reference']
        # these get a 2 because they are duplicates...
        df['RuDr'] = integrated_pairs['RuDr']['npairs']
        df['NRu2'] = integrated_pairs['RuDr']['n_unknown']
        df['NDr2'] = integrated_pairs['RuDr']['n_reference']
        df['DuRr'] = integrated_pairs['DuRr']['npairs']
        df['NDu2'] = integrated_pairs['DuRr']['n_unknown']
        df['NRr2'] = integrated_pairs['DuRr']['n_reference']

        # flatten index out
        df = df.reset_index()

        # now also get rid of any w that are clearly from `fake' jackknife regions.
        # for example, if we grouped by region on the sky (say SDSS north vs south)
        # then the northern pixels are not really jackknives for the southern
        # portion of the survey, and so on.
        subkeys = ['NDu', 'NDr', 'NRu', 'NRr']#, 'NDu2', 'NRu2', 'NRr2', 'NDr2']
        df = df.set_index([key for key in combine_groupby_keys if key != 'jackknife'])
        diff = df[subkeys] - df[df['jackknife'] == -1][subkeys]
        zeros = (diff != 0).any(axis=1)
        jackknives = df['jackknife'] - df[df['jackknife'] == -1]['jackknife'] - 1
        diff['jackknife'] = jackknives
        diff['zeros'] = zeros
        diff['keep'] = (diff['jackknife'] == -1) | (diff['zeros'])
        # heh this is a bit silly
        # either jackknife == -1, OR
        # diff != 0 for at least one index
        df_old = df.reset_index().set_index(combine_groupby_keys)[diff.reset_index().set_index(combine_groupby_keys)['keep']]


        # # TODO: something is wrong at this point!
        # # IE I am pretty sure I should be keeping points here that are not being kept!
        # conds_old = (((df['jackknife'] == -1).reset_index()['jackknife']) | ((diff != 0).any(axis=1).reset_index()[0]))
        # df_old = df.reset_index()[conds_old].reset_index()
        # dfo = dfi.reset_index()[conds_old].reset_index()
        # import ipdb; ipdb.set_trace()
        df, df_old = df_old.reset_index().copy(), df.reset_index().copy()

        # update groupby ignore
        groupby_ignore += ['w',
                           'DuDr', 'NDu', 'NDr',
                           'RuRr', 'NRu', 'NRr',
                           'RuDr', 'NRu2', 'NDr2',
                           'DuRr', 'NDu2', 'NRr2',]
        groupby_keys = [key for key in df.columns if key not in groupby_ignore]

    elif integrate == 'w':
        # make w, then integrate
        pass

    # drop all the w that are not finite at this point...
    # first convert nan to inf
    df.replace({'w': {np.inf: np.nan, -np.inf: np.nan}}, inplace=True)
    df.dropna(axis=0, subset=['w'], inplace=True)



    # also drop all redshifts outside our zmin and zmax
    df['Z'] = redshifts[df[redshift_str]]
    df = df[(df['Z'] > z_min) & (df['Z'] < z_max)]

    # maybe DON'T group?
    tomo_groupby = [key for key in groupby_keys if key != 'jackknife']
    # gotta handle what to do if there are no groupby keys
    if len(tomo_groupby) == 0:
        # put in a fake key so we can just keep using the pandas groupby
        tomo_groupby = ['_group']
        df['_group'] = 1
        combine_groupby_keys = tomo_groupby + combine_groupby_keys
        groupby_keys = tomo_groupby + groupby_keys

    ######
    # keys:
    # groupby_keys
    # combine_groupby_keys
    # tomo_groupby
    # size_keys
    # size_post_keys

    # left:
    # groupby_keys = size_keys
    # tomo_groupby = size_post_keys
    # combine_groupby = tomo_groupby + [jackknife, redshift_str] = size_keys + [redshift_str]

    # right:
    # groupby_keys = size_keys
    # tomo_groupby = size_post_keys
    # combine_groupby = tomo_groupby + [jackknife, redshift] = size + [redshift_str]

    # size_keys = [key for key in combine_groupby_keys if key not in [redshift_str]]
    # size_post_keys = [key for key in combine_groupby_keys if key not in save_keys]

    # print(groupby_keys)
    # print(combine_groupby_keys)
    # print(tomo_groupby)
    # print(size_keys)
    # print(size_post_keys)
    # import ipdb; ipdb.set_trace()

    # 
    # sizes = df.reset_index().groupby(size_keys).size()
    # sizes = sizes.reset_index().set_index(size_post_keys)
    # # TODO: now drop any that do not equal respective -1 index
    # conds = ((sizes['jackknife'] == -1) | (sizes[0] - sizes[sizes['jackknife'] == -1][0] == 0))
    # acceptable_jackknives = sizes[conds].reset_index()
    # # this is kinda nuts but we can use a right-join to ditch any samples that don't pass muster
    # df = df.set_index(size_keys).join(acceptable_jackknives.set_index(size_keys), how='right', rsuffix='DROP_').reset_index().drop([0], axis=1)

    # TODO: This definitely needs to be cleaned
    # TODO: This also doesn't seem to be correct.
    # for the jackknives we now need to drop any set which don't have a
    # complete Z in the range, where by 'complete Z' I mean it spans the same
    # number of points as the full sample
    sizes_in = df.reset_index().groupby(groupby_keys).size()
    # catch any size_post_keys
    sizes = sizes_in.reset_index().set_index(tomo_groupby)
    # TODO: now drop any that do not equal respective -1 index
    conds = ((sizes['jackknife'] == -1) | (sizes[0] - sizes[sizes['jackknife'] == -1][0] == 0))
    acceptable_jackknives = sizes[conds].reset_index()
    # this is kinda nuts but we can use a right-join to ditch any samples that don't pass muster
    df_dropped = df.set_index(groupby_keys).join(acceptable_jackknives.set_index(groupby_keys), how='right', rsuffix='DROP_').reset_index().drop([0], axis=1)
    # import ipdb; ipdb.set_trace()
    df, df_dropped = df_dropped.copy(), df.copy()

    # add a new basis that will be the key we sort by
    # create a new index based on allowed values of each groupby_key
    # this new index needs to have as the very last axis the redshift string
    # TODO: when matching things do I need to match the _basis now?
    full_groupby = tomo_groupby + [redshift_str]
    basis = np.cumprod([1] + [len(np.unique(df[key])) for key in full_groupby][::-1])[:-1]
    df['_basis'] = (df[full_groupby[::-1]] * basis).sum(axis=1)

    # calculate normalization
    # don't forget to exclude any points outside of z_min and z_max!
    df = df.groupby(groupby_keys).apply(trapezoid_groupby, key='w', outkey='norm')

    df['phi'] = df['w'] / df['norm']

    # TODO: stats and covariance go into a second function?

    # calculate mean redshift
    df['w_z'] = df['w'] * df['Z'] / df['norm']
    df = df.groupby(groupby_keys).apply(trapezoid_groupby, key='w_z', outkey='z_mean')

    # calculate width of distribution
    df['w_width**2'] = (df['Z'] - df['z_mean']) ** 2 * df['w'] / df['norm']
    df = df.groupby(groupby_keys).apply(trapezoid_groupby, key='w_width**2',
                                        outkey='z_width**2')
    df['z_width'] = np.sqrt(df['z_width**2'])

    # calculate covariance

    dfsansone = df[df['jackknife'] != -1]
    dfone = df[df['jackknife'] == -1]

    if ((cov_type == 'full') or (cov_type == 'diagonal')):
        # TODO: Need to confirm this doesn't break with multiple kwargs
        # basis = np.cumprod([1] + [len(np.unique(dfsansone[key])) for key in full_groupby][::-1])[:-1]
        # dfsansone['_basis'] = (dfsansone[full_groupby[::-1]] * basis).sum(axis=1)
        # asdf = dfsansone[dfsansone['jackknife'] == 462][full_groupby + ['_basis', 'jackknife']].sort('_basis')
        # import ipdb; ipdb.set_trace()
        # cov = make_cov_array(dfsansone, val_key='w', redshift_key='_basis', jackknife_key='jackknife')
        # using pivot is faster
        cov = _make_cov_array(dfsansone, val_key='w', redshift_key='_basis', jackknife_key='jackknife')

    elif cov_type == 'tomo':
        # TODO: I fixed 'full' by moving stuff around -- do I need to do the same here with the tomo?
        # TODO: Need to make sure the ordering is correct
        # cov_grouped = dfsansone.groupby(tomo_groupby).apply(make_cov_array, val_key='w', redshift_key=redshift_str, jackknife_key='jackknife')
        # using pivot is faster
        cov_grouped = dfsansone.groupby(tomo_groupby).apply(_make_cov_array, val_key='w', redshift_key=redshift_str, jackknife_key='jackknife')
        # blow up into full covariance matrix
        cov = np.zeros((len(dfone), len(dfone)), dtype=np.float)
        corner_indx = 0
        for cov_i in cov_grouped.values:
            shape = cov_i.shape[0]  # square
            cov[corner_indx:corner_indx + shape, corner_indx:corner_indx + shape] = cov_i
            corner_indx += shape

    if cov_type == 'diagonal':
        # convert to diagonal
        cov = np.eye(cov.shape[0]) * np.diag(cov)

    # correct cov normalization
    # TODO: Should I just have made the covariance from the normalized phi?
    norm_x, norm_y = np.meshgrid(dfone['norm'].values, dfone['norm'].values)
    cov = cov / (norm_x * norm_y)

    # while we're at it get the errors in norm, z_mean, z_width
    # first gotta get rid of all these extra redshift entries
    dfsansone_onez = dfsansone.groupby(groupby_keys).first().reset_index()
    N_jackknife = len(dfsansone['jackknife'].unique())
    err = (dfsansone_onez.groupby(tomo_groupby).agg(np.std)[['norm', 'z_mean', 'z_width']] * np.sqrt(N_jackknife - 1)).reset_index()

    # stupidly we have to redo the groupby to match
    dfone = dfone.set_index(tomo_groupby)
    err = err.set_index(tomo_groupby)
    for key in ['norm', 'z_mean', 'z_width']:
        dfone[key + '_error'] = err[key]
    dfone.reset_index()

    # add diagonal errorbars
    dfone['phi_err'] = np.sqrt(np.diag(cov))

    # re index
    dfone = dfone.reset_index()
    # print warning if npoints > njackknives
    if len(dfone) >= len(df['jackknife'].unique()) - 1:
        print('Warning! We have {0} data points but only {1} jackknife samples!'.format(len(dfone), len(df['jackknife'].unique()) - 1))
    # purge columns

    dfone = dfone.drop(['jackknife', 'w_width**2', 'w_z'], axis=1)

    stats = dfone.groupby(tomo_groupby).first()[['norm', 'norm_error', 'z_mean', 'z_mean_error', 'z_width', 'z_width_error']].reset_index()

    if len(dndz_path) > 0:
        # save
        # dfone
        dfone.to_hdf(dndz_path, 'results')

        # cov
        pd.DataFrame(cov).to_hdf(dndz_path, 'cov')

        # stats
        stats.to_hdf(dndz_path, 'stats')

        # groupby keys (for reference)
        pd.Series(tomo_groupby).to_hdf(dndz_path, 'tags')

        # for reference, save the jackknife samples
        dfsansone.to_hdf(dndz_path, 'jackknife')

    return dfone, cov, stats, tomo_groupby

def integrate_pairs(pairs, z, integrate_key='npairs',
                    r_min=0.2, r_max=10, power=-1, units='physical',
                    max_sep=180, min_sep=0.01, nbins=100, chunk_size=9999):
    # basic goal: take a recarray where one entry is (1000i8) and others are
    # what they are, and convert the 1000 to a single entry
    df = pd.DataFrame({key: pairs[key]
                       for key in pairs.dtype.names if key != integrate_key})
    npairs = pairs[integrate_key]

    # generate logr
    logr = generate_logr(max_sep, min_sep, nbins)
    radius, weights = return_weight(logr, z, r_min, r_max, power, units)

    # integrate
    npairs_integrated = trapezoid(npairs, weights, radius, chunk_size)

    # put back into df
    df[integrate_key] = npairs_integrated

    return df


def angular_correlation(pairs, w_estimator='LS'):
    if w_estimator == 'LS':
        DuDr = pairs['DuDr']['npairs'] / (
            pairs['DuDr']['n_reference'] *
            pairs['DuDr']['n_unknown'])
        RuRr = pairs['RuRr']['npairs'] / (
            pairs['RuRr']['n_reference'] *
            pairs['RuRr']['n_unknown'])
        RuDr = pairs['RuDr']['npairs'] / (
            pairs['RuDr']['n_reference'] *
            pairs['RuDr']['n_unknown'])
        DuRr = pairs['DuRr']['npairs'] / (
            pairs['DuRr']['n_reference'] *
            pairs['DuRr']['n_unknown'])
        results = (DuDr - DuRr - RuDr + RuRr) / (RuRr)
    elif w_estimator == 'Natural':
        DuDr = pairs['DuDr']['npairs'] / (
            pairs['DuDr']['n_reference'] *
            pairs['DuDr']['n_unknown'])
        RuRr = pairs['RuRr']['npairs'] / (
            pairs['RuRr']['n_reference'] *
            pairs['RuRr']['n_unknown'])
        results = (DuDr - RuRr) / (RuRr)
    return results

def angular_correlation_array(pairs, jackknife, w_estimator='LS',
                              max_sep=180, min_sep=0.01, nbins=100):
    # like angular correlation, but deal with the different shapes.
    # also return logr
    if w_estimator == 'LS':
        DuDr = pairs['DuDr'][jackknife]['npairs'] / (
            pairs['DuDr'][jackknife]['n_reference'][:, np.newaxis] *
            pairs['DuDr'][jackknife]['n_unknown'][:, np.newaxis])
        RuRr = pairs['RuRr'][jackknife]['npairs'] / (
            pairs['RuRr'][jackknife]['n_reference'][:, np.newaxis] *
            pairs['RuRr'][jackknife]['n_unknown'][:, np.newaxis])
        RuDr = pairs['RuDr'][jackknife]['npairs'] / (
            pairs['RuDr'][jackknife]['n_reference'][:, np.newaxis] *
            pairs['RuDr'][jackknife]['n_unknown'][:, np.newaxis])
        DuRr = pairs['DuRr'][jackknife]['npairs'] / (
            pairs['DuRr'][jackknife]['n_reference'][:, np.newaxis] *
            pairs['DuRr'][jackknife]['n_unknown'][:, np.newaxis])
        results = (DuDr - DuRr - RuDr + RuRr) / (RuRr)
    elif w_estimator == 'Natural':
        DuDr = pairs['DuDr'][jackknife]['npairs'] / (
            pairs['DuDr'][jackknife]['n_reference'][:, np.newaxis] *
            pairs['DuDr'][jackknife]['n_unknown'][:, np.newaxis])
        RuRr = pairs['RuRr'][jackknife]['npairs'] / (
            pairs['RuRr'][jackknife]['n_reference'][:, np.newaxis] *
            pairs['RuRr'][jackknife]['n_unknown'][:, np.newaxis])
        results = (DuDr - RuRr) / (RuRr)

    logr = generate_logr(max_sep, min_sep, nbins)
    return results, logr

def trapezoid(values, weights, radius, chunk_size=9999):
    """Trapezoidal integration

    Parameters
    ----------
    values : array(Nsamples, Nbins)
    weights : array(Nsamples, Nbins)
    radius : array(Nsamples, Nbins)

    chunk_size : int [Default: 9999] number of pairs we process at a time.
                 Useful if there are memory issues

    Returns
    -------
    integrated values (Nsamples)

    """
    # must do this in chunk_size because of memory issues
    n_iter = int(values.shape[0] / chunk_size) + 1
    integrated_counts = np.zeros(values.shape[0])

    for ith in xrange(n_iter):
        lower = ith * chunk_size
        upper = lower + chunk_size

        # NOTE: if there are memory issues, do the weight and radius creation here!
        values_ith = values[lower:upper]
        weights_ith = weights[lower:upper]
        radius_ith = radius[lower:upper]

        vw = values_ith * weights_ith
        integrated_counts_ith = np.sum(0.5 *
            (radius_ith[:, 1:] - radius_ith[:, :-1]) *
            (vw[:, 1:] + vw[:, :-1]),
            axis=1)
        integrated_counts[lower:upper] = integrated_counts_ith
    return integrated_counts

def trapezoid_groupby(df, key, outkey, zkey='Z'):
    df[outkey] = trapezoid_groupby_func(df, key, zkey)
    return df

def trapezoid_groupby_func(df, key, zkey='Z'):
    # trapezoidal integration with pandas dataframes
    # TODO: Assumes redshifts are sorted, but that is easy to add
    # TODO: am I going to need the z max and min for proper normalization?
    z = df[zkey]
    val = df[key]
    integrated = (0.5 * (z - z.shift()) * (val + val.shift())).dropna().sum()
    return integrated

def generate_logr(max_sep=180, min_sep=0.01, nbins=100):
    """Generate log(r) where r is separation in arcmin

    Parameters
    ----------
    min_sep : in arcmin, minimum separation we consider [default: 0.01]
    max_sep : in arcmax, maximum separation we consider [default: 180]
    nbins   : number of bins used [default: 100]

    Returns
    -------
    logr    : equally spaced bins in logspace. exp(logr) = r

    """
    sep_units = np.pi / 180 / 60  # arcmin to radians
    log_sep_units = np.log(sep_units)
    min_sep_radian = min_sep * sep_units  # arcmin to radians
    max_sep_radian = max_sep * sep_units  # arcmin to radians
    bin_size = np.log(max_sep_radian / min_sep_radian) / nbins
    # This makes nbins evenly spaced entries in log(r) starting with 0 with step bin_size
    logr = np.linspace(start=0, stop=nbins * bin_size, num=nbins,
                       endpoint=False)
    # Offset by the position of the center of the first bin.
    logr += np.log(min_sep_radian) + 0.5 * bin_size
    # And correct the units:
    logr -= log_sep_units

    return logr

def return_weight(logr, redshifts, r_min=0.2, r_max=10, power=-1, units='physical'):
    """Weight function f by some powerlaw for r

    Parameters
    ----------
    logr : array
        Array of floats for the log(theta) where theta is in arcmin

    redshifts : list, default empty
        If not empty, return len(redshifts) copies of the weight function
        and the radius in physical units

    r_min : float
        Minimum range we are interested in. Sets normalization.

    r_max : float
        Maximum range we are interested in. Sets normalization.

    power : float
        What power for the powerlaw weighting are we using?
        If power == 0, just return 1.

    units : ['physical'=default, 'angular', 'comoving'] Units we do integration
        angular: integrate in arcmin
        physical: integrate in physical distance ie Da(z) * theta
        comoving: (1 + z) * Da(z) * theta

    Returns
    -------
    radius : array
        Radii in either angle or physical units. Can gain depth if redshifts
        is not empty

    weight : array
        Weights normalized so that int_rmin^rmax dr w(r) = 1
    """
    if units == 'physical':
        conversion = angular_diameter_mpc_per_arcmin(redshifts)
    elif units == 'comoving':
        conversion = (1 + redshifts) * angular_diameter_mpc_per_arcmin(redshifts)

    # blow up logr for each redshift
    r = np.ones((len(redshifts), len(logr))) * np.exp(logr)[np.newaxis]
    if units in ['physical', 'comoving']:
        r *= conversion[:, np.newaxis]

    # calculate power and normalization for powerlaw weighting
    num = np.where((r > r_min) & (r < r_max), r ** power, 0)
    if power == -1:
        denom = np.log(r_max / r_min)
    else:
        denom = (r_max ** (power + 1) - r_min ** (power + 1)) / (power + 1)
    weights = num / denom

    return r, weights

def angular_diameter_mpc_per_arcmin(z, cosmo='WMAP9'):
    if cosmo == 'WMAP9':
        cosmo = WMAP9
    elif cosmo == 'Planck15':
        cosmo = Planck15
    d_a = cosmo.angular_diameter_distance(z).value
    mpc_per_arcmin = d_a * np.pi / 10800
    return mpc_per_arcmin

def _make_cov_array(df, val_key='w', redshift_key='reference__Z_region', jackknife_key='jackknife'):
    # too clever by half.
    # note that Cov_JK -> cov_ij = E((Xi - E(Xi))(Xj - E(Xj))) len(Nij - 1)
    N_jackknife = len(df[jackknife_key].unique())
    cov = df.pivot(index=jackknife_key, columns=redshift_key, values=val_key).cov().values * (N_jackknife - 1)

    return cov

def make_cov_array(df, val_key='w', redshift_key='reference__Z_region', jackknife_key='jackknife'):
    values = df[val_key].values
    redshifts = df[redshift_key].values
    jackknives = df[jackknife_key].values

    cov = np.zeros((len(np.unique(redshifts)), len(np.unique(redshifts))))
    for i, redshift_i in enumerate(np.unique(redshifts)):
        # select values
        values_i = values[redshifts == redshift_i]
        jackknives_i = jackknives[redshifts == redshift_i]
        for j, redshift_j in enumerate(np.unique(redshifts)):
            # select the samples
            values_j = values[redshifts == redshift_j]
            jackknives_j = jackknives[redshifts == redshift_j]
            # get the overlap
            conds_i = np.in1d(jackknives_i, jackknives_j)
            conds_j = np.in1d(jackknives_j, jackknives_i)
            values_ij = values_i[conds_i]
            values_ji = values_j[conds_j]
            mean_ij = np.mean(values_ij)
            mean_ji = np.mean(values_ji)
            N_ij = conds_i.sum()
            cov[i, j] = np.mean((values_ij - mean_ij) * (values_ji - mean_ji)) * (N_ij - 1)

    return cov

def combine_function(df,
                     sum_keys=[],
                     first_keys=[], func_keys={}):
    """
    func_keys = {'df_out_column': [function, df_in_key]}
    """
    dfs = {}
    for key in sum_keys:
        dfs[key] = df[key].sum()
    for key in first_keys:
        dfs[key] = df[key].iloc[0]
    for func_key in func_keys:
        function, key = func_keys[func_key]
        dfs[func_key] = function(df[key])
    series = pd.Series(dfs)
    return series

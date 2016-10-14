"""
Given sets of pairs, combine and create different jackknife assessments, too

.. module:: combine
"""

from __future__ import print_function, division

import numpy as np

from .dataset import index_to_radec, radec_to_index

from .paircounts import interpret_bins_integerized, generate_combinations

def combine(reference_catalogs, reference_bins, unknown_catalogs, unknown_bins, pairs, healpix_nside=32, jackknife_nside=16, combine_path='', load_combine=False,
            time0=0):
    """Given a list of paircounts, combine together with jackknife regions

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
    pairs : dictionary keyed by entry kind and healpixel pairs
    healpix_nside : int [default: 32] nside we divide the paircounts into for
                    later jackknifing. This therefore also controls the number
                    of files saved, and thus the harddisk space of the output
    jackknife_nside : int [default: 16] nside we use for jackknifing. Sets the
                      number of jackknives in our sample
    combine_path : string. If given, save the combined sets into this directory
    load_combine : bool [default: False] If True, just load combined counts

    Returns
    -------
    combined_pairs : keyed by entry kind, then jackknife sample (-1 == full)

    Notes
    -----
    TODO: also be able to do auto correlation stuff
    """

    if time0: from time import time

    combined_pairs = {}
    if load_combine:
        kinds = pairs.keys()
        for kind in kinds:
            combined_pairs[kind] = np.load('{0}{1}.npy'.format(combine_path, kind)).item()
    else:
        combined_pairs = {}
        kinds = pairs.keys()

        # do the combining
        for kind in kinds:
            # determine the jackknives available.
            # this SHOULD be the same across kinds, but if it isn't, we'll find
            # out here!!
            healpixels = np.unique([i[:2] for i in pairs[kind]])
            jackknives = np.unique(radec_to_index(*index_to_radec(healpixels,
                                                                  healpix_nside),
                                                  nside=jackknife_nside)).tolist()
            combined_pairs[kind] = {}

            if time0:
                print(kind, len(healpixels), len(jackknives), time() - time0)

            for row_i, row in enumerate(pairs[kind]):
                if time0:
                    print(kind, row_i, row_i / len(pairs[kind]), len(pairs[kind]), time() - time0)
                if type(row[-1]) == str:
                    r_hpix, u_hpix, pairs_i_path = row
                    pairs_i = np.load(pairs_i_path)
                else:
                    r_hpix, u_hpix, pairs_i = row
                # assign jackknife healpixel
                r_jackknife = radec_to_index(*index_to_radec(r_hpix, healpix_nside),
                                             nside=jackknife_nside)
                u_jackknife = radec_to_index(*index_to_radec(u_hpix, healpix_nside),
                                             nside=jackknife_nside)

                # add to all but jackknife healpixel
                if row_i == 0:
                    # do the n_reference, n_unknown
                    reference_bins_integerized = interpret_bins_integerized(reference_bins, 'reference__{0}_region')
                    unknown_bins_integerized = interpret_bins_integerized(unknown_bins, 'unknown__{0}_region')

                    # generate combinations of all the possible values
                    reference_bin_combinations = generate_combinations(reference_bins_integerized)
                    unknown_bin_combinations = generate_combinations(unknown_bins_integerized)
                    ref_len = len(reference_bin_combinations)
                    unk_len = len(unknown_bin_combinations)
                    # get the rcat and ucat
                    if 'Du' in kind:
                        ucat = unknown_catalogs[0]
                    elif 'Ru' in kind:
                        ucat = unknown_catalogs[1]
                    if 'Dr' in kind:
                        rcat = reference_catalogs[0]
                    elif 'Rr' in kind:
                        rcat = reference_catalogs[1]
                    for jackknife in [-1] + jackknives:
                        # make the array if it doesn't exist
                        combined_pairs[kind][jackknife] = pairs_i.copy()
                        # set the npairs to zero
                        combined_pairs[kind][jackknife]['npairs'] *= 0
                    # now go through the combos and assign
                    rcat['JACKKNIFE'] = radec_to_index(*index_to_radec(rcat['HPIX'], healpix_nside),
                                                    nside=jackknife_nside)
                    ucat['JACKKNIFE'] = radec_to_index(*index_to_radec(ucat['HPIX'], healpix_nside),
                                                    nside=jackknife_nside)
                    usize = ucat.groupby(unknown_bin_combinations[0].keys() + ['JACKKNIFE']).size().reset_index()
                    rsize = rcat.groupby(reference_bin_combinations[0].keys() + ['JACKKNIFE']).size().reset_index()
                    # now iterate
                    for entry_i, unknown_combos in enumerate(unknown_bin_combinations):
                        n_unknowns = usize[usize[unknown_combos.keys()].isin(unknown_combos).all(axis=1)]
                        for entry_j, reference_combos in enumerate(reference_bin_combinations):
                            n_references = rsize[rsize[reference_combos.keys()].isin(reference_combos).all(axis=1)]
                            # assign
                            entry_ij = entry_i * ref_len + entry_j
                            for jackknife_k, jackknife in enumerate([-1] + jackknives):
                                entry_ijk = entry_i * ref_len + entry_j + ref_len * unk_len * jackknife_k
                                if time0 and (entry_ijk % 1000 == 0):
                                    print('numbers', kind, entry_i, entry_j, jackknife, unk_len, ref_len, len(jackknives) + 1, time() - time0)

                                n_reference = n_references[n_references['JACKKNIFE'] != jackknife][0].sum()
                                n_unknown = n_unknowns[n_unknowns['JACKKNIFE'] != jackknife][0].sum()
                                combined_pairs[kind][jackknife]['n_reference'][entry_ij] = n_reference
                                combined_pairs[kind][jackknife]['n_unknown'][entry_ij] = n_unknown
                for jackknife in [-1] + jackknives:

                    # jackknives are all about the region you are excluding
                    # TODO: potentially only throw out based on r_jackknife...
                    if (jackknife != r_jackknife) & (jackknife != u_jackknife):
                        # only add to pairs, n_reference, n_unknown
                        # TODO: When I do the auto correlations, what other keys?
                        combined_pairs[kind][jackknife]['npairs'] += pairs_i['npairs']
                        # if r_hpix == u_hpix:
                        #     for key in ['n_reference', 'n_unknown']:
                        #         combined_pairs[kind][jackknife][key] += pairs_i[key]

            # save
            if len(combine_path) > 0:
                np.save('{0}{1}.npy'.format(combine_path, kind), combined_pairs[kind])

    return combined_pairs

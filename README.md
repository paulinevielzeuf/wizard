# wizard
Magically inferring redshift distributions from angular cross-correlations!

So you want to calculate redshift distributions, but you want to try something
other than photo-$z$ algorithms. This code is meant to make creating redshift
distributions from clustering measurements easy. Here is all you need:

  1.    A set of galaxies whose redshift distribution you would like to infer.
        These are the 'Du' objects.
  2.    A set of randoms associated with these galaxies. These randoms ideally
        plot out the survey area of where 'Du' galaxies could have been
        measured. These are the 'Ru' objects.
  3.    A set of reference galaxies whose redshifts you are pretty confident
        in. These are the 'Dr' objects.
  4.    A set of randoms associated with the reference galaxies. These are the
        'Rr' objects.

`wizard` (and clustering redshift algorithms generically) take advantage of the
fact that baryonic matter traces the underlying dark matter web. So, objects
at the same redshift are causally connected (== correlated on the sky) because
they follow the same underlying dark matter distribution, while objects at
different redshifts are not (== uncorrelated on the sky), modulo things like
magnification effects from gravitational lensing. Therefore, if we measure the
excess angular correlation between the unknown sample and reference galaxies
in a small redshift slice, the excess amount of correlation will be
proportional to the number of unknown objects at that redshift.

`wizard` will go from the four input catalogs and return a redshift distribution
with errorbars estimated from jackknife sampling on the sky. It also can
compare against other redshift distributions -- say from photo-$z$ -- and try
to come up with some sort of corrections for both sets in order to reach a
happy agreement.

Configuration is done via a yaml file. An example with reasonable values may be
found in the examples folder.

## The Pipeline

### Dataset

`wizard` first takes the input galaxies and adds appropriate columns for later
processing. The columns you specify here must be in both the data and randoms,
though it is not required that they be shared between reference and unknown.

### Paircounts

`wizard` next calculates the paircounts, organized by kind of pairing, and
healpix pairs. This is by far the slowest part, so see the Example for some
ways to take advantage of the embarrassingly parallel nature of this with your
batch farm.

In each healpixel pair, paircounts are calculated with treecorr for all the
different bin combinations presented in the dataset section. When there are no
pairs, then the 'npairs' column is all zeros.

### Combine

`wizard` needs to combine the paircounts from the healpixels (and other
potential combinations) for future analysis here. We use jackknife regions
based on healpix to calculate future covariances, so that set of sums is also
accomplished here.

### dndz

`wizard` needs to go from paircounts to an estimate of the redshift
distribution.

### analyze

`wizard` can also analyze the results, comparing against some other catalog
redshift distribution.

## Example -- COSMOS Galaxies


## Requirements

numpy
pandas
treecorr
fitsio
astropy
healpy

## The Backronym

We don't know yet :(

## License

MIT

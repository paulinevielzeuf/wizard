"""
Compare the redshift distribution of a clustering redshift and some other sample

.. module:: compare

"""

from __future__ import print_function, division
import emcee
import numpy as np

def _chi2(pz_other, pz, cov):
    from scipy import linalg
    cv_chol = linalg.cholesky(cov, lower=True)
    cv_sol = linalg.solve(cv_chol, pz - pz_other, lower=True)
    chi2_val = np.sum(cv_sol ** 2)
    return chi2_val

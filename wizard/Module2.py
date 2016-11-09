#this is the module that compute the observables shear/wtheta/gammart



import os, sys
import numpy as np
import matplotlib.pylab as plt
from astropy.table import Table, vstack, hstack
from  config_module1 import *
#%matplotlib inline


###############################################
################### functions #################
###############################################

def covariance_matrix_jck(cov,jk_r,len):
  #  Covariance estimation
  average=np.zeros((len,len))
  cov_jck=np.zeros((len,len))
  err_jck=np.zeros((len,len))
  for kk in range(jk_r):
    average+=cov[kk]
  average=average/(jk_r)
  for kk in range(jk_r):
    cov_jck+=(-average+cov[kk])*(-average+cov[kk])
  err_jck=np.sqrt(cov_jck*(jk_r-1)/(jk_r))
  average=average*(jk_r)/(jk_r-1)
  return {'cov' : cov_jck,
          'err' : err_jck,
          'mean': average}


###############################################
################### Cosmosis  runs#################
###############################################

####################################################################
#**** running on the reconstructed distribution ********



nz_vec=np.loadtxt('input/N_z.txt')
np.savetxt('n_of_z.txt',nz_vec)

'''
if not os.path.exists('./rec_observables/'):
        run_cosmosis_command(cmmd+" Obs_est.ini")
        shutil.move('./out/','./rec_observables/')
'''

####################################################################
#***** running on jackknife regions *******




if N_rec_jk=='true':
 for i_count in range(njck_rec):
  print i_count+1
  nz_vec=np.loadtxt('input/N_z_jck.txt',usecols=(i_count+1,))
  np.savetxt('n_of_z.txt',nz_vec)
  '''
  if not os.path.exists('./rec_observables/'):
        run_cosmosis_command(cmmd+" Obs_est.ini")
        shutil.move('./out/','./rec_observables/')




  nz_vec=np.vstack((N_true[:,0],N_z[:,i_count+1])).T
  np.savetxt('n_of_z.txt',nz_vec)
  if not os.path.exists('./output'+labb+'/'+str(i_count+1)):
        run_cosmosis_command(cmmd+" demo6.ini")
        delete_folders(del_fold)
        shutil.move('./out/','./output'+labb+'/'+str(i_count+1))
        shutil.move('n_of_z.txt','./output'+labb+'/n_of_z_'+str(i_count)+'.txt')
  '''

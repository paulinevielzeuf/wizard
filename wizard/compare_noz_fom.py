# this module aims to prepare the Nz for which one wants to process FoM analysis including the possibility to create biased Nz (shifted , spread,..)



from __future__ import print_function
import numpy as np
import matplotlib.pylab as plt
from astropy.table import Table, vstack, hstack
import pandas as pa
from .dataset import catalog_to_histogram, _modify_indices_redshift
from .compare import _chi2
#  config_module1 import *
#%matplotlib inline








def compare_noz_fom(pcs,tags,unknown_catalog,redshift_bins,compare_noz_fom_path,redshift_catalog_columns,
                     z_min, z_max,path_rec,
                     shift_value,spreading_factor,N_comp_jk,N_rec_jk,
                     njck_comp,njck_rec,
                     shift,spread,Kurtosis,Skewness,
                     apply_to_compared_nz,apply_to_rec_nz,
                     time0=0,**kwargs):
  '''
  #######################################
  ########loading the files #############
  #######################################
  '''
  if time0: from time import time


  print ('THE module from pauline is running')
  gg=0
  ggjck=0


  unique_pc_tags = pcs[tags].drop_duplicates().reset_index(drop=True)
  tags_unknown = [tag for tag in tags if tag in unknown_catalog.columns]
      # make the different plots based on the tags
  for index, key in unique_pc_tags.iterrows():

          for ith, redshift_catalog_column in enumerate(redshift_catalog_columns):
              # convert catalog to corresponding pdf
              conds_cat = (unknown_catalog[tags_unknown] == key[tags_unknown]).all(axis=1).values
              z_cat = unknown_catalog[conds_cat][redshift_catalog_column]
              if 'W' in unknown_catalog:
                weights = unknown_catalog[conds_cat]['W']
              else:
                weights = None
              pz_cat, z_vals ,edges= catalog_to_histogram(z_cat, redshift_bins,
                                                    weights=weights,
                                                    z_min=z_min, z_max=z_max)


              vars()['Nt'+str(key[0])]=Table()
              vars()['Nt'+str(key[0])]['z']=z_vals
              vars()['Nt'+str(key[0])]['N']=pz_cat
              vars()['Nt'+str(key[0])]['sample']='compared'
              vars()['Nt'+str(key[0])]['biasing']='original'
              vars()['Nt'+str(key[0])]['tomo']=key[0]
              if gg>0:

                    N_compare_full=vstack([N_compare_full,vars()['Nt'+str(key[0])]])
              if gg==0 :
                    N_compare_full=vars()['Nt'+str(key[0])]
                    gg=gg+1


              print ('The compare Nz have been prepared for tomo_bin'+ str(key[0]))


              ################### JCK ###################
              if N_comp_jk :
                        print ('CREATING JK')
                        for ith, redshift_catalog_column in enumerate(redshift_catalog_columns):

                                for jk_count in np.unique(unknown_catalog['HPIX']):
                                    # convert catalog to corresponding pdf
                                    unknown_catalog_jck=unknown_catalog[unknown_catalog['HPIX']==jk_count]
                                    conds_cat_jck = (unknown_catalog_jck[tags_unknown] == key[tags_unknown]).all(axis=1).values

                                    z_cat_jck = unknown_catalog_jck[conds_cat_jck][redshift_catalog_column]
                                    if 'W' in unknown_catalog_jck:
                                        weights_jck = unknown_catalog_jck[conds_cat_jck]['W']
                                    else:
                                        weights_jck = None


                                    pz_cat_jck, z_vals_jck ,edges_jck= catalog_to_histogram(z_cat_jck, redshift_bins,
                                                                      weights=weights_jck,
                                                                      z_min=z_min, z_max=z_max)


                                    vars()['Ntjck'+str(key[0])]=Table()
                                    vars()['Ntjck'+str(key[0])]['z']=z_vals_jck
                                    vars()['Ntjck'+str(key[0])]['N']=pz_cat_jck
                                    vars()['Ntjck'+str(key[0])]['sample']='compared'
                                    vars()['Ntjck'+str(key[0])]['biasing']='original'
                                    vars()['Ntjck'+str(key[0])]['tomo']=key[0]
                                    vars()['Ntjck'+str(key[0])]['jackknife']=jk_count
                                    if ggjck>0:

                                        N_compare_full_jck=vstack([N_compare_full_jck,vars()['Ntjck'+str(key[0])]])
                                    if ggjck==0 :
                                        N_compare_full_jck=vars()['Ntjck'+str(key[0])]
                                        ggjck=ggjck+1

              print ('The compare Nz jackknife have been prepared')

              ####### recovered nz ################
              Nr_f = pa.HDFStore(path_rec)
              Nr_f=Nr_f['results']
              conds_cat_rec = (Nr_f[tags_unknown] == key[tags_unknown]).all(axis=1).values
              vars()['Nr'+str(key[0])]=Table()
              vars()['Nr'+str(key[0])]['z']=Nr_f[conds_cat_rec]['Z']
              vars()['Nr'+str(key[0])]['N']=Nr_f[conds_cat_rec]['phi']
              vars()['Nr'+str(key[0])]['sample']='reconstructed'
              vars()['Nr'+str(key[0])]['biasing']='original'
              vars()['Nr'+str(key[0])]['tomo']=key[0]
              N_compare_full=vstack([N_compare_full,vars()['Nr'+str(key[0])]])
              if N_rec_jk :
                Nr_f_jck = pa.HDFStore(path_rec)
                Nr_f_jck = Nr_f_jck['jackknife']
                conds_cat_rec = (Nr_f_jck[tags_unknown] == key[tags_unknown]).all(axis=1).values
                Nr_f_jck=Nr_f_jck[conds_cat_rec]
                for jk_count in np.unique(Nr_f_jck['jackknife']):
                        Nr_f_jck_out=Nr_f_jck[Nr_f_jck['jackknife']==jk_count]
                        vars()['Nrjck'+str(key[0])]=Table()
                        vars()['Nrjck'+str(key[0])]['z']=Nr_f_jck_out['Z']
                        vars()['Nrjck'+str(key[0])]['N']=Nr_f_jck_out['phi']
                        vars()['Nrjck'+str(key[0])]['jackknife']=Nr_f_jck_out['jackknife']
                        vars()['Nrjck'+str(key[0])]['sample']='reconstructed'
                        vars()['Nrjck'+str(key[0])]['biasing']='original'
                        vars()['Nrjck'+str(key[0])]['tomo']=key[0]
                        #N_compare_full_jck=vars()['Nrjck'+str(key[0])]
                        N_compare_full_jck=vstack([N_compare_full_jck,vars()['Nrjck'+str(key[0])]])
                print ('The reconstructed Nz jackknife have been prepared')


  '''
  #########################################################
  ################### create the biaised Nz ###############
  #########################################################
  '''


  for index, key in unique_pc_tags.iterrows():
    biaised_nz=[]
    #shifted Nz
    print (str(key[0]))
    if shift :
       if apply_to_compared_nz :
            print ('Creating shifted Nz to the compared Nz')
            Nt_SHp=Table()
            Nt_SHp['N']=vars()['Nt'+str(key[0])]['N']
            Nt_SHp['z']=vars()['Nt'+str(key[0])]['z']+shift_value

            Nt_SHm=Table()
            Nt_SHm['N']=vars()['Nt'+str(key[0])]['N']
            Nt_SHm['z']=vars()['Nt'+str(key[0])]['z']-shift_value

            Nt_SHp['sample']='compared'
            Nt_SHm['sample']='compared'

            Nt_SHm['biasing']='shift_minus'
            Nt_SHp['biasing']='shift_plus'

            Nt_SHm['tomo']=key[0]
            Nt_SHp['tomo']=key[0]

            N_compare_full=vstack([N_compare_full,Nt_SHm])
            N_compare_full=vstack([N_compare_full,Nt_SHp])
            biaised_nz.append('t_SHp_')
            biaised_nz.append('t_SHm_')

       if apply_to_rec_nz :
            print ('Creating shifted Nz to the recovered Nz')
            Nr_SHp=Table()
            Nr_SHp['N']=vars()['Nr'+str(key[0])]['N']
            Nr_SHp['z']=vars()['Nr'+str(key[0])]['z']+shift_value

            Nr_SHm = Table()
            Nr_SHm['N']=vars()['Nr'+str(key[0])]['N']
            Nr_SHm['z']=vars()['Nr'+str(key[0])]['z']-shift_value

            Nr_SHp['sample']='reconstructed'
            Nr_SHm['sample']='reconstructed'

            Nr_SHm['biasing']='shift_minus'
            Nr_SHp['biasing']='shift_plus'

            Nr_SHm['tomo']=key[0]
            Nr_SHp['tomo']=key[0]

            N_compare_full=vstack([N_compare_full,Nr_SHp])
            N_compare_full=vstack([N_compare_full,Nr_SHm])

            biaised_nz.append('r_SHp_')
            biaised_nz.append('r_SHm_')
            '''
            for jk_count in np.unique(Nr_f_jck['jackknife']):
                        Nr_SHp_jck_out=Nr_f_jck[Nr_f_jck['jackknife']==jk_count]
                        vars()['Nrjck'+str(key[0])]=Table()
                        vars()['Nrjck'+str(key[0])]['z']=Nr_f_jck_out['Z']
                        vars()['Nrjck'+str(key[0])]['N']=Nr_f_jck_out['phi']
                        vars()['Nrjck'+str(key[0])]['jackknife']=Nr_f_jck_out['jackknife']
                        vars()['Nrjck'+str(key[0])]['sample']='reconstructed'
                        vars()['Nrjck'+str(key[0])]['biasing']='original'
                        vars()['Nrjck'+str(key[0])]['tomo']=key[0]
                        #N_compare_full_jck=vars()['Nrjck'+str(key[0])]
                        N_compare_full_jck=vstack([N_compare_full_jck,vars()['Nrjck'+str(key[0])]])
            '''
    # spread Nz
    if spread:
        if apply_to_compared_nz :
                print ('Creating spreading Nz to the compared Nz')
                Nt_SP=Table()
                Nt_SP['N']=vars()['Nt'+str(key[0])]['N']
                Nt_SP['z']=mean_dist(vars()['Nt'+str(key[0])]['N'],vars()['Nt'+str(key[0])]['z'])-(mean_dist(vars()['Nt'+str(key[0])]['N'],vars()['Nt'+str(key[0])]['z'])-vars()['Nt'+str(key[0])]['z'])*spreading_factor

                Nt_SP['sample']='compared'
                Nt_SP['biasing']='spread'
                Nt_SP['tomo']=key[0]
                N_compare_full=vstack([N_compare_full,Nt_SP])

                biaised_nz.append('t_SP_')
        if apply_to_rec_nz==True:
                print ('Creating spreading Nz to the recovered Nz')
                Nr_SP=Table()
                Nr_SP['N']=vars()['Nr'+str(key[0])]['N']
                Nr_SP['z']=mean_dist(vars()['Nr'+str(key[0])]['N'],vars()['Nr'+str(key[0])]['z'])-(mean_dist(vars()['Nr'+str(key[0])]['N'],vars()['Nr'+str(key[0])]['z'])-vars()['Nr'+str(key[0])]['z'])*spreading_factor

                Nr_SP['sample']='reconstructed'
                Nr_SP['biasing']='spread'
                Nr_SP['tomo']=key[0]

                N_compare_full=vstack([N_compare_full,Nr_SP])

                biaised_nz.append('r_SP_')


    #Kurtosis
    if Kurtosis :
                #Kurtosis is a measure of whether the data are heavy-tailed
                # or light-tailed relative to a normal distribution.
                #That is, data sets with high kurtosis tend
                #to have heavy tails, or outliers.
                #Data sets with low kurtosis tend to have light tails,
                # or lack of outliers.
                #A uniform distribution would be the extreme case.

                print (' ')
    #Skewness
    if Skewness:
                #Skewness is a measure of symmetry,
                # or more precisely, the lack of symmetry.
                #A distribution, or data set,
                #is symmetric if it looks the same to the left and right
                #of the center point.
                print (' ')

    '''
    #################################################################
    ##################plotting Nz that will enter the Fom analysis###
    ################################################################

    plt.figure()
    #plt.plot(z,Nt,'--')
    plt.plot(vars()['Nr'+str(key[0])]['z'],vars()['Nr'+str(key[0])]['N'],'--')
    print (biaised_nz)
    for i in range(len(biaised_nz)):
         nz=Table.read(compare_noz_fom_path+'N'+biaised_nz[i]+str(key[0])+'.fits')
         nzz=nz['z']
         nzn=nz['N']
         plt.plot(nzz,nzn,'r--',label=biaised_nz[i])
 	     #plt.plot(nz['z'],nz['N'],'r--',label=biaised_nz[i])
    plt.savefig(compare_noz_fom_path+'input_fom_nz'+str(key[0])+'.png')
    '''
  N_compare_full.write(compare_noz_fom_path, path='results',format='hdf5',overwrite=True)
  if N_rec_jk or N_comp_jk :
    N_compare_full_jck.write(compare_noz_fom_path, path='jackknifes',overwrite=True,append=True,format='hdf5')
  return   biaised_nz

###############################################
################### functions #################
###############################################
#mean of the distribution
def mean_dist(N,z):
  mean=0.
  norm_mean=0.
  for k in range(len(N)):
        mean += N[k]*z[k]
        norm_mean += N[k]
  mean=mean/norm_mean
  return mean

#this is the module that compute the observables shear/wtheta/gammart




import numpy as np
import matplotlib.pylab as plt
from astropy.table import Table, vstack, hstack
import pandas as pa
import os, sys
import shutil
import subprocess
import math

#%matplotlib inline



def FoM(Nz_path,FoM_path,load_FoM,cosmosis_run_path,cosmosis_cmd,theta_min,theta_max,number_of_ang_bin,time0=0):
 if load_FoM :
        FoM_path=Table.read(FoM_path+'FoM.fits')
        #FoM=pa.HDFStore(FoM_path)
        return FoM
###############################################
################### Cosmosis  runs#################
###############################################

####################################################################
#**** running for all Nz distribution ********

 if not load_FoM :
    Nzs_full=pa.HDFStore(Nz_path)
    Nzs=Nzs_full['results']
    FoM_results=Table()
    FoM_results_jck=Table()
    gg=0
    kk_it=0
    errr=0
    print (cosmosis_run_path)
    for samp in np.unique(Nzs['sample']):
        for biasing in np.unique(Nzs['biasing']):
            for tomo in np.unique(Nzs['tomo']):
              print (samp, biasing, tomo)
              sel=((Nzs['sample']==samp) & (Nzs['biasing']==biasing) & (Nzs['tomo']==tomo))
              n_of_z=np.array([Nzs['z'][sel],Nzs['N'][sel]]).T

              if len(n_of_z)!=0 :
                      n_of_z=n_of_z[np.argsort(n_of_z[:,0])]
                      ##############interpolation 0 points to avoid cosmosis negative problem#############
                      z=n_of_z[:,0]
                      n_of_z=n_of_z[:,1]
                      if all(v == 0 for v in n_of_z) :
                          'empty bin'
                          n_of_z[0]=1e-3
                      nule=20  #linear interpolation last point
                      n_of_z_new=np.zeros((1,len(z)+nule))
                      z_new=np.zeros(len(z)+nule)
                      z_new[nule:]=z
                      bin_width=(z[1]-z[0])
                      z_new[0]=z[0]-bin_width
                      for i_bin in range(nule):
                          z_new[i_bin]=(bin_width/nule)*i_bin+z_new[0]
                      n_of_z_new[:,nule:]=n_of_z
                      n_of_z_new[:,:nule]=n_of_z_new[:,nule]*(z_new[:nule]-z_new[0])/(z_new[nule]-z_new[0])
                      z,n_of_z=z_new,n_of_z_new
                      n_of_z=n_of_z[:,z>=0]
                      z=z[z>=0]
                      n_of_z=np.vstack([z,n_of_z]).T
                      ##################################################################################################
                      np.savetxt('n_of_z.txt',n_of_z)
                      if not os.path.exists(FoM_path+'output_'+samp+'_'+biasing+'_'+str(tomo)):
                          command=cosmosis_cmd+" Obs_est.ini"
                          path_to_cosmosis='.'
                          cmd = ['bash', '-c', 'source '+ cosmosis_run_path +' && {1}'.format(path_to_cosmosis, command)]
                          process = subprocess.Popen(cmd)
                          process.wait()

                          if gg >0 :
                            FoM_results_1=Table()
                            FoM_results_1['theta']=np.loadtxt('out/shear_xi/theta.txt')
                            FoM_results_1['xip']=np.loadtxt('out/shear_xi/xiplus_1_1.txt')
                            FoM_results_1['xim']=np.loadtxt('out/shear_xi/ximinus_1_1.txt')
                            FoM_results_1['ggl']=np.loadtxt('out/ggl_xi/tanshear_1_1.txt')
                            FoM_results_1['wmatter']=np.loadtxt('out/matter_xi/wmatter_1_1.txt')
                            FoM_results_1['sample']=samp
                            FoM_results_1['biasing']=biasing
                            FoM_results_1['tomo']=tomo

                            # rebining
                            select=np.array(FoM_results_1['theta']*(180/np.pi)*60)>theta_min
                            select1=np.array(FoM_results_1['theta']*(180/np.pi)*60)<theta_max*1.3
                            select=select & select1
                            FoM_results_1=FoM_results_1[select]
                            reduction=int((len(FoM_results_1))/number_of_ang_bin)* np.array(range(int(number_of_ang_bin)))
                            FoM_results_1=FoM_results_1[reduction]


                            FoM_results=vstack([FoM_results,FoM_results_1])

                          else :
                            FoM_results['theta']=np.loadtxt('out/shear_xi/theta.txt')
                            FoM_results['xip']=np.loadtxt('out/shear_xi/xiplus_1_1.txt')
                            FoM_results['xim']=np.loadtxt('out/shear_xi/ximinus_1_1.txt')
                            FoM_results['ggl']=np.loadtxt('out/ggl_xi/tanshear_1_1.txt')
                            FoM_results['wmatter']=np.loadtxt('out/matter_xi/wmatter_1_1.txt')
                            FoM_results['sample']=samp
                            FoM_results['biasing']=biasing
                            FoM_results['tomo']=tomo

                            # rebining
                            select=np.array(FoM_results['theta']*(180/np.pi)*60)>theta_min
                            select1=np.array(FoM_results['theta']*(180/np.pi)*60)<theta_max*1.3
                            select=select & select1
                            FoM_results=FoM_results[select]
                            print int((len(FoM_results))/number_of_ang_bin) * np.array(range(int(number_of_ang_bin)))
                            reduction=int((len(FoM_results))/number_of_ang_bin) * np.array(range(int(number_of_ang_bin)))
                            #reduction=int((len(FoM_results))/number_of_ang_bin)
                            print (reduction)
                            FoM_results=FoM_results[reduction]


                            gg=gg+1

              else:
                  continue
    print ('cosmosis run for full area')
    print FoM_path
    FoM_results.write(FoM_path+'FoM.fits',overwrite=True)
    #FoM_results.write(FoM_path, path='results',format='hdf5',overwrite=True)
    '''
    ####################################################################
             #***** running on jackknife regions *******
    ####################################################################
    '''
    FoM_results=Table.read(FoM_path+'FoM.fits')
    print ('computing FoM in every jk regions')

    Nzs_jck=Nzs_full['jackknifes']
    print (Nzs_jck)
    number_of_jk=len(np.unique(Nzs_jck['jackknife']))
    print number_of_jk
    for samp in np.unique(Nzs_jck['sample']):

      for biasing in np.unique(Nzs_jck['biasing']):

         for tomo in np.unique(Nzs_jck['tomo']):
            print (np.unique(Nzs_jck['jackknife']))
            itt=0
            for jk_count in np.unique(Nzs_jck['jackknife']):
              itt=itt+1
              print (itt)
              print ('JK LOOP',samp, biasing , tomo, jk_count)

              sel=((Nzs_jck['sample']==samp) & (Nzs_jck['biasing']==biasing) & (Nzs_jck['tomo']==tomo) & (Nzs_jck['jackknife']==jk_count))
              n_of_z=np.array([Nzs_jck['z'][sel],Nzs_jck['N'][sel]]).T

              print n_of_z
              if len(n_of_z)!=0 :
                      n_of_z=n_of_z[np.argsort(n_of_z[:,0])]
                      print ('NON zero lenght',samp,biasing,tomo,jk_count)
                      ##############interpolation 0 points to avoid cosmosis negative problem#############
                      z=n_of_z[:,0]
                      n_of_z=n_of_z[:,1]
                      if all(v == 0 for v in n_of_z) :
                          'empty bin'
                          n_of_z[0]=1e-3
                      nule=20  #linear interpolation last point
                      n_of_z_new=np.zeros((1,len(z)+nule))
                      z_new=np.zeros(len(z)+nule)
                      z_new[nule:]=z
                      bin_width=(z[1]-z[0])
                      z_new[0]=z[0]-bin_width
                      print (bin_width)
                      for i_bin in range(nule):
                          z_new[i_bin]=(bin_width/nule)*i_bin+z_new[0]
                      n_of_z_new[:,nule:]=n_of_z
                      n_of_z_new[:,:nule]=n_of_z_new[:,nule]*(z_new[:nule]-z_new[0])/(z_new[nule]-z_new[0])
                      z,n_of_z=z_new,n_of_z_new
                      n_of_z=n_of_z[:,z>=0]
                      z=z[z>=0]
                      n_of_z=np.vstack([z,n_of_z]).T
                      ##################################################################################################
                      np.savetxt('n_of_z.txt',n_of_z)
                      #if not os.path.exists(FoM_path+'output_'+samp+'_'+biasing+'_'+str(tomo)+str(jk_count)):
                      command=cosmosis_cmd+" Obs_est.ini"
                      path_to_cosmosis='.'
                      cmd = ['bash', '-c', 'source '+ cosmosis_run_path +' && {1}'.format(path_to_cosmosis, command)]
                      process = subprocess.Popen(cmd)
                      process.wait()
                      if kk_it>0:
                            FoM_results_1_jck=Table()
                            FoM_results_1_jck['theta']=np.loadtxt('out/shear_xi/theta.txt')
                            FoM_results_1_jck['xip']=np.loadtxt('out/shear_xi/xiplus_1_1.txt')
                            FoM_results_1_jck['xim']=np.loadtxt('out/shear_xi/ximinus_1_1.txt')
                            FoM_results_1_jck['ggl']=np.loadtxt('out/ggl_xi/tanshear_1_1.txt')
                            FoM_results_1_jck['wmatter']=np.loadtxt('out/matter_xi/wmatter_1_1.txt')
                            FoM_results_1_jck['sample']=samp
                            FoM_results_1_jck['biasing']=biasing
                            FoM_results_1_jck['tomo']=tomo
                            FoM_results_1_jck['jackknife']=jk_count
                            # rebining
                            select=np.array(FoM_results_1_jck['theta']*(180/np.pi)*60)>theta_min
                            select1=np.array(FoM_results_1_jck['theta']*(180/np.pi)*60)<theta_max*1.3
                            select=select & select1
                            FoM_results_1_jck=FoM_results_1_jck[select]
                            reduction=int((len(FoM_results_1_jck))/number_of_ang_bin)* np.array(range(int(number_of_ang_bin)))
                            FoM_results_1_jck=FoM_results_1_jck[reduction]

                            FoM_results_jck=vstack([FoM_results_jck,FoM_results_1_jck])
                      else :
                            FoM_results_jck['theta']=np.loadtxt('out/shear_xi/theta.txt')
                            FoM_results_jck['xip']=np.loadtxt('out/shear_xi/xiplus_1_1.txt')
                            FoM_results_jck['xim']=np.loadtxt('out/shear_xi/ximinus_1_1.txt')
                            FoM_results_jck['ggl']=np.loadtxt('out/ggl_xi/tanshear_1_1.txt')
                            FoM_results_jck['wmatter']=np.loadtxt('out/matter_xi/wmatter_1_1.txt')
                            FoM_results_jck['sample']=samp
                            FoM_results_jck['biasing']=biasing
                            FoM_results_jck['tomo']=tomo
                            FoM_results_jck['jackknife']=jk_count
                            # rebining
                            select=np.array(FoM_results_jck['theta']*(180/np.pi)*60)>theta_min
                            select1=np.array(FoM_results_jck['theta']*(180/np.pi)*60)<theta_max*1.3
                            select=select & select1
                            FoM_results_jck=FoM_results_jck[select]
                            reduction=int((len(FoM_results_jck))/number_of_ang_bin)* np.array(range(int(number_of_ang_bin)))
                            FoM_results_jck=FoM_results_jck[reduction]
                            kk_it=kk_it+1

              else :
                  print ('zero-size array')
            '''
            ####################################################################
            *************** compute errors and covariances *********************
            ####################################################################
            '''
            print 'computing error and covariances'
            xip_array=[]
            xim_array=[]
            #FoM_results_jck.write('samp.fits',overwrite=True)
            Full_jk_tomo=FoM_results_jck[(FoM_results_jck['sample']==samp) & (FoM_results_jck['biasing']==biasing) & (FoM_results_jck['tomo']==tomo)]
            print len(Full_jk_tomo)
            print (Full_jk_tomo)
            jkk=0
            if math.isnan(Full_jk_tomo['xip'][0])  :
                print ('empty bin no error computation')

            else :
             for j_count in np.unique(Full_jk_tomo['jackknife']):
                print np.unique(Full_jk_tomo['jackknife'])
                print (samp,biasing,tomo)
                if jkk>0 :
                    xip_array=np.column_stack((xip_array,np.array(Full_jk_tomo[Full_jk_tomo['jackknife']==j_count]['xip'])))
                    xim_array=np.column_stack((xim_array,np.array(Full_jk_tomo[Full_jk_tomo['jackknife']==j_count]['xim'])))
                    ggl_array=np.column_stack((ggl_array,np.array(Full_jk_tomo[Full_jk_tomo['jackknife']==j_count]['ggl'])))
                    wmatter_array=np.column_stack((wmatter_array,np.array(Full_jk_tomo[Full_jk_tomo['jackknife']==j_count]['wmatter'])))
                else:
                   xip_array=np.array(Full_jk_tomo[Full_jk_tomo['jackknife']==j_count]['xip'])
                   xim_array=np.array(Full_jk_tomo[Full_jk_tomo['jackknife']==j_count]['xim'])
                   ggl_array=np.array(Full_jk_tomo[Full_jk_tomo['jackknife']==j_count]['ggl'])
                   wmatter_array=np.array(Full_jk_tomo[Full_jk_tomo['jackknife']==j_count]['wmatter'])
                   jkk=jkk+1


             #xip
             print xip_array
             print xip_array.shape
             print len(xip_array)
             dict_xip=covariance_jck(xip_array.T,len(xip_array.T))
             cov_jck_xip=dict_xip['cov']
             print (dict_xip['err'].shape)
             #xim
             dict_xim=covariance_jck(xim_array.T,len(xip_array.T))
             cov_jck_xim=dict_xim['cov']
             #ggl
             dict_ggl=covariance_jck(ggl_array.T,len(xip_array.T))
             cov_jck_ggl=dict_ggl['cov']
             #wmatter
             dict_wmatter=covariance_jck(wmatter_array.T,len(xip_array.T))
             cov_jck_wmatter=dict_wmatter['cov']

             if errr>0 :
                FoM_results_err2=FoM_results[(FoM_results['sample']==samp) & (FoM_results['biasing']==biasing) & (FoM_results['tomo']==tomo)]
                FoM_results_err2['err_xip']=dict_xip['err']
                FoM_results_err2['err_xim']=dict_xim['err']
                FoM_results_err2['err_ggl']=dict_ggl['err']
                FoM_results_err2['err_wmatter']=dict_wmatter['err']

                FoM_results_err=vstack([FoM_results_err,FoM_results_err2])
             else :
               FoM_results_err=FoM_results[(FoM_results['sample']==samp) & (FoM_results['biasing']==biasing) & (FoM_results['tomo']==tomo)]
               print FoM_results_err
               FoM_results_err['err_xip']=dict_xip['err']
               FoM_results_err['err_xim']=dict_xim['err']
               FoM_results_err['err_ggl']=dict_ggl['err']
               FoM_results_err['err_wmatter']=dict_wmatter['err']

               errr=errr+1


    FoM_results_jck.write(FoM_path+'FoM_jck.fits',overwrite=True)
    FoM_results_err.write(FoM_path+'FoM_err.fits',overwrite=True)
    #FoM_results_jck.write(FoM_path, path='jackknifes',format='hdf5',overwrite=True)
#    FoM_results_err.write(FoM_path, path='errors',format='hdf5',overwrite=True)


    return FoM_results
###############################################
################### functions #################
###############################################
def run_cosmosis_command(command):
    cmd = [command]#['bash', '-c', 'source'+ cosmosis_run_path +' && {1}'.format(path_to_cosmosis, command)]
    process = subprocess.Popen(cmd)
    process.wait()
    return process


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

def covariance_jck(TOTAL_PHI,jk_r):
  #  Covariance estimation
  average=np.zeros(TOTAL_PHI.shape[1])
  cov_jck=np.zeros((TOTAL_PHI.shape[1],TOTAL_PHI.shape[1]))
  err_jck=np.zeros(TOTAL_PHI.shape[1])
  for kk in range(jk_r):
    average+=TOTAL_PHI[kk,:]
  average=average/(jk_r)
  for ii in range(TOTAL_PHI.shape[1]):
     for jj in range(ii+1):
          for kk in range(jk_r):
            cov_jck[ii,jj]+=TOTAL_PHI[kk,ii]*TOTAL_PHI[kk,jj]
          cov_jck[ii,jj]=(-average[ii]*average[jj]*jk_r+cov_jck[ii,jj])*(jk_r-1)/(jk_r)
          cov_jck[jj,ii]=cov_jck[ii,jj]
  for ii in range(TOTAL_PHI.shape[1]):
   err_jck[ii]=np.sqrt(cov_jck[ii,ii])
  #compute correlation
  corr=np.zeros((TOTAL_PHI.shape[1],TOTAL_PHI.shape[1]))
  for i in range(TOTAL_PHI.shape[1]):
      for j in range(TOTAL_PHI.shape[1]):
        corr[i,j]=cov_jck[i,j]/(np.sqrt(cov_jck[i,i]*cov_jck[j,j]))
  average=average
  return {'cov' : cov_jck,
          'err' : err_jck,
          'corr':corr,
          'mean':average}

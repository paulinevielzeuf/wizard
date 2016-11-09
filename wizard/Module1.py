# this module aims to prepare the Nz for which one wants to process FoM analysis including the possibility to create biased Nz (shifted , spread,..)




import numpy as np
import matplotlib.pylab as plt
from astropy.table import Table, vstack, hstack
from  config_module1 import *
#%matplotlib inline


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







#######################################
########loading the files #############
#######################################


Nz_to_compare_with='input/N_true'
Nz_rec='input/N_z'

Nt1=np.loadtxt(Nz_to_compare_with+'.txt',usecols=(0,1))
Nt=Table()
Nt['z']=Nt1[:,0]
Nt['N']=Nt1[:,1]
if N_comp_jk=='true' :
	Nt_jk=np.loadtxt(Nz_to_compare_with+'_jkc.txt',usecols=range(1,njck_comp))



Nr1=np.loadtxt(Nz_rec+'.txt',usecols=(0,1))
Nr=Table()
Nr['z']=Nr1[:,0]
Nr['N']=Nr1[:,1]
if N_rec_jk=='true' :
	        Nr_jk=np.loadtxt(Nz_rec+'_jck.txt',usecols=range(1,njck_rec))






#########################################################
################### create the biaised Nz ###############
#########################################################

biaised_nz=[]

#shifted Nz
if shift=='true':
 if apply_to_compared_nz=='true':
  print 'Creating shifted Nz to the compared Nz'
  Nt_SHp=Table()
  Nt_SHp['N']=Nt['N']
  Nt_SHp['z']=Nt['z']+shift_value

  Nt_SHm=Table()
  Nt_SHm['N']=Nt['N']
  Nt_SHm['z']=Nt['z']-shift_value
  Nt_SHp.write('input/Nt_SHp.fits',overwrite=True)
  Nt_SHm.write('input/Nt_SHm.fits',overwrite=True)
  biaised_nz.append('t_SHp')
  biaised_nz.append('t_SHm')
 if apply_to_rec_nz=='true':
   print 'Creating shifted Nz to the recovered Nz'
   Nr_SHp=Table()
   Nr_SHp['N']=Nr['N']
   Nr_SHp['z']=Nr['z']+shift_value

   Nr_SHm = Table()
   Nr_SHm['N']=Nr['N']
   Nr_SHm['z']=Nr['z']-shift_value


   Nr_SHp.write('input/Nr_SHp.fits',overwrite=True)
   Nr_SHm.write('input/Nr_SHm.fits',overwrite=True)
   biaised_nz.append('r_SHp')
   biaised_nz.append('r_SHm')

# spread Nz
if spread=='true':
 if apply_to_compared_nz=='true':
  print 'Creating spreading Nz to the compared Nz'

  Nt_SP=Table()
  Nt_SP['N']=Nt['N']
  Nt_SP['z']=mean_dist(Nt['N'],Nt['z'])-(mean_dist(Nt['N'],Nt['z'])-Nt['z'])*spreading_factor
  biaised_nz.append('t_SP')
  Nt_SP.write('input/Nt_SP.fits',overwrite=True)
 if apply_to_rec_nz=='true':
  print 'Creating spreading Nz to the recovered Nz'

  Nr_SP=Table()
  Nr_SP['N']=Nr['N']
  Nr_SP['z']=mean_dist(Nr['N'],Nr['z'])-(mean_dist(Nr['N'],Nr['z'])-Nr['z'])*spreading_factor
  Nr_SP.write('input/Nr_SP.fits',overwrite=True)
  biaised_nz.append('r_SP')


#Kurtosis
if Kurtosis=='true' :

 print ' '
# Skewness
if Skewness=='true':
 print ' '


#################################################################
##################plotting Nz that will enter the Fom analysis###
################################################################
plt.figure()
#plt.plot(z,Nt,'--')
plt.plot(Nr['z'],Nr['N'],'--')
for i in range(len(biaised_nz)):
	plt.plot(vars()['N'+biaised_nz[i]]['z'],vars()['N'+biaised_nz[i]]['N'],'r--')
plt.show()

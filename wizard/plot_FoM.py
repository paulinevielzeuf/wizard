import numpy as np
import matplotlib.pylab as plt
from astropy.table import Table, vstack, hstack
import pandas as pa
import os, sys
import shutil
import subprocess


def plot_FoM(fom_path,comp,bins,wiz_file,xip,xim,ggl,wmatter,plot_path,compared_plot,time0=0):
      bi=interpret_bins(bins)
      if wiz_file :
          FoM=Table.read(fom_path+'FoM.fits')
          FoM_err=Table.read(fom_path+'FoM_err.fits')
         # print FoM
          tomo_lab_ind=0
          for tomo in np.unique(FoM['tomo']):

               rec=FoM_err[(FoM_err['tomo']==tomo) & (FoM_err['sample']=='reconstructed')]
               compared=FoM_err[(FoM_err['tomo']==tomo) & (FoM_err['sample']=='compared')]


               if xip :
                   # make figure
                   fig, ax = make_figure(r'$\xi_+(\theta)$')
                   xi_rec_kwargs = {'x': rec['theta'], 'y': rec['xip'], 'yerr': rec['err_xip'],
                                        'marker': 'D',
                                        'linestyle': '-',
                                        'color': 'teal', 'alpha': 0.8, 'label' : 'crossX'}
                   ax.errorbar(**xi_rec_kwargs)

                   if len(compared)> 0:
                    xi_comp_kwargs = {'x': compared['theta'], 'y': compared['xip'], 'yerr': compared['err_xip'],
                                        'marker': 'D',
                                        'linestyle': '-',
                                        'color': 'purple', 'alpha': 0.8, 'label' : comp}
                    ax.errorbar(**xi_comp_kwargs)


                   for ith, sampbias in enumerate(compared_plot):
                                biasplot=FoM[(FoM['tomo']==tomo) & (FoM['sample']==sampbias[0]) &(FoM['biasing']==sampbias[1])]
                                if len(biasplot)>0 :
                                    if sampbias[0]=='reconstructed':
                                        col='teal'
                                    if sampbias[0]=='compared':
                                        col='purple'
                                    if sampbias[1]=='shift_minus' or sampbias[1]=='shift_plus' :
                                        lnst='--'
                                    if sampbias[1]=='spread':
                                        lnst='-.'
                                    ax.plot (biasplot['theta'],  biasplot['xip'],linestyle= lnst,color= col,alpha= 0.8, label= str(sampbias))




                   ax.legend(loc='lower right', fancybox=True).get_frame().set_alpha(0.5)
                   #ax.plot(**biasing_kwargs)
                   ax.set_xscale("log")

                   ax.set_title(str(bi[tomo_lab_ind])+r'$ < z < $'+str(bi[tomo_lab_ind+1]))
                   fig.savefig(plot_path+'xip_'+str(tomo))


               if xim :
                   # make figure
                   fig, ax = make_figure(r'$\xi_-(\theta)$')
                   xi_rec_kwargs = {'x': rec['theta'], 'y': rec['xim'], 'yerr': rec['err_xim'],
                                        'marker': 'D',
                                        'linestyle': '-',
                                        'color': 'teal', 'alpha': 0.8, 'label' : 'crossX'}
                   ax.errorbar(**xi_rec_kwargs)

                   if len(compared)> 0:
                     xi_comp_kwargs = {'x': compared['theta'], 'y': compared['xim'], 'yerr': compared['err_xim'],
                                        'marker': 'D',
                                        'linestyle': '-',
                                        'color': 'purple', 'alpha': 0.8, 'label' : comp}
                     ax.errorbar(**xi_comp_kwargs)

                   for ith, sampbias in enumerate(compared_plot):
                                biasplot=FoM[(FoM['tomo']==tomo) & (FoM['sample']==sampbias[0]) &(FoM['biasing']==sampbias[1])]
                                if len(biasplot)>0 :
                                    if sampbias[0]=='reconstructed':
                                        col='teal'
                                    if sampbias[0]=='compared':
                                        col='purple'
                                    if sampbias[1]=='shift_minus' or sampbias[1]=='shift_plus' :
                                        lnst='--'
                                    if sampbias[1]=='spread':
                                        lnst='-.'
                                    ax.plot (biasplot['theta'],  biasplot['xim'],linestyle= lnst,color= col,alpha= 0.8, label= sampbias)




                   ax.legend(loc='lower right', fancybox=True).get_frame().set_alpha(0.5)
                   ax.set_xscale("log")
                   ax.set_title(str(bi[tomo_lab_ind])+r'$ < z < $'+str(bi[tomo_lab_ind+1]))

                   plt.savefig(plot_path+'xim_'+str(tomo))
                   plt.close()


               if ggl :
                   # make figure
                   fig, ax = make_figure(r'$\gamma_t(\theta)$')
                   xi_rec_kwargs = {'x': rec['theta'], 'y': rec['ggl'], 'yerr': rec['err_ggl'],
                                        'marker': 'D',
                                        'linestyle': '-',
                                        'color': 'teal', 'alpha': 0.8, 'label' : 'crossX'}
                   ax.errorbar(**xi_rec_kwargs)

                   if len(compared)> 0:
                     xi_comp_kwargs = {'x': compared['theta'], 'y': compared['ggl'], 'yerr': compared['err_ggl'],
                                        'marker': 'D',
                                        'linestyle': '-',
                                        'color': 'purple', 'alpha': 0.8, 'label' : comp}
                     ax.errorbar(**xi_comp_kwargs)


                   for ith, sampbias in enumerate(compared_plot):
                                biasplot=FoM[(FoM['tomo']==tomo) & (FoM['sample']==sampbias[0]) &(FoM['biasing']==sampbias[1])]
                                if len(biasplot)>0 :
                                    if sampbias[0]=='reconstructed':
                                        col='teal'
                                    if sampbias[0]=='compared':
                                        col='purple'
                                    if sampbias[1]=='shift_minus' or sampbias[1]=='shift_plus' :
                                        lnst='--'
                                    if sampbias[1]=='spread':
                                        lnst='-.'
                                    ax.plot (biasplot['theta'],  biasplot['ggl'],linestyle= lnst,color= col,alpha= 0.8, label= sampbias)



                   ax.legend(loc='lower right', fancybox=True).get_frame().set_alpha(0.5)
                   ax.set_xscale("log")
                   ax.set_title(str(bi[tomo_lab_ind])+r'$ < z < $'+str(bi[tomo_lab_ind+1]))
                   plt.savefig(plot_path+'ggl_'+str(tomo))
                   plt.close()


               if wmatter :
                   # make figure
                   fig, ax = make_figure(r'$w(\theta)$')
                   xi_rec_kwargs = {'x': rec['theta'], 'y': rec['wmatter'], 'yerr': rec['err_wmatter'],
                                        'marker': 'D',
                                        'linestyle': '-',
                                        'color': 'teal', 'alpha': 0.8, 'label' : 'crossX'}

                   ax.errorbar(**xi_rec_kwargs)
                   if len(compared)> 0:

                     xi_comp_kwargs = {'x': compared['theta'], 'y': compared['ggl'], 'yerr': compared['err_wmatter'],
                                        'marker': 'D',
                                        'linestyle': '-',
                                        'color': 'purple', 'alpha': 0.8, 'label' : comp}
                     ax.errorbar(**xi_comp_kwargs)

                   for ith, sampbias in enumerate(compared_plot):
                                biasplot=FoM[(FoM['tomo']==tomo) & (FoM['sample']==sampbias[0]) &(FoM['biasing']==sampbias[1])]
                                if len(biasplot)>0 :
                                    if sampbias[0]=='reconstructed':
                                        col='teal'
                                    if sampbias[0]=='compared':
                                        col='purple'
                                    if sampbias[1]=='shift_minus' or sampbias[1]=='shift_plus' :
                                        lnst='--'
                                    if sampbias[1]=='spread':
                                        lnst='-.'
                                    ax.plot (biasplot['theta'],  biasplot['wmatter'],linestyle= lnst,color= col,alpha= 0.8, label= sampbias)




                   ax.legend(loc='lower right', fancybox=True).get_frame().set_alpha(0.5)
                   ax.set_xscale("log")
                   ax.set_title(str(bi[tomo_lab_ind])+r'$ < z < $'+str(bi[tomo_lab_ind+1]))
                   plt.savefig(plot_path+'wmatter_'+str(tomo))
                   plt.close()
               tomo_lab_ind=tomo_lab_ind+1

def make_figure(observable):
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(13, 8))
    ax.set_xlabel(r'$\theta$', fontsize=28)
    ax.set_ylabel(observable, fontsize=28, rotation=0)
    ax.yaxis.labelpad = 25
    return fig, ax


def interpret_bins(bins):
    results = {}
    for key in bins:
        if bins[key][0] == 'between':
            # if 3 entries, then turn into arange
            results[key] = np.arange(*bins[key][1])
            print results[key]
        elif bins[key][0] == 'equal':
            # else leave in
            results[key] = ['equal', bins[key][1]]

    return results[key]

import os
import numpy              as np
import matplotlib.pyplot  as plt
from   matplotlib         import rc
from   utils              import *

fname    = 'fig5.eps'
fontsize = 24
rc('font',**{'family':'serif','serif':['Computer Modern']})
rc('text',usetex=True)
rc('axes',labelsize=fontsize,titlesize=fontsize)
rc('xtick',labelsize=fontsize-3)
rc('ytick',labelsize=fontsize-3)

xfig = 9.0
yfig = 0.5*xfig

fig,ax = plt.subplots(2,2,figsize=(2*xfig,2*yfig))
fig.subplots_adjust(wspace=0.2,hspace=0.25)

# ==============================================================

import sys
sys.path.append('../code/')
from fig5a import make_figure_5a
from fig5b import make_figure_5b
from fig5c import make_figure_5c
from fig5d import make_figure_5d

color_dict = {'nmnm':c_list['black'],'romnm':c_list['mid-red'],'romrom':c_list['mid-green'],'qstrom':c_list['mid-blue']}
color_dict = {'nmnm':'o',            'romnm':'D',              'romrom':'s',                'qstrom':'*'}
make_figure_5a(ax[0,0])
make_figure_5b(ax[0,1])
make_figure_5c(ax[1,0])
make_figure_5d(ax[1,1])

for (x,y) in [(0,0),(0,1),(1,0)]:
    if((x,y)==(1,0)):
       fill_panel(ax[x,y],'',[0,2],[0,1,2],['PSPCz','2F-PSPCz','4F-PSPCz'],
                          r'$\Delta \Delta E_{ST} [eV]$',[0,28],[0,7,14,21,28],
                          ['0','7','14','21','28'],p=8.0,q=20.0)
    else:
       fill_panel(ax[x,y],'',[0,2],[0,1,2],['PSPCz','2F-PSPCz','4F-PSPCz'],
                          r'$\Delta \Delta E_{ST} [eV]$',[-7,21],[-7,0,7,14,21],
                          ['-7','0','7','14','21'],p=8.0,q=20.0)
    x0L,y0L = 0.10,0.60
    dxL,dyL = (1-2*x0L),0.20
    h,l = ax[x,y].get_legend_handles_labels()
    ax[x,y].axhline(0,color='black',linestyle=':')

fill_panel(ax[1,1],'calculated $\Delta E_{ST}$ [eV]',[0.0,0.9],[0.0,0.3,0.6,0.9],['0.0','0.3','0.6','0.9'],
                   r'experimental $\Delta E_{ST}$ [eV]',[-0.4,0.4],[-0.4,-0.2,0.0,0.2,0.4],['-0.4','-0.2','0.0','0.2','0.4'],
                   p=20.0,q=10.0)

x0L,y0L = 0.0,1.12
dxL,dyL = 2.0,0.3
ha,la = ax[0,1].get_legend_handles_labels()
hb,lb = ax[1,1].get_legend_handles_labels()
h,l = ha+hb,la+lb

ax[0,0].legend(h,l,fancybox=True,shadow=True,ncol=4,loc=3,
               bbox_to_anchor=(x0L,y0L,dxL,dyL),handlelength=2,fontsize=16)

fig.savefig(fname,format='eps')
os.system('ps2epsi '+fname)
os.system('mv '+fname+'i '+fname)


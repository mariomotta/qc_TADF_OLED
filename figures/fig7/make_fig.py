import os
import numpy              as np
import matplotlib.pyplot  as plt
from   matplotlib         import rc
from   utils              import *

fname    = 'fig7.eps'
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
from fig7a import make_figure_7a
from fig7b import make_figure_7b
from fig7c import make_figure_7c
from fig7d import make_figure_7d

color_dict = {'nmnm':c_list['black'],'romnm':c_list['mid-red'],'romrom':c_list['mid-green'],'qstrom':c_list['mid-blue']}
color_dict = {'nmnm':'o',            'romnm':'D',              'romrom':'s',                'qstrom':'*'}
make_figure_7a(ax[0,0])
make_figure_7b(ax[0,1])
make_figure_7c(ax[1,0])
make_figure_7d(ax[1,1])

for (x,y) in [(0,0),(0,1),(1,0)]:
    fill_panel(ax[x,y],'',[0,1],[0,1],['qst overlap + nm H','qst overlap + qst H'],
                       r'$\Delta \Delta E_{ST} [\mathrm{mHa}]$',[-6,24],[-6,0,6,12,18,24],
                       ['-6','0','6','12','18','24'],p=4.0,q=20.0)
    x0L,y0L = 0.80,0.70
    dxL,dyL = (1-2*x0L),0.20
    #h,l = ax[x,y].get_legend_handles_labels()
    #ax[x,y].legend(h,l,fancybox=True,shadow=True,ncol=1,loc=3,
    #               bbox_to_anchor=(x0L,y0L,dxL,dyL),handlelength=1.2,fontsize=fontsize-12) 
    ax[x,y].axhline(0,color='black',linestyle=':')

fill_panel(ax[1,1],'calculated $\Delta E_{ST}$ [eV]',[0.2,0.8],[0.2,0.4,0.6,0.8],['0.2','0.4','0.6','0.8'],
                   r'experimental $\Delta E_{ST}$ [eV]',[-0.4,0.4],[-0.4,-0.2,0.0,0.2,0.4],['-0.4','-0.2','0.0','0.2','0.4'],
                   p=8.0,q=10.0)

x0L,y0L = 0.02,1.12
dxL,dyL = 3,0.30
ha,la = ax[1,1].get_legend_handles_labels()
hb,lb = ax[0,0].get_legend_handles_labels()
h,l   = hb+ha,lb+la
ax[0,0].legend(h,l,fancybox=True,shadow=True,ncol=4,loc=3,
               bbox_to_anchor=(x0L,y0L,dxL,dyL),handlelength=2.5,fontsize=20,columnspacing=8)

# ==============================================================

fig.savefig(fname,format='eps')
os.system('ps2epsi '+fname)
os.system('mv '+fname+'i '+fname)


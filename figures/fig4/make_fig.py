import os
import numpy              as np
import matplotlib.pyplot  as plt
from   matplotlib         import rc
from   utils              import *

fname    = 'fig4.eps'
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

color_dict  = {'S_0':c_list['mid-green'],'T_1':c_list['mid-red'],'S_1':c_list['mid-blue']}
marker_dict = {'STO-3G':'D','6-31G(d)':'o'}
v_dict      = {}

for figure_label,(x,y),figure_title in zip(['a','b','c'],[(0,0),(0,1),(1,0)],['PSPCz','2F-PSPCz','4F-PSPCz']):
    fill_panel(ax[x,y],'',[0,2],[0,1,2],['exact','q-EOM','VQD'],
                       'Energy [Ha]',[-800,-400],[-800,-700,-600,-500,-400],['-800','-700','-600','-500','-400'],p=20.0,q=10.0)
    ax[x,y].set_title('(%s) %s' % (figure_label,figure_title))
    v = np.loadtxt('fig4%s.txt' % figure_label)
    for i,m in enumerate(['exact','q-EOM','VQD']):
        for k,b in enumerate(['STO-3G','6-31G(d)']):
            for j,s in enumerate(['S_0','T_1','S_1']):
                print(v.shape,i,j,k,j+3*k)
                v_dict[m+' '+figure_title+' '+s+' '+b] = v[i,j+3*k]
    
    for i,m in enumerate(['exact','q-EOM','VQD']):
        for k,b in enumerate(['STO-3G','6-31G(d)']):
            for j,s in enumerate(['S_0','T_1','S_1']):
                ax[x,y].scatter([i],[v_dict[m+' '+figure_title+' '+s+' '+b]],s=160,color=color_dict[s],
                                marker=marker_dict[b],linewidths=0.5,edgecolors='black',label=r'%s, %s' % (m,b))   

# ==============================================================

color_dict  = {'exact':c_list['mid-green'],'q-EOM':c_list['mid-red'],'VQD':c_list['mid-blue']}
style_dict  = {'exact':'--','q-EOM':'-.','VQD':':'}

fill_panel(ax[1,1],r'Calculated $\Delta E_{ST}$ [eV]',[0,2.5],[0.0,0.5,1.0,1.5,2.0,2.5],['0.0','0.5','1.0','1.5','2.0','2.5'],
                   r'Experimental $\Delta E_{ST}$ [eV]',[-0.4,0.4],[-0.4,-0.2,0.0,0.2,0.4],
                   ['-0.4','-0.2','0.0','0.2','0.4'],p=20.0,q=10.0)

v = np.loadtxt('fig4d.txt')
for a,f in enumerate(['PSPCz','2F-PSPCz','4F-PSPCz']):
    v_dict['GAP experiment '+f] = v[a,0]
    for i,m in enumerate(['exact','q-EOM','VQD']):
        for k,b in enumerate(['STO-3G','6-31G(d)']):
            v_dict['GAP '+m+' '+f+' '+b] = v[a,1+k*3+j]

for i,m in enumerate(['exact','q-EOM','VQD']):
    for k,b in enumerate(['STO-3G','6-31G(d)']):
        x,y = [],[]
        for a,f in enumerate(['PSPCz','2F-PSPCz','4F-PSPCz']):
            x.append(v_dict['GAP '+m+' '+f+' '+b])
            y.append(v_dict['GAP experiment '+f])
            ax[1,1].plot(x,y,color=color_dict[m],marker=marker_dict[b],markersize=13,markeredgecolor='black',
                         linestyle=style_dict[m],mew=0.5,lw=2,label=r'%s, %s' % (m,b))

x0L,y0L,dxL,dyL = 0.18,0.03,0.40,0.20
h,l = ax[1,1].get_legend_handles_labels()
h   = [h[x] for x in [0,6,12,3,9,15]]
l   = [l[x] for x in [0,6,12,3,9,15]]
ax[1,1].legend(h,l,fancybox=True,shadow=True,ncol=2,loc=3,
               bbox_to_anchor=(x0L,y0L,dxL,dyL),handlelength=1.4,handletextpad=0.5,fontsize=fontsize-6,columnspacing=0.5)

# ==============================================================

fig.savefig(fname,format='eps')
os.system('ps2epsi '+fname)
os.system('mv '+fname+'i '+fname)


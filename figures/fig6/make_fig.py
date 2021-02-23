import os
import numpy              as np
import matplotlib.pyplot  as plt
from   matplotlib         import rc
from   utils              import *

fname    = 'fig6.eps'
fontsize = 22
rc('font',**{'family':'serif','serif':['Computer Modern']})
rc('text',usetex=True)
rc('axes',labelsize=fontsize,titlesize=fontsize)
rc('xtick',labelsize=fontsize-3)
rc('ytick',labelsize=fontsize-3)

color_dict  = {'overlap':c_list['orange'],'energy':c_list['orchid'],r'$T_0$ (exact)':c_list['mid-blue'],'GS (exact)':c_list['mid_green'],r'$S_1$ (exact)':c_list['mid-red']}
v_dict      = {r'$T_0$ (exact)':-599.081000,'GS (exact)':-773.123458,r'$S_1$ (exact)':-580.846171}
style       = {r'$T_0$ (exact)':'-.','GS (exact)':'--',r'$S_1$ (exact)':':'}

xfig = 9.0
yfig = 0.66*xfig

fig,ax = plt.subplots(2,2,figsize=(2*xfig,2*yfig))
fig.subplots_adjust(wspace=0.30,hspace=0.45)
# ==============================================================

for m,(name,(x,y)) in enumerate(zip([r'$(a)$',r'$(b)$',r'$(c)$',r'$(d)$'],[(0,0),(0,1),(1,0),(1,1)])):
    v = np.loadtxt('plot%s.txt' % str(m+1))
    if(x==0):
       xlim,xticks,xlab=[0,250],[0,50,100,150,200,250],[0,50,100,150,200,250],
       fill_panel(ax[x,y],'VQD steps',[0,250],[0,50,100,150,200,250],[0,50,100,150,200,250],
                          'energy [mHa]',[-800,-400],[-800,-700,-600,-500,-400],[-800,-700,-600,-500,-400],p=20.0,q=20.0)
       ax[x,y].text(250,-410,name,horizontalalignment='center',verticalalignment='center',fontsize=fontsize)
    if(x==1):
       xlim,xticks,xlab=[0,300],[0,60,120,180,240,300],[0,60,120,180,240,300],
       fill_panel(ax[x,y],'VQD steps',[0,300],[0,60,120,180,240,300],[0,60,120,180,240,300],
                          'energy [mHa]',[-800,-400],[-800,-700,-600,-500,-400],[-800,-700,-600,-500,-400],p=20.0,q=20.0)
       ax[x,y].text(300,-410,name,horizontalalignment='center',verticalalignment='center',fontsize=fontsize)

    ax2 = ax[x,y].twinx()
    fill_panel(ax2,'',xlim,xticks,xlab,
                   'overlap',[0,1],[0,0.2,0.4,0.6,0.8,1.0],['0.0','0.2','0.4','0.6','0.8','1.0'],p=20.0,q=20.0)

    vi,ve,vo = v[:,0],v[:,1],-800+400*v[:,2]
    ax[x,y].plot(vi,ve,color=color_dict['energy'],label='energy',lw=2.5)
    ax[x,y].plot(vi,vo,color=color_dict['overlap'],label='overlap',lw=2.5) 
    for k in v_dict.keys():
        if(k==r'$S_1$ (exact)' and x==0): print(k)
        elif(k=='GS (exact)' and x==1): print(k)
        else:
           ax[x,y].axhline(v_dict[k],c=color_dict[k],ls=style[k],label=k,lw=2)

    x0L,y0L = 0.00,1.05
    dxL,dyL = 1,0.20
    h,l = ax[x,y].get_legend_handles_labels()
    ax[x,y].legend(h,l,fancybox=True,shadow=True,ncol=4,loc=3,
                   bbox_to_anchor=(x0L,y0L,dxL,dyL),handlelength=1.5,fontsize=18,labelspacing=0.5,handletextpad=0.5,columnspacing=0.5)

fig.savefig(fname,format='eps')
os.system('ps2epsi '+fname)
os.system('mv '+fname+'i '+fname)


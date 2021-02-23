import numpy as np
from matplotlib import pyplot as plt

def draw_n_agon(pan,x0,y0,R,n,theta0,everyother=False):
    theta0=np.pi*theta0/180.0
    for i in range(n):
        xi = R*np.cos(theta0+2*np.pi*i/float(n))
        yi = R*np.sin(theta0+2*np.pi*i/float(n))
        j = i+1
        xj = R*np.cos(theta0+2*np.pi*j/float(n))
        yj = R*np.sin(theta0+2*np.pi*j/float(n))
        print(xi,yi,xj,yj)
        if(everyother):
           if(i%2==0):
              pan.plot([x0+xi,x0+xj],[y0+yi,y0+yj],ls='-',c='black',lw=2,zorder=-1)
        else:
           pan.plot([x0+xi,x0+xj],[y0+yi,y0+yj],ls='-',c='black',lw=2)

def draw_atom(pan,x0,y0,Rd,Z,c='black'):
    circle1 = plt.Circle((x0,y0),Rd*3.0/4.0,color='w',zorder=45)
    pan.add_patch(circle1)
    pan.text(x0,y0,Z,zorder=46,horizontalalignment='center',verticalalignment='center',c=c,fontsize=16)

def draw_geometry(pan,i_list):
    L   = 1.0
    eps = 0.05
    
    R0 = L/2.0/np.sin(2*np.pi/10.0)
    V0 = L*np.cos(2*np.pi/10.0)+R0*np.cos(2*np.pi/12.0)
    draw_n_agon(pan,x0=0,y0=0,R=R0,n=5,theta0=0,everyother=False)
    draw_atom(pan,R0,0,0.4,'N')
    
    # -----
    
    R1 = R0
    pan.plot([R1,R1+L],[0,0],ls='-',c='black',lw=2)
    c,s=np.cos(108*np.pi/180.0),np.sin(108*np.pi/180.0)
    R0 = L
    draw_n_agon(pan,x0=V0*c,y0=V0*s,R=R0,n=6,theta0=108-90,everyother=False)
    draw_n_agon(pan,x0=V0*c,y0=V0*s,R=0.8*R0,n=6,theta0=108-90,everyother=True)
    
    draw_n_agon(pan,x0=V0*c,y0=-V0*s,R=R0,n=6,theta0=360-108-90,everyother=False)
    draw_n_agon(pan,x0=V0*c,y0=-V0*s,R=0.8*R0,n=6,theta0=360-108-90,everyother=True)
    
    draw_n_agon(pan,x0=R1+2*L,y0=0,R=R0,n=6,theta0=0,everyother=False)
    draw_n_agon(pan,x0=R1+2*L,y0=0,R=0.8*R0,n=6,theta0=180,everyother=True)
   
    x0,y0 = R1+2*L,0
    for i in i_list:
        xi = L*np.cos(2*np.pi*i/float(6))
        yi = L*np.sin(2*np.pi*i/float(6))
        pan.plot([x0+xi,x0+2*xi],[y0+yi,y0+2*yi],ls='-',c='red',lw=2)
        draw_atom(pan,x0+2*xi,y0+2*yi,0.4,'F',c='r')
   
    pan.plot([R1+3*L,R1+4*L-eps],[0,0],ls='-',c='black',lw=2)
    draw_atom(pan,R1+4*L,0,0.4,'S')
    x0,y0 = R1+4*L,0
    pan.plot([x0-eps,x0-eps],[y0,y0+L],ls='-',c='black',lw=2)
    pan.plot([x0+eps,x0+eps],[y0,y0+L],ls='-',c='black',lw=2)
    pan.plot([x0-eps,x0-eps+np.sqrt(3.0)*L/2.0],[y0,y0+L*0.5],ls='-',c='black',lw=2)
    draw_atom(pan,x0,y0+L,0.4,'O')
    pan.plot([x0+eps,x0+eps+np.sqrt(3.0)*L/2.0],[y0-eps,y0+L*0.5-eps],ls='-',c='black',lw=2)
    draw_atom(pan,x0+np.sqrt(3.0)*L/2.0,y0+L*0.5,0.4,'O')
    dx,dy = 0.5*L,-np.sqrt(3)*L/2.0
    pan.plot([x0,x0+dx],[y0,y0+dy],ls='-',c='black',lw=2)
    
    x0,y0 = x0+2*dx,y0+2*dy
    draw_n_agon(pan,x0=x0,y0=y0,R=R0,n=6,theta0=0,everyother=False)
    draw_n_agon(pan,x0=x0,y0=y0,R=0.8*R0,n=6,theta0=180,everyother=True)

import os
from   matplotlib         import rc

fontsize = 30
rc('font',**{'family':'serif','serif':['Computer Modern']})
rc('text',usetex=True)
rc('axes',labelsize=fontsize,titlesize=fontsize)
rc('xtick',labelsize=fontsize-3)
rc('ytick',labelsize=fontsize-3)

xfig = 5.0
yfig = 0.7*xfig

fig,ax = plt.subplots(1,3,figsize=(3*xfig,1*yfig))
fig.subplots_adjust(wspace=0,hspace=0)

draw_geometry(ax[0],[]); ax[0].axis('off'); ax[0].set_aspect('equal'); ax[0].text(3,-3.2,'PSPCz',fontsize=28,horizontalalignment='center',verticalalignment='center')
draw_geometry(ax[1],[2,4]); ax[1].axis('off'); ax[1].set_aspect('equal'); ax[1].text(3,-3.2,'2F-PSPCz',fontsize=28,horizontalalignment='center',verticalalignment='center')
draw_geometry(ax[2],[1,2,4,5]); ax[2].axis('off'); ax[2].set_aspect('equal'); ax[2].text(3,-3.2,'4F-PSPCz',fontsize=28,horizontalalignment='center',verticalalignment='center')

fname = 'fig1.eps'

fig.savefig(fname,format='eps')
os.system('ps2epsi '+fname)
os.system('mv '+fname+'i '+fname)


import pyx
from pyx.graph.axis import linear, texter
text = pyx.text
x = 110
y = 45
t = -5
z = 75
w = 140
text.set(mode="latex")
text.preamble(r"\usepackage{ulem}") # underlining text...

text.preamble(r"\usepackage{anyfontsize}") # one way of scaling fonts (except

c = pyx.canvas.canvas()
c.insert(pyx.epsfile.epsfile(0,0,"PSPCz-HOMO.eps", align="tl"))
c.insert(pyx.epsfile.epsfile(x,0,"2F-PSPCz-HOMO.eps", align="tc"))
c.insert(pyx.epsfile.epsfile(2*x,0,"4F-PSPCz-HOMO.eps", align="tr"))
c.insert(pyx.epsfile.epsfile(0,0,"PSPCz-LUMO.eps", align="bl"))
c.insert(pyx.epsfile.epsfile(x,0,"2F-PSPCz-LUMO.eps", align="bc"))
c.insert(pyx.epsfile.epsfile(2*x,0,"4F-PSPCz-LUMO.eps", align="br"))
c.text(0,y,r'\fontsize{180}{5}\selectfont $(a)$')
c.text(z,y,r'\fontsize{180}{5}\selectfont $(b)$')
c.text(w,y,r'\fontsize{180}{5}\selectfont $(c)$')
c.text(0,t,r'\fontsize{180}{5}\selectfont $(d)$')
c.text(z,t,r'\fontsize{180}{5}\selectfont $(e)$')
c.text(w,t,r'\fontsize{180}{5}\selectfont $(f)$')
c.writeEPSfile("combined.eps")



'''import matplotlib.pyplot as plt
from PIL import Image

#from matplotlib.cbook import get_sample_data

#im = plt.imread('2F-PSPCz-HOMO.eps',format='eps') #get_sample_data('grace_hopper.jpg'))

fig,ax = plt.subplots(2,3)
plt.subplots_adjust(hspace=0.0,wspace=0.0)

for ic,c in enumerate(['PSPCz','2F-PSPCz','4F-PSPCz']):
    im = Image.open('%s-HOMO.eps' % c) #,format='eps') #plt.imread('%s-HOMO.eps' % c,format='eps')
    ax[0,ic].imshow(im)
    ax[0,ic].text(0.5,0.9,'%s, HOMO' % c)
#    ax[0,ic].axis('off')

    im = Image.open('%s-LUMO.eps' % c) #plt.imread('%s-LUMO.eps' % c)
    ax[1,ic].imshow(im)
    ax[1,ic].text(0.5,0.9,'%s, LUMO' % c)
#    ax[1,ic].axis('off')

fig.savefig('out.eps',format='eps')
'''


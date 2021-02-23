rm fig3-figure0.*
rm fig3.auxlock
rm fig3.aux
rm fig3.dvi
rm fig3.log

latex --shell-escape fig3.tex

epstopdf fig3-figure0.eps

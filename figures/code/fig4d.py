import matplotlib
matplotlib.rcParams['text.usetex'] = True
import matplotlib.pyplot as plt  #to implement plots
import csv                       #to read and process CSV files
import numpy as np               #to format numbers
from matplotlib.ticker import (MultipleLocator,FormatStrFormatter,AutoMinorLocator)

expt = []    #y-axis data
eesto = []   #x-axis data 1
qeomsto = [] #x-axis data 2
vqdsto = []  #x-axis data 3
eedz = []    #x-axis data 4
qeomdz = []  #x-axis data 5
vqddz = []   #x-axis data 6

#open CSV file
with open('/Users/gvnjones/Projects/gavinjones/mccoled/docs/npjcompmat/fig4d.csv') as csvDataFile:
    csvReader = csv.reader(csvDataFile)    #file handler
    for row in csvReader:                  #loop through file handler
        expt.append(row[1])                #load up second column element
        eesto.append(row[2])               #load up third column element
        qeomsto.append(row[3])
        vqdsto.append(row[4])
        eedz.append(row[5])
        qeomdz.append(row[6])
        vqddz.append(row[7])
fig, ax = plt.subplots()

exptFloat=list(map(float, expt))
eestoFloat=list(map(float, eesto))
qeomstoFloat=list(map(float, qeomsto))
vqdstoFloat=list(map(float, vqdsto))
eedzFloat=list(map(float, eedz))
qeomdzFloat=list(map(float, qeomdz))
vqddzFloat=list(map(float, vqddz))

print(exptFloat)

ax.plot(eestoFloat,exptFloat, label="Exact (STO-3G)", color="green", linestyle="dotted", marker="d")
ax.plot(qeomstoFloat,exptFloat, label="qEOM-VQE (STO-3G)", color="red", linestyle=(0, (1,10)), marker="d")
ax.plot(vqdstoFloat,exptFloat, label="VQD (STO-3G)", color="blue", linestyle=(0, (5,10)), marker="d")
ax.plot(eedzFloat,exptFloat, label="Exact (6-31G(d))", color="green", linestyle="dotted", marker="o")
ax.plot(qeomdzFloat,exptFloat, label="qEOM-VQE (6-31G(d))", color="red", linestyle=(0, (1,10)), marker="o")
ax.plot(vqddzFloat,exptFloat, label="VQD (6-31G(d))", color="blue", linestyle=(0, (5,10)), marker="o")

ax.set_xlabel(r'Experimental $ \Delta E_{ST}$', fontsize=16)
ax.set_ylabel(r'Calculated $ \Delta E_{ST}$', fontsize=16)
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
ax.get_xaxis().set_major_formatter(FormatStrFormatter('% 1.1f'))
ax.get_yaxis().set_major_formatter(FormatStrFormatter('% 1.1f'))
plt.legend(fontsize=12)                #displays legend
plt.tight_layout()          #auto adjust plot for better viewing
plt.savefig('/Users/gvnjones/Projects/gavinjones/mccoled/docs/npjcompmat/fig4d.pdf', dpi=900)
plt.show()                  #show the plot

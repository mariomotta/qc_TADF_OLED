import matplotlib.pyplot as plt  #to implement plots
import csv                       #to read and process CSV files
import numpy as np
from matplotlib.ticker import (MultipleLocator,FormatStrFormatter,AutoMinorLocator)

method = [] #x-axis data
szsto = [] #y-axis data 1
tfsto = [] #y-axis data 2
sfsto = [] #y-axis data 3
szdz = [] #y-axis data 4
tfdz = [] #y-axis data 5
sfdz = [] #y-axis data 6

#open CSV file
with open('/Users/gvnjones/Projects/gavinjones/mccoled/docs/npjcompmat/fig4a.csv') as csvDataFile:
    csvReader = csv.reader(csvDataFile)    #file handler
    for row in csvReader:                  #loop through file handler
        method.append(row[0])              #load up first column element into
        szsto.append(row[1])               #load up second column element into
        tfsto.append(row[2])
        sfsto.append(row[3])
        szdz.append(row[4])
        tfdz.append(row[5])
        sfdz.append(row[6])
fig, ax = plt.subplots()
plt.rcParams['font.size'] = '14'

szstoFloat=list(map(float, szsto))
tfstoFloat=list(map(float, tfsto))
sfstoFloat=list(map(float, sfsto))
szdzFloat=list(map(float, szdz))
tfdzFloat=list(map(float, tfdz))
sfdzFloat=list(map(float, sfdz))

ax.scatter(method,szstoFloat, color="green", marker="|", linewidth=20, label="S0")   #plot the data with + as the marker for each point
ax.scatter(method,tfstoFloat, color="red", marker="|", linewidth=20, label="T1")     #plot the data with diamond as the marker for each point
ax.scatter(method,sfstoFloat, color="blue", marker="|", linewidth=20, label="S1")    #plot the data with X as the marker for each point
ax.scatter(method,szdzFloat, color="green", marker=".", linewidth=6)    #plot the data with triangle as the marker for each point
ax.scatter(method,tfdzFloat, color="red", marker=".", linewidth=6)      #plot the data with point as the marker for each point
ax.scatter(method,sfdzFloat, color="blue", marker=".", linewidth=6)     #plot the data with point as the marker for each point

ax.set_title('(a) PSPCz')
ax.set_ylim(-850,-300)
ax.set_ylabel('Energy [Ha]', fontstyle='italic', fontsize=16)
ax.get_yaxis().set_major_formatter(FormatStrFormatter('% 1.1f'))
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
#plt.legend(loc=0, fontsize=10)                #displays legend
plt.tight_layout()          #auto adjust plot for better viewing
plt.savefig('/Users/gvnjones/Projects/gavinjones/mccoled/docs/npjcompmat/fig4a.pdf', dpi=900)
plt.show()                  #show the plot

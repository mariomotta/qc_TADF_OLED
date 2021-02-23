c_list = {'mid-red'      : '#e61e1e',
          'mid-green'    : '#5cbf00',
          'orange'       : '#FF8856',
          'mid-blue'     : '#1e5ae6',
          'gray'         : '#505050'}

def make_figure_5d(ax):

	import matplotlib
	matplotlib.rcParams['text.usetex'] = True
	import matplotlib.pyplot as plt  #to implement plots
	import csv						 #to read and process CSV files
	import numpy as np				 #to format numbers
	from matplotlib.ticker import (MultipleLocator,FormatStrFormatter,AutoMinorLocator)
	from scipy.ndimage.filters import gaussian_filter1d
	from scipy.interpolate import make_interp_spline, BSpline

	expt = []	 #y-axis data
	nmnm = []	#x-axis data 1
	romnm = [] #x-axis data 2
	romrom = []  #x-axis data 3
	qstrom = []    #x-axis data 4

	#open CSV file
	with open('../fig5d.csv') as csvDataFile:
		csvReader = csv.reader(csvDataFile)    #file handler
		for row in csvReader:				   #loop through file handler
			expt.append(row[1])				   #load up second column element
			nmnm.append(row[2])				  #load up third column element
			romnm.append(row[3])
			romrom.append(row[4])
			qstrom.append(row[5])

	exptFloat=list(map(float, expt))
	nmnmFloat=list(map(float, nmnm))
	romnmFloat=list(map(float, romnm))
	romromFloat=list(map(float, romrom))
	qstromFloat=list(map(float, qstrom))

	ax.plot(nmnmFloat,exptFloat, label="VQE (nm), qEOM (nm)",     color=c_list['gray'],      linestyle="dotted",    ms=12, marker="o",mew=0.5,mec='black')
	ax.plot(romnmFloat,exptFloat, label="VQE (rom), qEOM (nm)",   color=c_list['mid-red'],   linestyle=(0, (1,10)), ms=12, marker="D",mew=0.5,mec='black')
	ax.plot(romromFloat,exptFloat, label="VQE (rom), qEOM (rom)", color=c_list['mid-green'], linestyle=(0, (5,10)), ms=12, marker="s",mew=0.5,mec='black')
	ax.plot(qstromFloat,exptFloat, label="VQE (qst), qEOM (rom)", color=c_list['mid-blue'],  linestyle="dotted",    ms=16, marker="*",mew=0.5,mec='black')


	ax.set_title('(d)', fontsize=24,pad=10)
	ax.set_xlabel(r'Calculated $ \Delta E_{ST}$ [eV]', fontsize=19)
	ax.set_ylabel(r'Experimental $ \Delta E_{ST}$ [eV]', fontsize=19)
	ax.get_xaxis().set_major_formatter(FormatStrFormatter('% 1.1f'))
	ax.get_yaxis().set_major_formatter(FormatStrFormatter('% 1.1f'))


c_list = {'mid-red' : '#e61e1e',
        'mid-green' : '#5cbf00',
        'orange' : '#FF8856',
        'mid-blue' : '#1e5ae6',
        'gray' : '#505050'}

def make_figure_7d(ax):
		import matplotlib
		matplotlib.rcParams['text.usetex'] = True
		import matplotlib.pyplot as plt  #to implement plots
		import csv						 #to read and process CSV files
		import numpy as np				 #to format numbers
		from matplotlib.ticker import (MultipleLocator,FormatStrFormatter,AutoMinorLocator)

		mols = []
		expt = []	 #y-axis data
		vqdnm = []	 #x-axis data 1
		vqdqst = [] #x-axis data 2

		#open CSV file
		with open('../fig7d.csv') as csvDataFile:
			csvReader = csv.reader(csvDataFile)    #file handler
			for row in csvReader:				   #loop through file handler
				mols.append(row[0])
				expt.append(row[1])				   #load up second column element
				vqdnm.append(row[2])			   #load up third column element
				vqdqst.append(row[3])
		plt.rcParams['font.size'] = '16'

		exptFloat=list(map(float, expt))
		vqdnmFloat=list(map(float, vqdnm))
		vqdqstFloat=list(map(float, vqdqst))

		ax.plot(vqdnmFloat,exptFloat, label="VQD (nm)", color=c_list['mid-red'], linestyle="solid", marker="D", ms=10, lw =2,mec='black',mew=0.5)
		ax.plot(vqdqstFloat,exptFloat, label="VQD (qst)", color=c_list['mid-blue'], linestyle="dotted", marker="o", ms =10, lw=2,mec='black',mew=0.5)

		for x_pos, y_pos, mols in zip(vqdnmFloat, exptFloat, mols):
			ax.annotate(mols,						# The label for this point
						xy=(x_pos, y_pos),			# Position of the corresponding point
						xytext=(10,0),				# Offset text by 7 points to the right
						textcoords='offset points', # use offset points
						ha='left',					# Horizontally aligned to the left
						va='center')				# Vertical alignment is centered

		ax.set_title('(d)',pad=10)
		ax.set_xlabel(r'Calculated $ \Delta E_{ST}$ [eV]', fontsize=20)
		ax.set_ylabel(r'Experimental $ \Delta E_{ST}$ [eV]', fontsize=20)
		ax.get_xaxis().set_major_formatter(FormatStrFormatter('% 1.1f'))
		ax.get_yaxis().set_major_formatter(FormatStrFormatter('% 1.1f'))
		#plt.xticks(fontsize=15)
		#plt.yticks(fontsize=15)
		#plt.legend()				#displays legend
		#plt.tight_layout()			#auto adjust plot for better viewing
		#plt.savefig('/Users/gvnjones/Projects/gavinjones/mccoled/docs/npjcompmat/fig7d.pdf', dpi=900)
		#plt.show()					#show the plot

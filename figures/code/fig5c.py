c_list = {'mid-red'      : '#e61e1e',
          'mid-green'    : '#5cbf00',
          'orange'       : '#FF8856',
          'mid-blue'     : '#1e5ae6',
          'gray'         : '#505050'}

def make_figure_5c(ax):

		#import matplotlib
		#matplotlib.rcParams['text.usetex'] = True
		#import matplotlib.pyplot as plt  #to implement plots
		import csv						 #to read and process CSV files
		import numpy as np				 #to format numbers
		from matplotlib.ticker import (MultipleLocator,FormatStrFormatter,AutoMinorLocator)

		mols = []	 #x-axis data
		nmnm = []	#y-axis data 1
		romnm = [] #y-axis data 2
		romrom = []  #y-axis data 3
		qstrom = []    #y-axis data 4

		#open CSV file
		with open('../fig5c_s1.csv') as csvDataFile:
				csvReader = csv.reader(csvDataFile)    #file handler
				for row in csvReader:				   #loop through file handler
						mols.append(row[0])				   #load up second column element
						nmnm.append(row[1])					 #load up third column element
						romnm.append(row[2])
						romrom.append(row[3])
						qstrom.append(row[4])

		nmnmFloat=list(map(float, nmnm))
		romnmFloat=list(map(float, romnm))
		romromFloat=list(map(float, romrom))
		qstromFloat=list(map(float, qstrom))

		x = np.arange(len(mols))				   #label locations
		width = 0.1								  #width of the bars
		#plt.rcParams['font.size'] = '16'

		rects1=ax.bar(x - 0.15,nmnmFloat, width, label="VQE (nm), qEOM (nm)", color=c_list['gray'])
		rects2=ax.bar(x - 0.05,romnmFloat, width, label="VQE (rom), qEOM (nm)", color=c_list['mid-red'])
		rects3=ax.bar(x + 0.05,romrom, width, label="VQE (rom), qEOM (rom)", color=c_list['mid-green'])
		rects4=ax.bar(x + 0.15,qstromFloat, width, label="VQE (qst), qEOM (rom)", color=c_list['mid-blue'])

		def autolabel(rects, xpos='center',jr=0):
					ha = {'center': 'center', 'right': 'left', 'left': 'right'}
					offset = {'center': 0, 'right': 1, 'left': -1}
			 
					for rect in rects:
						y_value = rect.get_height()
						x_value = rect.get_x() + rect.get_width() / 2
						space =  5
						va = 'bottom'
						if y_value < 0:
							space =1 #*= -1
							#va = 'top'
						label = "{:.1f}".format(y_value)
						ax.annotate(
							label,						# Use `label` as label
							(x_value,25-3*(jr//3)),		# Place label at end of the bar
							xytext=(0, space),			# Vertically shift label by `space`
							textcoords="offset points", # Interpret `xytext` as offset in points
							ha='center',				# Horizontally center label
							va=va,fontsize=20)			# Vertically align label differently for
														# positive and negative values.
						if(y_value>0): ax.plot([x_value,x_value],[25-3*(jr//3)+0.5,y_value+0.5],linestyle='--',color='black',linewidth=0.8)
						else:		   ax.plot([x_value,x_value],[25-3*(jr//3)+0.5,0.5],linestyle='--',color='black',linewidth=0.8)
						jr += 1
					return jr
			 
		jr = 0
		jr = autolabel(rects1, "center",jr)
		jr = autolabel(rects2, "center",jr)
		jr = autolabel(rects3, "center",jr)
		jr = autolabel(rects4, "center",jr)

		ax.set_xticks(x)			#axis setting of X ticks
		ax.set_xticklabels(mols)	#axis setting of X tick labels
		ax.set_title('(c) $S_1$ States', fontsize=24,pad=10)
		ax.set_ylabel(r' $ \Delta \Delta E_{ST}$ [eV]', fontsize=18)
		ax.get_yaxis().set_major_formatter(FormatStrFormatter('% 1.1f'))

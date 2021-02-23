c_list = {'mid-red'      : '#e61e1e',
          'mid-green'    : '#5cbf00',
          'orange'       : '#FF8856',
          'mid-blue'     : '#1e5ae6',
          'gray'         : '#505050'}

def make_figure_7b(ax):
		import matplotlib
		matplotlib.rcParams['text.usetex'] = True
		import matplotlib.pyplot as plt  #to implement plots
		import csv						 #to read and process CSV files
		import numpy as np				 #to format numbers
		from matplotlib.ticker import (MultipleLocator,FormatStrFormatter,AutoMinorLocator)

		es = []    #x-axis data
		qstnm = []	 #y-axis data 1
		qstqst = [] #y-axis data 2

		#open CSV file
		with open('../fig7b.csv') as csvDataFile:
			csvReader = csv.reader(csvDataFile)    #file handler
			for row in csvReader:				   #loop through file handler
				es.append(row[0])				   #load up second column element
				qstnm.append(row[1])			   #load up third column element
				qstqst.append(row[2])

		qstnmFloat=list(map(float, qstnm))
		qstqstFloat=list(map(float, qstqst))

		#fig, ax = plt.subplots()
		x = np.arange(len(es))					   #label locations
		width = 0.1								   #width of the bars
		plt.rcParams['font.size'] = '16'

		rects1=ax.bar(x - 0.15,qstnmFloat, width, label=r'$T_1$', color=c_list['mid-red'])
		rects2=ax.bar(x - 0.05,qstqstFloat, width, label=r'$S_1$', color=c_list['mid-blue'])

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
								label,											# Use `label` as label
								(x_value,20-3*(jr//2)),			# Place label at end of the bar
								xytext=(0, space),						# Vertically shift label by `space`
								textcoords="offset points", # Interpret `xytext` as offset in points
								ha='center',							# Horizontally center label
								va=va,fontsize=20)						# Vertically align label differently for
																						# positive and negative values.
						if(y_value>0): ax.plot([x_value,x_value],[20-3*(jr//2)+0.5,y_value+0.5],linestyle='--',color='black',linewidth=0.8)
						else:			   ax.plot([x_value,x_value],[20-3*(jr//2)+0.5,0.5],linestyle='--',color='black',linewidth=0.8)
						jr += 1
				return jr

		jr = 0
		jr = autolabel(rects1, "center",jr)
		jr = autolabel(rects2, "center",jr)

		ax.set_xticks(x)			#axis setting of X ticks
		ax.set_xticklabels(es)		#axis setting of X tick labels
		ax.set_title('(b) 2F-PSPCz',pad=10)
		ax.set_ylabel(r' $ \Delta \Delta E_{Exact}$',fontsize=20)
		ax.get_yaxis().set_major_formatter(FormatStrFormatter('% 1.1f'))

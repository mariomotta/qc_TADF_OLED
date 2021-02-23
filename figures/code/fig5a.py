c_list = {'mid-red'      : '#e61e1e',
          'mid-green'    : '#5cbf00',
          'orange'       : '#FF8856',
          'mid-blue'     : '#1e5ae6',
          'gray'         : '#505050'}

def make_figure_5a(ax):
    #import matplotlib
    #matplotlib.rcParams['text.usetex'] = True
    import matplotlib.pyplot as plt  #to implement plots
    import csv                       #to read and process CSV files
    import numpy as np               #to format numbers
    from matplotlib.ticker import (MultipleLocator,FormatStrFormatter,AutoMinorLocator)
   
    mols = [] #x-axis data
    nm   = [] #y-axis data 1
    rom  = [] #y-axis data 2
    nul  = [] #y-axis data 3
    qst  = [] #y-axis data 4
    ee   = [] #x-axis data 5
    
    #open CSV file
    with open('../fig5a_s0.csv') as csvDataFile:
        csvReader = csv.reader(csvDataFile)    #file handler
        for row in csvReader:                  #loop through file handler
            mols.append(row[0])                #load up second column element
            nm.append(row[1])                  #load up third column element
            rom.append(row[2])
            nul.append(row[3])
            qst.append(row[4])
   
    nmFloat=list(map(float, nm))
    romFloat=list(map(float, rom))
    qstFloat=list(map(float, qst))
    
    #fig, ax = plt.subplots()
    x = np.arange(len(mols))                   #label locations
    width = 0.1                               #width of the bars
    #plt.rcParams['font.size'] = '16'
   
    rects1=ax.bar(x - 0.15,nmFloat, width, label="VQE (nm)",   color = c_list['gray'])
    rects2=ax.bar(x - 0.05,romFloat, width, label="VQE (rom)", color = c_list['mid-red'])
    rects3=ax.bar(x + 0.05,nul, width)
    rects4=ax.bar(x + 0.15,qstFloat, width, label="VQE (qst)", color = c_list['mid-blue'])   

    def autolabel(rects, xpos='center',jr=0):
        ha = {'center': 'center', 'right': 'left', 'left': 'right'}
        offset = {'center': 0, 'right': 1, 'left': -1}

        for rect in rects:
            y_value = rect.get_height()
            x_value = rect.get_x() + rect.get_width() / 2
            space =  5 
            va = 'bottom'
            if y_value < 0:
                space *= -1
                va = 'top'
            label = "{:.1f}".format(y_value)
            ax.annotate(
                label,                      # Use `label` as label
                (x_value,17-4*(jr//3)),     # Place label at end of the bar
                xytext=(0, space),          # Vertically shift label by `space`
                textcoords="offset points", # Interpret `xytext` as offset in points
                ha='center',                # Horizontally center label
                va=va,fontsize=20)          # Vertically align label differently for
                                            # positive and negative values.
            ax.plot([x_value,x_value],[17-4*(jr//3)+0.5,y_value+0.5],linestyle='--',color='black',linewidth=0.8)
            jr += 1
        return jr
   
    jr = 0
    jr = autolabel(rects1, "center", jr)
    jr = autolabel(rects2, "center", jr)
    jr = autolabel(rects4, "center", jr)
    
    ax.set_xticks(x)            #axis setting of X ticks
    ax.set_xticklabels(mols)    #axis setting of X tick labels
    ax.set_title('(a) $S_0$ states', fontsize=24,pad=10)
    ax.set_ylabel(r' $ \Delta \Delta E_{ST}$ [eV]', fontsize=18)
    ax.get_yaxis().set_major_formatter(FormatStrFormatter('% 1.1f'))

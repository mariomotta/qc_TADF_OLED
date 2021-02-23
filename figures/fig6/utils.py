c_list = {'purple'       : '#B284BE',
          'jacaranda'    : '#888FC7',
          'light_yellow' : '#EEEA62',
          'gold'         : '#FFD300',
          'earwax'       : '#DAA520',
          'brown'        : '#7B3F00',
          'light_blue'   : '#bcd9ea',
          'cerulean'     : '#5ba4cf',
          'cobalt'       : '#0079bf',
          'dark_blue'    : '#055a8c',
          'light_green'  : '#acdf87',
          'yellow_green' : '#B9D146',
          'mid_green'    : '#68bb59',
          'dark_green'   : '#1e5631',
          'orchid'       : '#DA70D6',
          'orange'       : '#FFA500',
          'red'          : '#DC343B',
          'light-gray'   : '#C0C0C0',
          'palatinate'   : '#72246C',
          'mid-red'      : '#e61e1e',
          'mid-green'    : '#5cbf00',
          'mid-blue'     : '#1e5ae6'}

def fill_panel(pan,xlabel,xlim,xticks,xticklabels,ylabel,ylim,yticks,yticklabels,p=20.0,q=10.0):
    x0,x1 = xlim
    xlim  = [x0-(x1-x0)/p,x1+(x1-x0)/p]
    pan.set_xlabel(xlabel)
    pan.set_xlim(xlim)
    pan.set_xticks(xticks)
    pan.set_xticklabels(xticklabels)
    pan.set_ylabel(ylabel)
    y0,y1 = ylim
    ylim  = [y0-(y1-y0)/q,y1+(y1-y0)/q]
    pan.set_ylim(ylim)
    pan.set_yticks(yticks)
    pan.set_yticklabels(yticklabels)
    pan.tick_params(direction='in',which='both')



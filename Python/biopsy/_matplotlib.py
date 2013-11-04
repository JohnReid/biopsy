#
# Copyright John Reid 2006
#

def set_matplotlib_params_for_eps( fig_width_pt = 360.0 ):
    import pylab

    # fig_width_pt = 246.0  # Get this from LaTeX using \showthe\columnwidth

    inches_per_pt = 1.0/72.27               # Convert pt to inch
    golden_mean = (pylab.sqrt(5)-1.0)/2.0         # Aesthetic ratio
    fig_width = fig_width_pt*inches_per_pt  # width in inches
    fig_height = fig_width*golden_mean      # height in inches
    fig_size =  [fig_width,fig_height]

    params = {
            'backend': 'PS',
            'axes.labelsize': 11,
            'text.fontsize': 11,
            'xtick.labelsize': 8,
            'ytick.labelsize': 8,
            'text.usetex': False,
            'figure.figsize': fig_size
    }

    pylab.rcParams.update(params)

def _example_matplotlib_plot():
    import pylab
    from pylab import arange,pi,sin,cos,sqrt

    # Generate data
    x = arange(-2*pi,2*pi,0.01)
    y1 = sin(x)
    y2 = cos(x)

    # Plot data
    pylab.ioff()
    # print 'Plotting data'
    pylab.figure(1)
    pylab.clf()
    # print 'Setting axes'
    pylab.axes([0.125,0.2,0.95-0.125,0.95-0.2])
    # print 'Plotting'
    pylab.plot(x,y1,'g:',label='sin(x)')
    pylab.plot(x,y2,'-b',label='cos(x)')
    # print 'Labelling'
    pylab.xlabel('x (radians)')
    pylab.ylabel('y')
    pylab.legend()
    print 'Saving to fig1.eps'
    pylab.savefig('fig1.eps')
    pylab.ion()



if __name__ == '__main__':
    print 'setting matplotlib parameters'
    set_matplotlib_params_for_eps( fig_width_pt = 246.0 )
    print 'example plot'
    _example_matplotlib_plot()

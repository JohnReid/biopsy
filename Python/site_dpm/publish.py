#
# Copyright John Reid 2008, 2009, 2010, 2011
#

"""
Code to configure pylab for publish quality figures.
"""

def configure_pylab_for_publish(fig_width_pt=335.):
    """
    Configures pylab for publishable quality figures.
    
    Get fig_width_pt from LaTeX using \showthe\columnwidth
    """
    import pylab
    from math import sqrt
    inches_per_pt = 1/72.27                            # Convert pt to inches
    golden_ratio_conjugate = (sqrt(5)-1.0)/2.0         # Aesthetic ratio
    fig_width = fig_width_pt*inches_per_pt             # width in inches
    fig_height = fig_width*golden_ratio_conjugate      # height in inches
    fig_size = (fig_width, fig_height)
    params = {
      'backend' : 'eps',
      'axes.labelsize'  : 8,
      'axes.titlesize'  : 8,
      'text.fontsize'   : 8,
      'legend.fontsize' : 8,
      'xtick.labelsize' : 6,
      'ytick.labelsize' : 6,
      'text.usetex' : True,
      'figure.figsize' : fig_size,
      'figure.subplot.top'    : 0.93,
      'figure.subplot.bottom' : 0.1,
      'figure.subplot.left'   : 0.15,
      'figure.subplot.right'  : 0.98,
      'patch.facecolor'    : '.6',
      'lines.color'        : 'black',
      'axes.color_cycle' : ['gray', 'k', 'g', 'r', 'c', 'm', 'y', 'b'],
    }
    pylab.rcParams.update(params)
    return fig_size


# Settings to make LaTeX-compatible figures.
#
# Use as follows:
#
# with matplotlib.rc_context(rc = RCPARAMS_LATEX_DOUBLE_COLUMN):
#     fig, ax = plt.subplots()
#     ax.plot(...)
#     ...
#     save_figure(fig, 'very-plot-such-amazing-wow')
#
# Amit Moscovich, Tel Aviv University, 2022.

import matplotlib

FIGURES_PATH = 'figures/'
DPI = 600


_RCPARAMS_LATEX_SINGLE_COLUMN = {
    'font.family': 'serif',
    'text.usetex': True,

    'axes.labelsize': 12,
    'axes.titlesize': 12,
    'legend.fontsize': 12,
    'xtick.labelsize': 11,
    'ytick.labelsize': 12,
     
    #'axes.prop_cycle': matplotlib.pyplot.cycler('color', ['#006496', '#ff816b', '#fbca60', '#6d904f', '#8b8b8b']) + matplotlib.pyplot.cycler('marker', ['o', 'd', 's', '*', '>']),
    'axes.prop_cycle': matplotlib.pyplot.cycler('color', ['#ff7d66', '#ffdc30', '#40a0cc', '#529915', '#8b8b8b']) + matplotlib.pyplot.cycler('marker', ['d', 's', 'o', r'$\clubsuit$', '>']),

    'lines.markersize': 9,
    'lines.markeredgewidth': 0.75,
    'lines.markeredgecolor': 'k',
                               
    'grid.color': '#C0C0C0', # 25% black

    'legend.fancybox': True, # Rounded legend box
    'legend.framealpha': 0.8,

    'axes.linewidth': 1,
}

# This is the right width (in inches) for a 'letter' page LaTeX document that imports the geometry package with default parameters.
_PAGE_WIDTH_INCHES = 6.775
_GOLDEN_RATIO = (5**0.5 - 1)/2
_WIDTH = _PAGE_WIDTH_INCHES
_HEIGHT = _GOLDEN_RATIO*_WIDTH
RCPARAMS_LATEX_DOUBLE_COLUMN = {**_RCPARAMS_LATEX_SINGLE_COLUMN, 'figure.figsize': (_WIDTH/2, _HEIGHT/2)}
RCPARAMS_LATEX_SINGLE_COLUMN_WIDE = {**_RCPARAMS_LATEX_SINGLE_COLUMN, 'figure.figsize': (_WIDTH, _HEIGHT/2)}
RCPARAMS_LATEX_SINGLE_COLUMN_LARGE = {**_RCPARAMS_LATEX_SINGLE_COLUMN, 'figure.figsize': (_WIDTH, _HEIGHT)}
RCPARAMS_LATEX_SINGLE_COLUMN_LARGE_SHORT = {**_RCPARAMS_LATEX_SINGLE_COLUMN, 'figure.figsize': (_WIDTH, _HEIGHT*0.75)}
RCPARAMS_LATEX_SINGLE_COLUMN_LARGE_TALL = {**_RCPARAMS_LATEX_SINGLE_COLUMN, 'figure.figsize': (_WIDTH, 1.2*_WIDTH)}



def save_figure(fig, name):
    import os
    os.makedirs(FIGURES_PATH, exist_ok=True)
    filename = os.path.join(FIGURES_PATH, name).replace('.','_') + '.pdf'
    print(f'Saving figure to "{os.path.realpath(filename)}"')
    fig.savefig(filename, dpi=DPI, bbox_inches='tight')


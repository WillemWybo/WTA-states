import os
import matplotlib
if os.uname()[0] == 'Linux' and (os.uname()[1] == 'pc59' or os.uname()[1] == 'pc58'):
    matplotlib.use("Agg")

import matplotlib.pyplot as pl
import matplotlib.animation as manimation
from matplotlib.gridspec import GridSpec
from matplotlib.patches import Rectangle
import matplotlib.patches as mpatches
import matplotlib.markers as mmark
import matplotlib.lines as mlines
import matplotlib.colors as mcolors
import matplotlib.axes as maxes
import matplotlib.patheffects as mpatheffects
import matplotlib.ticker as mticker
import matplotlib.transforms as mtransforms
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib import rc, rcParams
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset
from matplotlib.offsetbox import AnchoredText
from mpl_toolkits.mplot3d import proj3d
from mpl_toolkits.axisartist.axislines import SubplotZero

import numpy as np

import string
import copy

# colours = ['DeepPink', 'Purple', 'MediumSlateBlue', 'Blue', 'Teal',
#                 'ForestGreen',  'DarkOliveGreen', 'DarkGoldenRod',
#                 'DarkOrange', 'Coral', 'Red', 'Sienna', 'Black', 'DarkGrey']
colours = list(pl.rcParams['axes.prop_cycle'].by_key()['color'])
hatchpatterns = [ "|" , "\\" , "/" , "+" , "-", ".", "*", "x", "o", "O" ]

# matplotlib settings
legendsize = 12
labelsize = 18
ticksize = 15
lwidth = 1.5
markersize = 6.
fontsize = 16
lettersize = 26.
#~ font = {'family' : 'serif',
        #~ 'weight' : 'normal',
        #~ 'size'   : fontsize}
        #'sans-serif':'Helvetica'}
#'family':'serif','serif':['Palatino']}
#~ rc('font', **font)
# rc("font", family='sans-serif')#, weight="normal", size="18")
rc('font',**{'family':'sans-serif','sans-serif':['stixsans', 'helvetica', 'arial']})
rc('mathtext',**{'fontset': 'stixsans'})
# rc('text', usetex=True)
# rcParams['text.latex.preamble'] += r"\usepackage{amsmath}\usepackage{xfrac}"
# rc('legend',**{'fontsize': 'medium'})
# rc('xtick',**{'labelsize': 'small'})
# rc('ytick',**{'labelsize': 'small'})
# rc('axes',**{'labelsize': 'large', 'labelweight': 'normal'})
rc('legend',**{'fontsize': labelsize})
rc('xtick',**{'labelsize': ticksize})
rc('ytick',**{'labelsize': ticksize})
rc('axes',**{'labelsize': labelsize, 'labelweight': 'normal'})

# transparent backgrounds for everything
rcParams.update({
    "figure.facecolor":  (1.0, 1.0, 1.0, 0.0),  # red   with alpha = 30%
    "axes.facecolor":    (1.0, 1.0, 1.0, 0.0),  # green with alpha = 50%
    "savefig.facecolor": (1.0, 1.0, 1.0, 0.0),  # blue  with alpha = 20%
})


cs = ['r', 'b', 'g', 'c', 'y']
css = ['r', 'b', 'g', 'c', 'y']
cfl = ['fuchsia', 'lime']
cll = [colours[4], colours[2]]
mfs = ['D', 'o', 'v', '^', 's', 'p']
mls = ['+', '*', 'x', '1', '2']
lss = ['-', '--', '-.', ':']
cmap = pl.get_cmap('jet')


def getXCoords(spacings):
    coords = np.cumsum(spacings)
    coords /= coords[-1]
    return coords

def getYCoords(spacings):
    coords = np.cumsum(spacings)
    coords /= coords[-1]
    coords = coords[:-1]
    return coords

def myAx(ax):
    # customize the ax
    ax.spines['top'].set_color('none')
    ax.spines['right'].set_color('none')
    ax.yaxis.set_ticks_position('left')
    ax.xaxis.set_ticks_position('bottom')

    return ax

def noFrameAx(ax):
    ax.spines['top'].set_color('none')
    ax.spines['bottom'].set_color('none')
    ax.spines['right'].set_color('none')
    ax.spines['left'].set_color('none')
    ax.draw_frame = False

    ax.set_xticks([])
    ax.set_yticks([])

    return ax

def drawScaleBars(ax,
        xlabel=None, ylabel=None,
        lx_offset=.1, ly_offset=.1,
        bx_offset=.05, by_offset=.05, bc_offset=.1,
        fstr_xlabel=r'%.2g ', fstr_ylabel=r'%.2g ',
        text_kwargs_x=dict(size=ticksize, rotation=0, va='center'),
        text_kwargs_y=dict(size=ticksize, rotation=90, ha='center'),
    ):
    """
    scalebar location in axes coordinates:
    lx(y)_offset: vertical (horizontal) offset of the x(y)-label
    bx(y)_offset: vertical (horizontal) offset of the x(y)-scalebar
    bc_offset: bottom right corner offset of the scalebars
    """
    ax.spines['top'].set_color('none')
    ax.spines['bottom'].set_color('none')
    ax.spines['right'].set_color('none')
    ax.spines['left'].set_color('none')
    ax.set_xticklabels([])
    ax.set_yticklabels([])

    if xlabel is not None:
        xticks = ax.get_xticks()
    if ylabel is not None:
        yticks = ax.get_yticks()
    xlim = ax.get_xlim()
    ylim = ax.get_ylim()

    transf = ax.transData.inverted() + ax.transAxes
    p0, p1 = transf.transform((0., 0.)), transf.transform((-bc_offset,bc_offset))
    bcx_offset = p1[0] - p0[0]
    bcy_offset = p1[1] - p0[1]
    p0, p1 = transf.transform((0., 0.)), transf.transform((by_offset,-bx_offset))
    bx_offset = p1[1] - p0[1]
    by_offset = p1[0] - p0[0]
    p0, p1 = transf.transform((0., 0.)), transf.transform((ly_offset,-lx_offset))
    lx_offset = p1[1] - p0[1]
    ly_offset = p1[0] - p0[0]

    if xlabel is not None:
        # position and length
        sblen = xticks[-1] - xticks[-2]
        xpos = xlim[1] + bcx_offset
        ypos = ylim[0] + bx_offset

        px = (xpos - sblen / 2., ypos)
        xbar = ((xpos - sblen, xpos), (ypos, ypos))

        # draw the scale bar
        ax.plot(*xbar, 'k-', lw=1.5*lwidth, clip_on=False)
        ax.annotate(fstr_xlabel%sblen + xlabel,
                        xy=px, xytext=(px[0], px[1]+lx_offset), annotation_clip=False, transform=ax.transData,
                        ha='center',
                        **text_kwargs_x)

    if ylabel is not None:
        # position and length
        sblen_ = yticks[1] - yticks[0]
        ypos_ = ylim[0] + bcy_offset
        xpos_ = xlim[1] + by_offset

        py = (xpos_, ypos_ + sblen_ / 2.)
        ybar = ((xpos_, xpos_), (ypos_, ypos_ + sblen_))

        # draw y scalebar
        ax.plot(*ybar, 'k-', lw=1.5*lwidth, clip_on=False)
        ax.annotate(fstr_ylabel%sblen_ + ylabel,
                        xy=py, xytext=(py[0]+ly_offset, py[1]), annotation_clip=False, transform=ax.transData,
                        va='center',
                        **text_kwargs_y)

    ax.tick_params(axis='both', which='both', length=0, color='none')

    ax.set_xlim(xlim)
    ax.set_ylim(ylim)


def plotSpikeRaster(ax, spike_arr, y0=0., margin=.2, cs='k', plotargs={},
                                   inv_order=True, tlim=None, no_frame=True):
    if no_frame:
        ax = noFrameAx(ax)

    ny = len(spike_arr)
    hy = 1.

    ypos = np.arange(ny) * (hy+margin)
    ycoo = [[y0+yp, y0+yp+hy] for yp in ypos]
    if inv_order:
        ycoo = ycoo[::-1]
    yloc = np.mean(ycoo, axis=1)

    if not isinstance(cs, list):
        cs = [cs for _ in range(ny)]
    else:
        assert len(cs) == ny

    if 'c' in plotargs:
        del plotargs['c']

    for ii in range(ny):
        for tsp in spike_arr[ii]:
            ax.plot([tsp, tsp], ycoo[ii], c=cs[ii], **plotargs)

    if tlim is not None:
        ax.set_xlim(tlim)

    return yloc


def myLegend(ax, add_frame=True, **kwarg):
    leg = ax.legend(**kwarg)
    if add_frame:
        frame = leg.get_frame()
        frame.set_color('white')
        frame.set_alpha(0.8)
    else:
        frame = leg.get_frame()
        frame.set_alpha(0.)
    return leg


def myColorbar(ax, im, **kwargs):
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", "5%", pad="3%")
    return pl.colorbar(im, cax=cax, **kwargs)


def getAnnotations(*args, size_factor=1.0):
    assert len(args) > 0
    if isinstance(args[0], int):
        n_panel = args[0]
        alphabet = string.ascii_uppercase
        labels = [AnchoredText(r''+letter, loc=2, prop=dict(size=lettersize*size_factor),
                               pad=0., borderpad=-1.5, frameon=False) \
                  for letter in alphabet[:n_panel]]
    else:
        labels = [AnchoredText(labelstr, loc=2, prop=dict(size=lettersize*size_factor),
                               pad=0., borderpad=-1.5, frameon=False) \
                  for labelstr in args]
    return labels


class TransformedCMap(mcolors.Colormap):
    def __init__(self, func, cmap):
        '''
        `func`: bijective function on interval [0,1]
        '''
        self.func = func
        self.cmap = cmap

        # copy cmap variables to this class
        orig_keys = set(cmap.__dict__.keys())
        for key in orig_keys:
            self.__dict__[key] = copy.deepcopy(cmap.__dict__[key])

    def __call__(self, x, alpha=None, bytes=False):
        return self.cmap(self.func(x), alpha=alpha, bytes=bytes)


def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=100):
    new_cmap = mcolors.LinearSegmentedColormap.from_list(
        'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval),
        cmap(np.linspace(minval, maxval, n)))
    return new_cmap


class Arrow3D(mpatches.FancyArrowPatch):
    """
    3D arrow
    """
    def __init__(self, xs, ys, zs, *args, **kwargs):
        mpatches.FancyArrowPatch.__init__(self, (0,0), (0,0), *args, **kwargs)
        self._verts3d = xs, ys, zs

    def draw(self, renderer):
        xs3d, ys3d, zs3d = self._verts3d
        xs, ys, zs = proj3d.proj_transform(xs3d, ys3d, zs3d, renderer.M)
        self.set_positions((xs[0],ys[0]),(xs[1],ys[1]))
        mpatches.FancyArrowPatch.draw(self, renderer)


def multicolor_label(ax, list_of_strings, list_of_colors, axis='x', anchorpad=0, **kw):
    """this function creates axes labels with multiple colors
    ax specifies the axes object where the labels should be drawn
    list_of_strings is a list of all of the text items
    list_if_colors is a corresponding list of colors for the strings
    axis='x', 'y', or 'both' and specifies which label(s) should be drawn"""
    from matplotlib.offsetbox import AnchoredOffsetbox, TextArea, HPacker, VPacker

    # x-axis label
    if axis=='x' or axis=='both':
        boxes = [TextArea(text, textprops=dict(color=color, ha='left',va='bottom',**kw))
                    for text,color in zip(list_of_strings,list_of_colors) ]
        xbox = HPacker(children=boxes,align="center",pad=0, sep=5)
        anchored_xbox = AnchoredOffsetbox(loc=3, child=xbox, pad=anchorpad,frameon=False,bbox_to_anchor=(0.2, -0.09),
                                          bbox_transform=ax.transAxes, borderpad=0.)
        ax.add_artist(anchored_xbox)

    # y-axis label
    if axis=='y' or axis=='both':
        boxes = [TextArea(text, textprops=dict(color=color, ha='left',va='bottom',rotation=90,**kw))
                     for text,color in zip(list_of_strings[::-1],list_of_colors[::-1]) ]
        ybox = VPacker(children=boxes,align="center", pad=0, sep=5)
        anchored_ybox = AnchoredOffsetbox(loc=3, child=ybox, pad=anchorpad, frameon=False, bbox_to_anchor=(-0.10, 0.2),
                                          bbox_transform=ax.transAxes, borderpad=0.)
        ax.add_artist(anchored_ybox)


def curlyBrace(*args, **kwargs):
    curlybrace.curlyBrace(*args, **kwargs)


def arrowed_spines(ax, arrow_style='-|>', arrow_size=10):
    # for direction in ["xzero", "yzero"]:
    #     # adds arrows at the ends of each axis
    #     ax.axis[direction].set_axisline_style(arrow_style)

    #     # adds X and Y-axis from the origin
    #     ax.axis[direction].set_visible(True)

    # for direction in ["left", "right", "bottom", "top"]:
    #     # hides borders
    #     ax.axis[direction].set_visible(False)


    ax.plot((1), (0), ls="-", marker=">", ms=arrow_size, color="k",
            transform=ax.get_yaxis_transform(), clip_on=False)
    ax.plot((0), (1), ls="", marker="^", ms=arrow_size, color="k",
            transform=ax.get_xaxis_transform(), clip_on=False)




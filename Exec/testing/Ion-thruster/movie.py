
import sys  # nopep8
cmd_folder = "../../../vis"  # nopep8
if cmd_folder not in sys.path:  # nopep8
    sys.path.insert(0, cmd_folder)

from tile_mov import tile_movie
from make_mov import make_all, get_particle_trajectories
import numpy as np
import pylab as plt
import matplotlib as mpl
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.image import NonUniformImage
from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle

import pdb 

# ==============================================================================
# MAKE MOVIES
# ==============================================================================

def get_velocity_magnitude(ds, c):
    x, u = ds.get("x_vel-electron")
    x, v = ds.get("y_vel-electron")

    return {"x":x[0], "y":x[1], "value":np.sqrt(u**2 + v**2)}

def get_boxes(ds, c):
    boxes = ds.get_boxes()
    return {"boxes":boxes}

def plot(frame, data, output_name):

    xc = data["vel_mag"]["x"]
    yc = data["vel_mag"]["y"]
    v = data["vel_mag"]["value"]

    boxes = data["boxes"]["boxes"][()]

    # print('no crash 0')


    vmax = frame["vel_mag"]["max"]
    vmin = frame["vel_mag"]["min"] #-vmax 

    limits = frame["q"]["xy_limits"]

    # print('no crash 1')
    fig = plt.figure()
    ax = fig.add_subplot(111)



    norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)
    im = NonUniformImage(ax, interpolation='bilinear', extent=(limits[0][0], limits[1][0], limits[0][1], limits[1][1]),
        norm=norm)
    # print('no crash 2')


    im.set_data(xc, yc,np.transpose(v)) #np.transpose(a)
    ax.images.append(im)
    # print('no crash 3')

    # cs = ax.contour(xc, yc, np.transpose(v), levels=[0.5], colors=['k'], linewidths=[0.25])
    # line_x = cs.allsegs[0][0][:,0].tolist()
    # line_y = cs.allsegs[0][0][:,1].tolist()
    # plt.fill(line_x+[limits[1][0]], line_y+[limits[0][1]], color='k')

    ax.text(0.05, 0.05, 'velocity', horizontalalignment='left', 
        verticalalignment='bottom', transform=ax.transAxes, fontsize=18, color='w')

    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    plt.colorbar(im, cax)
   


    # plot boxes, i believe this is the amr grids
    # grid = []
    # for box in boxes:
    #     sz = box[1] - box[0]
    #     rect = Rectangle((box[0][0], box[0][1]), sz[0], sz[1])
    #     grid.append(rect)

    # pc = PatchCollection(grid, facecolor='none', alpha=1.0, edgecolor='w', linewidth=0.25)
    # ax.add_collection(pc)
    #AMR grids
    # pc = PatchCollection(grid, facecolor='none', alpha=1.0, edgecolor='w', linewidth=0.25)
    # ax.add_collection(pc)

    ax.set_xlim(limits[0][0], limits[1][0])
    ax.set_ylim(limits[0][1], limits[1][1])

    ax.set_aspect(1)
    ax.axes.xaxis.set_visible(False)
    ax.axes.yaxis.set_visible(False)
    fig.tight_layout()
    fig.savefig(output_name, dpi=300, bbox_inches="tight")
    plt.close(fig)

    return



dt = 0.005

Q = []

q = {}
q["files_dir"] = "./plotfiles"
q["level"] = -1

# all the data we need to retrieve
q["get"] = [

    {"func":get_velocity_magnitude, "tag":"vel_mag"},
    {"func":get_boxes, "tag":"boxes"},
]

# how to make a frame
q["plot"] = plot
q["name"] = "movie"

##
q["framerate"] = 60
q["mov_save"] = q["files_dir"] + "/mov"
q["offset"] = [0.0, 0.0]
q["xy_limits"] = [[0.0, 0.0], [4.0, 2.0]]
q["file_include"] = ["plt"]
q["file_exclude"] = ["chk"]
q["cores"] = 1
q["time_span"] = [] #np.arange(0,0.1,0.005).tolist()
q["force_data"] = True
q["force_frames"] = True
q["only_frames"] = False
q["redo_streaks"] = False
q["dpi"] = 300

q["normalize"] = "all"

Q.append(q)

make_all(Q)

print("DONE")


import sys, pdb
cmd_folder = "../../../vis"
if cmd_folder not in sys.path:
    sys.path.insert(0, cmd_folder)
    
from get_boxlib import ReadBoxLib, get_files

import numpy as np
import pylab as plt
from matplotlib.image import NonUniformImage
import matplotlib.ticker as ticker
import h5py
import pdb;

#==============================================================================
# 
#==============================================================================

plt.rc("font", family="serif")
plt.rc("font", size=14)
plt.rc("mathtext", fontset="cm")
params= {'text.latex.preamble' : [r'\usepackage{amsmath}']}
plt.rcParams.update(params)


def sortPltChk(inputString):                                                                                                                                                                                  
  stepKeyword = "step"            
  if 'plt' in inputString and 'chk' in inputString:
    print("\nError in get_box_lib sortPltChk function\n")
    return inputString            
  if 'plt' not in inputString and 'chk' not in inputString and 'step' not in inputString:
    if "Step" in inputString:     
      stepKeyword = "Step_"       
    else:                         
      print(f"\nError in get_box_lib sortPltChk function: {inputString}\n")
      return inputString          
                                  
  if stepKeyword in inputString:     
    return int(inputString.split(stepKeyword)[1].split("_level")[0])
  elif 'plt' in inputString: splitOn = 'plt'
  elif 'chk' in inputString: splitOn = 'chk'
  return int(inputString.split(splitOn)[1])
#==============================================================================
# 
#==============================================================================

# get a list of all the files in this directory
files = get_files('.', include=['.plt'], exclude=["temp"], times=[], tol=1e-4, get_all=True)
files = sorted(files, key=sortPltChk) 

f = files[-1]
rh5 = ReadBoxLib(f)
t = rh5.time
    
rh5 = ReadBoxLib(f, max_level=-1)

print(f"\n\n#######\nHard coded for t = {t}\n\tfile:\t{f}")
#print("\n\n#######\nGamma assume 5./3.\n\n")
gamma = 5./3.

x, rho_e, ratio = rh5.get('rho-electrons', get_refinement=True)
x, ux_e  = rh5.get("x_vel-electrons")
x, uy_e  = rh5.get("y_vel-electrons")
x, uz_e  = rh5.get("z_vel-electrons")
x, p_e   =  rh5.get("p-electrons")
mass_e = rh5.get("mass-electrons")[0][0];
charge_e = rh5.get("charge-electrons")[0][0];

T_e = p_e / (rho_e/mass_e)

eden_e = p_e/(gamma-1) + 0.5*rho_e*(ux_e*ux_e + uy_e*uy_e + uz_e*uz_e);

x, rho_i = rh5.get("rho-ions")
x, ux_i  = rh5.get("x_vel-ions")
x, uy_i  = rh5.get("y_vel-ions")
x, uz_i  = rh5.get("z_vel-ions")
x, p_i   = rh5.get("p-ions")
mass_i = rh5.get("mass-ions")[0][0];
charge_i = rh5.get("charge-ions")[0][0];
T_i = p_i / (rho_i/mass_i)
eden_i = p_i/(gamma-1) + 0.5*rho_i*(ux_i*ux_i + uy_i*uy_i + uz_i*uz_i);

x, Bz    = rh5.get("z_B-field")
x, By    = rh5.get("y_B-field")
#x, Bx    = rh5.get("x_B-field")
x, Dx    = rh5.get("x_D-field")

J_z = uz_e * rho_e/mass_e * charge_e + uz_i * rho_i/mass_i * charge_i 

"""
x, rho_mhd = rh5.get("rho-mhd")
x, ux_mhd  = rh5.get("x_vel-mhd")
x, uy_mhd  = rh5.get("y_vel-mhd")
x, uz_mhd  = rh5.get("z_vel-mhd")
x, p_mhd   =  rh5.get("p-mhd")
x, Bz_mhd    = rh5.get("z_B-mhd")
x, By_mhd    = rh5.get("y_B-mhd")
eden_mhd = p_mhd/(gamma-1) + 0.5*rho_mhd*(ux_mhd*ux_mhd + uy_mhd*uy_mhd + uz_mhd*uz_mhd);
"""
# =============================================================================
# 
# =============================================================================

rho_ref = "Loverich-1-rho.hdf5"
ux_ref  = "Loverich-1-ux.hdf5"
uy_ref  = "Loverich-1-uy.hdf5"
uz_ref  = "Loverich-1-uz.hdf5"
p_ref   = "Loverich-1-p.hdf5"
By_ref  = "Loverich-1-By.hdf5"
Bz_ref  = "Loverich-1-Bz.hdf5"

def ref_data(ax, fname, label=False):
    h5 = h5py.File(fname,'r')
    
    for item in h5.items():
        if "Line" in item[0]:
            X = item[1]['X'][()]
            Y = item[1]['Y'][()]
            
            if X.size > 8:
                if label:
                    set_label = label
                else:
                     set_label = None
                     
                ax.plot(X, Y,'ko', lw=0.3 ,ms=2.0, mew=0.3, mfc='none', markevery=20, label = set_label)
                # ax.plot(X, Y,'r--', lw=0.5, label = set_label)

def sim_data(ax, x, y, label=None, lnStyle=None):
    if lnStyle == None:
      lnStyle = 'k-'

    ax.plot(x, y, lnStyle, lw=0.75, label=label)

def mhd_data(ax, x, y, label=None):
    ax.plot(x, y,'k-.',lw=0.5, label=label)

# =============================================================================
# 
# =============================================================================

axes = []
fig = plt.figure(figsize=(16, 16))

nr = 3
nc = 2

x = x[0]

ax = fig.add_subplot(nr,nc,1); axes.append(ax)
sim_data(ax, x, ux_e, r"$u_{x,e}$", "-")
sim_data(ax, x, ux_i, r"$u_{x,i}$", "--")
sim_data(ax, x, uz_e, r"$u_{z,e}$", "-.")
sim_data(ax, x, uz_i, r"$u_{z,i}$", ":")
ax.set_ylabel(r"$u$")
ax.legend()

ax = fig.add_subplot(nr,nc,2); axes.append(ax)
sim_data(ax, x, uy_e, r"$u_{y,e}$","-")
sim_data(ax, x, uy_i, r"$u_{y,i}$", "--")
ax.set_ylabel(r"$u_{y,i}$")
ax.legend()

ax = fig.add_subplot(nr,nc,3); axes.append(ax)
sim_data(ax, x, J_z, r"$J_z$")
ax.set_ylabel(r"$J_z$")
ax.legend()

ax = fig.add_subplot(nr,nc,4); axes.append(ax)
sim_data(ax, x, Dx, r"$D_x$")
ax.set_ylabel(r"$D_x$")
ax.legend()

ax = fig.add_subplot(nr,nc,5); axes.append(ax)
sim_data(ax, x, By, r"$B_y$")
#mhd_data(ax, x, By_mhd)
ax.set_ylabel(r"$B_y$", fontsize=20)

ax = fig.add_subplot(nr,nc,6); axes.append(ax)
sim_data(ax, x, Bz, r"$B_z$")
#mhd_data(ax, x, Bz_mhd)
ax.set_ylabel(r"$B_z$")

for ax in axes:
    ax.set_xlim(-0.5, 0.5)  

for ax in axes[0:4]:
    ax.set_xticklabels([])

axes[-1].set_xlabel(r"$x$")

fig.tight_layout()

fig.savefig(f"plot_t_{t:.3}_fullBraginskiiRefactor.pdf", dpi=300)
plt.close(fig)
    
    
print("DONE")

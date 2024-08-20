#!python

import numpy as np
from netCDF4 import Dataset
import sys,os
from shutil import copyfile

"""

setup=NUUP_200W

batfile=tests/bathy_FCOO-$setup.nc
perifery=tests/maskreg_$setup.txt

# 2) Save land and sea as constants (0 or 1)

cdo setmisstoc,0. -setrtoc,0.,9999.,1. $batfile tests/tmp0.nc

# 3) Mask the area outside the region

cdo maskregion,$perifery $batfile tests/tmp.nc
cdo setmisstoc,0. -setrtoc,0.,9999.,3. tests/tmp.nc tests/tmp1.nc

# 4) Determine points on the the mask edge that are sea points
#
# Better do this in python

# 4a) Determine edge between sea points in tmp0.nc and mask of tmp1.nc

run seapoints_bdy.py tests/tmp0.nc tests/tmp1.nc mapsta_$setup.nc mapsta.$setup

sea_land_nc='tests/tmp0.nc'
mask_nc='tests/tmp1.nc'
dest_nc='tests/mapsta_SETUP.nc'
dest_mask='tests/mapsta.SETUP'

"""

sea_land_nc=sys.argv[1]
mask_nc=sys.argv[2]
dest_nc=sys.argv[3]
dest_mask=sys.argv[4]


sea_land_set = Dataset(sea_land_nc,'r')
mask_set = Dataset(mask_nc,'r')

if os.path.exists(dest_nc):
    answer=input("Do you want to overwrite existing file " + dest_nc + "?[y/n]")
    if not answer[0] == 'y':
        raise RuntimeError, "Move away " + dest_nc + ". Then try again."
copyfile(mask_nc, dest_nc)

dest_set = Dataset(dest_nc,'r+')

variables=mask_set.variables.keys()
 
for var in variables:
    if 'lon' in var:
        vlon = var
    elif 'lat' in var:
        vlat = var
    elif 'bathy' in var or 'epth' in var:
        vbathy = var

mask=mask_set[vbathy][:]
sea_masked=sea_land_set[vbathy][:]

masked_sea = (mask == 0) * (sea_masked == 1)

sea_masked[masked_sea] = 3.

# Is there a 'land' border all around?
ei=False # Assume not
if np.all(sea_masked[:,0]==0) and np.all(sea_masked[:,-1]==0) \
   and np.all(sea_masked[0,:]==0) and np.all(sea_masked[-1,:]==0):
  ei=True
  # Fill with undefined all around
  sea_masked[:,0]=3
  sea_masked[:,-1]=3
  sea_masked[0,:]=3
  sea_masked[-1,:]=3

# Edges
#North:
north = np.zeros_like(sea_masked,dtype=np.bool)
north[:-1,:] = sea_masked[1:,:] - sea_masked[:-1,:] == 2

#East:
east = np.zeros_like(sea_masked,dtype=np.bool)
east[:,:-1] = sea_masked[:,1:] - sea_masked[:,:-1] == 2

# South:
south = np.zeros_like(sea_masked,dtype=np.bool)
south[1:,:] = sea_masked[:-1,:] - sea_masked[1:,:] == 2

# West:
west = np.zeros_like(sea_masked,dtype=np.bool)
west[:,1:] = sea_masked[:,:-1] - sea_masked[:,1:] == 2

if not ei:
  # Where water right to the edge, set it as a boundary
  north[-1,:] = sea_masked[-1,:] == 1
  east[:,-1] = sea_masked[:,-1] == 1
  south[0,:] = sea_masked[0,:] == 1
  west[:,0] = sea_masked[:,0] == 1

bdy = north+south+east+west

sea_masked[bdy] = 2

# Neighbours to bdy that have land on all other sides? (or have mask?) 

land = sea_masked == 0
open_sea = sea_masked == 1
open_sea_i =  open_sea[1:-1,1:-1]
# examples: j,i=(133,63)
# j,i=(133,76)
# east of bdy:
# bdy[j,i-1] * sea[j,i] * land[j,i+1] * land[j-1,i] * land[j+1,i]
isol = bdy[1:-1,:-2] * open_sea_i \
    * land[1:-1,2:] * land[:-2,1:-1] * land[2:,1:-1]
# north of bdy:
isol += bdy[:-2,1:-1] * open_sea_i \
    * land[2:,1:-1] * land[1:-1,:-2] * land[1:-1,2:]
# west of bdy:
isol += bdy[1:-1,2:] * open_sea_i \
    * land[1:-1,:-2] * land[:-2,1:-1] * land[2:,1:-1]
# south of bdy:
isol += bdy[2:,1:-1] * open_sea_i \
    * land[:-2,1:-1] * land[1:-1,:-2] * land[1:-1,2:]

sea_masked[1:-1,1:-1][isol] = 3    

open_sea_i[isol] = False

# Boundary cell has no water on any side?

where_no_bdy = np.where( bdy[1:-1,1:-1] * ~open_sea[:-2,1:-1] \
               * ~open_sea[2:,1:-1] * ~open_sea[1:-1,:-2] * ~open_sea[1:-1,2:] )

sea_masked[1:-1,1:-1][where_no_bdy] = 3

# Was there a 'land' border all around?
# Then restore to land
if ei:
  sea_masked[:,0]=0
  sea_masked[:,-1]=0
  sea_masked[0,:]=0
  sea_masked[-1,:]=0


dest_set.renameVariable(vbathy,'MAPSTA')
dest_set.variables['MAPSTA'][:] = sea_masked

dest_set.close()

np.savetxt(dest_mask, sea_masked, fmt="%d")


#!python


import numpy as np
from netCDF4 import Dataset
import sys,os
from shutil import copyfile

"""

# In Bash:
# Generate a MAPSTA file for grid NUUP_600W without bdy conditions from NUUP_200
cd $BUILDDIR/inp/NUUP/setup
ln -sf  ww3_grid_NUUP_600W_no2wayin.nml  ww3_grid.nml
$BUILDDIR/bin/ww3_grid > ww3_grid_NUUP_600W_no2wayin.log
mv mapsta.ww3 mapsta.NUUP_600W_no2wayin
cp -a mapsta.NUUP_600W_no2wayin $WDIR/grid/region_perifery/tests

# In Python:

setup='NUUP_200W'
outbound='NUUP_600W'

mask_nc='tests/mapsta_NUUP_200W.nc'

outdepths_nc='tests/bathy_FCOO-NUUP_600W.nc'
dest_mapsta='tests/mapsta.FCOO-NUUP_600W'

out_mapsta='tests/mapsta.NUUP_600W_no2wayin'

# Separately run:
run outbound_bdy.py $mask_nc $out_mapsta $outdepths_nc $dest_mapsta

"""

try:
    mask_nc = sys.argv[1]
    out_mapsta = sys.argv[2]
    outdepths_nc = sys.argv[3]
    dest_mapsta = sys.argv[4]
except:
    setup='NUUP_200W'
    outbound ='NUUP_600W'

    mask_nc='tests/mapsta_NUUP_200W.nc'

    outdepths_nc='tests/bathy_FCOO-NUUP_600W.nc'
    dest_mapsta='tests/mask_FCOO-NUUP_600W'

    out_mapsta='tests/mapsta.NUUP_600W_no2wayin'
    

mapsta_set = Dataset(mask_nc,'r')

variables=mapsta_set.variables.keys()
for var in variables:
    if 'lon' in var:
        vlon = var
    elif 'lat' in var:
        vlat = var
    elif 'MAPSTA' in var or 'bathy' in var:
        vmapsta = var

# Inner grid mapsta:

m_i = mapsta_set[vmapsta][:]


outdepths = Dataset(outdepths_nc,'r')

variables=outdepths.variables.keys()
for var in variables:
    if 'lon' in var:
        olon = var
    elif 'lat' in var:
        olat = var
    elif 'depth' in var or 'bathy' in var:
        odepth = var
        
olons=outdepths.variables[olon]
olats=outdepths.variables[olat]

vlon0=mapsta_set.variables[vlon][0]

# Read mapsta and fill 'MAPSTA' with the values

read_mapsta = np.loadtxt(out_mapsta)

dest2_mapsta = read_mapsta[::50]

"""
Why only every 50'th value?
See https://github.com/NOAA-EMC/WW3/issues/391
This bug was introduced on Sep 29, 2017 commit b38fdec l 3879ff. The change was relative to Parent commit 0ac55e6 along with a total of "68 changed files with 4,322 additions and 1,864 deletions". Maybe part of a transition from svn to git?
"""
assert len(olats)*len(olons) == dest2_mapsta.shape[0]
shape_out=dest2_mapsta.shape[0]

m_out = dest2_mapsta.reshape(len(olats),len(olons))[::-1,:]
ilon0,ilon1=mapsta_set.variables[vlon][[0,-1]]
ilat0,ilat1=mapsta_set.variables[vlat][[0,-1]]
dlat=(ilat1-ilat0)/(len(mapsta_set.variables[vlat])-1)
dlon=(ilon1-ilon0)/(len(mapsta_set.variables[vlon])-1)
sy=ilat0-olats[0]
sx=ilon0-olons[0]
oy=np.round(sy/dlat).astype(int)
ox=np.round(sx/dlon).astype(int)
assert abs(sy/dlat-oy) < 0.00001
assert abs(sx/dlon-ox) < 0.00001

dolat=(olats[-1]-olats[0])/(len(olats)-1)
dolon=(olons[-1]-olons[0])/(len(olons)-1)
muly = np.round(dolat/dlat,decimals=5)
mulx = np.round(dolon/dlon,decimals=5)

assert muly-np.round(muly) == 0.
assert mulx-np.round(mulx) == 0.



"""
# Add MAPSTA = 2 from inner mapsta

# First, find inner grid points at >= 2-cell distance from bdy

->                    -1  -2  -3  -4  -5
->                     4   3   2   1   0
<-                     0   1   2   3   4

# inner grid:  3   3   2   1   1   1   1   1   0   0   l
#              i   i   b   w   w   w   w   w   l   l   l
#
# outer grid:      1           2           1           0
#                  w           b           w           l

"""

i_d =   (m_i[4:  ,2:-2] == 2) \
      * (m_i[1:-3,2:-2] == 1) * (m_i[2:-2,2:-2] == 1) * (m_i[3:-1,2:-2] == 1) \
      * (m_i[2:-2,1:-3] == 1) * (m_i[2:-2,3:-1] == 1)

i_u =   (m_i[ :-4,2:-2] == 2)  \
      * (m_i[3:-1,2:-2] == 1) * (m_i[2:-2,2:-2] == 1) * (m_i[1:-3,2:-2] == 1) \
      * (m_i[2:-2,1:-3] == 1) * (m_i[2:-2,3:-1] == 1)

i_l =   (m_i[2:-2,4:  ] == 2)  \
      * (m_i[2:-2,1:-3] == 1) * (m_i[2:-2,2:-2] == 1) * (m_i[2:-2,3:-1] == 1) \
      * (m_i[1:-3,2:-2] == 1) * (m_i[3:-1,2:-2] == 1)

i_r =   (m_i[2:-2, :-4] == 2)  \
      * (m_i[2:-2,3:-1] == 1) * (m_i[2:-2,2:-2] == 1) * (m_i[2:-2,1:-3] == 1) \
      * (m_i[1:-3,2:-2] == 1) * (m_i[3:-1,2:-2] == 1)



#          0         1      2         3       4         5        6        7
shs =   [ [1,1],    [1,0], [1,-1],   [0,-1], [-1,-1],  [-1, 0], [-1, 1], [0, 1]]
combs = ['i_u*i_r', 'i_u', 'i_u*i_l', 'i_l', 'i_d*i_l', 'i_d', 'i_d*i_r', 'i_r']

# ox,oy = [ 0, 0 ]

# mulx, muly = [ 3, 3 ]

m_o = np.empty_like(m_i)
m_o[:] = m_i[:]

i_edge = [i_u, i_l, i_d, i_r ]
she = [[ 1,0], [ 0,-1], [ -1, 0], [ 0, 1] ]

# merge the corners [i_r, i_u], [i_u, i_l], [i_l, i_d], [i_d, i_r]
# to extend i_edge -> combs, i_e -> i_c and she -> shs:

ii=0
i_c=[]
shs = []
for jj in range(len(i_edge)):
    i_m = i_edge[jj-1]
    i_e = i_edge[jj]
    i_me = i_m[:]*i_e[:]
    
    shs += [ np.add(she[jj-1],she[jj]), she[jj] ]
    i_c += [ i_me[:] ]
    
    # Here the two sides have been used to create the corner:
    i_m[i_me] = False
    i_e[i_me] = False
    i_c += [ i_e[:] ]

m_o00=m_o[0,0]
assert m_o00 != 1

oy0=oy/muly
ox0=ox/mulx

show_proc = True
if show_proc:
    v_u = len(i_c)+4
    for ii in range(len(i_c)):    
        m_o[2:-2,2:-2][i_c[ii]] = v_u+ii

for ii in range(len(i_c)):    
  
  wo0,wo1 = np.where(i_c[ii])
  wo0=(wo0+oy+2.)/muly; wo1=(wo1+ox+2.)/mulx
  
  if shs[ii][0] == 1:
      wo0 = np.ceil(wo0)
  elif shs[ii][0] == -1:
      wo0 = np.floor(wo0)
  else:
      nomatch=~(np.remainder(wo0,1.)==0)
      wo0[nomatch]=0
      wo1[nomatch]=0
  
  if shs[ii][1] == 1:
      wo1 = np.ceil(wo1)
  elif shs[ii][1] == -1:
      wo1 = np.floor(wo1)
  else:
      nomatch=~(np.remainder(wo1,1.)==0)
      wo0[nomatch]=0
      wo1[nomatch]=0

  wo0[wo0==0] = oy0
  wo1[wo1==0] = ox0

  # In devel, test points in inner grid
  wo0=(wo0*muly).astype(int)-oy; wo1=(wo1*mulx).astype(int)-ox;
 
  m_o[ tuple( [ wo0, wo1 ] ) ] = ii + 4

muly=int(muly)
mulx=int(mulx)
for ii in range(len(i_c)):    
  if shs[ii][0]*shs[ii][1] ==0: continue

  # Skip 'hidden corners'
  wt0,wt1 = np.where(m_o==ii + 4)
  for kk in range(len(wt0)):
      if m_o[ wt0[kk]-shs[ii][0]*muly, wt1[kk] ] > 3 \
         and m_o[ wt0[kk], wt1[kk]-shs[ii][1]*mulx ] > 3:
          print ii,kk,'m_o[ ',wt0[kk],',', wt1[kk],' ] = ii + 12'
          m_o[ wt0[kk], wt1[kk] ] = ii + 12
      
m_o[0,0] = m_o00

# Should any neighbor to a bdy point in m_o be a land point in the outer grid?
for ii in range(len(i_c)):
    wt0,wt1 = np.where(m_o==ii + 4)
    for kk in range(len(wt0)):
        # left-parallel neighbor
        spary = -shs[ii][1]*mulx; sparx = shs[ii][0]*muly
        for sgn in [1,-1]:
            # sgn == -1 shifts to right-parallel neighbor
            spary*=sgn; sparx*=sgn;
            tpy=wt0[kk]+spary; tpx=wt1[kk]+sparx
            # Not an already skipped point
            if m_o[tpy, tpx] > 1:  continue
            # Must have more than 2 land and some water around
            if (m_o[tpy-1:tpy+2, tpx-1:tpx+2] == 0).sum() < 3: continue
            if np.all(m_o[tpy-1:tpy+2, tpx-1:tpx+2] != 1): continue

            print ii,kk,'m_o[ ',tpy,',', tpx,' ] = 20'
            m_o[ tpy, tpx ] = 20


if show_proc:
    # In devel, view points in inner grid
    copyfile(mask_nc, 'tests/test_mask.nc')
    testbound = Dataset('tests/test_mask.nc','r+')

    testbound.variables[vmapsta][:] = m_o[:]
    testbound.close()

where_bdyo = np.where((m_o > 3) * (m_o < 12))

# Make the bdy points to be sea points
where_bdyl = np.where(m_i[where_bdyo] == 0)

# For test set o_repl = 4
o_repl=1
n_repl=1

#o_repl=4
#n_repl=5

# m_i[where_bdyo][where_bdyl[0]] = o_repl
# This was wrong, see
# https://scipy-cookbook.readthedocs.io/items/ViewsVsCopies.html
# Instead, setitems for all where_bdyo directly to ocean:

# Skip point if more that one neighbour is land
sea_added=False
for ibdyo in range(len(where_bdyo[0])):
    
    iy=where_bdyo[0][ibdyo]
    ix=where_bdyo[1][ibdyo]
    # Are all or some of their neighbours sea points?
    numsea =np.sum(m_i[iy-1:iy+2,ix])+np.sum(m_i[iy,ix-1:ix+2])
    if numsea == 6: continue
    if numsea > 3:
        sea_added=True
        m_o[iy,ix] = 1
        # Also make their neighbours sea points
        m_i[iy-1:iy+2,ix] = 1; m_i[iy,ix-1:ix+2] = 1
    else:
        m_o[iy,ix] = 0
        # Remove the index by replicating the previous
        where_bdyo[0][ibdyo] = where_bdyo[0][ibdyo-1]
        where_bdyo[1][ibdyo] = where_bdyo[1][ibdyo-1]


m_i[where_bdyo] = o_repl

if o_repl > 1:
    # Mark up their neighbours
    m_i[(where_bdyo[0]-1,where_bdyo[1])] = n_repl
    m_i[(where_bdyo[0]+1,where_bdyo[1])] = n_repl
    m_i[(where_bdyo[0],where_bdyo[1]-1)] = n_repl
    m_i[(where_bdyo[0],where_bdyo[1]+1)] = n_repl

    copyfile(mask_nc, 'tests/o_repl.nc')
    testbound = Dataset('tests/o_repl.nc','r+')
    testbound.variables[vmapsta][:] = m_i[:]
    testbound.close()
    m_i[m_i > 3] = 1

if sea_added:
    print 'Sea points added. Consider to exchange '+mask_nc+' with '+new_mask_nc
    copyfile(mask_nc, 'tests/new_inner_mapsta.nc')
    testbound = Dataset('tests/new_inner_mapsta.nc','r+')
    testbound.renameVariable(vmapsta,'MAPSTA')
    testbound.variables['MAPSTA'][:] = m_i[:]
    testbound.close()

wo0,wo1 = where_bdyo

assert np.all((wo0+oy)/muly == 1./muly*(wo0+oy))
assert np.all((wo1+ox)/mulx == 1./mulx*(wo1+ox))

# Boundary points in outer
bdt=2
if o_repl > 1: bdt = 5 # temporarily '5':
    
wo0=(wo0+oy)/muly; wo1=(wo1+ox)/mulx

m_out[ tuple([wo0,wo1]) ] = bdt


# Side to bdy point in inner grid must be land point in outer
# Added land points
bds=0
if o_repl > 1: bds = 6

wos0,wos1 = np.where(m_o==20)
wos0=(wos0+oy)/muly; wos1=(wos1+ox)/mulx

print "Side bdy points in outer: ", m_out[ tuple([wos0,wos1]) ]

m_out[ tuple([wos0,wos1]) ] = bds

if o_repl > 1:
    print "Test MAPSTA to " + dest_mapsta+'_test.nc, not to ' +dest_mapsta+'.nc' 
    dest_mapsta += '_test'
else:
    print 'New MAPSTA for outer grid in '+dest_mapsta+'and '+dest_mapsta+'.nc'

copyfile(outdepths_nc, dest_mapsta+'.nc')
outbound = Dataset(dest_mapsta+'.nc','r+')

outbound.renameVariable(odepth,'MAPSTA')
outbound.variables['MAPSTA'][:] = m_out

outbound.close()

if o_repl <= 1: np.savetxt(dest_mapsta, m_out, fmt="%d")

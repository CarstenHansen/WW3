#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Carsten Hansen"
__contact__ = "cha@fcoo.dk"
__copyright__ = "Copyright 2016, Danish Defence"
__license__ = "GNU General Public License v3.0"

from netCDF4 import Dataset
from optparse import OptionParser
from numpy import savetxt as np_savetxt
from itertools import izip

usage = """usage: %prog [--format='%7.1f'] [--selvar='bathymetry'] file.nc [file.dat | --output=file.dat]
or simply (if you know that for bathymetry the missing_value=-9999.f):
run nc2dat.py -s '%7.1f' -v bathymetry bathy.nc"""

parser = OptionParser(usage=usage)

parser.add_option("-s", "--format",
                  dest="format",metavar="PRINTFORMAT",
                  action="store",type="string",default="%7.1f",
                  help="Print format. Default '%7.1f'")
parser.add_option("-o", "--output",
                  dest="file_out",metavar="OUTPUT_FILE",
                  action="store",type="string",default=".dat",
                  help="Print format. Default prefix(ncfile)+'.dat'")
parser.add_option("-v", "--selvar",
                  dest="selvar",metavar="VARIABLE",
                  action="store",type="string",default=None,
                  help="Select variable")
parser.add_option("-m", "--scale",
                  dest="var_scale",metavar="VARIABLE",
                  action="store",type="float",default=1.,
                  help="Variable multiplier")

(options, args) = parser.parse_args()

print 'options =',options

format = options.format
selvar = options.selvar
var_scale= options.var_scale


if len(args)>=1:
    file_in = args[0]
else:
    file_in = 'bathy.nc'

if len(args)>=2:
    file_out = args[1]
else:
    file_out = options.file_out
    if file_out[0]=='.':
        filesplit=file_in.split('.')
        try:
            ncindx = filesplit.index('nc')
            filesplit[ncindx] = file_out[1:]
        except:
            file_out = file_in+'.'+file_out[1:]
        else:
            file_out = '.'.join(filesplit)

print file_in, '->', file_out

fm = Dataset(file_in,'r')

for dim in fm.dimensions.keys():
    if dim.lower()[1:3] == 'on':
        lon=dim
    elif dim.lower()[1:3] == 'at':
        lat=dim

varname = None
for key in fm.variables.keys():
    if key==lon or key==lat or key=='time':
        continue
    else:
        if varname is not None and key != selvar:
            raise RuntimeError, 'More than one variable in file. You must specify --selvar=<variable>'
        if selvar is not None and key != selvar:
            continue
        varname = key
        if key == selvar:
            break
try:
    var=fm.variables[varname][:]
except:
    raise RuntimeError, 'Cannot read variable "'+varname+'"'
try:
    var_units=fm.variables[varname].units
except:
    print 'WARNING: No units, assume dimensionless'
    var_units='1'

# Programmer's presumptions:
latlen=len(fm.variables[lat])
lonlen=len(fm.variables[lon])
var_shape=var.shape
if len(var_shape) > 2:
    var_shape=(var_shape[-2],var_shape[-1])
assert var_shape == (latlen, lonlen)


vals_perline = 78/( len(format%0.0) + 1 ) # Presume format is a float '%NN.Nf'

lon_linerange=range(0,lonlen,vals_perline)
lon_endrange=range(vals_perline,lonlen+vals_perline,vals_perline)

if var_units == '%' or var_units == '0.1':
    print "Rescale from units of '"+var_units+"' by multiplier", 
    if var_units == '%': var_scale*=0.01
    if var_units == '0.1': var_scale*=0.1
    print var_scale

fo = open(file_out,'w')
for lati in range(latlen):
    varlati=var[lati:lati+1,:]*var_scale
    for loni,lonj in izip(lon_linerange, lon_endrange):
        fo.write(' ') # At beginning of each line
        np_savetxt(fo,varlati[0:1,loni:lonj], fmt=format)
fo.close()

lonse = fm.variables[lon][[0,(lonlen-1)]]
latse = fm.variables[lat][[0,(latlen-1)]]
delta_lon = (lonse[1]-lonse[0])/(lonlen-1)
delta_lat = (latse[1]-latse[0])/(latlen-1)

print ''
print ' For usage in e.g. ww3 grid setup:'
print ' Lengths of lon and lat dimensions: %5d %12d' %(lonlen, latlen)
print ' Delta_lon, delta_lat:  %12.8f %12.8f' %(delta_lon, delta_lat)
print ' Origo_lon, origo_lat: %10.5f %12.5f' %(lonse[0], latse[0])
print ' Name:', file_out

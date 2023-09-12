#!/bin/env bash
set -o errexit
ext=''
[[ "$1" ]] && ext="_$1" && wgrid=$1
[[ "$ext" ]] && ln -sf wind$ext.nc wind.nc && ln -sf mod_def.$wgrid mod_def.ww3
ln -sf ww3_prnc_wind$ext.inp ww3_prnc.inp; ./ww3_prnc; rm ww3_prnc.inp
[[ "$wgrid" ]] && mv wind.ww3 wind.$wgrid
if [[ "$2" ]] && [[ -f ww3_prnc_icec$ext.inp ]]; then
    [[ "$ext" ]] && ln -sf icec$ext.nc icec.nc
    ln -sf ww3_prnc_icec$ext.inp ww3_prnc.inp; ./ww3_prnc; rm ww3_prnc.inp
    [[ "$wgrid" ]] && mv ice.ww3 ice.$wgrid
fi
if [[ "$2" ]] && [[ -f ww3_prnc_icet$ext.inp ]]; then
    # icec.nc also contains ice thickness
    ln -sf ww3_prnc_icet$ext.inp ww3_prnc.inp; ./ww3_prnc; rm ww3_prnc.inp
    [[ "$wgrid" ]] && mv ice1.ww3 ice1.$wgrid
fi
[[ "$ext" ]] && rm wind.nc && rm mod_def.ww3
[[ "$2" ]] && [[ "$ext" ]] && [[ -f ww3_prnc_icec$ext.inp ]] && rm -f icec.nc


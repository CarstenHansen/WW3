#!/bin/env bash
set -o errexit
ln -sf ww3_prnc_wind.inp ww3_prnc.inp; ./ww3_prnc; rm ww3_prnc.inp
ln -sf ww3_prnc_icec.inp ww3_prnc.inp; ./ww3_prnc; rm ww3_prnc.inp
ln -sf ww3_prnc_icet.inp ww3_prnc.inp; ./ww3_prnc; rm ww3_prnc.inp


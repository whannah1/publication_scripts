#!/usr/bin/env python3
import os, glob, subprocess as sp

'''
SCRATCH=/global/cfs/cdirs/m3312/whannah/2023-CPL

# Generate ERA5 grid file
LATDIR=n2s; ncremap -G ttl='Equi-Angular grid 721x1440'#latlon=721,1440#lat_typ=uni#lat_drc=${LATDIR}#lon_typ=grn_ctr  -g ${SCRATCH}/scrip_ERA5_721x1440_${LATDIR}.nc
LATDIR=s2n; ncremap -G ttl='Equi-Angular grid 721x1440'#latlon=721,1440#lat_typ=uni#lat_drc=${LATDIR}#lon_typ=grn_ctr  -g ${SCRATCH}/scrip_ERA5_721x1440_${LATDIR}.nc

# Generate map file for ERA5 => ne30pg2
LATDIR=n2s; ncremap --alg_typ=tempest -a fv2fv --src_grd=${SCRATCH}/scrip_ERA5_721x1440_${LATDIR}.nc --dst_grd=${SCRATCH}/scrip_ne30pg2.nc --map_file=${SCRATCH}/map_721x1440_${LATDIR}_to_ne30pg2.nc
LATDIR=s2n; ncremap --alg_typ=tempest -a fv2fv --src_grd=${SCRATCH}/scrip_ERA5_721x1440_${LATDIR}.nc --dst_grd=${SCRATCH}/scrip_ne30pg2.nc --map_file=${SCRATCH}/map_721x1440_${LATDIR}_to_ne30pg2.nc

'''
map_file = '/global/cfs/cdirs/m3312/whannah/2023-CPL/map_721x1440_s2n_to_ne30pg2.nc'
src_file = '/global/cfs/cdirs/e3sm/diagnostics/observations/Atm/climatology_1985-2014/ERA5/ERA5_ANN_198501_201412_climo.nc'
dst_file = '/global/cfs/cdirs/m3312/whannah/2023-CPL/ERA5/ERA5_ANN_198501_201412_climo.ne30pg2.nc'

cmd = f'ncremap -m {map_file} -i {src_file} -o {dst_file}  '
print(cmd)
sp.check_output(cmd,shell=True,universal_newlines=True)

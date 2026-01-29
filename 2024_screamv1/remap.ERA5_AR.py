
#-------------------------------------------------------------------------------
'''
SCRATCH_PATH=/pscratch/sd/w/whannah/e3sm_scratch
SRC_NY=721; SRC_NX=1440; SRC_GRID=${SRC_NY}x${SRC_NX}; SRC_LATDIR=n2s; SRC_GRID_FILE=${SCRATCH_PATH}/files_grid/${SRC_GRID}_${SRC_LATDIR}_scrip.nc
DST_NY=256; DST_NX=512 ; DST_GRID=${DST_NY}x${DST_NX}; DST_LATDIR=s2n; DST_GRID_FILE=${SCRATCH_PATH}/files_grid/${DST_GRID}_${DST_LATDIR}_scrip.nc

MAP_FILE=${SCRATCH_PATH}/files_map/map_${SRC_GRID}_${SRC_LATDIR}_to_${DST_GRID}_${DST_LATDIR}_aave.nc

ncremap -G ttl='Equi-Angular grid ${SRC_GRID}'#latlon=${SRC_NY},${SRC_NX}#lat_typ=uni#lat_drc=${SRC_LATDIR}#lon_typ=grn_ctr -g ${SRC_GRID_FILE}
ncremap -G ttl='Equi-Angular grid ${DST_GRID}'#latlon=${DST_NY},${DST_NX}#lat_typ=uni#lat_drc=${DST_LATDIR}#lon_typ=grn_ctr -g ${DST_GRID_FILE}

ncremap -6 --alg_typ=aave --src_grd=${SRC_GRID_FILE} --dst_grd=${DST_GRID_FILE} --map=${MAP_FILE}

echo ${MAP_FILE}
'''
#-------------------------------------------------------------------------------
import os
class clr:END,RED,GREEN,MAGENTA,CYAN = '\033[0m','\033[31m','\033[32m','\033[35m','\033[36m'

map_file = '/pscratch/sd/w/whannah/e3sm_scratch/files_map/map_721x1440_n2s_to_256x512_s2n_aave.nc'

src_path = '/pscratch/sd/p/paullric/SCREAM_40day'
dst_path = os.getenv('HOME')+'/Research/E3SM/pub_figs/screamv1/data'

src_file_list = []
src_file_list.append(f'{src_path}/apr/e5.oper.an.ar_climo.2013040300_2013050923.nc')
src_file_list.append(f'{src_path}/dy1/e5.oper.an.ar_climo.2016080300_2016090923.nc')
src_file_list.append(f'{src_path}/dy2/e5.oper.an.ar_climo.2020012200_2020022823.nc')
src_file_list.append(f'{src_path}/oct/e5.oper.an.ar_climo.2013100300_2013110923.nc')

for src_file in src_file_list:

  # src_file_stripped = src_file.split('/')[-1]
  dst_file = src_file.split('/')[-1].replace('.nc','.remap.nc')
  dst_file = f'{dst_path}/{dst_file}'

  cmd = f'ncremap -m {map_file} -i {src_file} -o {dst_file}'

  msg = clr.GREEN + cmd + clr.END
  print(f'\n{msg}')
  os.system(cmd)

#-------------------------------------------------------------------------------
  
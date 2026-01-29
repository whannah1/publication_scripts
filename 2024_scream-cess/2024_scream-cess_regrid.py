#!/usr/bin/env python
#---------------------------------------------------------------------------------------------------
'''
SRC_GRID=/pscratch/sd/w/whannah/e3sm_scratch/files_grid/ne1024pg2_scrip.nc
DST_GRID=/pscratch/sd/w/whannah/e3sm_scratch/files_grid/ne30pg2_scrip.nc
MAP_FILE=/pscratch/sd/w/whannah/e3sm_scratch/files_map/map_ne1024pg2_to_ne30pg2_traave.nc
ncremap --alg_typ=traave --grd_src=${SRC_GRID} --grd_dst=${DST_GRID} --map=${MAP_FILE}

SRC_GRID=/pscratch/sd/w/whannah/e3sm_scratch/files_grid/ne256pg2_scrip.nc
DST_GRID=/pscratch/sd/w/whannah/e3sm_scratch/files_grid/ne30pg2_scrip.nc
MAP_FILE=/pscratch/sd/w/whannah/e3sm_scratch/files_map/map_ne256pg2_to_ne30pg2_traave.nc
ncremap --alg_typ=traave --grd_src=${SRC_GRID} --grd_dst=${DST_GRID} --map=${MAP_FILE}
'''
#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------
import sys,os,glob,subprocess as sp
do_ha, do_h0, do_h1, do_h2, do_elm, overwrite = False, False, False, False, False, False

src_root = '/global/cfs/cdirs/e3smdata/simulations/Cess'
dst_root = '/global/cfs/cdirs/m3312/whannah/2024-SCREAM-CESS'

case_list = []

case_list.append('cess-control.ne1024pg2_ne1024pg2.F2010-SCREAMv1.cess-oct2')
case_list.append('cess-plus4k.ne1024pg2_ne1024pg2.F2010-SCREAMv1.cess-oct2')
# case_list.append('cess-control.ne256pg2_ne256pg2.F2010-SCREAMv1.cess-oct2')
# case_list.append('cess-plus4k.ne256pg2_ne256pg2.F2010-SCREAMv1.cess-oct2')


htype = 'output.scream.Cess.monthly_ne1024.AVERAGE.nmonths_x1'

# src_sub = 'archive/atm/hist'
src_sub = 'run'
dst_sub = f'data_remap_ne30pg2'

log_file = './regrid.case.out'

class tcolor: ENDC,RED,GREEN,MAGENTA,CYAN = '\033[0m','\033[31m','\033[32m','\033[35m','\033[36m'

#---------------------------------------------------------------------------------------------------
print()
for case in case_list:
    #---------------------------------------------------------------------------
    if 'ne1024pg2' in case: map_file = f'/pscratch/sd/w/whannah/e3sm_scratch/files_map/map_ne1024pg2_to_ne30pg2_traave.nc'
    if  'ne256pg2' in case: map_file = f'/pscratch/sd/w/whannah/e3sm_scratch/files_map/map_ne256pg2_to_ne30pg2_traave.nc'
    #---------------------------------------------------------------------------
    src_dir = f'{src_root}/{case}/{src_sub}'     # Input directory
    dst_dir = f'{dst_root}/{case}/{dst_sub}'     # Output directory
    #---------------------------------------------------------------------------
    print(f'  case     : {case}')
    print(f'  src_dir  : {src_dir}')
    print(f'  dst_dir  : {dst_dir}')
    print(f'  map_file : {map_file}')
    print('')
    #---------------------------------------------------------------------------
    # if not os.path.exists(dst_dir): os.mkdir(dst_dir)
    if not os.path.exists(dst_dir): os.makedirs(dst_dir, exist_ok=True)
    #---------------------------------------------------------------------------
    file_list = sorted( glob.glob( f'{src_dir}/{htype}*nc') )
    
    for src_file in file_list : 

        dst_file = src_file.replace(src_dir,dst_dir)
        dst_file = dst_file.replace('.nc','.remap_ne30pg2.nc')

        cmd = f'ncremap -P eamxx -m {map_file} -i {src_file} -o {dst_file}'

        print(f'\n{tcolor.GREEN}{cmd}{tcolor.ENDC}\n')
        # sp.check_output(cmd,shell=True,universal_newlines=True)
        os.system(cmd)
#---------------------------------------------------------------------------------------------------
print()
print('done.')
print()
#===============================================================================
#===============================================================================

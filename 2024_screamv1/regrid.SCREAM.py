#!/usr/bin/env python
#---------------------------------------------------------------------------------------------------
### map file for CMIP remap method:
# GRIDIN=ne30pg2; GRIDOUT=180x360 ; ncremap -6 --alg_typ=aave --grd_src=grid_files/${GRIDIN}_scrip.nc --grd_dst=grid_files/cmip6_${GRIDOUT}_scrip.20181001.nc --map=map_files/map_${GRID}_to_cmip6_${GRIDOUT}_aave.nc
# GRIDIN=ne30pg2; GRIDOUT=90x180 ; ncremap -6 --alg_typ=aave --grd_src=grid_files/${GRIDIN}_scrip.nc --grd_dst=grid_files/cmip6_${GRIDOUT}_scrip.20181001.nc --map=map_files/map_${GRID}_to_cmip6_${GRIDOUT}_aave.nc

#---------------------------------------------------------------------------------------------------
# command to create CMIP6 HighResMIP vertical pressure coordinate file
# ncap2 -O -v -s 'defdim("plev",27);plev[$plev]={100000,97500,95000,92500,90000,87500,85000,82500,80000,77500,75000,70000,65000,60000,55000,50000,45000,40000,35000,30000,25000,22500,20000,17500,15000,12500,10000};' ~/E3SM/vert_grid_files/vrt_prs_CMIP6.nc

# ERA5 vertical pressure coordinate
# ncap2 -O -v -s 'defdim("plev",37);plev[$plev]={100000,97500,95000,92500,90000,87500,85000,82500,80000,77500,75000,70000,65000,60000,55000,50000,45000,40000,35000,30000,25000,22500,20000,17500,15000,12500,10000,7000,5000,3000,2000,1000,700,500,300,200,100};' grid_files/vrt_prs_ERA5.nc
#---------------------------------------------------------------------------------------------------
# time python -u code_data/regrid.E3SM.py > regrid.L72.out &
#---------------------------------------------------------------------------------------------------
import sys,os,subprocess as sp
from pathlib import Path

src_root = '/global/cfs/cdirs/e3sm/gsharing/EAMxx'
dst_root = '/pscratch/sd/w/whannah/scream_scratch'
atm_comp = 'eam'


cases = []
cases.append('DYAMOND2_SCREAMv1')
cases.append('Apr1_2013_SCREAMv1')
cases.append('DYAMOND1_SCREAMv1')
cases.append('Oct1_2013_SCREAMv1')

# omit flags from case list
for c in cases: 
    if c[0]=='-': cases.remove(c)


file_type = 'output.scream.SurfVars.INSTANT'

execute   = True
overwrite = True
print_cmd = True

dst_grid = 'ne120pg2'
map_file = '/pscratch/sd/w/whannah/e3sm_scratch/files_map/map_ne1024pg2_to_ne120pg2_nco_20221130.nc'

vert_pressure_remap = False
pressure_level_file = 'grid_files/vrt_prs_ERA5.nc'

# log_file = './regrid.case.out'

class tcolor: ENDC,RED,GREEN,MAGENTA,CYAN = '\033[0m','\033[31m','\033[32m','\033[35m','\033[36m'

#===============================================================================
#===============================================================================
for case in cases :
    #---------------------------------------------------------------------------
    src_path = f'{src_root}/{case}' # Input directory
    dst_path = f'{dst_root}/{case}/{dst_grid}'# Output directory
    # if not os.path.exists(dst_path): path.mkdir(dst_path,parents=True)
    path = Path(dst_path)
    path.mkdir(parents=True, exist_ok=True)
    #---------------------------------------------------------------------------
    print(f'case     : {case}')
    print(f'src path : {src_path}')
    print(f'dst path : {dst_path}')
    print(f'map file : {map_file}')
    print('')
    #---------------------------------------------------------------------------
    files = os.listdir(src_path)
    cnt = 0 
    for src_file in files : 
        remap_flag = False
        if file_type in src_file : remap_flag = True
        if remap_flag :
            cnt += 1
            dst_file = src_file.replace('.nc',f'.remap_{dst_grid}.nc')
            #-------------------------------------------------------------------
            src_file_path = f'{src_path}/{src_file}'
            dst_file_path = f'{dst_path}/{dst_file}'
            #-------------------------------------------------------------------
            # print(src_file_path)
            # print(dst_file_path)
            # exit()
            #-------------------------------------------------------------------
            if os.path.isfile(dst_file_path):
                if overwrite: os.remove(dst_file_path)
                else: continue
            #-------------------------------------------------------------------
            cmd = f'ncremap -m {map_file} -i {src_file_path} -o {dst_file_path}'
            if file_type=='output.scream.SurfVars.INSTANT':
                # cmd += ' --xcl --vars=horiz_winds@bot,surf_mom_flux'
                cmd += ' -v T_2m,qv_2m,wind_speed_10m,ps,precip_liq_surf_mass,precip_ice_surf_mass'
            if vert_pressure_remap: cmd += f' --vrt_fl={pressure_level_file}'
            if print_cmd: print(f'\n{tcolor.GREEN}{cmd}{tcolor.ENDC}\n')
            #-------------------------------------------------------------------
            if execute: 
                try:
                    sp.check_output(cmd,shell=True,universal_newlines=True)
                except sp.CalledProcessError as error:
                    print(error.output)
                    exit()
#===============================================================================
#===============================================================================
    if cnt==0:
        print('\nNo files found for remapping...?\n')
    else:
        print(); print('done.'); print()
#===============================================================================
#===============================================================================

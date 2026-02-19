#!/usr/bin/env python
#===============================================================================
'''
commands to create grid and map files:

NE=30
SRC_GRID=ne${NE}pg2
DST_NY=90
DST_NX=180
DST_GRID=${DST_NY}x${DST_NX}

TIMESTAMP=20230222
SRC_GRID_FILE=${SRC_GRID}_scrip.${TIMESTAMP}.nc
DST_GRID_FILE=${DST_GRID}_scrip.${TIMESTAMP}.nc
MAP_FILE=map_${SRC_GRID}_to_${DST_GRID}_aave.${TIMESTAMP}.nc

or without time stamp
SRC_GRID_FILE=~/E3SM/data_grid/${SRC_GRID}_scrip.nc
DST_GRID_FILE=~/E3SM/data_grid/${DST_GRID}_scrip.nc
MAP_FILE=~/maps/map_${SRC_GRID}_to_${DST_GRID}_aave.nc

# generate model grid file
GenerateCSMesh --alt --res ${NE} --file ne${NE}.g
GenerateVolumetricMesh --in ne${NE}.g --out ne${NE}pg2.g --np 2 --uniform
ConvertMeshToSCRIP --in ne${NE}pg2.g --out ne${NE}pg2_scrip.nc
ConvertExodusToSCRIP --in ne${NE}pg2.g --out ne${NE}pg2_scrip.nc

# generate lat/lon grid file
ncremap -g ${DST_GRID_FILE} -G ttl="Equi-Angular grid, dimensions ${DST_GRID}, cell edges on Poles/Equator and Prime Meridian/Date Line"#latlon=${DST_NY},${DST_NX}#lat_typ=uni#lon_typ=grn_wst

# generate map file
ncremap -6 --alg_typ=aave --grd_src=$SRC_GRID_FILE --grd_dst=$DST_GRID_FILE --map=$MAP_FILE
'''
#===============================================================================
#===============================================================================
import sys,os,glob
import subprocess as sp
#===============================================================================
#-------------------------------------------------------------------------------
case_list = []
src_top_dir_list,src_sub_dir_list = [],[]
dst_top_dir_list,dst_sub_dir_list = [],[]
def add_case(case_in,src_dir_in,src_sub_in,dst_dir_in,dst_sub_in):
   global case_list,src_top_dir_list,src_sub_dir_list
   case_list.append(case_in)
   src_top_dir_list.append(src_dir_in); src_sub_dir_list.append(src_sub_in)
   dst_top_dir_list.append(dst_dir_in); dst_sub_dir_list.append(dst_sub_in)
#-------------------------------------------------------------------------------
#===============================================================================
do_atm_h0, do_atm_h1, do_atm_h2, do_lnd_h0, overwrite = False, False, False, False, False

# nlat_dst,nlon_dst,map_file =  90,180,os.getenv('HOME')+'/maps/map_ne30pg2_to_90x180_aave.nc'
nlat_dst,nlon_dst,map_file = 73,144,'/global/homes/w/whannah/maps/map_ne30pg2_to_73x144_traave.nc'


tmp_path,tmp_sub = '/global/cfs/cdirs/m3312/whannah/e3smv3_amip','archive/atm/hist'

# add_case('v3.LR.amip_0101',tmp_path,tmp_sub,tmp_path,f'data_remap_{nlat_dst}x{nlon_dst}')
# add_case('v3.LR.amip_0151',tmp_path,tmp_sub,tmp_path,f'data_remap_{nlat_dst}x{nlon_dst}')
# add_case('v3.LR.amip_0201',tmp_path,tmp_sub,tmp_path,f'data_remap_{nlat_dst}x{nlon_dst}')

src_root,src_sub = '/global/cfs/cdirs/m4310/data/sims','archive/atm/hist'
dst_root = '/pscratch/sd/w/whannah/e3sm_scratch/pm-cpu'
add_case('v3.LR.amip_0101.QBObenchmark.20241008',src_root,src_sub,dst_root,f'data_remap_{nlat_dst}x{nlon_dst}')

htype = 'eam.h1'; var_list_arg = '-v area,PRECT,FLUT,U850'
# htype = 'eam.h3'; var_list_arg = '-v U850'

execute   = True
overwrite = False
print_cmd = True

# yr1,yr2 = 1994,2004
yr1,yr2 = 2005,2005

#===============================================================================
class tcolor: ENDC,RED,GREEN,MAGENTA,CYAN = '\033[0m','\033[31m','\033[32m','\033[35m','\033[36m'
def run_cmd(cmd,print_cmd=True,suppress_output=False):
    if suppress_output : cmd = cmd + ' > /dev/null'
    if print_cmd: print(tcolor.GREEN + cmd + tcolor.ENDC)
    if execute:
        try:
            sp.check_output(cmd,shell=True,universal_newlines=True)
        except sp.CalledProcessError as error:
            print(error.output)
            exit()
    return
#===============================================================================
def check_path_and_create(path_in):
    if not os.path.exists(path_in): os.mkdir(path_in)
    return       
#===============================================================================
for c,case in enumerate(case_list) :
    src_top_dir, src_sub_dir = src_top_dir_list[c], src_sub_dir_list[c]
    dst_top_dir, dst_sub_dir = dst_top_dir_list[c], dst_sub_dir_list[c]
    #---------------------------------------------------------------------------
    # Set in/out paths and create output directory if it doesn't exist
    src_dir = f'{src_top_dir}/{case}/{src_sub_dir}'
    dst_dir = f'{dst_top_dir}/{case}/{dst_sub_dir}'
    if not os.path.exists(dst_dir):
        print(); print(f'Creating output directory: {dst_dir}'); print()
        check_path_and_create(f'{dst_top_dir}')
        check_path_and_create(f'{dst_top_dir}/{case}')
        check_path_and_create(f'{dst_top_dir}/{case}/{dst_sub_dir}')
    #---------------------------------------------------------------------------
    print(f'case    : {case}')
    print(f'src_dir : {src_dir}')
    print(f'dst_dir : {dst_dir}')
    print('')
    atm_comp,lnd_comp = 'eam','elm'
    #---------------------------------------------------------------------------
    # get list of all files to loop over
    # files = sorted( os.listdir(src_dir) )
    #---------------------------------------------------------------------------
    # if 'yr1' in locals() \
    # or 'yr2' in locals():
    #     yr_list = []
    #     for f,f_in in enumerate(files):
    #         yr = f_in
    #         # for htype in ['h0','h1','h2','h3','h4','h5']:
    #             # yr = yr.replace(f'{case}.{atm_comp}.{htype}.','')
    #             # yr = yr.replace(f'{case}.{lnd_comp}.{htype}.','')
    #         yr = yr.replace(f'.nc','')
    #         yr = int(yr.split('-')[0])
    #         yr_list.append(yr)
    #---------------------------------------------------------------------------
    cnt = 0
    file_list = sorted( glob.glob(f'{src_dir}/{case}.{htype}.*') )
    for f,f_in in enumerate(file_list):
        remap_flag = True

        # don't regrid data that's already been regridded
        if '.remap_' in f_in : remap_flag = False

        yr = int( f_in.replace(f'{src_dir}/{case}.{htype}.','').replace(f'.nc','').split('-')[0] )

        if 'yr1' in locals() and yr<yr1: remap_flag = False
        if 'yr2' in locals() and yr>yr2: remap_flag = False

        if remap_flag :
            cnt += 1

            src_file_name = f_in #f'{src_dir}/{f_in}'
            dst_file_name = src_file_name.replace('.nc',f'.remap_{nlat_dst}x{nlon_dst}.nc')

            if os.path.isfile(dst_file_name) :
                if overwrite : os.remove(dst_file_name)
                else : continue

            # Change path of output file
            dst_file_name = dst_file_name.replace(src_dir,dst_dir)

            cmd = f'ncremap -m {map_file} -i {src_file_name} -o {dst_file_name} {var_list_arg}'

            run_cmd(cmd)
#===============================================================================
    if cnt==0:
        print(files)
        print()
        print(f'  dst_grid : {nlat_dst} x {nlon_dst}')
        print(f'  atm_comp : {atm_comp}')
        print('\nNo files found for remapping...?\n')
#===============================================================================
print('Done.')


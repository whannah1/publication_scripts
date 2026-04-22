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

nlat_dst,nlon_dst,map_file =  90,180,os.getenv('HOME')+'/maps/map_ne30pg2_to_90x180_aave.nc'
# nlat_dst,nlon_dst,map_file =  90,180,os.getenv('HOME')+'/maps/map_ne30pg2_to_cmip6_90x180_aave.nc'

# cases = []
# cases.append('20230629.v3alpha02.amip.chrysalis.L72')
# cases.append('20230629.v3alpha02.amip.chrysalis.L80')
# cases.append('E3SM.INCITE2022-LOW-CLD-CESS.ne30pg2.F2010-MMF1.SSTP_0K')
# cases.append('E3SM.INCITE2022-CPL.ne30pg2_EC30to60E2r2.WCYCL20TR')


# tmp_path,tmp_sub = '/global/cfs/cdirs/m3312/whannah/2023-CPL/','archive/atm/hist'
# add_case('E3SM.INCITE2023-CPL.ne30pg2_EC30to60E2r2.WCYCL20TR-MMF1',tmp_path,tmp_sub,tmp_path,f'data_remap_{nlat_dst}x{nlon_dst}')

tmp_path,tmp_sub = '/global/cfs/cdirs/m3312/whannah/e3smv2_historical','archive/atm/hist'
add_case('v2.LR.historical_0101',tmp_path,tmp_sub,tmp_path,f'data_remap_{nlat_dst}x{nlon_dst}')

# tmp_path,tmp_sub = '/global/cfs/cdirs/m3312/whannah/2023-MJO','archive/atm/hist'
# add_case('E3SM.2023-MJO-test.GNUGPU.ne30pg2_EC30to60E2r2.F2010-MMF1.L60.NX_32.RX_4.DX_4000',tmp_path,tmp_sub,tmp_path,f'data_remap_{nlat_dst}x{nlon_dst}')
# add_case('E3SM.2023-MJO-test.GNUGPU.ne30pg2_EC30to60E2r2.F2010-MMF1.L60.NX_64.RX_4.DX_2000',tmp_path,tmp_sub,tmp_path,f'data_remap_{nlat_dst}x{nlon_dst}')
# add_case('E3SM.2023-MJO-test.GNUGPU.ne30pg2_EC30to60E2r2.F2010-MMF1.L72.NX_32.RX_4.DX_4000',tmp_path,tmp_sub,tmp_path,f'data_remap_{nlat_dst}x{nlon_dst}')
# add_case('E3SM.2023-MJO-test.GNUGPU.ne30pg2_EC30to60E2r2.F2010-MMF1.L72.NX_64.RX_4.DX_2000',tmp_path,tmp_sub,tmp_path,f'data_remap_{nlat_dst}x{nlon_dst}')


# tmp_path,tmp_sub = '/pscratch/sd/w/whannah/e3sm_scratch/pm-gpu','run'
# /global/cfs/cdirs/m4310/whannah/E3SM/
# add_case('E3SM.2023-SCIDAC-MMF.ne30pg2_EC30to60E2r2.F2010-MMF1.L60',tmp_path,tmp_sub,tmp_path,f'data_remap_{nlat_dst}x{nlon_dst}')
# add_case('E3SM.2023-SCIDAC-MMF.ne30pg2_EC30to60E2r2.F2010-MMF1.L64',tmp_path,tmp_sub,tmp_path,f'data_remap_{nlat_dst}x{nlon_dst}')
# add_case('E3SM.2023-SCIDAC-MMF.ne30pg2_EC30to60E2r2.F2010-MMF1.L72',tmp_path,tmp_sub,tmp_path,f'data_remap_{nlat_dst}x{nlon_dst}')

### 4xCO2 tests on Summit
# tmp_path,tmp_sub = '/gpfs/alpine/cli115/proj-shared/hannah6/e3sm_scratch','run'
# add_case('E3SM.2023-CO2-TEST-00.GNUGPU.ne30pg2_oECv3.WCYCL1850-MMF1.1xCO2',tmp_path,tmp_sub,tmp_path,f'data_remap_{nlat_dst}x{nlon_dst}')
# add_case('E3SM.2023-CO2-TEST-00.GNUGPU.ne30pg2_oECv3.WCYCL1850-MMF1.4xCO2',tmp_path,tmp_sub,tmp_path,f'data_remap_{nlat_dst}x{nlon_dst}')
# add_case('E3SM.INCITE2022-CPL.ne30pg2_EC30to60E2r2.WCYCL20TR',tmp_path,tmp_sub,tmp_path,f'data_remap_{nlat_dst}x{nlon_dst}')

### comment/uncomment to disable/enable
# do_atm_h0 = True
do_atm_h1 = True
# do_atm_h2 = True
# do_lnd_h0 = True

execute   = True
overwrite = True
print_cmd = True

yr1,yr2 = 1950,1999
# yr1,yr2 = 2000,2009
# yr1,yr2 = 2010,2014
# yr1,yr2 = 1,10


# src_top_dir,src_sub_dir = os.getenv('MEMBERWORK')+'/cli115/e3sm_scratch','run'
# dst_top_dir,dst_sub_dir = '/gpfs/alpine/cli115/proj-shared/hannah6/INCITE2022/LOW-CLD-CESS',f'data_remap_{nlat_dst}x{nlon_dst}'

### chrysalis v3 dev
# src_top_dir,src_sub_dir = '/lcrc/group/e3sm/ac.whannah/E3SMv3_dev','archive/atm/hist'
# dst_top_dir,dst_sub_dir = '/lcrc/group/e3sm/ac.whannah/E3SMv3_dev',f'data_remap_{nlat_dst}x{nlon_dst}'

# src_top_dir,src_sub_dir = '/global/cfs/cdirs/m3312/whannah/2022-CPL/','data_native'
# dst_top_dir,dst_sub_dir = '/global/cfs/cdirs/m3312/whannah/2022-CPL/',f'data_remap_{nlat_dst}x{nlon_dst}'

# src_top_dir,src_sub_dir = '/global/cfs/cdirs/m3312/whannah/2023-CPL/','archive/atm/hist'
# dst_top_dir,dst_sub_dir = '/global/cfs/cdirs/m3312/whannah/2023-CPL/',f'data_remap_{nlat_dst}x{nlon_dst}'

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
    file_list = sorted( glob.glob(f'{src_dir}/{case}.eam.h1.*') )
    for f,f_in in enumerate(file_list):
        # remap_flag = False
        # if do_atm_h0 and f'{atm_comp}.h0.' in f_in : remap_flag = True
        # if do_atm_h1 and f'{atm_comp}.h1.' in f_in : remap_flag = True
        # if do_atm_h2 and f'{atm_comp}.h2.' in f_in : remap_flag = True
        # if do_lnd_h0 and f'{lnd_comp}.h0.' in f_in : remap_flag = True

        remap_flag = True

        # don't regrid data that's already been regridded
        if '.remap_' in f_in : remap_flag = False

        yr = int( f_in.replace(f'{src_dir}/{case}.eam.h1.','').replace(f'.nc','').split('-')[0] )

        if 'yr1' in locals() and yr<yr1: remap_flag = False
        if 'yr2' in locals() and yr>yr2: remap_flag = False

        # if 'yr1' in locals() and yr_list[f]<yr1: remap_flag = False
        # if 'yr2' in locals() and yr_list[f]>yr2: remap_flag = False

        if remap_flag :
            cnt += 1

            src_file_name = f_in #f'{src_dir}/{f_in}'
            dst_file_name = src_file_name.replace('.nc',f'.remap_{nlat_dst}x{nlon_dst}.nc')

            if os.path.isfile(dst_file_name) :
                if overwrite : os.remove(dst_file_name)
                else : continue

            # Change path of output file
            dst_file_name = dst_file_name.replace(src_dir,dst_dir)

            cmd = f'ncremap -m {map_file} -i {src_file_name} -o {dst_file_name}'

            run_cmd(cmd)
#===============================================================================
    if cnt==0:
        print(files)
        print()
        print(f'  do_atm_h0: {do_atm_h0}')
        print(f'  do_atm_h1: {do_atm_h1}')
        print(f'  do_atm_h2: {do_atm_h2}')
        print(f'  do_lnd_h0: {do_lnd_h0}')
        print(f'  dst_grid : {nlat_dst} x {nlon_dst}')
        print(f'  atm_comp : {atm_comp}')
        print('\nNo files found for remapping...?\n')
#===============================================================================
print('Done.')


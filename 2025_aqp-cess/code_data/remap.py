#!/usr/bin/env python
#---------------------------------------------------------------------------------------------------
''' map file for CMIP remap method:

MAP_ROOT=/pscratch/sd/w/whannah/2024-AQP-CESS/map_files
GRID_ROOT=/pscratch/sd/w/whannah/2024-AQP-CESS/grid_files
GRID_DST=${GRID_ROOT}/cmip6_90x180_scrip.20181001.nc

# generate model grid file
NE=90 # 30 / 45 / 60 / 90
GenerateCSMesh --alt --res ${NE} --file ${GRID_ROOT}/ne${NE}.g
GenerateVolumetricMesh --in ${GRID_ROOT}/ne${NE}.g --out ${GRID_ROOT}/ne${NE}pg2.g --np 2 --uniform
ConvertMeshToSCRIP --in ${GRID_ROOT}/ne${NE}pg2.g --out ${GRID_ROOT}/ne${NE}pg2_scrip.nc

GRIDIN=ne30pg2; ncremap -6 --alg_typ=traave --grd_src=${GRID_ROOT}/${GRIDIN}_scrip.nc --grd_dst=${GRID_DST} --map=${MAP_ROOT}/map_${GRIDIN}_to_90x180_traave.nc
GRIDIN=ne45pg2; ncremap -6 --alg_typ=traave --grd_src=${GRID_ROOT}/${GRIDIN}_scrip.nc --grd_dst=${GRID_DST} --map=${MAP_ROOT}/map_${GRIDIN}_to_90x180_traave.nc
GRIDIN=ne60pg2; ncremap -6 --alg_typ=traave --grd_src=${GRID_ROOT}/${GRIDIN}_scrip.nc --grd_dst=${GRID_DST} --map=${MAP_ROOT}/map_${GRIDIN}_to_90x180_traave.nc
GRIDIN=ne90pg2; ncremap -6 --alg_typ=traave --grd_src=${GRID_ROOT}/${GRIDIN}_scrip.nc --grd_dst=${GRID_DST} --map=${MAP_ROOT}/map_${GRIDIN}_to_90x180_traave.nc
'''
#---------------------------------------------------------------------------------------------------
import sys,os,glob,subprocess as sp
do_h0, do_h1, overwrite = False, False, False

data_root = '/pscratch/sd/w/whannah/2024-AQP-CESS'

cases = []

# cases.append('E3SM.2024-AQP-CESS-00.FAQP.GNUCPU.ne30pg2_ne30pg2.NN_32.SSTP_0K')
cases.append('E3SM.2024-AQP-CESS-00.FAQP.GNUCPU.ne45pg2_ne45pg2.NN_64.SSTP_0K')
cases.append('E3SM.2024-AQP-CESS-00.FAQP.GNUCPU.ne60pg2_ne60pg2.NN_128.SSTP_0K')
cases.append('E3SM.2024-AQP-CESS-00.FAQP.GNUCPU.ne90pg2_ne90pg2.NN_256.SSTP_0K')

cases.append('E3SM.2024-AQP-CESS-00.FAQP-MMF1.GNUGPU.ne30pg2_ne30pg2.NN_32.SSTP_0K')
cases.append('E3SM.2024-AQP-CESS-00.FAQP-MMF1.GNUGPU.ne45pg2_ne45pg2.NN_64.SSTP_0K')
cases.append('E3SM.2024-AQP-CESS-00.FAQP-MMF1.GNUGPU.ne60pg2_ne60pg2.NN_128.SSTP_0K')
cases.append('E3SM.2024-AQP-CESS-00.FAQP-MMF1.GNUGPU.ne90pg2_ne90pg2.NN_256.SSTP_0K')


cases.append('E3SM.2024-AQP-CESS-00.FAQP.GNUCPU.ne30pg2_ne30pg2.NN_32.SSTP_0K.ALT-NCPL_72')
cases.append('E3SM.2024-AQP-CESS-00.FAQP.GNUCPU.ne45pg2_ne45pg2.NN_64.SSTP_0K.ALT-NCPL_72')
cases.append('E3SM.2024-AQP-CESS-00.FAQP.GNUCPU.ne60pg2_ne60pg2.NN_128.SSTP_0K.ALT-NCPL_72')
cases.append('E3SM.2024-AQP-CESS-00.FAQP.GNUCPU.ne90pg2_ne90pg2.NN_256.SSTP_0K.ALT-NCPL_72')


# omit flags from case list
for c in cases: 
    if c[0]=='-': cases.remove(c)

### comment/uncomment to disable/enable
# do_h0 = True
do_h1 = True

execute   = True
overwrite = False
print_cmd = True
write_log = False


nlat_dst,nlon_dst =  90,180
# nlat_dst,nlon_dst = 180,360



idir_sub = 'archive/atm/hist'
odir_sub = f'data_remap_{nlat_dst}x{nlon_dst}'

# vert_pressure_remap,pressure_level_file = False,'grid_files/vrt_prs_ERA5.nc'
# if vert_pressure_remap: odir_sub = odir_sub+'_prs'

log_file = './remap.out'

class tcolor: ENDC,RED,GREEN,MAGENTA,CYAN = '\033[0m','\033[31m','\033[32m','\033[35m','\033[36m'

#===============================================================================
#===============================================================================
for case in cases :
    atm_comp = 'eam'
    
    idir = f'{data_root}/{case}/{idir_sub}'     # Input directory
    odir = f'{data_root}/{case}/{odir_sub}'     # Output directory
    #---------------------------------------------------------------------------
    map_root = '/pscratch/sd/w/whannah/2024-AQP-CESS/map_files'
    map_file = None
    if 'ne30pg2' in case: map_file = f'{map_root}/map_ne30pg2_to_{nlat_dst}x{nlon_dst}_traave.nc'
    if 'ne45pg2' in case: map_file = f'{map_root}/map_ne45pg2_to_{nlat_dst}x{nlon_dst}_traave.nc'
    if 'ne60pg2' in case: map_file = f'{map_root}/map_ne60pg2_to_{nlat_dst}x{nlon_dst}_traave.nc'
    if 'ne90pg2' in case: map_file = f'{map_root}/map_ne90pg2_to_{nlat_dst}x{nlon_dst}_traave.nc'
    if map_file is None: raise ValueError(f'map_file not set correctly!\ncase: {case}')
    #---------------------------------------------------------------------------
    print(f'case     : {case}')
    print(f'data dir : {idir}')
    print(f'out dir  : {odir}')
    print(f'map file : {map_file}')
    # print(f'log file : {log_file}')
    print('')
    #---------------------------------------------------------------------------
    if not os.path.exists(odir): os.mkdir(odir)

    files = sorted( glob.glob(f'{idir}/{case}.{atm_comp}.h*nc') )
    
    cnt = 0
    for f_in in files : 
        remap_flag = False
        if do_h0  and f'{atm_comp}.h0' in f_in : remap_flag = True
        if do_h1  and f'{atm_comp}.h1' in f_in : remap_flag = True

        # don't remap already remapped data
        if '.remap_' in f_in : remap_flag = False
        
        if remap_flag :
            cnt += 1
            f_out = f_in.replace('.nc','.remap.nc')

            if os.path.isfile(f_out) :
                if overwrite : os.remove(f_out)
                else : continue

            src_file_name = f_in
            dst_file_name = src_file_name.replace('.nc',f'.remap_{nlat_dst}x{nlon_dst}.nc')

            # Change directory
            dst_file_name = dst_file_name.replace(f'/{idir_sub}/',f'/{odir_sub}/')

            cmd = f'ncremap --var_lst=PRECT,FLNT,U850 -m {map_file} -i {src_file_name} -o {dst_file_name} '


            # if vert_pressure_remap: cmd += f' --vrt_fl={pressure_level_file}'

            if     print_cmd: print(f'\n{tcolor.GREEN}{cmd}{tcolor.ENDC}\n')
            if not print_cmd: print('    '+f_in+'  >  '+idir+f_out)
            if write_log: cmd += ' > '+log_file

            if execute: 
                # os.system(cmd)
                try:
                    sp.check_output(cmd,shell=True,universal_newlines=True)
                except sp.CalledProcessError as error:
                    print(error.output)
                    exit()
#===============================================================================
#===============================================================================
    if cnt==0:
        print(files)
        print()
        print(f'  do_h0    : {do_h0}')
        print(f'  do_h1    : {do_h1}')
        print(f'  dst_grid : {nlat_dst} x {nlon_dst}')
        print(f'  atm_comp : {atm_comp}')
        print(f'  odir     : {odir}')

        print('\nNo files found for remapping...?\n')

print()
print('done.')
print()
#===============================================================================
#===============================================================================

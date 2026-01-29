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
do_h0, do_h1, do_h2, do_clm, overwrite = False, False, False, False, False

# scratch_dir = '/global/cfs/cdirs/m3312/whannah/E3SM'
# scratch_dir = '/pscratch/sd/w/whannah/e3sm_scratch/pm-cpu'
scratch_dir = '/pscratch/sd/w/whannah/e3sm_scratch/pm-gpu'
data_sub = 'archive/atm/hist'
atm_comp = 'eam'

cases = []
cases.append('E3SM.2023-RAD-SENS-00.GNUGPU.ne30pg2_oECv3.F2010-MMF1.NXY_64x1.RNX_1.MCICA_OFF')
# cases.append('E3SM.2023-RAD-SENS-00.GNUGPU.ne30pg2_oECv3.F2010-MMF1.NXY_64x1.RNX_8.MCICA_OFF')
# cases.append('E3SM.2023-RAD-SENS-00.GNUGPU.ne30pg2_oECv3.F2010-MMF1.NXY_64x1.RNX_64.MCICA_OFF')
# cases.append('E3SM.2023-RAD-SENS-00.GNUGPU.ne30pg2_oECv3.F2010-MMF1.NXY_64x1.RNX_1.MCICA_ON')
# cases.append('E3SM.2023-RAD-SENS-00.GNUGPU.ne30pg2_oECv3.F2010-MMF1.NXY_64x1.RNX_8.MCICA_ON')
# cases.append('E3SM.2023-RAD-SENS-00.GNUGPU.ne30pg2_oECv3.F2010-MMF1.NXY_64x1.RNX_64.MCICA_ON')

# omit flags from case list
for c in cases: 
    if c[0]=='-': cases.remove(c)

### comment/uncomment to disable/enable
# do_h0 = True
do_h1 = True
# do_h2 = True
# do_clm = True

execute   = True
overwrite = True
print_cmd = True
# write_log = False

nlat_dst,nlon_dst =  90,180
# nlat_dst,nlon_dst = 180,360


vert_pressure_remap = False
pressure_level_file = 'grid_files/vrt_prs_ERA5.nc'

# log_file = './regrid.case.out'

class tcolor: ENDC,RED,GREEN,MAGENTA,CYAN = '\033[0m','\033[31m','\033[32m','\033[35m','\033[36m'

#===============================================================================
#===============================================================================
for case in cases :
    #---------------------------------------------------------------------------
    # Input directory
    data_dir = f'{scratch_dir}/{case}/{data_sub}/'
    #---------------------------------------------------------------------------
    # Output directory
    if vert_pressure_remap: 
        odir_sub = f'data_remap_{nlat_dst}x{nlon_dst}_prs'
    else:
        odir_sub = f'data_remap_{nlat_dst}x{nlon_dst}'
    odir = f'{scratch_dir}/{case}/{odir_sub}'
    if not os.path.exists(odir): os.mkdir(odir)
    #---------------------------------------------------------------------------
    # Map file
    src_grid_name = 'ne30pg2'
    map_file = f'map_files/map_{src_grid_name}_to_cmip6_{nlat_dst}x{nlon_dst}_aave.nc'
    #---------------------------------------------------------------------------
    print(f'case     : {case}')
    print(f'data dir : {data_dir}')
    print(f'out dir  : {odir}')
    print(f'map file : {map_file}')
    # print(f'log file : {log_file}')
    print('')
    # exit()
    #---------------------------------------------------------------------------
    files = os.listdir(data_dir)
    cnt = 0 
    for f_in in files : 
        remap_flag = False
        if do_h0  and f'{atm_comp}.h0' in f_in : remap_flag = True
        if do_h1  and f'{atm_comp}.h1' in f_in : remap_flag = True
        if do_h2  and f'{atm_comp}.h2' in f_in : remap_flag = True
        # don't remap already remapped data
        if '.remap_' in f_in : remap_flag = False
        if remap_flag :
            cnt += 1
            f_out = f_in.replace('.nc','.remap.nc')
            if os.path.isfile(f_out) :
                if overwrite : os.remove(f_out)
                else : continue
            src_file_name = data_dir+f_in
            #-------------------------------------------------------------------
            dst_file_name = src_file_name.replace('.nc',f'.remap_{nlat_dst}x{nlon_dst}.nc')
            # Change directory
            dst_file_name = dst_file_name.replace(f'/{data_sub}/',f'/{odir_sub}/')
            #-------------------------------------------------------------------
            cmd = f'ncremap -m {map_file} -i {src_file_name} -o {dst_file_name}'
            if vert_pressure_remap: cmd += f' --vrt_fl={pressure_level_file}'
            if print_cmd:
                print(f'\n{tcolor.GREEN}{cmd}{tcolor.ENDC}\n')
            else:
                print('    '+f_in+'  >  '+data_dir+f_out)
            # if write_log: cmd += ' > '+log_file
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
        print(f'  do_h2    : {do_h2}')
        print(f'  do_clm   : {do_clm}')
        print(f'  dst_grid : {nlat_dst} x {nlon_dst}')
        print(f'  atm_comp : {atm_comp}')
        print(f'  odir     : {odir}')
        print('\nNo files found for remapping...?\n')
    else:
        print(); print('done.'); print()
#===============================================================================
#===============================================================================

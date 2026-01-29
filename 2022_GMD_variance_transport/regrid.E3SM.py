#!/usr/bin/env python

### commands for making grid files
# ncremap -G ttl='FV grid 73x144'#latlon=73,144#lat_typ=fv#lon_typ=grn_ctr -g $HOME/E3SM/data_grid/73x144_scrip.nc

### commands for making map files
# SRC=ne30pg2; ncremap --src_grd=${HOME}/E3SM/data_grid/${SRC}_scrip.nc --dst_grd=$HOME/HICCUP/files_grid/scrip_90x180_s2n.nc -m $HOME/maps/map_${SRC}_to_90x180_aave.nc

# GRID=ne30np4;NLAT=73;NLON=144;ncremap --alg_typ=aave --grd_src=$HOME/E3SM/data_grid/${GRID}_scrip.nc  --grd_dst=$HOME/E3SM/data_grid/${NLAT}x${NLON}_scrip.nc --map=$HOME/maps/map_${GRID}_to_${NLAT}x${NLON}_aave.nc
# GRID=ne30pg2;NLAT=73;NLON=144;ncremap --alg_typ=aave --grd_src=$HOME/E3SM/data_grid/${GRID}_scrip.nc  --grd_dst=$HOME/E3SM/data_grid/${NLAT}x${NLON}_scrip.nc --map=$HOME/maps/map_${GRID}_to_${NLAT}x${NLON}_aave.nc

### map file for CMIP remap method:
# GRIDIN=ne45pg2; GRIDOUT=90x180 ; ncremap -6 --alg_typ=aave --grd_src=$HOME/E3SM/data_grid/${GRIDIN}_scrip.nc --grd_dst=$HOME/E3SM/data_grid/cmip6_${GRIDOUT}_scrip.20181001.nc --map=$HOME/maps/map_${GRID}_to_cmip6_180x360_aave.nc

#===============================================================================
#===============================================================================
import sys,os,subprocess as sp
from optparse import OptionParser
parser = OptionParser()
parser.add_option('-n',dest='num_file',default=None,help='sets number of files to print')
# parser.add_option('--comp',dest='component',default=None,help='model component history file to search for')
(opts, args) = parser.parse_args()
#===============================================================================
#===============================================================================
home = os.getenv('HOME')

# scratch_dir = os.getenv('SCRATCH')+'/e3sm_scratch/cori-knl/'        ### Cori scratch space
scratch_dir = os.getenv('MEMBERWORK')+'/cli115/e3sm_scratch/'       ### OLCF scratch space

if scratch_dir=='' :exit('ERROR: scratch directory not set!')
#===============================================================================
#===============================================================================
do_h0, do_h1, do_h2, do_clm, overwrite = False, False, False, False, False

if len(sys.argv) < 2 :
    print('ERROR: no case name provided!')
    exit()
else :
    cases = sys.argv[1:]

# if cases==['RGMA']:
#     cases = []
#     cases.append('E3SM.RGMA.ne120pg2_r05_oECv3.FC5AV1C-H01A.00')
#     cases.append('E3SM.RGMA.ne30pg2_r05_oECv3.F-MMF1.CRMNX_64.CRMDX_2000.RADNX_4.00')
#     cases.append('E3SM.RGMA.ne30pg2_r05_oECv3.FC5AV1C-L.00')


### comment/uncomment to disable/enable
# do_h0 = True
do_h1 = True
# do_h2 = True
# do_clm = True

execute   = True
overwrite = True
print_cmd = True
write_log = False

nlat_dst,nlon_dst =  73, 144
# nlat_dst,nlon_dst =  90,180
# nlat_dst,nlon_dst = 180,360
# nlat_dst,nlon_dst = 121, 240

log_file = './regrid.case.out'

class tcolor: ENDC,RED,GREEN,MAGENTA,CYAN = '\033[0m','\033[31m','\033[32m','\033[35m','\033[36m'

#===============================================================================
#===============================================================================
for case in cases :

    #---------------------------------------------------------------------------
    # Input directory
    #---------------------------------------------------------------------------
    data_sub = 'run'
    data_dir = f'{scratch_dir}/{case}/{data_sub}/'
    # if 'INCITE2019' in case: data_dir = scratch_dir+case+'/data_native/'
    # if 'E3SM.RGMA.' in case: 
    #     data_sub = 'data_native'
    #     data_dir = f'/global/project/projectdirs/m3305/RGMA-2020/{case}/{data_sub}/'
    #---------------------------------------------------------------------------
    # Output directory
    #---------------------------------------------------------------------------
    odir_sub = f'data_remap_{nlat_dst}x{nlon_dst}'
    odir = f'{scratch_dir}/{case}/{odir_sub}'
    # odir = f'{scratch_dir}/{case}/data_remap_pmp_{nlat_dst}x{nlon_dst}'
    # if 'E3SM.RGMA.' in case: odir = f'/global/project/projectdirs/m3305/RGMA-2020/{case}/{odir_sub}'
    #---------------------------------------------------------------------------
    #---------------------------------------------------------------------------
    print('case     : '+case)
    print('data dir : '+data_dir)
    print('out dir  : '+odir)
    print('log file : '+log_file)
    print('')
    # exit()

    files = os.listdir(data_dir)

    if 'ne30_'     in case: src_grid_name = 'ne30np4'
    if 'ne30pg2_'  in case: src_grid_name = 'ne30pg2'
    if 'ne30pg2.'  in case: src_grid_name = 'ne30pg2'

    #---------------------------------------------------------------------------
    # specify map file
    #---------------------------------------------------------------------------

    # if 'map_file' not in locals() and nlat_dst==180 and nlon_dst==360:
    #     map_file = f'{home}/maps/map_{src_grid_name}_to_cmip6_{nlat_dst}x{nlon_dst}_aave.20200624.nc'

    if 'map_file' not in locals():
        map_file = f'{home}/maps/map_{src_grid_name}_to_{nlat_dst}x{nlon_dst}_aave.nc'

    #---------------------------------------------------------------------------
    #---------------------------------------------------------------------------
    if not os.path.exists(odir): os.mkdir(odir)

    atm_comp,lnd_comp = 'eam','elm'
    # if 'E3SM.PI-CPL.v1.' in case: atm_comp,lnd_comp = 'cam','clm'
        
    cnt = 0
    for f_in in files : 
        remap_flag = False
        if do_h0  and f'{atm_comp}.h0' in f_in : remap_flag = True
        if do_h1  and f'{atm_comp}.h1' in f_in : remap_flag = True
        if do_h2  and f'{atm_comp}.h2' in f_in : remap_flag = True
        if do_clm and f'{lnd_comp}2.h' in f_in : remap_flag = True

        # don't remap already remapped data
        if '.remap_' in f_in : remap_flag = False

        if remap_flag:
            if do_h1  and f'{atm_comp}.h1' in f_in : tmp = f'{atm_comp}.h1'
            if do_h2  and f'{atm_comp}.h2' in f_in : tmp = f'{atm_comp}.h2'

            ### special logic for KE spectra - use bilin for wind data but conservative for PS
            # if do_h1  and f'{atm_comp}.h1' in f_in : map_file = f'{home}/maps/map_{src_grid_name}_to_{nlat_dst}x{nlon_dst}.20201019.nc'
            # if do_h2  and f'{atm_comp}.h2' in f_in : map_file = f'{home}/maps/map_{src_grid_name}_to_{nlat_dst}x{nlon_dst}_bilin.20201019.nc'

        if '.remap_' in f_in : remap_flag = False
        
        if remap_flag :
            cnt += 1
            f_out = f_in.replace('.nc','.remap.nc')

            if os.path.isfile(f_out) :
                if overwrite : os.remove(f_out)
                else : continue

            src_file_name = data_dir+f_in
            dst_file_name = src_file_name.replace('.nc',f'.remap_{nlat_dst}x{nlon_dst}.nc')


            # Change directory
            dst_file_name = dst_file_name.replace(f'/{data_sub}/',f'/{odir_sub}/')

            cmd = f'ncremap -m {map_file} -i {src_file_name} -o {dst_file_name}'
            if print_cmd:
                msg = tcolor.GREEN + cmd + tcolor.ENDC
                print('\n'+msg+'\n')
            else:
                print('    '+f_in+'  >  '+data_dir+f_out)
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
        print(f'  do_h2    : {do_h2}')
        print(f'  do_clm   : {do_clm}')
        print(f'  dst_grid : {nlat_dst} x {nlon_dst}')
        print(f'  atm_comp : {atm_comp}')
        print(f'  odir     : {odir}')
        

        print('\nNo files found for remapping...?\n')

#===============================================================================
#===============================================================================

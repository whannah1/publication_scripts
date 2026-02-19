#!/usr/bin/env python
#===============================================================================
'''
commands to create grid and map files:

GRID_ROOT=/global/homes/w/whannah/Research/E3SM/pub_figs/2025_scream_decadal/files_grid
MAPS_ROOT=/global/homes/w/whannah/Research/E3SM/pub_figs/2025_scream_decadal/files_maps

NE=30
SRC_GRID=ne${NE}pg2
# DST_NY=90
# DST_NX=180
DST_NY=73
DST_NX=144
DST_GRID=${DST_NY}x${DST_NX}

# TIMESTAMP=20230222
# SRC_GRID=${GRID_ROOT}/${SRC_GRID}_scrip.${TIMESTAMP}.nc
# DST_GRID=${GRID_ROOT}/${DST_GRID}_scrip.${TIMESTAMP}.nc
# MAP_FILE=${MAPS_ROOT}/map_${SRC_GRID}_to_${DST_GRID}_aave.${TIMESTAMP}.nc

# or without time stamp
SRC_GRID=${GRID_ROOT}/${SRC_GRID}_scrip.nc
DST_GRID=${GRID_ROOT}/${DST_GRID}_scrip.nc
MAP_FILE=${MAPS_ROOT}/map_${SRC_GRID}_to_${DST_GRID}_aave.nc

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
htype_list = []
var_list_arg_list = []
def add_case(case_in,src_dir_in,src_sub_in,dst_dir_in,dst_sub_in,htype,var_list_arg):
   global case_list,src_top_dir_list,src_sub_dir_list
   case_list.append(case_in)
   src_top_dir_list.append(src_dir_in); src_sub_dir_list.append(src_sub_in)
   dst_top_dir_list.append(dst_dir_in); dst_sub_dir_list.append(dst_sub_in)
   htype_list.append(htype); var_list_arg_list.append(var_list_arg)
#-------------------------------------------------------------------------------
#===============================================================================

# nlat_dst,nlon_dst,map_file = 90,180,os.getenv('HOME')+'/maps/map_ne30pg2_to_90x180_aave.nc'
# nlat_dst,nlon_dst,map_file = 90,180,os.getenv('HOME')+'/Research/E3SM/pub_figs/2025_scream_decadal/files_maps/map_ne30pg2_to_90x180_traave.20250624.nc'
nlat_dst,nlon_dst,map_file = 73,144,'/global/homes/w/whannah/maps/map_ne30pg2_to_73x144_traave.nc'

# data_root,tmp_sub = '/global/cfs/cdirs/m3312/whannah/e3smv3_amip','archive/atm/hist'
# add_case('v3.LR.amip_0101',data_root,tmp_sub,data_root,f'data_remap_{nlat_dst}x{nlon_dst}')
# add_case('v3.LR.amip_0151',data_root,tmp_sub,data_root,f'data_remap_{nlat_dst}x{nlon_dst}')
# add_case('v3.LR.amip_0201',data_root,tmp_sub,data_root,f'data_remap_{nlat_dst}x{nlon_dst}')
# htype = 'eam.h3'


#===============================================================================
# ne1024 decadal cases

ne1024_data_root = '/global/cfs/cdirs/e3smdata/simulations/scream-decadal'

# tmp_src_sub  = 
# tmp_dst_sub  = f'data_remap_{nlat_dst}x{nlon_dst}'#/output.scream.decadal.3hourlyINST_ne30pg2'
# htype        = 
# var_list_arg = 

# tmp_src_sub  = 'run/output.scream.decadal.6hourlyINST_ne30pg2'
# tmp_dst_sub  = f'data_remap_{nlat_dst}x{nlon_dst}'#/output.scream.decadal.6hourlyINST_ne30pg2'
# htype        = 
# var_list_arg = 

# add_case('decadal-production-run6', 
#          src_dir_in=ne1024_data_root,
#          src_sub_in='run/output.scream.decadal.3hourlyINST_ne30pg2',
#          dst_dir_in=ne1024_data_root,
#          dst_sub_in=f'data_remap_{nlat_dst}x{nlon_dst}/output.scream.decadal.3hourlyINST_ne30pg2',
#          htype='output.scream.decadal.3hourlyINST_ne30pg2.INSTANT.nhours_x3',
#          var_list_arg='-v area,U_at_850hPa,LiqWaterPath,IceWaterPath' )

# add_case('decadal-production-run6', 
#          src_dir_in=ne1024_data_root,
#          src_sub_in='run/output.scream.decadal.6hourlyINST_ne30pg2',
#          dst_dir_in=ne1024_data_root,
#          dst_sub_in=f'data_remap_{nlat_dst}x{nlon_dst}/output.scream.decadal.6hourlyINST_ne30pg2',
#          htype='output.scream.decadal.6hourlyINST_ne30pg2.INSTANT.nhours_x6',
#          var_list_arg='-v area,ps,T_mid,qv,qc,qi,U,omega' )

# add_case('decadal-production-run6', 
#          src_dir_in=ne1024_data_root,
#          src_sub_in='run/output.scream.decadal.6hourlyAVG_ne30pg2',
#          dst_dir_in=ne1024_data_root,
#          dst_sub_in=f'data_remap_{nlat_dst}x{nlon_dst}/output.scream.decadal.6hourlyAVG_ne30pg2',
#          htype='output.scream.decadal.6hourlyAVG_ne30pg2.AVERAGE.nhours_x6',
#          var_list_arg='-v area,ps,LW_flux_up_at_model_top,precip_ice_surf_mass_flux,precip_liq_surf_mass_flux' ) 

#===============================================================================
# ne256 decadal cases

# htype        = '6ha_ne30pg2.AVERAGE.nhours_x6'
# var_list_arg = '-v area,ps,LW_flux_up_at_model_top,precip_ice_surf_mass_flux,precip_liq_surf_mass_flux'

# add_case(case_in='ne256pg2_ne256pg2.F20TR-SCREAMv1.May-12.with.rain.frac.n0128',
#          src_dir_in='/global/cfs/cdirs/e3sm/beydoun/ne256pg2_ne256pg2.F20TR-SCREAMv1.May-12.with.rain.frac.n0128.completed',
#          src_sub_in='run',
#          dst_dir_in='/global/cfs/cdirs/e3smdata/simulations/scream-decadal-ne256',
#          dst_sub_in=f'data_remap_{nlat_dst}x{nlon_dst}/6ha_ne30pg2.AVERAGE.nhours_x6',
#          htype='6ha_ne30pg2.AVERAGE.nhours_x6',
#          var_list_arg='-v area,ps,LW_flux_up_at_model_top,precip_ice_surf_mass_flux,precip_liq_surf_mass_flux' ) 

# add_case(case_in='ne256pg2_ne256pg2.F20TR-SCREAMv1.July-1.spanc800.2xauto.acc150.n0032.test2.1', 
#          src_dir_in='/global/cfs/cdirs/e3sm/terai/EAMxx/ne256_Fcase_sim',
#          src_sub_in='run',
#          dst_dir_in='/global/cfs/cdirs/e3smdata/simulations/scream-decadal-ne256',
#          dst_sub_in=f'data_remap_{nlat_dst}x{nlon_dst}/6ha_ne30pg2.AVERAGE.nhours_x6',
#          htype='6ha_ne30pg2.AVERAGE.nhours_x6',
#          var_list_arg='-v area,ps,LW_flux_up_at_model_top,precip_ice_surf_mass_flux,precip_liq_surf_mass_flux' ) 


# add_case(case_in='ne256pg2_ne256pg2.F20TR-SCREAMv1.May-12.with.rain.frac.n0128',
#          src_dir_in='/global/cfs/cdirs/e3sm/beydoun/ne256pg2_ne256pg2.F20TR-SCREAMv1.May-12.with.rain.frac.n0128.completed',
#          src_sub_in='run',
#          dst_dir_in='/global/cfs/cdirs/e3smdata/simulations/scream-decadal-ne256',
#          dst_sub_in=f'data_remap_{nlat_dst}x{nlon_dst}/3hi_ne30pg2.INSTANT.nhours_x3',
#          htype='3hi_ne30pg2.INSTANT.nhours_x3',
#          var_list_arg='-v area,U_at_850hPa,LiqWaterPath,IceWaterPath' )

# add_case(case_in='ne256pg2_ne256pg2.F20TR-SCREAMv1.July-1.spanc800.2xauto.acc150.n0032.test2.1', 
#          src_dir_in='/global/cfs/cdirs/e3sm/terai/EAMxx/ne256_Fcase_sim',
#          src_sub_in='run',
#          dst_dir_in='/global/cfs/cdirs/e3smdata/simulations/scream-decadal-ne256',
#          dst_sub_in=f'data_remap_{nlat_dst}x{nlon_dst}/3hi_ne30pg2.INSTANT.nhours_x3',
#          htype='3hi_ne30pg2.INSTANT.nhours_x3',
#          var_list_arg='-v area,U_at_850hPa,LiqWaterPath,IceWaterPath' )

#===============================================================================
# DYNAMO MJO Hindcasts


# SCREAM.2025-DYNAMO-00.ne1024pg2_ICOS10.2011-11-10.rfrac_fix_0
# SCREAM.2025-DYNAMO-00.ne1024pg2_ICOS10.2011-11-10.rfrac_fix_1

# SCREAM.2025-DYNAMO-01.ne1024pg2_ICOS10.2011-11-15.rfrac_fix_0
# SCREAM.2025-DYNAMO-01.ne1024pg2_ICOS10.2011-11-15.rfrac_fix_1

# SCREAM.2025-DYNAMO-01.ne1024pg2_ICOS10.2011-11-20.rfrac_fix_0
# SCREAM.2025-DYNAMO-01.ne1024pg2_ICOS10.2011-11-20.rfrac_fix_1

hx_dir_in = '/global/cfs/cdirs/e3smdata/simulations/scream-decadal/DYNAMO_MJO_hindcasts'
hx_1hr_var_list = '-v ps,precip_total_surf_mass_flux,VapWaterPath,LiqWaterPath,IceWaterPath,RainWaterPath,surf_sens_flux,surface_upward_latent_heat_flux,U_at_850hPa,U_at_200hPa,SW_flux_up_at_model_top,SW_flux_dn_at_model_top,LW_flux_up_at_model_top'
hx_6hr_var_list = '-v ps,omega,horiz_winds,qv,qc,qi,qr,qm,RelativeHumidity,T_mid,z_mid,rad_heating_pdel,cldfrac_tot_for_analysis,cldfrac_liq,cldfrac_ice_for_analysis'

# 1-hr 2D output
add_case(case_in='SCREAM.2025-DYNAMO-01.ne1024pg2_ICOS10.2011-11-10.rfrac_fix_0', 
         src_dir_in=hx_dir_in, src_sub_in='run', dst_dir_in=hx_dir_in, dst_sub_in=f'data_remap_{nlat_dst}x{nlon_dst}',
         htype='output.scream.2D.1hr.ne30pg2.AVERAGE.nhours_x1', var_list_arg=hx_1hr_var_list )
add_case(case_in='SCREAM.2025-DYNAMO-01.ne1024pg2_ICOS10.2011-11-10.rfrac_fix_1', 
         src_dir_in=hx_dir_in, src_sub_in='run', dst_dir_in=hx_dir_in, dst_sub_in=f'data_remap_{nlat_dst}x{nlon_dst}',
         htype='output.scream.2D.1hr.ne30pg2.AVERAGE.nhours_x1', var_list_arg=hx_1hr_var_list )
# add_case(case_in='SCREAM.2025-DYNAMO-01.ne1024pg2_ICOS10.2011-11-15.rfrac_fix_0', 
#          src_dir_in=hx_dir_in, src_sub_in='run', dst_dir_in=hx_dir_in, dst_sub_in=f'data_remap_{nlat_dst}x{nlon_dst}',
#          htype='output.scream.2D.1hr.ne30pg2.AVERAGE.nhours_x1', var_list_arg=hx_1hr_var_list )
# add_case(case_in='SCREAM.2025-DYNAMO-01.ne1024pg2_ICOS10.2011-11-15.rfrac_fix_1', 
#          src_dir_in=hx_dir_in, src_sub_in='run', dst_dir_in=hx_dir_in, dst_sub_in=f'data_remap_{nlat_dst}x{nlon_dst}',
#          htype='output.scream.2D.1hr.ne30pg2.AVERAGE.nhours_x1', var_list_arg=hx_1hr_var_list )
# add_case(case_in='SCREAM.2025-DYNAMO-01.ne1024pg2_ICOS10.2011-11-20.rfrac_fix_0', 
#          src_dir_in=hx_dir_in, src_sub_in='run', dst_dir_in=hx_dir_in, dst_sub_in=f'data_remap_{nlat_dst}x{nlon_dst}',
#          htype='output.scream.2D.1hr.ne30pg2.AVERAGE.nhours_x1', var_list_arg=hx_1hr_var_list )
# add_case(case_in='SCREAM.2025-DYNAMO-01.ne1024pg2_ICOS10.2011-11-20.rfrac_fix_1', 
#          src_dir_in=hx_dir_in, src_sub_in='run', dst_dir_in=hx_dir_in, dst_sub_in=f'data_remap_{nlat_dst}x{nlon_dst}',
#          htype='output.scream.2D.1hr.ne30pg2.AVERAGE.nhours_x1', var_list_arg=hx_1hr_var_list )

# # 6-hr 3D output - don't need this yet
# add_case(case_in='SCREAM.2025-DYNAMO-01.ne1024pg2_ICOS10.2011-11-15.rfrac_fix_0', 
#          src_dir_in=hx_dir_in, src_sub_in='run', dst_dir_in=hx_dir_in, dst_sub_in=f'data_remap_{nlat_dst}x{nlon_dst}',
#          htype='output.scream.3D.6hr.ne30pg2.AVERAGE.nhours_x6', var_list_arg=hx_6hr_var_list )
# add_case(case_in='SCREAM.2025-DYNAMO-01.ne1024pg2_ICOS10.2011-11-15.rfrac_fix_1', 
#          src_dir_in=hx_dir_in, src_sub_in='run', dst_dir_in=hx_dir_in, dst_sub_in=f'data_remap_{nlat_dst}x{nlon_dst}',
#          htype='output.scream.3D.6hr.ne30pg2.AVERAGE.nhours_x6', var_list_arg=hx_6hr_var_list )
# add_case(case_in='SCREAM.2025-DYNAMO-01.ne1024pg2_ICOS10.2011-11-20.rfrac_fix_0', 
#          src_dir_in=hx_dir_in, src_sub_in='run', dst_dir_in=hx_dir_in, dst_sub_in=f'data_remap_{nlat_dst}x{nlon_dst}',
#          htype='output.scream.3D.6hr.ne30pg2.AVERAGE.nhours_x6', var_list_arg=hx_6hr_var_list )
# add_case(case_in='SCREAM.2025-DYNAMO-01.ne1024pg2_ICOS10.2011-11-20.rfrac_fix_1', 
#          src_dir_in=hx_dir_in, src_sub_in='run', dst_dir_in=hx_dir_in, dst_sub_in=f'data_remap_{nlat_dst}x{nlon_dst}',
#          htype='output.scream.3D.6hr.ne30pg2.AVERAGE.nhours_x6', var_list_arg=hx_6hr_var_list )

#===============================================================================

execute   = True
overwrite = False
print_cmd = True

# yr1,yr2 = 1994,2005
# yr1,yr2 = 1995,1995
# yr1,yr2 = 1995,2004
# yr1,yr2 = 2000,2005
# yr1,yr2 = 2004,2005
yr1,yr2 = 2011,2011


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
    if not os.path.exists(path_in):
        os.makedirs(path_in)
    return       
#===============================================================================
for c,case in enumerate(case_list) :
    src_top_dir, src_sub_dir = src_top_dir_list[c], src_sub_dir_list[c]
    dst_top_dir, dst_sub_dir = dst_top_dir_list[c], dst_sub_dir_list[c]
    htype, var_list_arg = htype_list[c], var_list_arg_list[c]
    #---------------------------------------------------------------------------
    # Set in/out paths and create output directory if it doesn't exist
    src_dir = f'{src_top_dir}/{case}/{src_sub_dir}'
    dst_dir = f'{dst_top_dir}/{case}/{dst_sub_dir}'
    if not os.path.exists(dst_dir):
        print(); print(f'Creating output directory: {dst_dir}'); print()
        # check_path_and_create(f'{dst_top_dir}')
        # check_path_and_create(f'{dst_top_dir}/{case}')
        check_path_and_create(f'{dst_top_dir}/{case}/{dst_sub_dir}')
    #---------------------------------------------------------------------------
    print(f'case    : {case}')
    print(f'src_dir : {src_dir}')
    print(f'dst_dir : {dst_dir}')
    print('')
    #---------------------------------------------------------------------------
    # get list of all files to loop over
    # files = sorted( os.listdir(src_dir) )
    #---------------------------------------------------------------------------
    # if 'scream' in htype:
    #     file_list = sorted( glob.glob(f'{src_dir}/{htype}.*.nc') )
    # else:
    #     file_list = sorted( glob.glob(f'{src_dir}/{case}.{htype}.*.nc') )
    #---------------------------------------------------------------------------
    file_path = f'{src_dir}/*{htype}.*.nc'
    file_list = sorted( glob.glob(file_path) )
    #---------------------------------------------------------------------------
    if file_list==[]:
        print('\nNo files found...?\n')
        print(f'  file_path: {file_path}')
        exit()
    #---------------------------------------------------------------------------
    # cnt = 0
    # for f in file_list:
    #     print(f)
    #     cnt+=1
    #     if cnt>=5: break
    # exit()
    #---------------------------------------------------------------------------
    cnt = 0
    for f,f_in in enumerate(file_list):
        remap_flag = True

        # don't regrid data that's already been regridded
        if '.remap_' in f_in : remap_flag = False

        if 'scream' in htype:
            f_tmp =  f_in.replace(f'{src_dir}/{htype}.','')
        else:
            f_tmp =  f_in.replace(f'{src_dir}/{case}.{htype}.','')

        yr = int( f_tmp.replace(f'.nc','').split('.')[-1].split('-')[0] )

        if 'yr1' in locals() and yr<yr1: remap_flag = False
        if 'yr2' in locals() and yr>yr2: remap_flag = False

        # skip files with tilde in the name
        if '.nc~' in f_in: remap_flag = False

        src_file_name = f_in #f'{src_dir}/{f_in}'
        dst_file_name = src_file_name.replace('.nc',f'.remap_{nlat_dst}x{nlon_dst}.nc')

        # Change path of output file
        dst_file_name = dst_file_name.replace(src_dir,dst_dir)

        if os.path.isfile(dst_file_name):
            if overwrite : os.remove(dst_file_name)
            else : print(f'skipping {src_file_name}') ; remap_flag = False

        if remap_flag :
            cnt += 1

            tmp_file_root = dst_dir

            cmd = f'ncremap --pdq_opt=time,lev,ncol --tmp_dir={tmp_file_root} -m {map_file} -i {src_file_name} -o {dst_file_name} {var_list_arg}'

            run_cmd(cmd)

            # exit()

#===============================================================================
    if cnt==0:
        # print(files)
        print()
        print(f'  dst_grid : {nlat_dst} x {nlon_dst}')
        print('\nNo files found for remapping...?\n')
#===============================================================================
print('Done.')


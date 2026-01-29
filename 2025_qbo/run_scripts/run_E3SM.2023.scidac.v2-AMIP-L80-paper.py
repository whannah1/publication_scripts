#!/usr/bin/env python3
#---------------------------------------------------------------------------------------------------
class clr:END,RED,GREEN,MAGENTA,CYAN = '\033[0m','\033[31m','\033[32m','\033[35m','\033[36m'
def run_cmd(cmd): print('\n'+clr.GREEN+cmd+clr.END); os.system(cmd); return
#---------------------------------------------------------------------------------------------------
# NOTE: The branch (whannah/scidac-2023) started 
# with master @ Aug 7 w/ v3atm/eam/master_MAM5_wetaero_chemdyg:
#     commit 915e929b5118243bc2f15b90c4ca911352985524
#     Merge: 855bd01e21 b3454bb277
#     Author: Wuyin Lin <wlin@bnl.gov>
#     Date:   Mon Aug 7 07:47:35 2023 -0700
# but also includes commits to add the hdepth_scaling_factor namelist parameter
#---------------------------------------------------------------------------------------------------
import os, subprocess as sp
newcase,config,build,clean,submit,continue_run,reset_resub,st_archive = False,False,False,False,False,False,False,False

acct = 'm4310'
src_dir  = os.getenv('HOME')+'/E3SM/E3SM_SRC3' # branch => maint-2.1

# clean        = True
# newcase      = True
# config       = True
# build        = True
submit       = True
continue_run = True
### reset_resub  = True

disable_bfb = True

# queue = 'batch' # batch / debug
num_nodes   = 22
grid        = 'ne30pg2_EC30to60E2r2'
compset     = 'F20TR' # F20TR / F20TR_chemUCI-Linozv3-mam5

# stop_opt,stop_n,resub,walltime = 'nsteps',10,0,'0:30:00' 
# stop_opt,stop_n,resub,walltime = 'ndays',32,0,'1:00:00'
# stop_opt,stop_n,resub,walltime = 'ndays',365,0,'4:00:00'
# stop_opt,stop_n,resub,walltime = 'ndays',365*2,5-1,'6:00:00' # 10 years
# stop_opt,stop_n,resub,walltime = 'ndays',365*2,15-1,'6:00:00' # 30 years


nlev_list = []
# nlev_list.append( 72)
# nlev_list.append( 80)
nlev_list.append(128)

din_loc_root = '/global/cfs/cdirs/e3sm/inputdata'
init_scratch = '/global/cfs/cdirs/m4310/whannah/HICCUP/data/'

# land_init_path = '/pscratch/sd/w/whannah/e3sm_scratch/init_scratch/'
land_data_path = f'{din_loc_root}/lnd/clm2/surfdata_map'
if grid=='ne30pg2_EC30to60E2r2': 
   # land_init_file = 'ELM_spinup.ICRUELM.ne30pg2_oECv3.20-yr.2010-01-01.elm.r.2010-01-01-00000.nc'
   if 'F2010' in compset: land_data_file = 'surfdata_ne30pg2_simyr2010_c210402.nc'
   if 'F20TR' in compset: land_data_file = 'surfdata_ne30pg2_simyr1850_c210402.nc'

RUN_REFDATE = '1985-01-01'
REF_SCRATCH = '/global/cfs/cdirs/m4310/whannah/E3SM/'
RUN_REFCASE = '20220504.v2.LR.bi-grid.amip.chemMZT.chrysalis'
#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------
def main(nlev=None):
   if nlev is None: exit(' one or more arguments not provided?')

   case = '.'.join(['E3SM','2023-SCIDAC-v2-AMIP',grid,f'L{nlev}'])

   # print(case); return

   case_root = f'/pscratch/sd/w/whannah/e3sm_scratch/pm-cpu/{case}'

   if nlev== 72: init_file_atm = f'{init_scratch}/20220504.v2.LR.bi-grid.amip.chemMZT.chrysalis.eam.i.1985-01-01-00000.nc'
   if nlev== 80: init_file_atm = f'{init_scratch}/20220504.v2.LR.bi-grid.amip.chemMZT.chrysalis.eam.i.1985-01-01-00000.L80_c20230629.nc'
   if nlev==128: init_file_atm = f'{init_scratch}/20220504.v2.LR.bi-grid.amip.chemMZT.chrysalis.eam.i.1985-01-01-00000.L128_c20230629.nc'

   #------------------------------------------------------------------------------------------------
   #------------------------------------------------------------------------------------------------
   print('\n  case : '+case+'\n')
   atm_ntasks,atm_nthrds = num_nodes*128,1
   #------------------------------------------------------------------------------------------------
   # Create new case
   if newcase :
      # Check if directory already exists   
      if os.path.isdir(case_root): exit(f'\n{clr.RED}This case already exists!{clr.END}\n')
      cmd = f'{src_dir}/cime/scripts/create_newcase'
      cmd += f' --case {case}'
      cmd += f' --output-root {case_root} '
      cmd += f' --script-root {case_root}/case_scripts '
      cmd += f' --handle-preexisting-dirs u '
      cmd += f' --compset {compset}'
      cmd += f' --res {grid} '
      cmd += f' --machine pm-cpu '
      cmd += f' --pecount {atm_ntasks}x{atm_nthrds} '
      cmd += f' --project {acct} '
      cmd += f' --walltime {walltime} '
      run_cmd(cmd)
   #------------------------------------------------------------------------------------------------
   os.chdir(f'{case_root}/case_scripts')
   #------------------------------------------------------------------------------------------------
   # Configure
   if config :
      run_cmd(f'./xmlchange EXEROOT={case_root}/bld ')
      run_cmd(f'./xmlchange RUNDIR={case_root}/run ')
      #-------------------------------------------------------
      # when specifying ncdata, do it here to avoid an error message
      file=open('user_nl_eam','w');file.write(get_atm_nl_opts(init_file_atm));file.close()
      run_cmd(f'./xmlchange --append --id CAM_CONFIG_OPTS --val \" -nlev {nlev} \" ')
      run_cmd(f'./xmlchange --append --id CAM_CONFIG_OPTS --val=\'-cosp\'')
      #-------------------------------------------------------
      run_cmd(f'./xmlchange --file env_mach_pes.xml NTASKS_OCN={(16*128)}')
      run_cmd(f'./xmlchange --file env_mach_pes.xml NTASKS_ICE={(16*128)}')
      #-------------------------------------------------------
      # run_cmd(f'./xmlchange DOUT_S_ROOT={case_root}/archive ')
      #-------------------------------------------------------
      if clean : run_cmd('./case.setup --clean')
      run_cmd('./case.setup --reset')
   #------------------------------------------------------------------------------------------------
   # Build
   if build : 
      run_cmd(f'./xmlchange RUN_TYPE=hybrid,GET_REFCASE=TRUE ')
      run_cmd(f'./xmlchange RUN_REFCASE={RUN_REFCASE},RUN_REFDATE={RUN_REFDATE} ')
      run_cmd(f'./xmlchange RUN_REFDIR={REF_SCRATCH}/{RUN_REFCASE}/archive/rest/{RUN_REFDATE}-00000 ')
      #-------------------------------------------------------
      if 'debug-on' in case : run_cmd('./xmlchange --file env_build.xml --id DEBUG --val TRUE ')
      if clean : run_cmd('./case.build --clean')
      run_cmd('./case.build')
   #------------------------------------------------------------------------------------------------
   # Write the namelist options and submit the run
   if submit : 
      # run_cmd(f'./xmlchange DIN_LOC_ROOT=/global/cfs/cdirs/e3sm/inputdata')
      #-------------------------------------------------------
      file=open('user_nl_eam','w');file.write(get_atm_nl_opts(init_file_atm));file.close()
      file=open('user_nl_elm','w');file.write(get_lnd_nl_opts());file.close()
      #-------------------------------------------------------
      run_cmd(f'./xmlchange RUN_STARTDATE={RUN_REFDATE}')
      if 'queue'    in globals(): run_cmd(f'./xmlchange JOB_QUEUE={queue}')
      if 'stop_opt' in globals(): run_cmd(f'./xmlchange STOP_OPTION={stop_opt}')
      if 'stop_n'   in globals(): run_cmd(f'./xmlchange STOP_N={stop_n}')
      # if 'resub' in globals() and not continue_run: run_cmd(f'./xmlchange RESUBMIT={resub}')
      # if 'resub' in globals() and reset_resub: run_cmd(f'./xmlchange RESUBMIT={resub}')
      if 'resub' in globals(): run_cmd(f'./xmlchange RESUBMIT={resub}')
      #-------------------------------------------------------
      if     disable_bfb: run_cmd('./xmlchange BFBFLAG=FALSE')
      if not disable_bfb: run_cmd('./xmlchange BFBFLAG=TRUE')
      #-------------------------------------------------------
      if     continue_run: run_cmd('./xmlchange --file env_run.xml CONTINUE_RUN=TRUE ')   
      if not continue_run: run_cmd('./xmlchange --file env_run.xml CONTINUE_RUN=FALSE ')
      #-------------------------------------------------------
      run_cmd('./case.submit')
   #------------------------------------------------------------------------------------------------
   # Print the case name again
   print(f'\n  case : {case}\n')
#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------
def get_atm_nl_opts(ncdata_file):
   return f'''
 ncdata = '{ncdata_file}'
 
 cosp_lite = .true.

 inithist = 'NONE'

 empty_htapes = .false.
 avgflag_pertape = 'A','A','I'
 nhtfrq = 0,-3,-6
 mfilt  = 1,40,20
 fincl1 = 'Z3','UTGWSPEC','BUTGWSPEC','UTGWORO','BTAUE','BTAUW'
          ,'BUTEND1','BUTEND2','BUTEND3','BUTEND4','BUTEND5'
          ,'BTAUXSp00'
          ,'BTAUXSp01','BTAUXSp02','BTAUXSp03','BTAUXSp04'
          ,'BTAUXSp05','BTAUXSp06','BTAUXSp07','BTAUXSp08'
          ,'BTAUXSp09','BTAUXSp10','BTAUXSp11','BTAUXSp12'
          ,'BTAUXSp13','BTAUXSp14','BTAUXSp15','BTAUXSp16'
          ,'BTAUXSp17','BTAUXSp18','BTAUXSp19','BTAUXSp20'
          ,'BTAUXSp21','BTAUXSp22','BTAUXSp23','BTAUXSp24'
          ,'BTAUXSp25','BTAUXSp26','BTAUXSp27','BTAUXSp28'
          ,'BTAUXSp29','BTAUXSp30','BTAUXSp31','BTAUXSp32'
          ,'BTAUXSn01','BTAUXSn02','BTAUXSn03','BTAUXSn04'
          ,'BTAUXSn05','BTAUXSn06','BTAUXSn07','BTAUXSn08'
          ,'BTAUXSn09','BTAUXSn10','BTAUXSn11','BTAUXSn12'
          ,'BTAUXSn13','BTAUXSn14','BTAUXSn15','BTAUXSn16'
          ,'BTAUXSn17','BTAUXSn18','BTAUXSn19','BTAUXSn20'
          ,'BTAUXSn21','BTAUXSn22','BTAUXSn23','BTAUXSn24'
          ,'BTAUXSn25','BTAUXSn26','BTAUXSn27','BTAUXSn28'
          ,'BTAUXSn29','BTAUXSn30','BTAUXSn31','BTAUXSn32'
 fincl2 = 'PS','TS','PSL'
          ,'PRECT','TMQ'
          ,'PRECC','PRECL'
          ,'LHFLX','SHFLX'
          ,'FSNT','FLNT','FLUT'
          ,'FLNS','FSNS'
          ,'FSNTC','FLNTC'
          ,'TGCLDLWP','TGCLDIWP'
          ,'TUQ','TVQ'
          ,'TBOT:I','QBOT:I','UBOT:I','VBOT:I'
          ,'T900:I','Q900:I','U900:I','V900:I'
          ,'T850:I','Q850:I','U850:I','V850:I'
          ,'Z300:I','Z500:I'
          ,'OMEGA850:I','OMEGA500:I'
          ,'U200:I','V200:I'
 fincl3 = 'PS','TS','PSL'
          ,'T','Q','Z3'
          ,'U','V','OMEGA'
          ,'QRL','QRS'
          ,'CLDLIQ','CLDICE'
          ,'UTGWORO','UTGWSPEC','BUTGWSPEC'

'''

# vars = "FSNTOA,FLUT,FSNT,FLNT,FSNS,FLNS,SHFLX,QFLX,TAUX,TAUY,PRECC,PRECL,PRECSC,PRECSL,TS,TREFHT,CLDTOT_CAL,CLDHGH_CAL,CLDMED_CAL,CLDLOW_CAL,U,ICEFRAC,LANDFRAC,OCNFRAC,AODALL,AODDUST,AODVIS,PS"

# vars = "FSNTOA,FLUT,FSNT,FLNT,FSNS,FLNS,SHFLX,QFLX,TAUX,TAUY,PRECC,
# PRECL,PRECSC,PRECSL,TS,TREFHT,CLDTOT,CLDHGH,CLDMED,CLDLOW,U,ICEFRAC,
# LANDFRAC,OCNFRAC,Mass_bc,Mass_dst,Mass_mom,Mass_ncl,Mass_pom,Mass_so4,Mass_soa,
# AODALL,AODBC,AODDUST,AODPOM,AODSO4,AODSOA,AODSS,AODVIS,PS"
#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------
def get_lnd_nl_opts():
   return f'''
hist_dov2xy = .true.,.true.
hist_fincl1 = 'SNOWDP'
hist_mfilt = 1
hist_nhtfrq = 0
hist_avgflag_pertape = 'A'
check_finidat_year_consistency = .false.
check_dynpft_consistency = .false.
flanduse_timeseries = '{din_loc_root}/lnd/clm2/surfdata_map/landuse.timeseries_ne30np4.pg2_hist_simyr1850-2015_c210113.nc'
fsurdat = \'{land_data_path}/{land_data_file}\'
'''
# finidat = \'{land_init_path}/{land_init_file}\'
#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------
# def get_din_loc_root():
#    (din_loc_root, err) = sp.Popen('./xmlquery DIN_LOC_ROOT    -value', \
#                                  stdout=sp.PIPE, shell=True, \
#                                  universal_newlines=True).communicate()
#    return din_loc_root
#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------
if __name__ == '__main__':

   for n in range(len(nlev_list)):
      # print('-'*80)
      main(nlev=nlev_list[n])
#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------

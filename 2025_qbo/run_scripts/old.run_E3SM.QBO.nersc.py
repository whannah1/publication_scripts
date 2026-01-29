#!/usr/bin/env python
#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------
import os, datetime, subprocess as sp, numpy as np
from shutil import copy2
newcase,config,build,clean,submit,continue_run = False,False,False,False,False,False

acct = 'm3312'    # m3312 / m3305

top_dir  = os.getenv('HOME')+'/E3SM/'
case_dir = top_dir+'Cases/'
src_dir  = top_dir+'E3SM_SRC3/'

# clean        = True
# newcase      = True
# config       = True
# build        = True
submit       = True
# continue_run = True

queue = 'regular'  # regular / debug 

if queue=='debug'  : stop_opt,stop_n,resub,walltime = 'ndays',5, 0,'0:30:00'
if queue=='regular': stop_opt,stop_n,resub,walltime = 'ndays',73*5,10-1,'5:00:00'
# if queue=='regular': stop_opt,stop_n,resub,walltime = 'ndays',73,0,'2:00:00'

compset,ne,npg,grid = 'F2010',30,2,'ne30pg2_EC30to60E2r2'

# nlev= 72;         case='.'.join(['E3SM','QBO-TEST',compset,f'ne{ne}pg{npg}',f'L{nlev}'])
# nlev= 72;         case='.'.join(['E3SM','QBO-TEST',compset,f'ne{ne}pg{npg}',f'L{nlev}-alt'])
# nlev= 72; nsm= 5; case='.'.join(['E3SM','QBO-TEST',compset,f'ne{ne}pg{npg}',f'L{nlev}-nsm{nsm:02d}'])
# nlev= 72; nsm=40; case='.'.join(['E3SM','QBO-TEST',compset,f'ne{ne}pg{npg}',f'L{nlev}-nsm{nsm:02d}'])
# nlev= 72; nsm=40; case='.'.join(['E3SM','QBO-TEST',compset,f'ne{ne}pg{npg}',f'L{nlev}-nsu{nsm:02d}'])
# nlev= 60;         case='.'.join(['E3SM','QBO-TEST',compset,f'ne{ne}pg{npg}',f'L{nlev}'])
# nlev=100;         case='.'.join(['E3SM','QBO-TEST',compset,f'ne{ne}pg{npg}',f'L{nlev}'])
# case += '.00' # initial version - 2021/12/17

# nlev= 72;         case='.'.join(['E3SM','QBO-TEST',compset,f'ne{ne}pg{npg}',f'L{nlev}'])
# nlev= 72; nsm=40; case='.'.join(['E3SM','QBO-TEST',compset,f'ne{ne}pg{npg}',f'L{nlev}-nsu{nsm:02d}']) # smoothing in upper levels
nlev= 80; nsm=40; case='.'.join(['E3SM','QBO-TEST',compset,f'ne{ne}pg{npg}',f'L{nlev}-rsu{nsm:02d}']) # smoothing + refinement

case += '.01' # new version + more output and refined grid - 2021/12/17

### debugging options
# case += '.debug-on'
# case += '.checks-on'

init_scratch = '/global/cscratch1/sd/whannah/HICCUP/data'
if nlev== 60: init_file_atm = f'{init_scratch}/cami_mam3_Linoz_ne30np4_L60_c20211217.nc'
if nlev==100: init_file_atm = f'{init_scratch}/cami_mam3_Linoz_ne30np4_L100_c20211217.nc'

# if 'L72-alt' in case: init_file_atm = f'{init_scratch}/cami_mam3_Linoz_ne30np4_L72-alt_c20211217.nc'
# if 'L72-nsm' in case: init_file_atm = f'{init_scratch}/cami_mam3_Linoz_ne30np4_L72_nsmooth_{nsm:02d}_c20211217.nc'
# if 'L72-nsu' in case: init_file_atm = f'{init_scratch}/cami_mam3_Linoz_ne30np4_L72_nsmooth_{nsm:02d}_upper_lev_only_c20211217.nc'
# if 'L72-rsu' in case: init_file_atm = f'{init_scratch}/cami_mam3_Linoz_ne30np4_L72_nsmooth_40_upper_lev_only_w-refinement_c20211217.nc'
if 'L72-nsu' in case: init_file_atm = f'{init_scratch}/20180716.DECKv1b_A3.ne30_oEC.edison.cam.i.2010-01-01-00000.L72_nsmooth_{nsm:02d}_upper_lev_only_c20220315.nc'
if 'L80-rsu' in case: init_file_atm = f'{init_scratch}/20180716.DECKv1b_A3.ne30_oEC.edison.cam.i.2010-01-01-00000.L80_nsmooth_{nsm:02d}_upper_lev_only_w-refinement_c20220315.nc'


# init_path_lnd = '/global/cscratch1/sd/whannah/e3sm_scratch/init_files/PI_control_v2/archive/rest/0501-01-01-00000'
# init_file_lnd = f'{init_path_lnd}/v2.LR.piControl.elm.r.0501-01-01-00000.nc'

# init_path_lnd = '/global/cscratch1/sd/whannah/e3sm_scratch/init_files/v2.LR.amip_0101/init'
# init_file_lnd = f'{init_path_lnd}/v2.LR.historical_0101.elm.r.1870-01-01-00000.nc'

init_path_lnd = '/global/cscratch1/sd/whannah/e3sm_scratch/init_files/'
init_case_lnd = 'ELM_spinup.ICRUELM.ne30pg2_EC30to60E2r2.20-yr.2010-01-01'
init_file_lnd = f'{init_path_lnd}/{init_case_lnd}.elm.r.2010-01-01-00000.nc'
#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------
print('\n  case : '+case+'\n')

#-------------------------------------------------------------------------------
# Define run command
#-------------------------------------------------------------------------------
# Set up terminal colors
class tcolor:
   ENDC,RED,GREEN,MAGENTA,CYAN = '\033[0m','\033[31m','\033[32m','\033[35m','\033[36m'
def run_cmd(cmd,suppress_output=False,execute=True):
   if suppress_output : cmd = cmd + ' > /dev/null'
   msg = tcolor.GREEN + cmd + tcolor.ENDC
   print(f'\n{msg}')
   if execute: os.system(cmd)
   return
#---------------------------------------------------------------------------------------------------
# Create new case
#---------------------------------------------------------------------------------------------------
if newcase :
   if os.path.isdir(f'{case_dir}/{case}'): exit(tcolor.RED+"\nThis case already exists!\n"+tcolor.ENDC)

   cmd = src_dir+f'cime/scripts/create_newcase -case {case_dir}/{case}'
   cmd += f' -compset {compset} -res {grid} '
   # cmd += ' --pecount 5400x1'
   run_cmd(cmd)

   # Copy this run script into the case directory
   timestamp = datetime.datetime.utcnow().strftime('%Y-%m-%d.%H%M%S')
   run_cmd(f'cp {os.path.realpath(__file__)} {case_dir}/{case}/run_script.{timestamp}.py')
#---------------------------------------------------------------------------------------------------
# Configure
#---------------------------------------------------------------------------------------------------
os.chdir(case_dir+case+'/')
if config : 
   #-------------------------------------------------------
   # if specifying ncdata, do it here to avoid an error message
   if 'init_file_atm' in locals():
      file = open('user_nl_eam','w')
      file.write(f' ncdata = \'{init_file_atm}\' \n')
      file.close()
   #-------------------------------------------------------
   run_cmd(f'./xmlchange --append -id CAM_CONFIG_OPTS -val \" -nlev {nlev} \" ')
   #-------------------------------------------------------
   # Run case setup
   if clean : run_cmd('./case.setup --clean')
   run_cmd('./case.setup --reset')
#---------------------------------------------------------------------------------------------------
# Build
#---------------------------------------------------------------------------------------------------
if build : 
   #-------------------------------------------------------
   # Build the case
   if 'debug-on' in case : run_cmd('./xmlchange -file env_build.xml -id DEBUG -val TRUE ')
   if clean : run_cmd('./case.build --clean')
   run_cmd('./case.build')
#---------------------------------------------------------------------------------------------------
# Write the namelist options and submit the run
#---------------------------------------------------------------------------------------------------
if submit : 
   #-------------------------------------------------------
   # First query some stuff about the case
   #-------------------------------------------------------
   # (din_loc_root   , err) = sp.Popen('./xmlquery DIN_LOC_ROOT    -value', \
   #                                   stdout=sp.PIPE, shell=True).communicate()
   # (config_opts, err) = sp.Popen('./xmlquery CAM_CONFIG_OPTS -value', \
   #                                   stdout=sp.PIPE, shell=True, universal_newlines=True).communicate()
   # config_opts = ' '.join(config_opts.split())   # remove extra spaces to simplify string query
   #-------------------------------------------------------
   # Namelist options
   #-------------------------------------------------------
   nfile = 'user_nl_eam'
   file = open(nfile,'w') 
   #------------------------------
   # Specify history output frequency and variables
   #------------------------------
   file.write(' nhtfrq    = 0,-3,-6 \n')
   file.write(' mfilt     = 1,8,4 \n')
   file.write(" fincl1 = 'Z3'") # this is for easier use of height axis on profile plots
   file.write(          ",'UTGWSPEC','BUTGWSPEC','UTGWORO','BTAUE','BTAUW'")   # fields sugested by Yaga
   file.write(          ",'BUTEND1','BUTEND2','BUTEND3','BUTEND4','BUTEND5'") # Beres U-tendency
   # Beres tau at each phase speed
   file.write(          ",'BTAUXSp00','BTAUXSp01','BTAUXSp02','BTAUXSp03','BTAUXSp04'")
   file.write(                      ",'BTAUXSp05','BTAUXSp06','BTAUXSp07','BTAUXSp08'") 
   file.write(                      ",'BTAUXSp09','BTAUXSp10','BTAUXSp11','BTAUXSp12'")
   file.write(                      ",'BTAUXSp13','BTAUXSp14','BTAUXSp15','BTAUXSp16'")
   file.write(                      ",'BTAUXSp17','BTAUXSp18','BTAUXSp19','BTAUXSp20'")
   file.write(                      ",'BTAUXSp21','BTAUXSp22','BTAUXSp23','BTAUXSp24'")
   file.write(                      ",'BTAUXSp25','BTAUXSp26','BTAUXSp27','BTAUXSp28'")
   file.write(                      ",'BTAUXSp29','BTAUXSp30','BTAUXSp31','BTAUXSp32'")
   file.write(                      ",'BTAUXSn01','BTAUXSn02','BTAUXSn03','BTAUXSn04'")
   file.write(                      ",'BTAUXSn05','BTAUXSn06','BTAUXSn07','BTAUXSn08'")
   file.write(                      ",'BTAUXSn09','BTAUXSn10','BTAUXSn11','BTAUXSn12'")
   file.write(                      ",'BTAUXSn13','BTAUXSn14','BTAUXSn15','BTAUXSn16'")
   file.write(                      ",'BTAUXSn17','BTAUXSn18','BTAUXSn19','BTAUXSn20'")
   file.write(                      ",'BTAUXSn21','BTAUXSn22','BTAUXSn23','BTAUXSn24'")
   file.write(                      ",'BTAUXSn25','BTAUXSn26','BTAUXSn27','BTAUXSn28'")
   file.write(                      ",'BTAUXSn29','BTAUXSn30','BTAUXSn31','BTAUXSn32'")
   file.write('\n')
   file.write(" fincl2 = 'PS','TS','PSL'")
   file.write(          ",'PRECT','TMQ'")
   file.write(          ",'PRECC','PRECL'")
   file.write(          ",'LHFLX','SHFLX'")             # surface fluxes
   file.write(          ",'FSNT','FLNT','FLUT'")        # Net TOM heating rates
   file.write(          ",'FLNS','FSNS'")               # Surface rad for total column heating
   file.write(          ",'FSNTC','FLNTC'")             # clear sky heating rates for CRE
   file.write(          ",'TGCLDLWP','TGCLDIWP'")       # liq & ice water path
   file.write(          ",'TUQ','TVQ'")                 # vapor transport for AR tracking
   # variables for tracking stuff like hurricanes
   file.write(          ",'TBOT:I','QBOT:I','UBOT:I','VBOT:I'") # lowest model leve
   file.write(          ",'T900:I','Q900:I','U900:I','V900:I'") # 900mb data
   file.write(          ",'T850:I','Q850:I','U850:I','V850:I'") # 850mb data
   file.write(          ",'Z300:I','Z500:I'")
   file.write(          ",'OMEGA850:I','OMEGA500:I'")
   file.write(          ",'U200:I','V200:I'")
   file.write('\n')
   # h2 is mainly for calculating TEM 
   file.write(" fincl3 = 'PS','TS','PSL'")
   file.write(          ",'T','Q','Z3'")                      # 3D thermodynamic budget components
   file.write(          ",'U','V','OMEGA'")                    # 3D velocity components
   file.write(          ",'QRL','QRS'")                        # 3D radiative heating profiles
   file.write(          ",'CLDLIQ','CLDICE'")                  # 3D cloud fields
   file.write(          ",'UTGWORO','UTGWSPEC','BUTGWSPEC'")   # gravity wave U tendencies
   file.write('\n')

   # ,'MSKtem','VTH2d','UV2d','UW2d','U2d','V2d','TH2d','W2d' # FV TE terms
   
   #------------------------------
   # Other namelist stuff
   #------------------------------   
   if 'init_file_atm' in locals(): file.write(f' ncdata = \'{init_file_atm}\' \n')

   # file.write(" inithist = \'ENDOFRUN\' \n")

   if 'checks-on' in case : file.write(' state_debug_checks = .true. \n')

   # close atm namelist file
   file.close()

   #-------------------------------------------------------
   # LND namelist
   #-------------------------------------------------------
   if 'init_file_lnd' in locals():
      nfile = 'user_nl_elm'
      file = open(nfile,'w')
      file.write(f' finidat = \'{init_file_lnd}\' \n')
      # file.write(f' check_finidat_fsurdat_consistency = .false. \n')
      file.write(f' fsurdat = \'/global/cfs/cdirs/e3sm/inputdata/lnd/clm2/surfdata_map/surfdata_ne30pg2_simyr2000_c210402.nc\' \n')
      # file.write(f' fsurdat = \'/global/cfs/cdirs/e3sm/inputdata/lnd/clm2/surfdata_map/surfdata_ne30pg2_simyr2010_c210402.nc\' \n')
      file.close()
   #-------------------------------------------------------
   # Set some run-time stuff
   #-------------------------------------------------------
   run_cmd(f'./xmlchange STOP_OPTION={stop_opt},STOP_N={stop_n},RESUBMIT={resub}')
   run_cmd(f'./xmlchange JOB_QUEUE={queue},JOB_WALLCLOCK_TIME={walltime}')
   run_cmd(f'./xmlchange CHARGE_ACCOUNT={acct},PROJECT={acct}')

   if continue_run :
      run_cmd('./xmlchange CONTINUE_RUN=TRUE')   
   else:
      run_cmd('./xmlchange CONTINUE_RUN=FALSE')
   #-------------------------------------------------------
   # Submit the run
   #-------------------------------------------------------
   run_cmd('./case.submit')

#---------------------------------------------------------------------------------------------------
# Print the case name again
#---------------------------------------------------------------------------------------------------
print(f'\n  case : {case}\n') 
#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------

#!/usr/bin/env python3
#---------------------------------------------------------------------------------------------------
class clr:END,RED,GREEN,MAGENTA,CYAN = '\033[0m','\033[31m','\033[32m','\033[35m','\033[36m'
def run_cmd(cmd): print('\n'+clr.GREEN+cmd+clr.END) ; os.system(cmd); return
#---------------------------------------------------------------------------------------------------
import os, subprocess as sp
from optparse import OptionParser
parser = OptionParser()
parser.add_option('--radnx',dest='radnx',default=None,help='sets number of rad columns')
(opts, args) = parser.parse_args()
newcase,config,build,clean,submit,continue_run = False,False,False,False,False,False

acct = 'm3312'    # m3312 / m3305
case_dir = os.getenv('HOME')+'/E3SM/Cases'
src_dir  = os.getenv('HOME')+'/E3SM/E3SM_SRC2' # branch => whannah/mmf/reduced-rad-sensitivity

# clean        = True
# newcase      = True
# config       = True
# build        = True
submit       = True
continue_run = True

# queue = 'batch' # batch / debug

disable_bfb = False

compset,num_nodes,arch = 'F2010-MMF1',64,'GNUGPU'

# stop_opt,stop_n,resub,walltime = 'ndays',11,0,'2:00:00' 

rad_arch_list,rad_nx_list = [],[]
rad_sort_flag = []
gcm_dt_list = []
#rad_nx_list.append(128) ; stop_n,resub,walltime = 73,5*4-1,'6:00:00'; rad_sort_flag.append(False); gcm_dt_list.append(None)
#rad_nx_list.append( 64) ; stop_n,resub,walltime = 73,5*4-1,'6:00:00'; rad_sort_flag.append(False); gcm_dt_list.append(None)
#rad_nx_list.append( 32) ; stop_n,resub,walltime = 73,5*4-1,'6:00:00'; rad_sort_flag.append(False); gcm_dt_list.append(None)
#rad_nx_list.append( 16) ; stop_n,resub,walltime = 73,5*4-1,'6:00:00'; rad_sort_flag.append(False); gcm_dt_list.append(None)
#rad_nx_list.append(  8) ; stop_n,resub,walltime = 73,5*4-1,'6:00:00'; rad_sort_flag.append(False); gcm_dt_list.append(None)
#rad_nx_list.append(  4) ; stop_n,resub,walltime = 73,5*4-1,'6:00:00'; rad_sort_flag.append(False); gcm_dt_list.append(None)
#rad_nx_list.append(  2) ; stop_n,resub,walltime = 73,5*4-1,'6:00:00'; rad_sort_flag.append(False); gcm_dt_list.append(None)
#rad_nx_list.append(  1) ; stop_n,resub,walltime = 73,5*4-1,'6:00:00'; rad_sort_flag.append(False); gcm_dt_list.append(None)

rad_nx_list.append(128) ; rad_sort_flag.append(False); gcm_dt_list.append(None)
rad_nx_list.append( 64) ; rad_sort_flag.append(False); gcm_dt_list.append(None)
rad_nx_list.append( 32) ; rad_sort_flag.append(False); gcm_dt_list.append(None)
rad_nx_list.append( 16) ; rad_sort_flag.append(False); gcm_dt_list.append(None)
rad_nx_list.append(  8) ; rad_sort_flag.append(False); gcm_dt_list.append(None)
rad_nx_list.append(  4) ; rad_sort_flag.append(False); gcm_dt_list.append(None)
rad_nx_list.append(  2) ; rad_sort_flag.append(False); gcm_dt_list.append(None)
rad_nx_list.append(  1) ; rad_sort_flag.append(False); gcm_dt_list.append(None)

# rad_nx_list.append(128) ; stop_n,resub,walltime = 32,5,'3:00:00'; rad_sort_flag.append(True); gcm_dt_list.append(None)
# rad_nx_list.append(  8) ; stop_n,resub,walltime = 32,5,'3:00:00'; rad_sort_flag.append(True); gcm_dt_list.append(None)

# rad_nx_list.append(  8) ; stop_n,resub,walltime = 32,0,'5:00:00'; rad_sort_flag.append(False); gcm_dt_list.append(2*60)

#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------
def main(rad_nx=None,rad_arch='gpu',use_rad_sort=False,gcm_dt=None):

   if rad_nx is None: exit('rad_nx argument not provided?')

   ne,npg         = 30,2
   use_momfb      = True
   use_vt         = True
   crm_dx         = 1000
   crm_nx,crm_ny  = 128,1

   res=f'ne{ne}' if npg==0 else f'ne{ne}pg{npg}'; grid = f'{res}_oECv3' # f'{res}_r05_oECv3' / f'{res}_{res}'

   case_list = ['E3SM','RAD-SENS',arch,grid,compset,f'NXY_{crm_nx}x{crm_ny}',f'RNX_{rad_nx}']
   if use_rad_sort: case_list.append('RAD_SORT')
   if gcm_dt is not None: case_list.append(f'GDT_{gcm_dt}')

   ### batch/version number
   batch_num = '00' # initial runs
   # batch_num = '01' # ???

   case_list.append(batch_num)
   case = '.'.join(case_list)

   # case = case+'.debug-on'

   ### specify land initial condition file
   land_init_path = '/pscratch/sd/w/whannah/e3sm_scratch/init_scratch/'
   if grid=='ne30pg2_r05_oECv3':land_init_file = 'ELM_spinup.ICRUELM.ne30pg2_r05_oECv3.20-yr.2010-01-01.elm.r.2010-01-01-00000.nc'
   if grid=='ne30pg2_oECv3':    land_init_file = 'ELM_spinup.ICRUELM.ne30pg2_oECv3.20-yr.2010-01-01.elm.r.2010-01-01-00000.nc'
   #------------------------------------------------------------------------------------------------
   #------------------------------------------------------------------------------------------------
   print('\n  case : '+case+'\n')

   if 'CPU' in arch : max_mpi_per_node,atm_nthrds  = 64,1
   if 'GPU' in arch : max_mpi_per_node,atm_nthrds  =  4,8
   atm_ntasks = max_mpi_per_node*num_nodes

   if gcm_dt is not None: 
      dtime = gcm_dt; ncpl = int( 86400 / dtime )
   #------------------------------------------------------------------------------------------------
   # Create new case
   #------------------------------------------------------------------------------------------------
   if newcase :
      # Check if directory already exists   
      if os.path.isdir(f'{case_dir}/{case}'): exit(f'\n{clr.RED}This case already exists!{clr.END}\n')
      cmd = f'{src_dir}/cime/scripts/create_newcase -case {case_dir}/{case}'
      cmd = cmd + f' -compset {compset} -res {grid}  '
      if arch=='GNUCPU' : cmd += f' -mach pm-cpu -compiler gnu    -pecount {atm_ntasks}x{atm_nthrds} '
      if arch=='GNUGPU' : cmd += f' -mach pm-gpu -compiler gnugpu -pecount {atm_ntasks}x{atm_nthrds} '
      run_cmd(cmd)
      os.chdir(f'{case_dir}/{case}/')
      run_cmd(f'./xmlchange --file env_mach_pes.xml --id MAX_MPITASKS_PER_NODE --val {max_mpi_per_node} ')
   else:
      os.chdir(f'{case_dir}/{case}/')
   #------------------------------------------------------------------------------------------------
   # Configure
   #------------------------------------------------------------------------------------------------
   if config :
      #-------------------------------------------------------
      # Set some non-default stuff for all MMF experiments here
      if 'MMF' in compset:
         rad_ny = rad_nx if crm_ny>1 else 1
         # run_cmd(f'./xmlchange --append --id CAM_CONFIG_OPTS --val \" -crm_dt {crm_dt} \" ')
         run_cmd(f'./xmlchange --append --id CAM_CONFIG_OPTS --val \" -crm_dx {crm_dx}  \" ')
         run_cmd(f'./xmlchange --append --id CAM_CONFIG_OPTS --val \" -crm_nx {crm_nx} -crm_ny {crm_ny} \" ')
         run_cmd(f'./xmlchange --append --id CAM_CONFIG_OPTS --val \" -crm_nx_rad {rad_nx} -crm_ny_rad {rad_ny} \" ')
         if use_vt: run_cmd(f'./xmlchange --append --id CAM_CONFIG_OPTS --val \" -use_MMF_VT \" ')
      #-------------------------------------------------------
      # Add special MMF options based on case name
      cpp_opt = ''
      if 'debug-on' in case : cpp_opt += ' -DYAKL_DEBUG'

      # if  crm_ny==1: cpp_opt += ' -DMMF_DIR_NS' # no longer needed
      if crm_ny==1 and use_momfb: cpp_opt += ' -DMMF_ESMT -DMMF_USE_ESMT'
      if crm_ny>1  and use_momfb: cpp_opt += ' -DMMF_MOMENTUM_FEEDBACK'

      if use_rad_sort: cpp_opt += ' -DMMF_RAD_SORT'

      if cpp_opt != '' :
         cmd  = f'./xmlchange --append --file env_build.xml --id CAM_CONFIG_OPTS'
         cmd += f' --val \" -cppdefs \' {cpp_opt} \'  \" '
         run_cmd(cmd)
      #-------------------------------------------------------
      # Set tasks for all components
      if ne>=30:
         cmd = './xmlchange --file env_mach_pes.xml '
         alt_ntask = 512; cmd += f'NTASKS_OCN={alt_ntask},NTASKS_ICE={alt_ntask}'
         alt_ntask = 512; cmd += f',NTASKS_LND={alt_ntask},NTASKS_CPL={alt_ntask}'
         alt_ntask = max_mpi_per_node
         cmd += f',NTASKS_ROF={alt_ntask},NTASKS_WAV={alt_ntask},NTASKS_GLC={alt_ntask}'
         cmd += f',NTASKS_ESP=1,NTASKS_IAC=1'
         run_cmd(cmd)
      #-------------------------------------------------------
      # 64_data format is needed for large numbers of columns (GCM or CRM)
      # run_cmd('./xmlchange PIO_NETCDF_FORMAT=\"64bit_data\" ')
      #-------------------------------------------------------
      # Run case setup
      if clean : run_cmd('./case.setup --clean')
      run_cmd('./case.setup --reset')
   #------------------------------------------------------------------------------------------------
   # Build
   #------------------------------------------------------------------------------------------------
   if build : 
      if 'debug-on' in case : run_cmd('./xmlchange --file env_build.xml --id DEBUG --val TRUE ')
      if clean : run_cmd('./case.build --clean')
      run_cmd('./case.build')
   #------------------------------------------------------------------------------------------------
   # Write the namelist options and submit the run
   #------------------------------------------------------------------------------------------------
   if submit : 
      init_scratch = '/global/cfs/cdirs/e3sm/inputdata'
      run_cmd(f'./xmlchange DIN_LOC_ROOT={init_scratch}')
      #-------------------------------------------------------
      # ATM namelist
      #-------------------------------------------------------
      nfile = 'user_nl_eam'
      file = open(nfile,'w') 
      #------------------------------
      # history output variables
      #------------------------------
      if batch_num=='00':
         if grid=='ne4pg2_ne4pg2':
            file.write(' nhtfrq = 0,-1,-24 \n')
            file.write(' mfilt  = 1,24,1 \n')
         else:
            file.write(' nhtfrq = 0,-3,-24 \n')
            file.write(' mfilt  = 1, 8,1 \n')
         file.write(" fincl1 = 'Z3'") # this is for easier use of height axis on profile plots   
         file.write(         ",'CLOUD','CLDLIQ','CLDICE'")
         file.write('\n')
         file.write(" fincl2 = 'PS','TS','PSL'")
         file.write(          ",'PRECT','TMQ'")
         file.write(          ",'LHFLX','SHFLX'")             # surface fluxes
         file.write(          ",'FSNT','FLNT'")               # Net TOM heating rates
         file.write(          ",'FLNS','FSNS'")               # Surface rad for total column heating
         file.write(          ",'FSNTC','FLNTC'")             # clear sky heating rates for CRE
         file.write(          ",'TGCLDLWP','TGCLDIWP'")
         file.write(          ",'TAUX','TAUY'")               # surface stress
         file.write(          ",'TUQ','TVQ'")                         # vapor transport
         file.write(          ",'TBOT:I','QBOT:I','UBOT:I','VBOT:I'") # lowest model leve
         file.write(          ",'T900:I','Q900:I','U900:I','V900:I'") # 900mb data
         file.write(          ",'T850:I','Q850:I','U850:I','V850:I'") # 850mb data
         file.write(          ",'Z300:I','Z500:I'")
         file.write(          ",'OMEGA850:I','OMEGA500:I'")
         file.write('\n')
         if grid=='ne4pg2_ne4pg2':
            file.write(" fincl3 =  'PS','PSL'")
            file.write(          ",'T','Q','Z3' ")               # 3D thermodynamic budget components
            file.write(          ",'U','V','OMEGA'")             # 3D velocity components
            file.write(          ",'CLOUD','CLDLIQ','CLDICE'")   # 3D cloud fields
            file.write(          ",'QRL','QRS'")                 # 3D radiative heating profiles
            file.write('\n')
      #------------------------------
      # Other namelist stuff
      #------------------------------
      if gcm_dt is not None:
         if gcm_dt < 60 :
            file.write(f'dt_tracer_factor = 1 \n')
            file.write(f'dt_remap_factor = 1 \n')
            file.write(f'se_tstep = {dtime} \n')
            file.write(f'hypervis_subcycle_q = 1 \n')
         if gcm_dt == 1*60 :
            file.write(f'dt_tracer_factor = 1 \n')
            file.write(f'dt_remap_factor = 1 \n')
            file.write(f'se_tstep = 60 \n')
            file.write(f'hypervis_subcycle_q = 1 \n')
         if gcm_dt == 2*60 :
            file.write(f'dt_tracer_factor = 1 \n')
            file.write(f'dt_remap_factor = 1 \n')
            file.write(f'se_tstep = 60 \n')
            file.write(f'hypervis_subcycle_q = 1 \n')
         if gcm_dt == 5*60 :
            file.write(f'dt_tracer_factor = 5 \n')
            file.write(f'dt_remap_factor = 1 \n')
            file.write(f'se_tstep = 60 \n')
            file.write(f'hypervis_subcycle_q = 5 \n')
         if gcm_dt == 10*60 :
            file.write(f'dt_tracer_factor = 5 \n')
            file.write(f'dt_remap_factor = 1 \n')
            file.write(f'se_tstep = 60 \n')
            file.write(f'hypervis_subcycle_q = 5 \n')

      ### limit dynamics tasks
      # num_dyn = ne*ne*6
      # ntask_dyn = int( num_dyn / atm_nthrds )
      # file.write(f' dyn_npes = {ntask_dyn} \n')
      # file.write(" inithist = \'ENDOFRUN\' \n") # ENDOFRUN / NONE
      file.close()
      #-------------------------------------------------------
      # ELM namelist
      #-------------------------------------------------------
      nfile = 'user_nl_elm'
      file = open(nfile,'w')
      if ne==30 and npg==2: file.write(f' fsurdat = \'{init_scratch}/lnd/clm2/surfdata_map/surfdata_ne30pg2_simyr2010_c210402.nc\' \n')
      if 'land_init_file' in locals(): file.write(f' finidat = \'{land_init_path}/{land_init_file}\' \n')
      file.close()
      #-------------------------------------------------------
      # Set some run-time stuff
      #-------------------------------------------------------
      if gcm_dt is not None: run_cmd(f'./xmlchange ATM_NCPL={ncpl}')

      if 'queue'    in globals(): run_cmd(f'./xmlchange JOB_QUEUE={queue}')
      if 'stop_opt' in globals(): run_cmd(f'./xmlchange STOP_OPTION={stop_opt}')
      if 'stop_n'   in globals(): run_cmd(f'./xmlchange STOP_N={stop_n}')
      if 'resub'    in globals(): run_cmd(f'./xmlchange RESUBMIT={resub}')
      if 'walltime' in globals(): run_cmd(f'./xmlchange JOB_WALLCLOCK_TIME={walltime}')
      if 'GPU' in arch: run_cmd(f'./xmlchange CHARGE_ACCOUNT={acct}_g,PROJECT={acct}_g')
      if 'CPU' in arch: run_cmd(f'./xmlchange CHARGE_ACCOUNT={acct},PROJECT={acct}')
      
      if     disable_bfb: run_cmd('./xmlchange BFBFLAG=FALSE')
      if not disable_bfb: run_cmd('./xmlchange BFBFLAG=TRUE')

      if     continue_run: run_cmd('./xmlchange --file env_run.xml CONTINUE_RUN=TRUE ')   
      if not continue_run: run_cmd('./xmlchange --file env_run.xml CONTINUE_RUN=FALSE ')
      #-------------------------------------------------------
      # Submit the run
      #-------------------------------------------------------
      run_cmd('./case.submit')
   #------------------------------------------------------------------------------------------------
   # Print the case name again
   #------------------------------------------------------------------------------------------------
   print(f'\n  case : {case}\n')
#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------
if __name__ == '__main__':

   for n in range(len(rad_nx_list)):
      print('-'*80)
      main(rad_nx=rad_nx_list[n],use_rad_sort=rad_sort_flag[n],gcm_dt=gcm_dt_list[n])


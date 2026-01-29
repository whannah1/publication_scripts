import os

cases = []

git_hash = 'cbe53b'
cases.append(f'E3SM.PGVAL.ne30_r05_oECv3.F2010SC5-CMIP6.master-{git_hash}')        
cases.append(f'E3SM.PGVAL.ne30pg2_r05_oECv3.F2010SC5-CMIP6.master-{git_hash}')     
cases.append(f'E3SM.PGVAL.ne30pg3_r05_oECv3.F2010SC5-CMIP6.master-{git_hash}')     
cases.append(f'E3SM.PGVAL.ne30pg4_r05_oECv3.F2010SC5-CMIP6.master-{git_hash}')     
cases.append(f'E3SM.PGVAL.conusx4v1_r05_oECv3.F2010SC5-CMIP6.master-{git_hash}')   
cases.append(f'E3SM.PGVAL.conusx4v1pg2_r05_oECv3.F2010SC5-CMIP6.master-{git_hash}')

git_hash = 'efa06a'
ntask = 1350
cases.append(f'E3SM.PGVAL-TIMING.ne30_ne30.F-EAM-AQP1.master-{git_hash}.ntask_{ntask}.nthrds_2')        
cases.append(f'E3SM.PGVAL-TIMING.ne30pg2_ne30pg2.F-EAM-AQP1.master-{git_hash}.ntask_{ntask}.nthrds_2')  
cases.append(f'E3SM.PGVAL-TIMING.ne30pg3_ne30pg3.F-EAM-AQP1.master-{git_hash}.ntask_{ntask}.nthrds_2')  
cases.append(f'E3SM.PGVAL-TIMING.ne30pg4_ne30pg4.F-EAM-AQP1.master-{git_hash}.ntask_{ntask}.nthrds_2')  
ntask = 2700
cases.append(f'E3SM.PGVAL-TIMING.ne30_ne30.F-EAM-AQP1.master-{git_hash}.ntask_{ntask}.nthrds_2')        
cases.append(f'E3SM.PGVAL-TIMING.ne30pg2_ne30pg2.F-EAM-AQP1.master-{git_hash}.ntask_{ntask}.nthrds_2')  
cases.append(f'E3SM.PGVAL-TIMING.ne30pg3_ne30pg3.F-EAM-AQP1.master-{git_hash}.ntask_{ntask}.nthrds_2')  
cases.append(f'E3SM.PGVAL-TIMING.ne30pg4_ne30pg4.F-EAM-AQP1.master-{git_hash}.ntask_{ntask}.nthrds_2')  
ntask = 5400
cases.append(f'E3SM.PGVAL-TIMING.ne30_ne30.F-EAM-AQP1.master-{git_hash}.ntask_{ntask}.nthrds_2')       
cases.append(f'E3SM.PGVAL-TIMING.ne30pg2_ne30pg2.F-EAM-AQP1.master-{git_hash}.ntask_{ntask}.nthrds_2') 
cases.append(f'E3SM.PGVAL-TIMING.ne30pg3_ne30pg3.F-EAM-AQP1.master-{git_hash}.ntask_{ntask}.nthrds_2') 
cases.append(f'E3SM.PGVAL-TIMING.ne30pg4_ne30pg4.F-EAM-AQP1.master-{git_hash}.ntask_{ntask}.nthrds_2') 


for case in cases:

  timing_dir = os.getenv('HOME')+f'/E3SM/Cases/{case}/timing'
  
  new_dir = f'/global/homes/w/whannah/Research/E3SM/pub_figs/2021_JAMES_physgrid/timing_data/{case}'

  if not os.path.isdir(new_dir): os.mkdir(new_dir)

  cmd = f'cp {timing_dir}/*  {new_dir}/ '
  
  print()
  print(cmd)
  os.system(cmd)
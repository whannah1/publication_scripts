# plot the fractional occurence of all possible patterns
# v1 - calculate the occurrence of each pattern
# v2 - identify patterns while retaining time dimension
# v3 - special case for high short term simulation w/ high-freq output
import os, ngl, sys, numba, copy, xarray as xr, numpy as np, string
import hapy_common as hc, hapy_E3SM as he, hapy_setres as hs
import pg_checkerboard_utilities as pg

case,name,clr,dsh,pat = [],[],[],[],[]
def add_case(case_in,n='',c='black',d=0,p=0):
   global name,case
   case.append(case_in); name.append(n); clr.append(c); dsh.append(d); pat.append(p)
var,lev_list = [],[]
def add_var(var_name,lev=-1): var.append(var); lev_list.append(lev)
#-------------------------------------------------------------------------------
# add_case('E3SM.VTVAL-HC.GNUGPU.ne30pg2_r05_oECv3.F-MMFXX.MOMFB.00',           n='E3SM-MMF')
# add_case('E3SM.VTVAL-HC.GNUGPU.ne30pg2_r05_oECv3.F-MMFXX.MOMFB.MOMVT.VT_0.00',n='E3SM-MMF+BVT')
# add_case('E3SM.VTVAL-HC.GNUGPU.ne30pg2_r05_oECv3.F-MMFXX.MOMFB.MOMVT.VT_1.00',n='E3SM-MMF+FVT')
# add_case('E3SM.VTVAL-HC.GNUGPU.ne30pg2_r05_oECv3.F-MMFXX.MOMFB.VT_0.00',n='E3SM-MMF + BVT')
# add_case('E3SM.VTVAL-HC.GNUGPU.ne30pg2_r05_oECv3.F-MMFXX.MOMFB.VT_1.00',n='E3SM-MMF+FVT1')
# add_case('E3SM.VTVAL-HC.GNUGPU.ne30pg2_r05_oECv3.F-MMFXX.MOMFB.VT_2.00',n='E3SM-MMF+FVT2')
# add_case('E3SM.VTVAL-HC.GNUGPU.ne30pg2_r05_oECv3.F-MMFXX.MOMFB.VT_4.00',n='E3SM-MMF+FVT4')
# add_case('E3SM.VTVAL-HC.GNUGPU.ne30pg2_r05_oECv3.F-MMFXX.MOMFB.VT_8.00',n='E3SM-MMF+FVT8')

add_case('E3SM.VTVAL-HC.GNUGPU.ne30pg2_r05_oECv3.F-MMFXX.MOMFB.MOMVT.VT_0.TEST_QU.00',n='MMF+BVT only Q+U')
add_case('E3SM.VTVAL-HC.GNUGPU.ne30pg2_r05_oECv3.F-MMFXX.MOMFB.MOMVT.VT_0.TEST_TQ.00',n='MMF+BVT only T+Q')
# add_case('E3SM.VTVAL-HC.GNUGPU.ne30pg2_r05_oECv3.F-MMFXX.MOMFB.MOMVT.VT_0.TEST_TU.00',n='MMF+BVT only T+U')
# add_case('E3SM.VTVAL-HC.GNUGPU.ne30pg2_r05_oECv3.F-MMFXX.MOMFB.MOMVT.VT_0.TEST_T.00', n='MMF+BVT only T'  )
# add_case('E3SM.VTVAL-HC.GNUGPU.ne30pg2_r05_oECv3.F-MMFXX.MOMFB.MOMVT.VT_0.TEST_Q.00', n='MMF+BVT only Q'  )
# add_case('E3SM.VTVAL-HC.GNUGPU.ne30pg2_r05_oECv3.F-MMFXX.MOMFB.MOMVT.VT_0.TEST_U.00', n='MMF+BVT only U'  )
#-------------------------------------------------------------------------------
var = []
var.append('TGCLDLWP')
var.append('PRECT')
# var.append('TGCLDIWP')
# var.append('Q850')
# var.append('T850')
# var.append('U850')
# var.append('TMQ')
# var.append('LHFLX')
# var.append('FLNT')
# var.append('FSNT')

htype = 'h1' ; num_files = 10 # use all data from short term cases

print_stats       = False   # print stats for debugging/sanity check

subset_min_length = 4       # threshold for partial checkerboard

#---------------------------------------------------------------------------------------------------
# generate sets of possible neighbor states
#---------------------------------------------------------------------------------------------------
num_neighbors = 8

### only pure checkboard
# tmp_file='data/occurrence-chx-only'; rotate_sets=False; sets=pg.chx_only_sets

### full set of all possible patterns
# scratch = '/global/cscratch1/sd/whannah/e3sm_scratch/cori-knl/chx-detection-data'
scratch = '/gpfs/alpine/scratch/hannah6/cli115/analysis_files/vt_validation/chx-detection-data'
tmp_file=f'{scratch}/occurence-partial-chx-over-time'; rotate_sets=True; sets=pg.all_possible_sets

(num_set,set_len) = sets.shape
sets.assign_coords(coords={'set':np.arange(num_set),'neighbors':np.arange(set_len)})

num_case,num_var,num_set = len(case),len(var),len(sets)


#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------
cnt_ds_list = []
for v in range(num_var) :
   print('  var: '+hc.tcolor.GREEN+f'{var[v]:10}'+hc.tcolor.ENDC)
   for c in range(num_case):
      
      #-------------------------------------------------------------------------
      #-------------------------------------------------------------------------

      case_tmp_file = f'{tmp_file}.high-freq.{case[c]}.{var[v]}.sml_{subset_min_length}.nc'

      print('\n    case: '+hc.tcolor.CYAN+case[c]+hc.tcolor.ENDC+'    '+case_tmp_file)
      #-------------------------------------------------------------------------
      #-------------------------------------------------------------------------
      case_obj = he.Case( name=case[c] )
      if 'lev' not in vars() : lev = np.array([-1])
      comp = 'eam'
      use_remap = False
      remap_str=f'remap_ne30pg2'
      #-------------------------------------------------------------------------
      # Load data
      #-------------------------------------------------------------------------

      print('      loading data...')
      data = case_obj.load_data(var[v],htype=htype,lev=lev,component=comp,
                                num_files=num_files,first_file=0,
                                use_remap=use_remap,remap_str=remap_str)

      ### Convert to daily mean
      # data = data.resample(time='D').mean(dim='time')

      if print_stats: hc.print_stat(data,name=var[v],stat='naxsh',indent='    ')
      #-------------------------------------------------------------------------
      #-------------------------------------------------------------------------
      if 'FV' in case[c]: 
         scripfile_path = 'scrip_files/cmip6_180x360_scrip.20210722.nc'
      else:
         scripfile_path = 'scrip_files/ne30pg2_scrip.nc'
      
      # exit(data)

      ### find neighbors and count the occurrence of each set
      cnt_ds = pg.find_neighbors_and_count_sets(data, scripfile_path, sets, rotate_sets, 
                                                keep_time=True, verbose=True)
      
      # print()
      # print(cnt_ds)
      # print()

      #-------------------------------------------------------------------------
      #-------------------------------------------------------------------------
      ### combine sets that contain partial checkerboard

   
      chk_idx,nox_idx = [],[]
      for s in range(num_set):
         if pg.is_partial_checkerboard(sets[s,:],subset_length=subset_min_length):
            chk_idx.append(s)
         else:
            nox_idx.append(s)

      ds_sum_nox = cnt_ds.isel(set=nox_idx).sum(dim='set')
      ds_sum_chk = cnt_ds.isel(set=chk_idx).sum(dim='set')

      ds_sum_nox = ds_sum_nox.expand_dims('set').assign_coords(coords={'set':np.array([0])})
      ds_sum_chk = ds_sum_chk.expand_dims('set').assign_coords(coords={'set':np.array([1])})

      cnt_ds = xr.concat( [ ds_sum_nox, ds_sum_chk ], 'set' )

      # print()
      # print(cnt_ds)
      # print()
      # exit()

      # cnt_ds_list = cnt_ds_list_tmp
      # set_labels = ['no checkerboard','partial checkerboard']
      # num_set = len(set_ladbels)

      #-------------------------------------------------------------------------
      #-------------------------------------------------------------------------
      cnt_ds.to_netcdf(path=case_tmp_file,mode='w')

#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------
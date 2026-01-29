# plot the fractional occurence of all possible patterns
import os, ngl, sys, numba, copy, xarray as xr, numpy as np, re
import hapy_common as hc, hapy_E3SM as he, hapy_setres as hs
import pg_checkerboard_utilities as pg

case,name,clr,dsh,pat = [],[],[],[],[]
def add_case(case_in,n='',c='black',d=0,p=0):
   global name,case
   case.append(case_in); name.append(n); clr.append(c); dsh.append(d); pat.append(p)
var,lev_list = [],[]
def add_var(var_name,lev=-1): var.append(var); lev_list.append(lev)
#-------------------------------------------------------------------------------
# add_case('ERA5-FV', n='ERA5 1x1 degree'                                      ,c='gray', p=0)
# add_case('ERA5-PG', n='ERA5 ne30pg2'                                         ,c='gray', p=0)
add_case('MAC-FV',  n='MAC 1x1 deg'                                          ,c='black',p=0)
add_case('MAC-PG',  n='MAC ne30pg2'                                          ,c='gray',p=0)
add_case('E3SM.CHX.ne30pg2_ne30pg2.FC5AV1C-L.no_dcape.00',n='E3SM'           ,c='green',p=0)
# add_case('E3SM.CHX.ne30pg2_ne30pg2.FC5AV1C-L.00',         n='E3SM'           ,c='red'  ,p=0)
add_case('E3SM.CHX.ne30pg2_ne30pg2.F-MMF1.00',            n='E3SM-MMF'       ,c='blue' ,p=0)

# Obs data is here: /global/cscratch1/sd/whannah/Obs/

#-------------------------------------------------------------------------------

var = []
# var.append('PRECT')
var.append('TGCLDLWP')
# var.append('TGCLDIWP')
# var.append('TMQ')
# var.append('LHFLX')
# var.append('FLNT')
# var.append('FSNT')
# var.append('Q850')
# var.append('T850')
# var.append('U850')


# use_regional_subset = False
# lat1,lat2,lon1,lon2 = -20,40,80,180+30


fig_type = 'png'
fig_file = 'figs/F04-chx-occurrence'
tmp_file = 'data/chx-occurrence'

# htype,first_file,num_files = 'h1',0,31
htype,first_file,num_files = 'h1',0,365*5

# ERA5 and MAC daily data comes in monthly files, 
# so we need a different way to specify the number of files
num_files_obs = 12*5


recalculate       = False


use_daily         = True
convert_to_freq   = True   # convert count to fractional occurrence

show_as_diff      = False   # calculate difference from first case
show_as_ratio     = False
add_diff_plot     = True

ocean_only        = False   # ocean only to compare with obs

subset_min_length = 4      # threshold for partial checkerboard when combine_mode=1

combine_sets = True
combine_mode = 2     # 1 = chk vs no chk / 2 = continuity metric
alt_chk_label = True

var_x_case        = False   # controls plot panel arrangement
print_stats       = False   #

sort_sets = True

#---------------------------------------------------------------------------------------------------
# generate sets of possible neighbor states
#---------------------------------------------------------------------------------------------------
num_neighbors = 8

if combine_sets: alt_chk_label = False

# rotate_sets = True
# sets = xr.DataArray(np.array([[0,1,0,1,0,1,0,1]]),dims=['set','neighbors'] )

# rotate_sets = False
# sets = xr.DataArray(np.array([[0,1,0,1,0,1,0,1], [1,0,1,0,1,0,1,0]]),dims=['set','neighbors'] )

### full set of patterns
rotate_sets = True; sets = pg.all_possible_sets

# use for plotting sum of sets
dummy_set = xr.DataArray(np.array([[-1]*num_neighbors]),dims=['set','neighbors'] )

(num_set,set_len) = sets.shape
sets.assign_coords(coords={'set':np.arange(num_set),'neighbors':np.arange(set_len)})

set_labels = pg.get_set_labels(sets)

#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------
def get_comp(case):
   comp = 'eam'
   if 'ERA5' in case: comp = None
   if 'MAC'  in case: comp = None
   if 'CESM' in case: comp = 'cam'
   return comp
#---------------------------------------------------------------------------------------------------
# sort sets according to partial checkerboard
#---------------------------------------------------------------------------------------------------

sort_idx = [-1]*num_set
if sort_sets:
   ### reorder sets to put partial checkerboards on one side
   nox_sets = []
   chk_sets = []
   for s in range(num_set): 
      if pg.is_partial_checkerboard(sets[s,:].values,subset_length=subset_min_length):
         chk_sets.append(sets[s,:].values)
      else:
         nox_sets.append(sets[s,:].values)
   # for s in range(len(nox_sets)): print(f'  {nox_sets[s]}')
   sets_sorted = xr.DataArray(np.array( nox_sets + chk_sets), dims=['set','neighbors'] )
   for s in range(num_set):
      for ss in range(num_set):
         if np.all(sets[s,:].values==sets_sorted[ss,:].values): 
            sort_idx[ss] = s
   ### useful for checking sort_idx
   # for s in range(num_set): print(f'  {set_labels[s]}    {set_labels[sort_idx[s]]}')
   set_labels_sorted = []
   for s in range(num_set):
      set_labels_sorted.append( set_labels[sort_idx[s]] )
   set_labels = set_labels_sorted
else:
   for s in range(num_set): sort_idx[s] = s

if alt_chk_label:
   found_first = False; pchk_start = -1
   for s in range(len(sets)):
      ss = sort_idx[s]
      if pg.is_partial_checkerboard(sets[ss,:],subset_length=subset_min_length): 
         if not found_first: 
            pchk_start = s
            found_first = True; break
else:
   pchk_start = num_set 


#---------------------------------------------------------------------------------------------------
# Set up plot resources
#---------------------------------------------------------------------------------------------------
num_case,num_var,num_set = len(case),len(var),len(sets)

wkres = ngl.Resources()
# npix=2**13; wkres.wkWidth,wkres.wkHeight=npix,npix # use this for plotting all patterns w/ rotation
wks = ngl.open_wks(fig_type,fig_file,wkres)
# plot = [None]*(num_var)
res = hs.res_xy()
res.vpHeightF = 0.4
if convert_to_freq: res.tiYAxisString = 'Fractional Occurrence'

if combine_sets==True and combine_mode==2: res.tiXAxisString = 'Number of Local Extrema in Neighborhood'

pgres = ngl.Resources()
# pgres.nglDraw,pgres.nglFrame = True,False
pgres.nglDraw,pgres.nglFrame = False,False

#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------
cnt_ds_list = []
for c in range(num_case):
   print('\n  case: '+hc.tcolor.CYAN+case[c]+hc.tcolor.ENDC)

   tname = case[c]
   if 'MAC'  in case[c]: tname = 'MAC'
   if 'ERA5' in case[c]: tname = 'ERA5'
   case_obj = he.Case( name=tname, time_freq='daily' )
   if 'lev' not in vars() : lev = np.array([-1])

   comp = 'eam'
   if case[c]=='EAR5': comp = None
   if case[c]=='E3SM.RGMA.ne30pg2_r05_oECv3.FC5AV1C-L.00': comp = 'cam'

   use_remap = False
   remap_str=f'remap_ne30pg2'
   if case[c]=='ERA5-FV': use_remap = False
   if case[c]=='ERA5-PG': use_remap = True
   if case[c]=='MAC-FV' : use_remap = False
   if case[c]=='MAC-PG' : use_remap = True

   #----------------------------------------------------------------------------
   #----------------------------------------------------------------------------
   for v in range(num_var) :
      if use_daily:
         case_tmp_file = f'{tmp_file}.daily.{case[c]}.{var[v]}.nc'
      else:
         case_tmp_file = f'{tmp_file}.{case[c]}.{var[v]}.nc'
      print('    var: '+hc.tcolor.GREEN+f'{var[v]:10}'+hc.tcolor.ENDC+'    '+case_tmp_file)
      if recalculate :
         if v==0: print(hc.tcolor.YELLOW+'    recalculating pattern occurrence...'+hc.tcolor.ENDC)
         #----------------------------------------------------------------------
         # Load data
         #----------------------------------------------------------------------
         num_files_tmp = num_files

         if 'MAC' in case[c]: num_files_tmp  = num_files_obs
         if 'ERA5' in case[c]: num_files_tmp = num_files_obs
            
         if case[c] in ['MAC-FV']: case_obj.file_path_daily = case_obj.data_dir+case_obj.name+'/daily/*'
         if case[c] in ['MAC-PG']: case_obj.file_path_daily = case_obj.data_dir+case_obj.name+'/daily_ne30pg2/*'
      
         if case[c] in ['ERA5-FV']: case_obj.file_path_daily = case_obj.data_dir+case_obj.name+'/daily/*'
         if case[c] in ['ERA5-PG']: case_obj.file_path_daily = case_obj.data_dir+case_obj.name+'/daily_ne30pg2/*'

         print('      loading data...')
         data = case_obj.load_data(var[v],htype=htype,lev=lev,component=comp,
                                   num_files=num_files_tmp,first_file=first_file,
                                   use_remap=use_remap,remap_str=remap_str)

         # reshape 2D data to be 1D
         if 'FV' in case[c]: 
            # lat,lon = data['lat'],data['lon']
            data = data.stack(ncol=("lat", "lon"))
            data['ncol'] = np.arange(len(data['ncol']))

         # Convert to daily mean
         if use_daily: data = data.resample(time='D').mean(dim='time')

         # print(hc.tcolor.RED+'\nWARNING - truncating the data for debugging purposes!!!! - WARNING\n'+hc.tcolor.ENDC)
         # data = data.isel(time=slice(0,4))

         if print_stats: hc.print_stat(data,name=var[v],stat='naxsh',indent='    ')

         # # for debugging - mask out everything except areas where ice is a problem
         # print();print(data);print()
         # if 'FV' in case[c]: 
         #    lat2D = xr.DataArray(np.repeat( lat.values[...,None],len(lon),axis=1),dims=['lat','lon'])
         #    lon2D = xr.DataArray(np.repeat( lon.values[...,None],len(lat),axis=1),dims=['lat','lon'])
         #    lat = lat2D.stack(ncol=("lat", "lon"))
         #    lon = lon2D.stack(ncol=("lat", "lon"))
         # data.load()
         # mask = xr.DataArray( np.ones(len(data['ncol']),dtype=bool), coords=[data['ncol']] )
         # mask = mask & (lat.values>=88) #& (lon>=lon1) & (lon<=lon2)
         # data = data.where(mask,drop=True)

         # print();print();print()
         # print();print(data);print()
         # hc.print_stat(data,name=var[v],stat='naxsh',indent='    ')
         # exit()

         #----------------------------------------------------------------------
         #----------------------------------------------------------------------
         if 'FV' in case[c]: 
            scripfile_path = 'scrip_files/cmip6_180x360_scrip.20210722.nc'
         else:
            scripfile_path = 'scrip_files/ne30pg2_scrip.nc'
         
         cnt_ds = pg.find_neighbors_and_count_sets(data, scripfile_path, sets, rotate_sets, verbose=True)

         print('      writing to file: '+case_tmp_file)
         cnt_ds.to_netcdf(path=case_tmp_file,mode='w')
      else:
         cnt_ds = xr.open_dataset( case_tmp_file )

      #-------------------------------------------------------------------------
      # print the time length
      #-------------------------------------------------------------------------
      ntime = cnt_ds.num_time.values
      
      msg = f'    num_time: {ntime}'
      if use_daily:
         msg += ' days'
         nyear = ntime/365
         msg += f' ({nyear} years)'
      else:
         if htype=='h1': msg += ' x1hr'
         if htype=='h2': msg += ' x3hr'
      print(msg)

      #-------------------------------------------------------------------------
      # Ocean mask
      #-------------------------------------------------------------------------

      if ocean_only: 
         if 'FV' in case[c]: 
            landfrac_ds = xr.open_dataset('data/ERA5_land-sea-mask_180x360.nc').isel(time=0)
            landfrac_ds = landfrac_ds.stack(ncol=('latitude', 'longitude'))
            landfrac_ds['ncol'] = np.arange(len(landfrac_ds['ncol']))
            landfrac = landfrac_ds['lsm']

            # costal_cnt = 0
            # for i in range(len(cnt_ds['ncol'])):
            #    if landfrac[i]<0.5:
            #       if any( landfrac[cnt_ds['neighbors'][i,:]]>0.5 ) : costal_cnt += 1

            # num_col = len(cnt_ds['ncol'])
            # coastal_frac = costal_cnt/num_col
            # print(f'total ncol   : {num_col}')
            # print(f'coastal count: {costal_cnt}')
            # print(f'coastal frac : {coastal_frac}')
            # exit()
            
         else:
            landfrac_ds = xr.open_dataset('data/E3SM_landfrac_ne30pg2.nc').isel(time=0)
            landfrac = landfrac_ds['LANDFRAC']
         mask = xr.DataArray( np.ones(len(cnt_ds['ncol']),dtype=bool), coords=[cnt_ds['ncol']] )
         mask = mask & (landfrac.values<0.5)
         num_time = cnt_ds['num_time']
         cnt_ds = cnt_ds['cnt'].where(mask,drop=True).to_dataset()
         cnt_ds['num_time'] = num_time
      #-------------------------------------------------------------------------
      # apply regional subset
      #-------------------------------------------------------------------------
      # if use_regional_subset: 
      #    lat = case_obj.load_data('lat',htype=htype,num_files=1,component=comp,use_remap=use_remap,remap_str=remap_str)
      #    lon = case_obj.load_data('lon',htype=htype,num_files=1,component=comp,use_remap=use_remap,remap_str=remap_str)

      #    mask = xr.DataArray( np.ones(len(cnt_ds['ncol']),dtype=bool), coords=[cnt_ds['ncol']] )
      #    mask = mask & (lat>=lat1) & (lat<=lat2)
      #    mask = mask & (lon>=lon1) & (lon<=lon2)
         
      #    num_time = cnt_ds['num_time']
      #    cnt_ds = cnt_ds['cnt'].where(mask,drop=True).to_dataset()
      #    cnt_ds['num_time'] = num_time

      #-------------------------------------------------------------------------
      ### debugging stuff - print some diagnostic stuff to explore problem with NaN detection
      #-------------------------------------------------------------------------
      # ratio = xr.DataArray(np.ma.masked_invalid( ( cnt_ds['cnt'].sum(dim='set') / cnt_ds['num_valid'] ).values ))
      
      # print()
      # # print( ( cnt_ds['cnt']                /     cnt_ds['num_valid']  / len(cnt_ds['ncol']) ).sum().values )
      # # print( ( cnt_ds['cnt'].sum(dim='set') /     cnt_ds['num_valid']  / len(cnt_ds['ncol']) ).sum().values )
      # tot_ncol = len(cnt_ds['ncol'])
      # cnt_sum  = cnt_ds['cnt'].sum().values
      # tot_val  = np.sum(np.ma.masked_invalid(cnt_ds['num_valid'].values))
      # num_time = cnt_ds['num_time'].values
      # print(f'num_set : {num_set}' )
      # print(f'num_time: {num_time}' )
      # print(f'ncol    : {tot_ncol}' )
      # print(f'cnt_sum : {cnt_sum}' )
      # print(f'tot_val : {tot_val}'  )

      # ratio1 = cnt_sum / tot_ncol / cnt_ds['num_time'].values
      # ratio2 = ( cnt_ds['cnt'] / tot_val ).sum().values
      # print(f'ratio1 : {ratio1}')
      # print(f'ratio2 : {ratio2}')
      # exit()

      # ratio = ( cnt_ds['cnt'] / cnt_ds['num_valid'] ).sum(dim='set')
      # ratio = cnt_ds['cnt'].sum(dim='set') / cnt_ds['num_valid']
      # ratio = np.ma.masked_invalid( ratio.values )
      # print()
      # cnt = 0
      # for n in range(len(ratio)):
      #    if ratio[n]>2.:
      #       tmp_count = cnt_ds['cnt'][n].sum(dim='set').values
      #       tmp_valid = cnt_ds['num_valid'][n].values
      #       print(f'  ratio: {ratio[n].values:8.1f} cnt: {tmp_count:8.1f}  num valid: {tmp_valid:8.1f}')
      #       cnt += 1
      #       if cnt>10: break
      # print()
      # hc.print_stat( cnt_ds['cnt'].sum(dim='set') , name='count sum (per column)', compact=True )
      # hc.print_stat( cnt_ds['num_valid'] ,          name='num_valid (per column)', compact=True )
      # hc.print_stat( ratio ,                        name='ratio     (per column)', compact=True )
      # print()
      # exit()

      
      #-------------------------------------------------------------------------
      # Plot data on map for debugging
      #-------------------------------------------------------------------------
      # tmp_data = cnt_ds['cnt'].sum(dim='set') #/ cnt_ds['num_valid']
      # tmp_data = xr.DataArray(np.ma.masked_invalid(tmp_data.values))

      # landfrac_ds = xr.open_dataset('data/ERA5_land-sea-mask_180x360.nc').isel(time=0)
      # landfrac_ds = landfrac_ds.stack(ncol=('latitude', 'longitude'))
      # landfrac_ds['ncol'] = np.arange(len(landfrac_ds['ncol']))
      # landfrac = landfrac_ds['lsm']

      # scripfile_tmp = xr.open_dataset('scrip_files/cmip6_180x360_scrip.20210722.nc')

      # # no_land_neighbors = xr.DataArray( np.ones(len(cnt_ds['ncol']),dtype=bool), coords=[cnt_ds['ncol']] )
      # # for i in range(len(cnt_ds['ncol'])):
      # #    if landfrac[i]<0.1:
      # #       if any( landfrac[cnt_ds['neighbors'][i,:]]>0.1 ) : no_land_neighbors[i] = 0
      
      # # mask = xr.DataArray( np.ones(len(cnt_ds['ncol']),dtype=bool), coords=[cnt_ds['ncol']] )
      # # mask = mask & (landfrac.values<0.1) & no_land_neighbors
      # # tmp_data_masked = tmp_data.where(mask,drop=True)#.to_dataset()
      # # mask = mask.rename({'ncol'  :'grid_size'})
      # # scripfile_masked = scripfile_tmp.where(mask,drop=True)

      # fig_file = 'F04_tmp_map'
      # wkres = ngl.Resources()
      # npix = 2048; wkres.wkWidth,wkres.wkHeight=npix,npix
      # wks_tmp = ngl.open_wks('png',fig_file,wkres)
      # plot_tmp = []; #[None]*2
      # res_tmp = hs.res_contour_fill_map()
      # res_tmp.cnFillMode    = "CellFill"
      # res_tmp.sfXArray      = scripfile_tmp['grid_center_lon'].values
      # res_tmp.sfYArray      = scripfile_tmp['grid_center_lat'].values
      # res_tmp.sfXCellBounds = scripfile_tmp['grid_corner_lon'].values
      # res_tmp.sfYCellBounds = scripfile_tmp['grid_corner_lat'].values

      # ### raw or normalized counts data
      # plot_tmp.append( ngl.contour_map(wks_tmp,tmp_data.values,res_tmp) )

      # ### data with land mask applied
      # # res_tmp_masked = copy.deepcopy(res_tmp)
      # # res_tmp_masked.sfXArray      = scripfile_masked['grid_center_lon'].values
      # # res_tmp_masked.sfYArray      = scripfile_masked['grid_center_lat'].values
      # # res_tmp_masked.sfXCellBounds = scripfile_masked['grid_corner_lon'].values
      # # res_tmp_masked.sfYCellBounds = scripfile_masked['grid_corner_lat'].values
      # # plot_tmp[1] = ngl.contour_map(wks_tmp,tmp_data_masked.values,res_tmp_masked) 

      # ngl.panel(wks_tmp,plot_tmp,[len(plot_tmp),1],hs.setres_panel())
      # ngl.end()
      # hc.trim_png(fig_file)
      # exit()
      #-------------------------------------------------------------------------
      #-------------------------------------------------------------------------

      ### sum across all columns
      cnt_ds['cnt'] = cnt_ds['cnt'].sum(dim='ncol')

      ### convert to frequency
      if convert_to_freq:
         # cnt_ds['cnt'] = cnt_ds['cnt'] / ( len(cnt_ds['num_valid']) * len(cnt_ds['ncol']) )
         num_valid_sum  = np.sum(np.ma.masked_invalid(cnt_ds['num_valid'].values))
         cnt_ds['cnt'] = cnt_ds['cnt'] / num_valid_sum

      if print_stats: hc.print_stat(cnt_ds['cnt'],name='final count dataset',stat='naxs',indent='    ')

      cnt_ds_list.append(cnt_ds)

#---------------------------------------------------------------------------------------------------
# combine sets
#---------------------------------------------------------------------------------------------------

### combine sets that contain partial checkerboard
if combine_sets and combine_mode==1:
   cnt_ds_list_tmp = []
   for cnt_ds in cnt_ds_list:
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
      # cnt_ds = xr.concat( [ ds_sum_nox, ds_sum_chk ], 'set' )
      cnt_ds = xr.concat( [ ds_sum_chk ], 'set' ) ### only show partial chx

      cnt_ds_list_tmp.append(cnt_ds)

   cnt_ds_list = cnt_ds_list_tmp
   # set_labels = ['no checkerboard','partial checkerboard']
   set_labels = ['partial checkerboard'] ### only show partial chx
   num_set = len(set_labels)
   sort_idx = [c for c in range(num_set)]

   
   

### combine sets based on a measure of discontinuities (# of local min/max)
if combine_sets and combine_mode==2:
   chk_idx,nox_idx = [],[]
   

   ### method for counting "transitions"
   # continuity_idx = [[] for _ in range(5)]
   # for s in range(num_set):
   #    tset = sets[s,:].values
   #    cont = np.sum(np.abs(tset - np.roll(tset,1)))
   #    continuity_idx[int(cont/2)].append(s)
   #    # print(f'set                                    : {tset}')
   #    # print(f'np.roll(tset,1)                        : {np.roll(tset,1)}' )
   #    # print(f'np.abs(tset - np.roll(tset,1))         : {np.abs(tset - np.roll(tset,1))}' )
   #    # print(f'np.sum(np.abs(tset - np.roll(tset,1))) : {np.sum(np.abs(tset - np.roll(tset,1)))}' )
   #    # print()

   # cnt_ds_list_tmp = []
   # for cnt_ds in cnt_ds_list:
   #    ds_tmp_0 = cnt_ds.isel(set=continuity_idx[0]).sum(dim='set')
   #    ds_tmp_1 = cnt_ds.isel(set=continuity_idx[1]).sum(dim='set')
   #    ds_tmp_2 = cnt_ds.isel(set=continuity_idx[2]).sum(dim='set')
   #    ds_tmp_3 = cnt_ds.isel(set=continuity_idx[3]).sum(dim='set')
   #    ds_tmp_4 = cnt_ds.isel(set=continuity_idx[4]).sum(dim='set')
   #    cnt_ds = xr.concat( [ ds_tmp_0, ds_tmp_1, ds_tmp_2, ds_tmp_3, ds_tmp_4 ], 'set' )
   #    cnt_ds_list_tmp.append(cnt_ds)


   ### method for counting number of local min/max
   continuity_idx = [[] for _ in range(9)]
   for s in range(num_set):
      tset = sets[s,:].values
      x = len(tset)
      lmax,lmin = [0]*x,[0]*x
      # lmax,lmin = np.zeros(x),np.zeros(x)
      for i in range(x):
         il,ir = i-1,i+1
         if i==0  : il = x-1
         if i==x-1: ir = 0
         if all(tset[i]>[tset[il],tset[ir]]): lmax[i] = 1
         if all(tset[i]<[tset[il],tset[ir]]): lmin[i] = 1
      cnt = np.sum(lmax) + np.sum(lmin)
      continuity_idx[cnt].append(s)
      # print(f'set  : {tset}')
      # print(f'lmax : {lmax}'); print(f'lmin : {lmin}')
      # print(f'cnt  : {cnt}'); print()

   cnt_ds_list_tmp = []
   for cnt_ds in cnt_ds_list:
      ds_tmp_0 = cnt_ds.isel(set=continuity_idx[0]).sum(dim='set')
      ds_tmp_1 = cnt_ds.isel(set=continuity_idx[1]).sum(dim='set')
      ds_tmp_2 = cnt_ds.isel(set=continuity_idx[2]).sum(dim='set')
      ds_tmp_3 = cnt_ds.isel(set=continuity_idx[3]).sum(dim='set')
      ds_tmp_4 = cnt_ds.isel(set=continuity_idx[4]).sum(dim='set')
      ds_tmp_5 = cnt_ds.isel(set=continuity_idx[5]).sum(dim='set')
      # ds_tmp_6 = cnt_ds.isel(set=continuity_idx[6]).sum(dim='set')
      # ds_tmp_7 = cnt_ds.isel(set=continuity_idx[7]).sum(dim='set')
      ds_tmp_8 = cnt_ds.isel(set=continuity_idx[8]).sum(dim='set')
      # cnt_ds = xr.concat( [ ds_tmp_0, ds_tmp_1, ds_tmp_2, ds_tmp_3, ds_tmp_4, ds_tmp_5, ds_tmp_6, ds_tmp_7, ds_tmp_8 ], 'set' )
      cnt_ds = xr.concat( [ ds_tmp_0, ds_tmp_1, ds_tmp_2, ds_tmp_3, ds_tmp_4, ds_tmp_5, ds_tmp_8 ], 'set' )
      cnt_ds_list_tmp.append(cnt_ds)

   cnt_ds_list = cnt_ds_list_tmp
   
   # sets = xr.DataArray(np.array([[0,1,0,1,0,1,0,1], [1,0,1,0,1,0,1,0]]),dims=['set','neighbors'] )
   
   # set_labels = [f'{c*2}' for c in range(5)]
   # sort_idx = [c for c in range(5)]

   # set_labels = [f'{c}' for c in range(9)]
   # sort_idx = [c for c in range(9)]
   sort_idx = [0,1,2,3,4,5,6]
   set_labels = [f'{c}' for c in [0,1,2,3,4,5,8]]

   num_set = len(set_labels)

#---------------------------------------------------------------------------------------------------
# Difference from first case
#---------------------------------------------------------------------------------------------------


#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------

num_plot = 1
if add_diff_plot: 
   num_plot += 1
   diff_plot_idx = 1

plot = [None]*(num_plot)

show_as_diff_save = show_as_diff
cnt_ds_list_save  = copy.deepcopy(cnt_ds_list)

for ip in range(num_plot):

   cnt_ds_list = copy.deepcopy(cnt_ds_list_save)

   if (add_diff_plot and ip==diff_plot_idx):
      show_as_diff = True
   else:
      show_as_diff = show_as_diff_save

   if show_as_diff or show_as_ratio or add_diff_plot:
      for v in range(num_var):
         c = 0
         baseline = cnt_ds_list[c*num_var+v]['cnt'].values
         for c in range(0,num_case):
            if show_as_diff:  cnt_ds_list[c*num_var+v]['cnt'] = cnt_ds_list[c*num_var+v]['cnt'] - baseline
            if show_as_ratio: cnt_ds_list[c*num_var+v]['cnt'] = cnt_ds_list[c*num_var+v]['cnt'] / baseline

   data_min,data_max = 1e10,0
   for cnt_ds in cnt_ds_list:
      data_min = np.min([ data_min, np.min(cnt_ds['cnt'].values) ])
      data_max = np.max([ data_max, np.max(cnt_ds['cnt'].values) ])
   if not show_as_diff and not show_as_ratio: data_min = 0

   x_values = np.arange(0,num_set,1,dtype=float)

   x_min = np.min(x_values)-1.5
   x_max = np.max(x_values)+1.5

   res.trYMinF,res.trYMaxF = data_min, data_max
   res.trXMinF,res.trXMaxF = x_min, x_max

   res.tmXBLabelAngleF  = -90.
   res.tmXBMode         = 'Explicit'
   res.tmXBValues       = x_values[:pchk_start]
   res.tmXBLabels       = set_labels[:pchk_start]
   res.xyMarkLineMode   = 'Lines'
   res.tmXBLabelFontHeightF = 0.0005
   res.tmXBLabelFontColor = 'black'

   if combine_sets and combine_mode==1:
      res.tmXBLabelAngleF  = -45.
      res.tmXBLabelFontHeightF = 0.0020

   if combine_sets and combine_mode==2:
      res.tmXBLabelAngleF  = 0.
      res.tmXBLabelFontHeightF = 0.0020

   # resources for multi-colored axis lables
   # tres = copy.deepcopy(res)
   if alt_chk_label:
      res.tmXBLabelFontColor = 'gray'
      tres = hs.res_default()
      tres.vpHeightF = 0.4
      tres.tmYLOn = False
      tres.trYMinF,tres.trYMaxF = data_min, data_max
      tres.trXMinF,tres.trXMaxF = x_min, x_max
      tres.tmXBLabelAngleF       = -90.
      tres.tmXBLabelFontHeightF = 0.0005
      tres.tmXBMode             = 'Explicit'
      tres.tmXBLabelFontColor   = 'black'
      tres.tmXBValues           = x_values[pchk_start:]
      tres.tmXBLabels           = set_labels[pchk_start:]

      # print(f'  data_min: {data_min}  data_max: {data_max}')
      # ### special option to only show partial checkerboard and reset vertical extent
      # data_max = 0
      # for c in range(num_case):
      #    print()
      #    for s in range(num_set):
      #       ss = sort_idx[s]
      #       cnt_ds = cnt_ds_list[c]
      #       if ss>=pchk_start:
      #          data_max = np.max([ data_max, cnt_ds['cnt'][ss].values ])
      #          val = cnt_ds['cnt'][ss].values
      #          print(f'  val: {val:8.4f}  data_max: {data_max:8.4f}')
      # print(f'  data_min: {data_min}  data_max: {data_max}')
      # data_max = 0.03
      # res.trYMinF,res.trYMaxF = data_min, data_max
      # res.trXMinF,res.trXMaxF = x_values[pchk_start]-0.5, x_max
      # tres.trYMinF,tres.trYMaxF = res.trYMinF,res.trYMaxF
      # tres.trXMinF,tres.trXMaxF = res.trXMinF,res.trXMaxF
    
   dx = 0.3
   ymin = 1.0 if show_as_ratio else 0.0
   bar_width_perc=0.6
   dxp = (dx * bar_width_perc)/2.

   for v in range(num_var):
      # c = 0
      # yy_baseline = cnt_ds_list[c*num_var+v]['cnt'].values
      for c in range(num_case):
         yy_all_sets = cnt_ds_list[c*num_var+v]['cnt'].values
         for s in range(num_set):

            ss = sort_idx[s]

            xx = x_values[s]
            yy = yy_all_sets[ss]

            xbar = np.array([ xx-dxp, xx+dxp, xx+dxp, xx-dxp, xx-dxp])
            ybar = np.array([ ymin, ymin,  yy,  yy, ymin])

            ### Shift to accomadate multiple cases
            for b in range(len(xbar)) : xbar[b] = xbar[b] - dxp*(num_case-1) + dxp*2*c

            ### plot polygon outline
            tplot = ngl.xy(wks,xbar,ybar,res)
            if c==0 and s==0:
               plot[ip] = tplot
            else:
               ngl.overlay( plot[ip], tplot )

            ### Add filled polygon
            pgres.gsFillColor = clr[c]
            pgres.gsFillIndex = 0
            polygon_dummy = ngl.add_polygon(wks,plot[ip],xbar,ybar,pgres)

            

            ### overlay pattern
            if pat[c]>0:
               pgres.gsFillColor = 'black'
               pgres.gsFillIndex = pat[c]
               polygon_dummy = ngl.add_polygon(wks,plot[ip],xbar,ybar,pgres)

         # overlay for colored axis lables       
         if alt_chk_label: 
            ngl.overlay( plot[ip], ngl.blank_plot(wks, tres) )
      
      subtitle_font_height = 0.01
      var_name = var[v]
      if var_name=='TGCLDLWP': var_name = 'Liq Water Path'
      # lstr,cstr,rstr = '','',var[v]
      lstr,cstr,rstr = var_name,'',''
      if show_as_diff: lstr,cstr,rstr = var_name,'',f'diff from {name[0]}'
      hs.set_subtitles(wks, plot[ip], lstr, cstr, rstr, font_height=subtitle_font_height)

#---------------------------------------------------------------------------------------------------
# Finalize plot
#---------------------------------------------------------------------------------------------------

#----------------------------------------------------------------------------
# Ad legend
#----------------------------------------------------------------------------
lgres = ngl.Resources()
lgres.vpWidthF, lgres.vpHeightF  = 0.09, 0.09
lgres.lgLabelFontHeightF = 0.01
lgres.lgLineThicknessF   = 20
lgres.lgMonoLineColor    = False
lgres.lgMonoDashIndex    = True
lgres.lgLineColors       = clr
lgres.lgDashIndexes      = 0#dsh
lgres.lgLabelJust    = 'CenterLeft'
lname = [f' {n}' for n in name]
pid = ngl.legend_ndc(wks, len(name), lname, 0.3, 0.6, lgres)  # 1x2


hs.set_plot_labels(wks, plot, font_height=0.01, justify='left')


# layout = [1,num_var]
# layout = [int(np.ceil(num_var/3.)),3]
layout = [int(np.ceil(num_plot/2.)),2]

pnl_res = hs.setres_panel()
pnl_res.nglPanelYWhiteSpacePercent = 5

# pnl_res.nglPanelFigureStrings            = list(string.ascii_lowercase)
# pnl_res.nglPanelFigureStringsJust        = "TopLeft"
# pnl_res.nglPanelFigureStringsFontHeightF = 0.01

ngl.panel(wks,plot,layout,pnl_res)
ngl.end()

hc.trim_png(fig_file)
#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------
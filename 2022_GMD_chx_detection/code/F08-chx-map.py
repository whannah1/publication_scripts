# plot the fractional occurence of all possible patterns
import os, ngl, sys, numba, copy, xarray as xr, numpy as np, re, string
import hapy_common as hc, hapy_E3SM as he, hapy_setres as hs
import pg_checkerboard_utilities as pg

case,name,clr,dsh,pat = [],[],[],[],[]
def add_case(case_in,n='',c='black',d=0,p=0):
   case.append(case_in); name.append(n); clr.append(c); dsh.append(d); pat.append(p)

var,var_str,lev_list = [],[],[]
def add_var(var_name,n=None,lev=-1): 
   if n is None: n = var_name
   var.append(var_name); var_str.append(n); lev_list.append(lev)
#-------------------------------------------------------------------------------
add_case('MAC-PG',  n='MAC')
# add_case('OBS-PG',  n='OBS')
add_case('E3SM.CHX.ne30pg2_ne30pg2.F-MMF1.00',            n='E3SM-MMF'       )
add_case('E3SM.CHX.ne30pg2_ne30pg2.FC5AV1C-L.no_dcape.00',n='E3SM (no DCAPE)')
add_case('E3SM.CHX.ne30pg2_ne30pg2.FC5AV1C-L.00',         n='E3SM'           )
#-------------------------------------------------------------------------------

# add_var('PRECT',n='Precipitation')
add_var('TGCLDLWP',n='Liq Water Path')


# use_regional_subset = False
lat1,lat2 = -60,60
# lat1,lat2,lon1,lon2 = -20,40,80,180+30


fig_type = 'png'
fig_file = 'figs/F08-chx-map'


use_daily         = True
convert_to_freq   = True   # convert count to fractional occurrence
common_colorbar   = True   # use common set of color levels for all plots

ocean_only = False

subset_min_length = 4
sort_sets    = True
combine_sets = True
combine_mode = 1     # 1 = chk vs no chk / 2 = continuity metric

var_x_case           = True   # controls plot panel arrangement
num_plot_col         = 2      # controls plot arrangement for single variable or single case
print_stats          = False   #
verbose              = True    #

#---------------------------------------------------------------------------------------------------
# generate sets of possible neighbor states
#---------------------------------------------------------------------------------------------------

### only pure checkboard
# tmp_file='data/occurrence-chx-only'; rotate_sets=False; sets=pg.chx_only_sets

### full set of all possible patterns
tmp_file='data/chx-occurrence'; rotate_sets=True; sets=pg.all_possible_sets
# tmp_file='data/occurrence-all-sets'; rotate_sets=True; sets=pg.all_possible_sets

### full set of patterns
rotate_sets = True; sets = pg.all_possible_sets

(num_set,set_len) = sets.shape
set_coord,nn_coord = np.arange(num_set),np.arange(set_len)
sets.assign_coords(coords={'set':set_coord,'neighbors':nn_coord})

set_labels = pg.get_set_labels(sets)

#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------
def get_comp(case):
   comp = 'eam'
   if case=='MAC': comp = None
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
   # for s in range(len(nox_sets)):
   #    print(f'  {nox_sets[s]}')
   # exit()
   sets_sorted = xr.DataArray(np.array( nox_sets + chk_sets), dims=['set','neighbors'] )
   for s in range(num_set):
      for ss in range(num_set):
         if np.all(sets[s,:].values==sets_sorted[ss,:].values): 
            sort_idx[ss] = s
   ### useful for checking sort_idx
   # for s in range(num_set): print(f'  {set_labels[s]}    {set_labels[sort_idx[s]]}')
   # exit()
   set_labels_sorted = []
   for s in range(num_set):
      set_labels_sorted.append( set_labels[sort_idx[s]] )
   set_labels = set_labels_sorted
else:
   for s in range(num_set): sort_idx[s] = s


#---------------------------------------------------------------------------------------------------
# Set up plot resources
#---------------------------------------------------------------------------------------------------
num_case,num_var,num_set = len(case),len(var),len(sets)

wkres = ngl.Resources()
# npix=2**13; wkres.wkWidth,wkres.wkHeight=npix,npix # use this for plotting all patterns w/ rotation
wks = ngl.open_wks(fig_type,fig_file,wkres)
plot = [None]*(num_var)

res = hs.res_contour_fill_map()
res.tmXBOn = False
res.tmYLOn = False
# res.lbTitlePosition = 'Bottom'
# res.lbTitleFontHeightF = 0.01
# res.lbTitleString = 'Count'
if common_colorbar:
   res.lbLabelBarOn = False
else:
   res.lbLabelBarOn = True

if 'lat1' in vars() : res.mpMinLatF = lat1
if 'lat2' in vars() : res.mpMaxLatF = lat2
if 'lon1' in vars() : res.mpMinLonF = lon1
if 'lon2' in vars() : res.mpMaxLonF = lon2

if 'lev' not in vars() : lev = np.array([-1])

#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------
cnt_ds_list = []
for c in range(num_case):
   print('\n  case: '+hc.tcolor.CYAN+case[c]+hc.tcolor.ENDC)

   comp = 'eam'
   if case[c]=='ERA5': comp = None
   if case[c]=='E3SM.RGMA.ne30pg2_r05_oECv3.FC5AV1C-L.00': comp = 'cam'

   use_remap = False
   remap_str=f'remap_ne30pg2'
   if case[c]=='MAC-FV' : use_remap = False
   if case[c]=='MAC-PG' : use_remap = True

   #----------------------------------------------------------------------------
   #----------------------------------------------------------------------------
   for v in range(num_var) :
      tname = case[c]
      if 'MAC' in case[c]: tname = 'MAC'
      if 'OBS' in case[c] and var[v]=='TGCLDLWP': tname = 'MAC'
      if 'OBS' in case[c] and var[v]=='PRECT'   : tname = 'GPM'
      if 'OBS' in case: alt_name = alt_name.replace('OBS',tname)
      case_obj = he.Case( name=tname, time_freq='daily' )

      alt_case = case[c]
      if 'OBS' in case[c]: alt_case = alt_case.replace('OBS',tname)

      case_tmp_file = f'{tmp_file}.daily.{alt_case}.{var[v]}.nc'
      
      print('    var: '+hc.tcolor.GREEN+f'{var[v]:10}'+hc.tcolor.ENDC+'    '+case_tmp_file)
      
      cnt_ds = xr.open_dataset( case_tmp_file )

      #-------------------------------------------------------------------------
      # print the time length
      #-------------------------------------------------------------------------
      msg = '    num_time: '
      msg += str(cnt_ds['num_time'].values)
      msg += ' days'
      print(msg)

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
      #-------------------------------------------------------------------------

      ### sum across all columsn
      # cnt_ds['cnt'] = cnt_ds['cnt'].sum(dim='ncol')

      ### convert to frequency
      # cnt_ds['cnt'] = cnt_ds['cnt'] / np.ma.masked_invalid(cnt_ds['num_valid'].values)
      
      # print()
      # print(cnt_ds['cnt'])
      # print()
      # print(cnt_ds['num_valid'])
      # print()
      
      cnt_ds['cnt'] = cnt_ds['cnt'] / cnt_ds['num_valid'] 
      # cnt_ds['cnt'] = cnt_ds['cnt'] / cnt_ds['num_time'] #/ len(cnt_ds['ncol']) #* 100.

      if print_stats: hc.print_stat(cnt_ds['cnt'],name='final count dataset',stat='naxs',indent='    ')

      cnt_ds_list.append(cnt_ds)


#---------------------------------------------------------------------------------------------------
# combine sets
#---------------------------------------------------------------------------------------------------

### combine sets that contain partial checkerboard
if combine_sets and combine_mode==1:
   only_show_partial = True
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
      if only_show_partial:
         cnt_ds = xr.concat( [ ds_sum_chk ], 'set' ) ### only show partial chx
      else:
         cnt_ds = xr.concat( [ ds_sum_nox, ds_sum_chk ], 'set' )

      cnt_ds_list_tmp.append(cnt_ds)

   cnt_ds_list = cnt_ds_list_tmp
   if only_show_partial:
      set_labels = ['partial checkerboard'] ### only show partial chx
   else:
      set_labels = ['no checkerboard','partial checkerboard']
   num_set = len(set_labels)
   sort_idx = [c for c in range(num_set)]

   
   

### combine sets based on a measure of discontinuities (# of local min/max)
if combine_sets and combine_mode==2:
   chk_idx,nox_idx = [],[]

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
   # exit()

   cnt_ds_list_tmp = []
   for cnt_ds in cnt_ds_list:
      ds_tmp_0 = cnt_ds.isel(set=continuity_idx[0]).sum(dim='set')
      ds_tmp_1 = cnt_ds.isel(set=continuity_idx[1]).sum(dim='set')
      ds_tmp_2 = cnt_ds.isel(set=continuity_idx[2]).sum(dim='set')
      ds_tmp_3 = cnt_ds.isel(set=continuity_idx[3]).sum(dim='set')
      ds_tmp_4 = cnt_ds.isel(set=continuity_idx[4]).sum(dim='set')
      ds_tmp_5 = cnt_ds.isel(set=continuity_idx[5]).sum(dim='set')
      ds_tmp_6 = cnt_ds.isel(set=continuity_idx[6]).sum(dim='set')
      ds_tmp_7 = cnt_ds.isel(set=continuity_idx[7]).sum(dim='set')
      ds_tmp_8 = cnt_ds.isel(set=continuity_idx[8]).sum(dim='set')
      cnt_ds = xr.concat( [ ds_tmp_0, ds_tmp_1, ds_tmp_2, ds_tmp_3, ds_tmp_4, ds_tmp_5, ds_tmp_6, ds_tmp_7, ds_tmp_8 ], 'set' )
      cnt_ds_list_tmp.append(cnt_ds)

   cnt_ds_list = cnt_ds_list_tmp
   
   # sets = xr.DataArray(np.array([[0,1,0,1,0,1,0,1], [1,0,1,0,1,0,1,0]]),dims=['set','neighbors'] )
   set_labels = [f'{c*2}' for c in range(5)]
   sort_idx = [c for c in range(5)]
   num_set = len(set_labels)

#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------

# if not make_plot: exit('\nmake_plot=False, so aborting the script now.\n')

res.cnFillPalette = "MPL_viridis"

if common_colorbar:
   data_min = np.zeros(num_set)
   data_max = np.zeros(num_set)
   for cnt_ds in cnt_ds_list:
      for s in range(num_set):
         data_min[s] = np.min([ data_min[s], np.min(cnt_ds['cnt'].isel(set=s).values) ])
         data_max[s] = np.max([ data_max[s], np.max(cnt_ds['cnt'].isel(set=s).values) ])
   
   # if combine_sets and num_set==2:
   #    data_min[0],data_min[0] = 0.,0.
   #    data_max[1],data_max[1] = 1.,0.5

   # if combine_sets and num_set==1:
   #    data_min[0] = 0.
   #    data_max[0] = 0.1

   # if convert_to_freq:
   #    data_min[:] = 0.
   #    data_max[:] = 1.
   #    # nlev = 20

   # print()
   # print(f'  data_min: {data_min}')
   # print(f'  data_max: {data_max}')
   # print()

#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------
list_cnt = 0
plot = [None]*(num_var*num_case*num_set)

for c in range(num_case):
   
   tname = case[c]
   if 'MAC' in case[c]: tname = 'MAC'
   if 'OBS' in case[c] and var[v]=='TGCLDLWP': tname = 'MAC'
   if 'OBS' in case[c] and var[v]=='PRECT'   : tname = 'GPM'
   if 'OBS' in case: alt_name = alt_name.replace('OBS',tname)
   case_obj = he.Case( name=tname, time_freq='daily' )

   scrip_file_path = 'scrip_files/ne30pg2_scrip.nc'
   if case[c]=='MAC-FV': scrip_file_path = 'scrip_files/cmip6_180x360_scrip.20210722.nc'

   scripfile = xr.open_dataset(scrip_file_path)

   if ocean_only: 
      if '-FV' in case[c] : 
         print('ocean_only not implemented for FV cases!')
      else:
         landfrac_ds = xr.open_dataset('data/E3SM_landfrac_ne30pg2.nc').isel(time=0)
         mask = xr.DataArray( np.ones(len(cnt_ds['ncol']),dtype=bool), coords=[cnt_ds['ncol']] )
         mask = mask & (landfrac_ds['LANDFRAC'].values<0.5)
         scripfile = scripfile.where(mask.rename({'ncol':'grid_size'}),drop=True)

   res.cnFillMode    = "CellFill"
   res.sfXArray      = scripfile['grid_center_lon'].values
   res.sfYArray      = scripfile['grid_center_lat'].values
   res.sfXCellBounds = scripfile['grid_corner_lon'].values
   res.sfYCellBounds = scripfile['grid_corner_lat'].values

   for v in range(num_var) :

      tres = copy.deepcopy(res)
   
      cnt_ds = cnt_ds_list[list_cnt]
      list_cnt += 1

      for s in range(num_set):

         if common_colorbar:
            # if 'nlev' not in locals(): nlev = 11
            # aboutZero = False
            # (cmin,cmax,cint) = ngl.nice_cntr_levels(data_min[s], data_max[s], 
            #                                        cint=None, max_steps=nlev, 
            #                                        returnLevels=False, aboutZero=aboutZero )
            # tres.cnLevels = np.linspace(cmin,cmax,num=nlev)
            tres.cnLevelSelectionMode = 'ExplicitLevels'
            tres.cnLevels = np.arange(5,60+5,5)/1e2
            

         # if num_var>1:
         #    ip = s*num_var+v if var_x_case else v*num_set+s
         # else:
         #    ip = s*num_case+c if var_x_case else c*num_set+s
         ip = c*num_var+v if var_x_case else v*num_case+c

         tmp_data = cnt_ds['cnt'].isel(set=s)

         if ocean_only: tmp_data = tmp_data.where(mask,drop=True)

         plot[ip] = ngl.contour_map(wks, tmp_data.values, tres) 

         set_str = set_labels[s]
         subtitle_font_height = 0.015
         # if num_set>3: subtitle_font_height = 0.002
         # if (num_var*num_case*num_set)>4: subtitle_font_height = 0.008
         # if num_set>12: subtitle_font_height = 0.002
         
         tname = name[c]
         alt_name = name[c]
         if 'OBS' in case[c] and var=='TGCLDLWP': tname = 'MAC'
         if 'OBS' in case[c] and var=='PRECT'   : tname = 'GPM'
         if 'OBS' in case[c]: alt_name = alt_name.replace('OBS',tname)

         hs.set_subtitles(wks, plot[ip], alt_name, '', var_str[v], font_height=subtitle_font_height)

#---------------------------------------------------------------------------------------------------
# Finalize plot
#---------------------------------------------------------------------------------------------------
# hs.set_plot_labels(wks, plot, font_height=0.01, justify='left')

# if num_var>1:
#    layout = [num_set,num_var] if var_x_case else [num_var,num_set]
# else:
#    layout = [num_set,num_case] if var_x_case else [num_case,num_set]

layout = [num_case,num_var] if var_x_case else [num_var,num_case]

# if num_case==1 and num_set>3: 
#    num_plot_col = 6
#    layout = [int(np.ceil(len(plot)/float(num_plot_col))),num_plot_col]

if num_case==1 or num_var==1:
   layout = [int(np.ceil(len(plot)/float(num_plot_col))),num_plot_col]

pnl_res = hs.setres_panel()
pnl_res.nglPanelYWhiteSpacePercent = 5
pnl_res.nglPanelFigureStrings            = list(string.ascii_lowercase)
pnl_res.nglPanelFigureStringsJust        = "TopLeft"
pnl_res.nglPanelFigureStringsFontHeightF = 0.01

if common_colorbar: 
   pnl_res.nglPanelLabelBar = True
   pnl_res.lbTitleFontHeightF = 0.015
   # pnl_res.lbLabelFontHeightF = 0.01
   pnl_res.nglPanelLabelBarLabelFontHeightF = 0.01
   pnl_res.lbTitlePosition = 'Bottom'
   pnl_res.lbTitleString = 'Fractional Occurrence of Partial Checkerboard'

ngl.panel(wks,plot,layout,pnl_res)
ngl.end()

hc.trim_png(fig_file)
#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------
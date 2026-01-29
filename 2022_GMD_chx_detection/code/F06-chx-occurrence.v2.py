# plot the fractional occurence of all possible patterns
# v2 - simplified from v1, use combine_mode=2 for neighborhood continuity index
import os, ngl, sys, numba, copy, xarray as xr, numpy as np, string
import hapy_common as hc, hapy_E3SM as he, hapy_setres as hs
import pg_checkerboard_utilities as pg

case,name,clr,dsh,pat = [],[],[],[],[]
def add_case(case_in,n='',c='black',d=0,p=0):
   global name,case
   case.append(case_in); name.append(n); clr.append(c); dsh.append(d); pat.append(p)
var,var_str,lev_list = [],[],[]
def add_var(var_name,n=None,lev=-1): 
   if n is None: n = var_name
   var.append(var_name); var_str.append(n); lev_list.append(lev)
#-------------------------------------------------------------------------------
# add_case('MAC-FV',    n='MAC 1x1 deg'                                        ,c='black',p=0)
# add_case('MAC-PG',    n='MAC ne30pg2'                                        ,c='gray', p=0)
# add_case('GPM-FV',    n='GPM 1x1 deg'                                        ,c='black',p=0)
# add_case('GPM-PG',    n='GPM ne30pg2'                                        ,c='gray', p=0)
add_case('OBS-PG',    n='OBS ne30pg2'                                        ,c='black',p=0)
add_case('OBS-FV',    n='OBS 1x1 deg'                                        ,c='gray', p=0)
add_case('E3SM.CHX.ne30pg2_ne30pg2.FC5AV1C-L.no_dcape.00',n='E3SM (no DCAPE)',c='green',p=0)
add_case('E3SM.CHX.ne30pg2_ne30pg2.FC5AV1C-L.00',         n='E3SM'           ,c='red'  ,p=0)
add_case('E3SM.CHX.ne30pg2_ne30pg2.F-MMF1.00',            n='E3SM-MMF'       ,c='blue' ,p=0)
add_case('E3SM.CHX.RAND',                                 n='RANDOM'         ,c='magenta',p=0)

#-------------------------------------------------------------------------------
var = []
var.append('TGCLDLWP')
# var.append('PRECT')

fig_type = 'png'
fig_file = 'figs/F06-chx-occurrence'
tmp_file = 'data/chx-occurrence'

ocean_only        = True   # ocean only to compare with obs

use_regional_subset = True
lat1,lat2 = 0,30
lon1,lon2 = 140,180+40

add_diff_plot     = True
print_stats       = False   # print stats for debugging/sanity check

#---------------------------------------------------------------------------------------------------
# generate sets of possible neighbor states
#---------------------------------------------------------------------------------------------------
num_neighbors = 8

### full set of all possible patterns
rotate_sets = True; sets = pg.all_possible_sets

(num_set,set_len) = sets.shape
sets.assign_coords(coords={'set':np.arange(num_set),'neighbors':np.arange(set_len)})

set_labels = pg.get_set_labels(sets)

#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------
def get_comp(case):
   comp = 'eam'
   if 'OBS'  in case: comp = None
   if 'MAC'  in case: comp = None
   if 'GPM'  in case: comp = None
   if 'CESM' in case: comp = 'cam'
   return comp
#---------------------------------------------------------------------------------------------------
# sort sets according to partial checkerboard
#---------------------------------------------------------------------------------------------------

sort_idx = [-1]*num_set

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
res.tiXAxisFontHeightF     = 0.015
res.tiYAxisFontHeightF     = 0.015
res.tmYLLabelFontHeightF   = 0.025
res.tmXBLabelFontHeightF   = 0.025
res.tiYAxisString = 'Fractional Occurrence'
res.tiXAxisString = 'Number of Extrema in Local Neighborhood'

pgres = ngl.Resources()
# pgres.nglDraw,pgres.nglFrame = True,False
pgres.nglDraw,pgres.nglFrame = False,False

lres = hs.res_xy()
lres.xyLineThicknessF = 1
lres.xyDashPattern    = 0
lres.xyLineColor      = 'black'

#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------
cnt_ds_list = []
for c in range(num_case):
   print('\n  case: '+hc.tcolor.CYAN+case[c]+hc.tcolor.ENDC)
   for v in range(num_var) :

      tname = case[c]
      if 'MAC'  in case[c]: tname = 'MAC'
      if 'GPM'  in case[c]: tname = 'GPM'
      if 'OBS'  in case[c]: 
         if var[v]=='TGCLDLWP': tname = 'MAC'
         if var[v]=='PRECT'   : tname = 'GPM'
      
      case_obj = he.Case( name=tname, time_freq='daily' )
      if 'lev' not in vars() : lev = np.array([-1])

      comp = 'eam'

      use_remap = False
      remap_str=f'remap_ne30pg2'
      if case[c]=='MAC-FV' : use_remap = False
      if case[c]=='MAC-PG' : use_remap = True
      if case[c]=='GPM-FV' : use_remap = True
      if case[c]=='GPM-PG' : use_remap = True

      if case[c]=='OBS-FV' and tname=='MAC' : use_remap = False
      if case[c]=='OBS-PG' and tname=='MAC' : use_remap = True
      if case[c]=='OBS-FV' and tname=='GPM' : use_remap = True
      if case[c]=='OBS-PG' and tname=='GPM' : use_remap = True
      #-------------------------------------------------------------------------
      #-------------------------------------------------------------------------
      case_tmp = case[c]
      if 'OBS' in case[c]: case_tmp = case_tmp.replace('OBS',tname)

      case_tmp_file = f'{tmp_file}.daily.{case_tmp}.{var[v]}.nc'
      
      # print('    var: '+hc.tcolor.GREEN+f'{var[v]:10}'+hc.tcolor.ENDC+'    '+case_tmp_file)
      
      cnt_ds = xr.open_dataset( case_tmp_file )

      ### print the time length of data
      ntime = cnt_ds.num_time.values; nyear = (ntime/365)
      # print(f'    num_time: {ntime} days ({nyear} years)')

      vstr = hc.tcolor.GREEN+f'{var[v]:10}'+hc.tcolor.ENDC
      print(f'    var: {vstr}  {case_tmp_file:90}  ({ntime:8.1f} days / {nyear:6.1f} years)')

      #-------------------------------------------------------------------------
      # Ocean mask
      #-------------------------------------------------------------------------
      if ocean_only: 
         if 'FV' in case[c]: 
            landfrac_ds = xr.open_dataset('data/land-sea-mask_180x360.nc').isel(time=0)
            landfrac_ds = landfrac_ds.stack(ncol=('latitude', 'longitude'))
            landfrac_ds['ncol'] = np.arange(len(landfrac_ds['ncol']))
            landfrac = landfrac_ds['lsm']
         else:
            landfrac_ds = xr.open_dataset('data/E3SM_landfrac_ne30pg2.nc').isel(time=0)
            landfrac = landfrac_ds['LANDFRAC']
         ocn_mask = xr.DataArray( np.ones(len(cnt_ds['ncol']),dtype=bool), coords=[cnt_ds['ncol']] )
         ocn_mask = ocn_mask & (landfrac.values<0.5)
         num_time = cnt_ds['num_time']
         # num_valid = cnt_ds['num_valid']
         # cnt_ds = cnt_ds['cnt'].where(ocn_mask,drop=True).to_dataset()
         cnt_ds['cnt'      ] = cnt_ds['cnt'      ].where(ocn_mask,drop=True)
         cnt_ds['num_valid'] = cnt_ds['num_valid'].where(ocn_mask,drop=True)
         cnt_ds['num_time'] = num_time
         # cnt_ds['num_valid'] = num_valid
      #-------------------------------------------------------------------------
      # apply regional subset
      #-------------------------------------------------------------------------
      if use_regional_subset: 
         ### load coordinate data from scrip file
         scripfile_path = 'scrip_files/ne30pg2_scrip.nc'
         if '-FV' in case[c]: scripfile_path = 'scrip_files/cmip6_180x360_scrip.20210722.nc'
         scrip_ds = xr.open_dataset(scripfile_path).rename({'grid_size':'ncol'})

         ### apply ocean mask
         if ocean_only: scrip_ds = scrip_ds.where(ocn_mask,drop=True)

         ### create new mask for regional subset
         mask = xr.DataArray( np.ones(len(cnt_ds['ncol']),dtype=bool), coords=[cnt_ds['ncol']] )
         if 'lat1' in locals(): mask = mask & (scrip_ds['grid_center_lat']>=lat1) & (scrip_ds['grid_center_lat']<=lat2)
         if 'lon1' in locals(): mask = mask & (scrip_ds['grid_center_lon']>=lon1) & (scrip_ds['grid_center_lon']<=lon2)
         
         num_time = cnt_ds['num_time']
         cnt_ds['cnt'      ] = cnt_ds['cnt'      ].where(mask,drop=True)
         cnt_ds['num_valid'] = cnt_ds['num_valid'].where(mask,drop=True)
         cnt_ds['num_time'] = num_time
      #-------------------------------------------------------------------------
      #-------------------------------------------------------------------------

      ### sum across all columns
      cnt_ds['cnt'] = cnt_ds['cnt'].sum(dim='ncol')

      ### convert to frequency
      num_valid_sum  = np.sum(np.ma.masked_invalid(cnt_ds['num_valid'].values))
      cnt_ds['cnt'] = cnt_ds['cnt'] / num_valid_sum

      if print_stats: hc.print_stat(cnt_ds['cnt'],name='final count dataset',stat='naxs',indent='    ')

      cnt_ds_list.append(cnt_ds)

#---------------------------------------------------------------------------------------------------
# combine sets based on a measure of discontinuities (# of local min/max)
#---------------------------------------------------------------------------------------------------
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

cnt_ds_list_tmp = []
for cnt_ds in cnt_ds_list:
   ds_tmp_0 = cnt_ds.isel(set=continuity_idx[0]).sum(dim='set')
   ds_tmp_1 = cnt_ds.isel(set=continuity_idx[1]).sum(dim='set')
   ds_tmp_2 = cnt_ds.isel(set=continuity_idx[2]).sum(dim='set')
   ds_tmp_3 = cnt_ds.isel(set=continuity_idx[3]).sum(dim='set')
   ds_tmp_4 = cnt_ds.isel(set=continuity_idx[4]).sum(dim='set')
   ds_tmp_5 = cnt_ds.isel(set=continuity_idx[5]).sum(dim='set')
   ds_tmp_8 = cnt_ds.isel(set=continuity_idx[8]).sum(dim='set')
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
#---------------------------------------------------------------------------------------------------

num_plot = 1
if add_diff_plot: num_plot = 2

plot = [None]*(num_plot*num_var)

cnt_ds_list_save  = copy.deepcopy(cnt_ds_list)

for n in range(num_plot):

   cnt_ds_list = copy.deepcopy(cnt_ds_list_save)

   show_as_diff = True if (add_diff_plot and n==1) else False

   if show_as_diff:
      for v in range(num_var):
         c = 0
         baseline = cnt_ds_list[c*num_var+v]['cnt'].values
         for c in range(0,num_case):
            cnt_ds_list[c*num_var+v]['cnt'] = cnt_ds_list[c*num_var+v]['cnt'] - baseline

   data_min,data_max = 1e10,0
   for cnt_ds in cnt_ds_list:
      data_min = np.min([ data_min, np.min(cnt_ds['cnt'].values) ])
      data_max = np.max([ data_max, np.max(cnt_ds['cnt'].values) ])
   if not show_as_diff: data_min = 0

   x_values = np.arange(0,num_set,1,dtype=float)

   x_min = np.min(x_values)-1.5
   x_max = np.max(x_values)+1.5

   pad = ( data_max - data_min ) * 0.02
   if show_as_diff:
      data_min, data_max = data_min-pad, data_max+pad
   else:
      data_min, data_max = data_min, data_max+pad
   res.trYMinF,res.trYMaxF = data_min, data_max
   res.trXMinF,res.trXMaxF = x_min, x_max
   res.tmXBLabelAngleF  = -90.
   res.tmXBMode         = 'Explicit'
   res.tmXBValues       = x_values[:pchk_start]
   res.tmXBLabels       = set_labels[:pchk_start]
   res.xyMarkLineMode   = 'Lines'
   res.tmXBLabelFontHeightF = 0.0005
   res.tmXBLabelFontColor = 'black'
   res.tmXBLabelAngleF,res.tmXBLabelFontHeightF  =   0.,0.0020
    
   dx = 0.3
   ymin = 0.0
   bar_width_perc=0.6
   dxp = (dx * bar_width_perc)/2.

   for v in range(num_var):
      ip = v*num_plot+n
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
            pgres.gsFillColor,pgres.gsFillIndex = clr[c],0
            polygon_dummy = ngl.add_polygon(wks,plot[ip],xbar,ybar,pgres)

            ### overlay pattern
            if pat[c]>0:
               pgres.gsFillColor = 'black'
               pgres.gsFillIndex = pat[c]
               polygon_dummy = ngl.add_polygon(wks,plot[ip],xbar,ybar,pgres)

      subtitle_font_height = 0.015
      var_name = var[v]
      if var_name=='TGCLDLWP': var_name = 'Liq Water Path'
      if var_name=='PRECT':    var_name = 'Precipitation'
      lstr,cstr,rstr = var_name,'',''
      tmp_base_name = name[0]
      if 'OBS' in tmp_base_name: 
         if var[v]=='TGCLDLWP': tmp_base_name = tmp_base_name.replace('OBS','MAC')
         if var[v]=='PRECT'   : tmp_base_name = tmp_base_name.replace('OBS','GPM')
      if show_as_diff: rstr = f'diff from {tmp_base_name}'
      hs.set_subtitles(wks, plot[ip], lstr, cstr, rstr, font_height=subtitle_font_height)


   ### add horizontal line at zero
   ngl.overlay( plot[ip], ngl.xy(wks,np.array([-1e3,1e3]),np.array([0,0]),lres) )

#---------------------------------------------------------------------------------------------------
# Add legend
#---------------------------------------------------------------------------------------------------
lgres = ngl.Resources()
lgres.lgLabelFontHeightF = 0.01
lgres.lgLineThicknessF   = 20
lgres.lgMonoLineColor,lgres.lgLineColors  = False, clr
lgres.lgMonoDashIndex,lgres.lgDashIndexes = True, 0
lgres.lgLabelJust    = 'CenterLeft'

if num_var==1:
   lgres.vpWidthF, lgres.vpHeightF  = 0.09, 0.13
   lname = [f' {n}' for n in name]
   xpos = 0.3 if add_diff_plot else 0.5
   ypos = 0.65 + (num_var-1)*0.1
   pid = ngl.legend_ndc(wks, len(name), lname, xpos, ypos, lgres)

if num_var==2:
   lgres.vpWidthF, lgres.vpHeightF  = 0.09, 0.11

   lname = [f' {n}' for n in name]
   for (n,l) in enumerate(lname):
      if 'OBS' in lname[n]: 
         lname[n] = lname[n].replace('OBS','MAC')
   xpos = 0.3 if add_diff_plot else 0.5
   ypos = 0.7 + (num_var-1)*0.1
   pid = ngl.legend_ndc(wks, len(name), lname, xpos, ypos, lgres)

   lname = [f' {n}' for n in name]
   for (n,l) in enumerate(lname):
      if 'OBS' in lname[n]: 
         lname[n] = lname[n].replace('OBS','IMERG')
   xpos = 0.3 if add_diff_plot else 0.5
   ypos = 0.35 + (num_var-1)*0.1
   pid = ngl.legend_ndc(wks, len(name), lname, xpos, ypos, lgres)

#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------
if use_regional_subset:

   ires = ngl.Resources()
   ires.nglDraw         = False
   ires.nglFrame        = False
   ires.nglMaximize     = False
   ires.tmXBOn          = False
   ires.tmYLOn          = False
   ires.vpHeightF       = 0.15
   ires.vpWidthF        = 0.15
   ires.mpGridAndLimbOn = False
   ires.mpCenterLonF    = 180.
   
   iplot = ngl.map(wks, ires) 

   ### draw box around map region used for inset
   bx = np.array([lon1,lon2,lon2,lon1,lon1])
   by = np.array([lat1,lat1,lat2,lat2,lat1])

   pgres = ngl.Resources()
   pgres.nglDraw,pgres.nglFrame = False,False
   pgres.gsLineColor = 'red'
   pgres.gsLineThicknessF = 6

   pdum = ngl.add_polyline(wks,iplot,bx,by,pgres)

   ### attach plot via annotation
   ares = ngl.Resources()
   ares.amZone = 1
   ares.amOrthogonalPosF = 0.35
   ares.amParallelPosF   = 0.75
   anno = ngl.add_annotation(plot[0], iplot ,ares)


#---------------------------------------------------------------------------------------------------
# Finalize plot
#---------------------------------------------------------------------------------------------------
if num_plot==1:
   num_plot_col = 2
   layout = [int(np.ceil(len(plot)/float(num_plot_col))),num_plot_col]
else:
   layout = [num_var,num_plot]

pnl_res = hs.setres_panel()
pnl_res.nglPanelYWhiteSpacePercent = 5
pnl_res.nglPanelFigureStrings            = list(string.ascii_lowercase)
pnl_res.nglPanelFigureStringsJust        = "TopLeft"
pnl_res.nglPanelFigureStringsFontHeightF = 0.01

ngl.panel(wks,plot,layout,pnl_res)
ngl.end()

hc.trim_png(fig_file)
#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------
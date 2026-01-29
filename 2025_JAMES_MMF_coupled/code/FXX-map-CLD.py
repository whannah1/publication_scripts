# v2 is a cleaned up version of v1 - functionality is mostly the same
import os, glob, ngl, subprocess as sp, numpy as np, xarray as xr, dask, copy, string, cmocean
import hapy_common as hc, hapy_E3SM   as he, hapy_setres as hs
host = hc.get_host()
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
name,case,case_dir,case_sub = [],[],[],[]
def add_case(case_in,n=None,p=None,s=None,g=None,c=None):
   global name,case,case_dir,case_sub
   tmp_name = case_in if n is None else n
   case.append(case_in); name.append(tmp_name); 
   case_dir.append(p); case_sub.append(s); 
#-------------------------------------------------------------------------------
var,lev_list,var_str = [],[],[]
def add_var(var_in,lev=-1,name=None): 
   var.append(var_in); lev_list.append(lev)
   var_str.append(var_in if name is None else name)
#-------------------------------------------------------------------------------
scrip_file_path  = os.getenv('HOME')+'/E3SM/data_grid/ne30pg2_scrip.nc'
tmp_path_hst_mmf = '/global/cfs/cdirs/m3312/whannah/2023-CPL'
tmp_path_hst_v2  = '/global/cfs/cdirs/m3312/whannah/e3smv2_historical'
tmp_path_co2_mmf = '/global/cfs/cdirs/m3312/whannah/2023-CPL'
tmp_sub = 'archive/atm/hist'

# add_case('CERES-EBAF',n='CERES-EBAF')

# add_case('v2.LR.historical_0101',                                  n='E3SMv2',  c='cyan', p=tmp_path_hst_v2, s='archive/atm/hist')
add_case('E3SM.INCITE2023-CPL.ne30pg2_EC30to60E2r2.WCYCL20TR-MMF1',n='E3SM-MMF',c='blue', p=tmp_path_hst_mmf,s='archive/atm/hist')

add_case('E3SM.2023-CO2-TEST-01.GNUGPU.ne30pg2_EC30to60E2r2.WCYCL1850-MMF1.1xCO2',n='1xCO2',c='blue' ,p=tmp_path_co2_mmf,s=tmp_sub)
add_case('E3SM.2023-CO2-TEST-01.GNUGPU.ne30pg2_EC30to60E2r2.WCYCL1850-MMF1.2xCO2',n='2xCO2',c='green',p=tmp_path_co2_mmf,s=tmp_sub)
add_case('E3SM.2023-CO2-TEST-01.GNUGPU.ne30pg2_EC30to60E2r2.WCYCL1850-MMF1.4xCO2',n='4xCO2',c='red'  ,p=tmp_path_co2_mmf,s=tmp_sub)
#-------------------------------------------------------------------------------

fig_file,fig_type = 'figs/F12-map-CLD','png'
tmp_file_head = 'F12-map-CLD'

add_var('TGCLDLWP',name='LWP')
add_var('TGCLDIWP',name='IWP')
add_var('SWCF'    ,name='SWCF')
add_var('LWCF'    ,name='LWCF')

#-------------------------------------------------------------------------------
# lat1,lat2 = -60,60

htype = 'ha'
# htype,yr1,yr2 = 'ha', 0, 120
# htype,first_file,num_files = 'ha',0,1
# htype,first_file,num_files = 'ha',100,20

use_remap,remap_grid = False,'90x180'

plot_diff,add_diff = True,False

recalculate = False

print_stats          = True
var_x_case           = False
num_plot_col         = 1 # len(case)
use_common_label_bar = False

#---------------------------------------------------------------------------------------------------
num_var,num_case = len(var),len(case)
diff_base = 0
if 'first_file'      not in locals(): first_file = 0
if 'num_files'       not in locals(): num_files = 0
if 'scrip_file_path' not in locals(): scrip_file_path = None
#---------------------------------------------------------------------------------------------------
wkres = ngl.Resources()
npix = 2048; wkres.wkWidth,wkres.wkHeight=npix,npix
wks = ngl.open_wks(fig_type,fig_file,wkres)
if plot_diff and add_diff: 
   plot = [None]*(num_var*(num_case*2-1))
else:
   plot = [None]*(num_var*num_case)

res = hs.res_contour_fill_map()
# res.tmYLLabelFontHeightF         = 0.01
# res.tmXBLabelFontHeightF         = 0.01
res.lbLabelFontHeightF           = 0.014
res.tmXBOn                       = False
res.tmYLOn                       = False
res.pmTickMarkDisplayMode        = 'Never'
res.mpProjection                 = 'Robinson'

if 'lat1' in locals() : res.mpMinLatF,res.mpMaxLatF = lat1, lat2
if 'lon1' in locals() : res.mpMinLonF,res.mpMaxLonF = lon1, lon2

subtitle_font_height = 0.01
#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------
for v in range(num_var):
   hc.printline()
   print('  var: '+hc.tcolor.MAGENTA+var[v]+hc.tcolor.ENDC)
   data_list = []
   glb_avg_list = []
   lat_list,lon_list = [],[]
   if 'lev_list' in locals(): lev = lev_list[v]
   for c in range(num_case):
      print('    case: '+hc.tcolor.GREEN+case[c]+hc.tcolor.ENDC)

      
      #-------------------------------------------------------------------------
      tmp_file = f'data/{tmp_file_head}.{case[c]}.{var[v]}.nc'
      scrip_ds = xr.open_dataset('grid_files/ne30pg2_scrip.nc')
      area = scrip_ds['grid_area'].rename({'grid_size':'ncol'})
      if recalculate:
      
         file_path = f'{case_dir[c]}/{case[c]}/{case_sub[c]}/*.eam.{htype}.*'
         # file_list_all = sorted(glob.glob(file_path))
         file_list = sorted(glob.glob(file_path))
         # if first_file!=0: file_list = file_list[first_file:]
         # if num_files !=0: file_list = file_list[:num_files]
         # file_list = []
         # for f in range(len(file_list_all)):
         #    if htype=='ha':
         #       yr = int(file_list_all[f][-7:-7+4])
         #       if yr>=yr1 and yr<=yr2: file_list.append(file_list_all[f])
         #    if htype=='h0':
         #       yrmn1,yrmn2 = yr1*1e2+mn1,yr2*1e2+mn2
         #       yr = int(file_list_all[f][-10:-10+4])
         #       mn = int(file_list_all[f][-5:-5+2])
         #       yrmn = yr*1e2 + mn
         #       if yrmn>=yrmn1 and yrmn<=yrmn2:
         #          file_list.append(file_list_all[f])
         # if file_list==[]: exit(f'\nERROR: no files found for file_path:\n{file_path}\n')
         ds = xr.open_mfdataset( file_list )
         data = ds[var[v]]
         #-------------------------------------------------------------------------
         # average over time dimension
         hc.print_time_length(data.time,indent=' '*6)
         data = data.mean(dim='time')
         #----------------------------------------------------------------------
         # Write to file 
         if os.path.isfile(tmp_file) : os.remove(tmp_file)
         tmp_ds = xr.Dataset()
         tmp_ds[var[v]] = data
         tmp_ds['area'] = area
         tmp_ds.to_netcdf(path=tmp_file,mode='w')
      else:
         tmp_ds = xr.open_dataset( tmp_file )
         data = tmp_ds[var[v]]
         # area = tmp_ds['area']
      #-------------------------------------------------------------------------
      if print_stats: hc.print_stat(data,name=var[v],stat='naxsh',indent=' '*6,compact=True)
      #-------------------------------------------------------------------------
      # Calculate area weighted global mean
      if 'area' in locals() :
         gbl_mean = ( (data*area).sum() / area.sum() ).values 
         print(hc.tcolor.CYAN+' '*6+f'Area Weighted Global Mean : {gbl_mean:6.4}'+hc.tcolor.ENDC)
         glb_avg_list.append(gbl_mean)
      #-------------------------------------------------------------------------
      # append to data lists
      data_list.append( data.values )
      # lat_list.append(ds['lat'].isel(time=0,missing_dims='ignore').values)
      # lon_list.append(ds['lon'].isel(time=0,missing_dims='ignore').values)
      lat_list.append(scrip_ds['grid_center_lat'].values)
      lon_list.append(scrip_ds['grid_center_lon'].values)
      #-------------------------------------------------------------------------
      # save baseline for diff map
      if plot_diff :
         if c==diff_base:
            data_baseline = data.copy().values
   #----------------------------------------------------------------------------
   # calculate common limits for consistent contour levels
   data_min = np.min([np.nanmin(d) for d in data_list])
   data_max = np.max([np.nanmax(d) for d in data_list])

   if plot_diff:
      rmse_list = []
      tmp_data_list = []
      for c in range(num_case): 
         rmse_tmp = np.sqrt( np.mean( np.square( data_list[c] - data_baseline )))
         rmse_list.append(rmse_tmp)
         tmp_data_list.append( data_list[c] - data_baseline )
      diff_data_min = np.min([np.nanmin(d) for d in tmp_data_list])
      diff_data_max = np.max([np.nanmax(d) for d in tmp_data_list])
   #----------------------------------------------------------------------------
   # Plot averaged data
   for c in range(num_case):

      data_dir_tmp,data_sub_tmp = None, None
      if use_remap: data_sub_tmp = f'data_remap_{remap_grid}/'
      if case_dir[c] is not None: data_dir_tmp = case_dir[c]
      if case_sub[c] is not None: data_sub_tmp = case_sub[c]
      
      #-------------------------------------------------------------------------
      # Set color palette
      tres = copy.deepcopy(res)
      tres.cnFillPalette = np.array( cmocean.cm.balance(np.linspace(0,1,256)) )
      #-------------------------------------------------------------------------
      # Set contour levels
      # tres.cnLevelSelectionMode = 'ExplicitLevels'
      # tres.cnLevels = np.arange(-110,110+20,20)
      #-------------------------------------------------------------------------
      # tres.cnFillMode = 'RasterFill'
      if not np.all(lat_list[c]==None):
         tres.sfXArray = lon_list[c]
         tres.sfYArray = lat_list[c]

      tres.lbLabelBarOn = False if use_common_label_bar else True      

      num_case_alt = num_case*2-1 if (plot_diff and add_diff) else num_case
      ip = v*num_case_alt+c if var_x_case else c*num_var+v

      if not plot_diff  or (plot_diff and add_diff) or (plot_diff and c==diff_base) : 

         plot[ip] = ngl.contour_map(wks,np.ma.masked_invalid(data_list[c]),tres) 
         
         # set plot subtitles
         # ctr_str = f'{glb_avg_list[c]:6.4}' if glb_avg_list!=[] else ''
         ctr_str = ''
         # if c==0: ctr_str = date_str

         # hs.set_subtitles( wks, plot[ip], 
         #                   left_string=name[c],
         #                   center_string=ctr_str,
         #                   right_string=var_str[v], 
         #                   center_sub_string=f'mean: {glb_avg_list[c]:5.2f} W/m2', 
         #                   font_height=subtitle_font_height)

         hs.set_subtitles( wks, plot[ip], 
                           left_string=name[c],
                           center_string=f'mean: {glb_avg_list[c]:5.2f} W/m2',
                           right_string=var_str[v], 
                           right_sub_string='',#f'{yr1}-{yr2}',
                           font_height=subtitle_font_height)

      #-------------------------------------------------------------------------
      # create difference plot
      if plot_diff and c!=diff_base :
         
         data_list[c] = data_list[c] - data_baseline
         
         # tres.cnFillPalette = 'BlueWhiteOrangeRed'
         tres.cnFillPalette = np.array( cmocean.cm.diff(np.linspace(0,1,256)) )
         # tres.cnFillPalette = np.array( cmocean.cm.delta(np.linspace(0,1,256)) )
         # tres.cnFillPalette = np.array( cmocean.cm.balance(np.linspace(0,1,256)) )

         if hasattr(tres,'cnLevels') : del tres.cnLevels
         # tres.cnLevels = np.arange(-55,55+10,10)
         # if var[v]=='SWCF'       : tres.cnLevels = np.arange(-50,50+10,10)
         # if var[v]=='LWCF'       : tres.cnLevels = np.arange(-50,50+10,10)
         if not hasattr(tres,'cnLevels') : 
            if np.min(data_list[c])==np.max(data_list[c]) : 
               print(hc.tcolor.RED+'WARNING: Difference is zero!'+hc.tcolor.ENDC)
            else:
               cmin,cmax,cint,clev = ngl.nice_cntr_levels(diff_data_min, diff_data_max, \
                                                          cint=None, max_steps=21, \
                                                          returnLevels=True, aboutZero=True )
               tres.cnLevels = np.linspace(cmin,cmax,num=21)
         
         tres.cnLevelSelectionMode = "ExplicitLevels"
         # tres.cnLevelSelectionMode = "AutomaticLevels" # override the level settings and just use auto

         ipd = ip
         if add_diff and     var_x_case: ipd = ip+(num_case-1)
         if add_diff and not var_x_case: ipd = ip+num_var*(num_case-1)

         tres.lbLabelBarOn = True

         plot[ipd] = ngl.contour_map(wks,np.ma.masked_invalid(data_list[c]),tres)
         
         ctr_str = ''
         if c>0: 
            ctr_str = f'RMSE: {rmse_list[c]:5.2f} W/m2'
         glb_diff = glb_avg_list[c] - glb_avg_list[diff_base]
         
         # hs.set_subtitles( wks, plot[ipd], 
         #                   left_string=name[c], 
         #                   center_string=ctr_str,
         #                   right_string=f'{var_str[v]} (diff)',
         #                   center_sub_string=f'mean: {glb_diff:5.2f} W/m2',
         #                   font_height=subtitle_font_height)

         hs.set_subtitles( wks, plot[ipd], 
                           left_string=name[c], 
                           center_string=f'mean: {glb_diff:5.2f} W/m2',
                           right_string=f'{var_str[v]} Bias',
                           right_sub_string='',#f'{yr1}-{yr2}',
                           font_height=subtitle_font_height)

         # lstr = f'{name[c]} - {name[0]}'
         # hs.set_subtitles(wks, plot[ipd], lstr, '', var_str[v], font_height=subtitle_font_height)

#---------------------------------------------------------------------------------------------------
# Finalize plot

num_case_alt = num_case*2-1 if (plot_diff and add_diff) else num_case
layout = [num_var,num_case_alt] if var_x_case else [num_case_alt,num_var]

if not (plot_diff and add_diff):
   if num_case==1 or num_var==1:
      layout = [int(np.ceil(len(plot)/float(num_plot_col))),num_plot_col]
   
pnl_res = hs.setres_panel()
if use_common_label_bar: pnl_res.nglPanelLabelBar = True
pnl_res.nglPanelYWhiteSpacePercent = 5
### add panel labels
pnl_res.nglPanelFigureStrings            = list(string.ascii_lowercase)
pnl_res.nglPanelFigureStringsJust        = "TopLeft"
pnl_res.nglPanelFigureStringsFontHeightF = 0.012

ngl.panel(wks,plot,layout,pnl_res)
ngl.end()

hc.trim_png(fig_file)
#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------

# v2 is a cleaned up version of v1 - functionality is mostly the same
import os, ngl, subprocess as sp, numpy as np, xarray as xr, dask, copy, string, cmocean, numba, glob
import hapy_common as hc, hapy_E3SM   as he, hapy_setres as hs
host = hc.get_host()
#-------------------------------------------------------------------------------
case,case_name,case_dir,case_sub = [],[],[],[]
clr,dsh = [],[]
def add_case(case_in,n=None,p=None,s=None,g=None,c='black',d=0):
   global name,case,case_dir,case_sub
   tmp_name = case_in if n is None else n
   case.append(case_in); case_name.append(tmp_name); 
   case_dir.append(p); case_sub.append(s); 
   clr.append(c); dsh.append(d)
#-------------------------------------------------------------------------------
var,ovar,var_str,var_unit = [],[],[],[]
def add_var(var_in,ovar_in,name=None,unit=None): 
   var.append(var_in)
   ovar.append(ovar_in)
   var_str.append(var_in if name is None else name)
   var_unit.append(unit)
#-------------------------------------------------------------------------------
scrip_file_path  = os.getenv('HOME')+'/E3SM/data_grid/ne30pg2_scrip.nc'
tmp_path_hst_mmf = '/global/cfs/cdirs/m3312/whannah/2023-CPL'
tmp_path_hst_v2  = '/global/cfs/cdirs/m3312/whannah/e3smv2_historical'
tmp_path_co2_mmf = '/global/cfs/cdirs/m3312/whannah/2023-CPL'
# obs_file_path = '/global/cfs/cdirs/e3sm/diagnostics/observations/Atm/climatology_1985-2014/ERA5/ERA5_ANN_198501_201412_climo.nc'
obs_file_path = '/global/cfs/cdirs/m3312/whannah/2023-CPL/ERA5/ERA5_ANN_198501_201412_climo.ne30pg2.nc'
tmp_sub = 'archive/atm/hist'
# add_case('ERA5',                                                     n='ERA5')
add_case('CERES-EBAF',n='CERES-EBAF')
add_case('E3SM.INCITE2023-CPL.ne30pg2_EC30to60E2r2.WCYCL20TR-MMF1',n='E3SM-MMF',c='blue', p=tmp_path_hst_mmf,s='archive/atm/hist')
add_case('v2.LR.historical_0101',                                  n='E3SMv2',  c='red',  p=tmp_path_hst_v2, s='archive/atm/hist')
# add_case('v2.LR.historical_0151',                                  n='E3SMv2',  c='red', p=tmp_path_hst_v2, s='archive/atm/hist')
# add_case('v2.LR.historical_0201',                                  n='E3SMv2',  c='red', p=tmp_path_hst_v2, s='archive/atm/hist')
# add_case('v2.LR.historical_0251',                                  n='E3SMv2',  c='red', p=tmp_path_hst_v2, s='archive/atm/hist')
# add_case('v2.LR.historical_0301',                                  n='E3SMv2',  c='red', p=tmp_path_hst_v2, s='archive/atm/hist')
# add_case('v2.LR.historical',                                       n='E3SMv2 ens',  c='red',  p=None, s=None)
#-------------------------------------------------------------------------------

fig_file,fig_type = 'figs/F05-zonal-mean-CRF','png'
tmp_file_head = 'data/zonal-mean'

# add_var('SWCF',    None,   name='SW CRE', unit='W/m2')
# add_var('LWCF',    None,   name='LW CRE', unit='W/m2')
# add_var('TGCLDLWP','clwvi',name='Liq Water Path', unit='kg/m2')
# add_var('TGCLDIWP','clivi',name='Ice Water Path', unit='kg/m2')

add_var('SWCF',    None, name='SW CRE',              unit='W/m2')
add_var('LWCF',    None, name='LW CRE',              unit='W/m2')
add_var('TGCLDLWP',None, name='Liq Water Path',      unit='kg/m2')
add_var('TGCLDIWP',None, name='Ice Water Path',      unit='kg/m2')
add_var('CLDLOW',  None, name='Low Cloud Fraction',  unit='fraction')
add_var('CLDHGH',  None, name='High Cloud Fraction', unit='fraction')

#-------------------------------------------------------------------------------

# htype,first_file,num_files = 'ha',35,30
# htype,yr1,yr2 = 'ha',1950,2014

htype,yr1,mn1,yr2,mn2 = 'h0',2000,3,2014,12
yrmn1,yrmn2 = yr1*1e2+mn1,yr2*1e2+mn2
# date_str = f'{yr1}/{mn1} - {yr2}/{mn2}'
date_str = f'{yr1}-{yr2}'


recalculate = False

# plot_diff            = False
print_stats          = True
var_x_case           = False

num_plot_col         = 2#len(var)

dlat = 2

#---------------------------------------------------------------------------------------------------
# Set up plot resources
if case==[]: raise ValueError('ERROR - case list is empty!')
num_var,num_case = len(var),len(case)

subtitle_font_height = 0.013


if 'scrip_file_path' not in locals(): scrip_file_path = None

wkres = ngl.Resources()
npix = 2048; wkres.wkWidth,wkres.wkHeight=npix,npix

wks = ngl.open_wks(fig_type,fig_file,wkres)
plot = [None]*(num_var)
# plot = [None]*(num_var*(1+int(plot_diff)))
   
res = hs.res_xy()
res.vpHeightF = 0.2
res.xyLineThicknessF = 10
res.tiXAxisString = 'Latitude'
res.tmYLAutoPrecision = False
res.tmYLPrecision = 2


lres = hs.res_xy()
lres.xyLineThicknessF = 1
lres.xyLineColor      = 'black'
lres.xyDashPattern    = 0

#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------
for v in range(num_var):
   hc.printline()
   print('  var: '+hc.tcolor.MAGENTA+var[v]+hc.tcolor.ENDC)
   lat_list = []
   data_list = []
   glb_avg_list = []
   for c in range(num_case):
      if case[c]=='CERES-EBAF' and var[v] not in ['SWCF','LWCF']: continue
      tmp_file = f'{tmp_file_head}.{case[c]}.{var[v]}.nc'
      print(' '*4+f'case: {hc.tclr.CYAN}{case[c]}{hc.tclr.END}  =>  {tmp_file}')
      if recalculate:
         scrip_ds = xr.open_mfdataset(scrip_file_path).rename({'grid_size':'ncol'})
         area = scrip_ds['grid_area']
         lat  = scrip_ds['grid_center_lat']
         #----------------------------------------------------------------------
         # if case[c]=='ERA5':
         #    scrip_ds = xr.open_dataset(scrip_file_path)
         #    area = scrip_ds['grid_area'].rename({'grid_size':'ncol'})
         #    ds = xr.open_dataset(obs_file_path)
         #    lat = ds['lat']
         #    data = ds[ovar[v]].isel(time=0)
         if case[c]=='CERES-EBAF':
            # obs_root = '/global/cfs/cdirs/e3sm/diagnostics/observations/Atm/climatology'
            # obs_file = f'{obs_root}/ceres_ebaf_toa_v4.1/ceres_ebaf_toa_v4.1_ANN_200101_201812_climo.nc'
            obs_root = '/global/cfs/cdirs/m3312/whannah/obs_data/CERES'
            obs_file = f'{obs_root}/CERES_EBAF-TOA_Ed4.2_Subset_200003-202406.remap_ne30pg2.nc'

            ds = xr.open_dataset(obs_file)
            time_list = []
            for t in range(len(ds['time'])):   
               yr,mn = ds['time.year'][t].values, ds['time.month'][t].values
               yrmn = yr*1e2 + mn
               if yrmn>=yrmn1 and yrmn<=yrmn2: time_list.append(t)
            ds = ds.isel(time=time_list)
            data = None
            if var[v]=='SWCF': data = -1 * ( ds['toa_sw_all_mon'] - ds['toa_sw_clr_c_mon'] )
            if var[v]=='LWCF': data = -1 * ( ds['toa_lw_all_mon'] - ds['toa_lw_clr_c_mon'] )
            data = data.mean(dim='time',skipna=True)
         else:
            # data_dir_tmp,data_sub_tmp = None, None
            # if case_dir[c] is not None: data_dir_tmp = case_dir[c]
            # if case_sub[c] is not None: data_sub_tmp = case_sub[c]
            # case_obj = he.Case( name=case[c], data_dir=data_dir_tmp, data_sub=data_sub_tmp )
            # case_obj.set_coord_names(var[v])
            # #-------------------------------------------------------------------
            # # read the data
            # with dask.config.set(**{'array.slicing.split_large_chunks': True}):
            #    tmp_first_file = first_file
            #    if 'v2.LR.historical' in case[c]: tmp_first_file = 50 + first_file
            #    lat  = case_obj.load_data('lat',  htype=htype)
            #    area = case_obj.load_data('area', htype=htype).astype(np.double)
            #    data = case_obj.load_data(var[v],htype=htype,ps_htype=htype,
            #                                    first_file=tmp_first_file,num_files=num_files,
            #                                    extrap_flag=False)
            #    hc.print_time_length(data.time,indent=' '*6)
            #    data = data.mean(dim='time',skipna=True)
            #-------------------------------------------------------------------
            file_path = f'{case_dir[c]}/{case[c]}/{case_sub[c]}/*.eam.{htype}.*'
            file_list = sorted(glob.glob(file_path))
            ds = xr.open_mfdataset( file_list )
            ds = ds.where( ds['time.year']*1e2 + ds['time.month']>=yrmn1, drop=True)
            ds = ds.where( ds['time.year']*1e2 + ds['time.month']<=yrmn2, drop=True)
            data = ds[var[v]].mean(dim='time',skipna=True)
         #----------------------------------------------------------------------
         # print stats after time averaging
         if print_stats: 
            hc.print_stat(data,name=var[v],stat='naxsh',indent='    ',compact=True)
         #----------------------------------------------------------------------
         # Calculate time and zonal mean
         bin_ds = hc.bin_YbyX( data, lat, bin_min=-90, bin_max=90, bin_spc=dlat, wgt=area )
         #----------------------------------------------------------------------
         bin_ds.to_netcdf(path=tmp_file,mode='w')
      else:
         bin_ds = xr.open_dataset( tmp_file, use_cftime=True )
      #----------------------------------------------------------------------
      # adjust units
      # if var[v]=='Q': data = data*1e3
      if var[v]=='TGCLDLWP': bin_ds['bin_val'] = bin_ds['bin_val']*1e3; var_unit[v] = 'g/m2'
      if var[v]=='TGCLDIWP': bin_ds['bin_val'] = bin_ds['bin_val']*1e3; var_unit[v] = 'g/m2'
      if var[v]=='CLDLOW'  : bin_ds['bin_val'] = bin_ds['bin_val']*1e2; var_unit[v] = '%'
      if var[v]=='CLDHGH'  : bin_ds['bin_val'] = bin_ds['bin_val']*1e2; var_unit[v] = '%'
      #-------------------------------------------------------------------------
      hc.print_stat(bin_ds['bin_val'],name=var[v],stat='naxsh',indent='    ',compact=True)
      #-------------------------------------------------------------------------
      data_list.append( bin_ds['bin_val'].values )
      lat_list.append( bin_ds['bins'].values )
   #----------------------------------------------------------------------------
   # Plot averaged data

   tres = copy.deepcopy(res)
   
   tres.trXMinF = np.min([np.nanmin(d) for d in  lat_list])
   tres.trXMaxF = np.max([np.nanmax(d) for d in  lat_list])
   tres.trYMinF = np.min([np.nanmin(d) for d in data_list])
   tres.trYMaxF = np.max([np.nanmax(d) for d in data_list])
   
   lat_tick = np.array([-90,-60,-30,0,30,60,90])
   tres.tmXBMode = "Explicit"
   tres.tmXBLabels = lat_tick

   tres.tiYAxisString = f'[{var_unit[v]}]'

   pcnt = 0
   for c in range(num_case):
      if case[c]=='CERES-EBAF' and var[v] not in ['SWCF','LWCF']: continue
      
      tres.xyLineColor   = clr[c]
      tres.xyDashPattern = dsh[c]
      tplot = ngl.xy(wks, lat_list[pcnt], np.ma.masked_invalid( data_list[pcnt] ), tres)
      if pcnt==0:
         plot[v] = tplot
      else:
         ngl.overlay( plot[v], tplot )
      pcnt += 1

   hs.set_subtitles(wks, plot[v], '', '', var_str[v], font_height=subtitle_font_height)

   #----------------------------------------------------------------------------
   # if plot_diff:
   #    baseline = np.ma.masked_invalid( data_list[0] )
   #    for c in range(1,num_case): data_list[c] = data_list[c] - baseline
   #    tres.trYMinF = np.min([np.nanmin(d) for d in data_list[1:]])
   #    tres.trYMaxF = np.max([np.nanmax(d) for d in data_list[1:]])
   #    # ip = num_var+v
   #    ip = v
   #    for c in range(1,num_case):
   #       tres.xyLineColor   = clr[c]
   #       tres.xyDashPattern = dsh[c]
   #       tmp_data = np.ma.masked_invalid( data_list[c] ) #- baseline
   #       tplot = ngl.xy(wks, lat_list[c], tmp_data, tres)
   #       if c==1:
   #          plot[ip] = tplot
   #       else:
   #          ngl.overlay( plot[ip], tplot )
   #    ngl.overlay( plot[ip], ngl.xy(wks, np.array([-1e8,1e8]), np.array([0,0]), lres) )
   #    hs.set_subtitles(wks, plot[ip], '', f'Bias wrt {name[0]}', var_str[v], font_height=subtitle_font_height)
#-------------------------------------------------------------------------------
# Add legend

lgres = ngl.Resources()
lgres.vpWidthF           = 0.03
lgres.vpHeightF          = 0.06
lgres.lgLabelFontHeightF = 0.008
lgres.lgLabelFont        = "courier"
lgres.lgMonoDashIndex    = False
lgres.lgLineLabelsOn     = False
lgres.lgLineThicknessF   = 20
lgres.lgLabelJust        = 'CenterLeft'
lgres.lgLineColors       = clr
lgres.lgDashIndexes      = dsh

indent = ' '*2
labels = case_name
for i in range(len(labels)): labels[i] = indent+labels[i] 

# pid = ngl.legend_ndc(wks, len(labels), labels, 0.2, 0.75, lgres) # 2x2
pid = ngl.legend_ndc(wks, len(labels), labels, 0.2, 0.785, lgres) # 3x2
#---------------------------------------------------------------------------------------------------
# Finalize plot

# layout = [num_var,1]
layout = [int(np.ceil(len(plot)/float(num_plot_col))),num_plot_col]
   
pnl_res = hs.setres_panel()

### add panel labels
pnl_res.nglPanelFigureStrings            = list(string.ascii_lowercase)
pnl_res.nglPanelFigureStringsJust        = "TopLeft"
pnl_res.nglPanelFigureStringsFontHeightF = 0.01

pnl_res.nglPanelYWhiteSpacePercent = 5

ngl.panel(wks,plot,layout,pnl_res)
ngl.end()

hc.trim_png(fig_file)
#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------

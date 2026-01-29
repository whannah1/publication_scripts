# v2 is a cleaned up version of v1 - functionality is mostly the same
import os, glob, ngl, subprocess as sp, numpy as np, xarray as xr, dask, copy, string, cmocean, pandas as pd
import hapy_common as hc, hapy_E3SM   as he, hapy_setres as hs
import scipy
from statsmodels.tsa.arima.model import ARIMA
host = hc.get_host()
#-------------------------------------------------------------------------------
''' mapping commands

ncremap -G ttl='BEST grid'#latlon==180,360#lat_typ=uni#lon_typ=180_wst -g grid_files/BEST_180x360_scrip.20241001.nc
SRC_GRID_FILE=grid_files/ne30pg2_scrip.nc
DST_GRID_FILE=grid_files/BEST_180x360_scrip.20241001.nc
MAP_FILE=map_files/map_ne30pg2_to_180x360_BEST_traave_20241001.nc
ncremap -a traave --src_grd=${SRC_GRID_FILE} --dst_grd=${DST_GRID_FILE} --map_file=${MAP_FILE}

SRC_GRID_FILE=grid_files/BEST_180x360_scrip.20241001.nc
DST_GRID_FILE=grid_files/ne30pg2_scrip.nc
MAP_FILE=map_files/map_180x360_BEST_to_ne30pg2_traave_20241001.nc
ncremap -a traave --src_grd=${SRC_GRID_FILE} --dst_grd=${DST_GRID_FILE} --map_file=${MAP_FILE}

SRC_FILE=/global/cfs/cdirs/m3312/whannah/obs_data/BEST/Land_and_Ocean_LatLong1.nc
DST_FILE=/global/cfs/cdirs/m3312/whannah/obs_data/BEST/Land_and_Ocean_LatLong1.remap_ne30pg2.nc
MAP_FILE=map_files/map_180x360_BEST_to_ne30pg2_traave_20241001.nc
ncremap -m ${MAP_FILE} -i ${SRC_FILE} -o ${DST_FILE}

ncremap -G ttl='HadCRU grid'#latlon==36,72#lat_typ=uni#lon_typ=grn_ctr -g grid_files/HadCRU_36x72_scrip.20241001.nc
SRC_GRID_FILE=grid_files/ne30pg2_scrip.nc
DST_GRID_FILE=grid_files/HadCRU_36x72_scrip.20241001.nc
MAP_FILE=map_files/map_ne30pg2_to_36x72_traave_20241001.nc
ncremap -a traave --src_grd=${SRC_GRID_FILE} --dst_grd=${DST_GRID_FILE} --map_file=${MAP_FILE}

DATA_ROOT=/global/cfs/cdirs/m3312/whannah/obs_data/HadCRU
SRC_DATA=${DATA_ROOT}/HadCRUT.5.0.1.0.analysis.anomalies.ensemble_mean.nc
DST_DATA=${DATA_ROOT}/HadCRUT.5.0.1.0.analysis.anomalies.ensemble_mean.remap_ne30pg2.nc
MAP_FILE=map_files/map_36x72_to_ne30pg2_traave_20241001.nc
ncremap -m ${MAP_FILE} -i ${SRC_DATA} -o ${DST_DATA}
'''
#-------------------------------------------------------------------------------
def run_cmd(cmd):
   msg = hc.tcolor.GREEN + cmd + hc.tcolor.ENDC ; print(f'\n{msg}')
   os.system(cmd); return
#-------------------------------------------------------------------------------
case_name,case,case_dir,case_sub = [],[],[],[]
clr,dsh = [],[]
obs_flag = []
def add_case(case_in,n=None,p=None,s=None,g=None,c='black',d=0,obs=False):
   global name,case,case_dir,case_sub
   tmp_name = case_in if n is None else n
   case.append(case_in); case_name.append(tmp_name); 
   case_dir.append(p); case_sub.append(s);
   obs_flag.append(obs)
   dsh.append(d) ; clr.append(c)
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

add_case('BEST', c='black',obs=True)

# add_case('v2.LR.historical_0101',                                  n='E3SMv2',  c='cyan', p=tmp_path_hst_v2, s='archive/atm/hist')
# add_case('v2.LR.historical',                                       n='E3SMv2 ens',  c='red',  p=None, s=None)
# add_case('E3SM.INCITE2023-CPL.ne30pg2_EC30to60E2r2.WCYCL20TR-MMF1',n='E3SM-MMF',c='blue', p=tmp_path_hst_mmf,s='archive/atm/hist')

#-------------------------------------------------------------------------------

fig_file,fig_type = 'figs/F08-ENSO-wavelet-spectra','png'
tmp_file_head = 'ENSO-power-spectra'

ens_spread_color = 'coral'

add_var('TS',name='Tsfc')

# region,lat1,lat2,lon1,lon2 = 'Nino3'  ,-5,5,360-150,360-90
# region,lat1,lat2,lon1,lon2 = 'Nino4'  ,-5,5,160,360-150
region,lat1,lat2,lon1,lon2 = 'Nino3.4',-5,5,190,240

#-------------------------------------------------------------------------------
# htype,yr1,yr2 = 'h0',1950,2014
htype,yr1,yr2 = 'h0',1980,2014
# date_str1 = f'{yr1}-{yr2} mean'
# date_str2 = f'{yr2}-{yr1} diff'

# plot_diff = True

recalculate_timeseries = False

print_stats          = True
var_x_case           = False
num_plot_col         = 1
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
# plot = [None]*(num_var*2*num_case)
plot = [None]*num_var

# res = hs.res_contour_fill_map()
# res.lbLabelFontHeightF           = 0.012
# res.tmXBOn                       = False
# res.tmYLOn                       = False
# res.mpCenterLonF                 = 0
# res.pmTickMarkDisplayMode        = 'Never'
# res.mpProjection                 = 'Robinson'

res = hs.res_xy()
res.vpHeightF = 0.5
res.xyMarkerSizeF = 0.008
res.xyMarker = 16
res.xyLineThicknessF = 12
# res.tmYLLabelFontHeightF         = 0.008
# res.tmXBLabelFontHeightF         = 0.008
# res.tmXBAutoPrecision = False
# res.tmXBPrecision = 2
res.tiXAxisString = 'Period [years]'
res.tiYAxisString = 'Variance [~F34~0~F~C~S~2~N~]'
tm_log = np.array([1,2,4,8,16,32,64])
res.tmXBMode      = 'Explicit'
res.tmXBValues    = tm_log
res.tmXBLabels    = tm_log

res.xyXStyle ='Log'


subtitle_font_height = 0.01
#---------------------------------------------------------------------------------------------------
def deseason(xraw):
    # Calculates the deseasonalized data
    months_per_year = 12
    # Create array to hold climatological values and deseasonalized data
    # Create months_per_year x 1 array of zeros
    xclim = np.zeros((months_per_year, 1))
    # Create array with same shape as xraw
    x_deseasoned = np.zeros(xraw.shape)
    # Iterate through all 12 months.
    for month in np.arange(months_per_year):
        # `xraw[month::12]` will return the data for this month every year (12 months)
        # (i.e., from month until the end of xraw, get every 12th month)
        # Get the mean of this month, using data from every year, ignoring NaNs
        xclim[month] = np.nanmean(xraw[month::months_per_year])
    num_years = int(np.floor(len(x_deseasoned) / months_per_year))
    # Iterate through all years in x_deseasoned (same number as in xraw)
    for year in np.arange(num_years):
        year_index = year * months_per_year
        # Iterate through all months of the year
        for month in np.arange(months_per_year):
            month_index = year_index + month
            print()
            print(x_deseasoned.shape)
            print()
            print(xraw.shape)
            print()
            print(xclim.shape)
            print()
            # Subtract the month's mean over num_years from xraw's data for this month in this year
            # i.e., get the difference between this month's value and it's "usual" value
            x_deseasoned[month_index] = xraw[month_index] - xclim[month]
    return x_deseasoned
#---------------------------------------------------------------------------------------------------
# def get_psd_from_wavelet(data):
#    deg = 6
#    period = np.arange(6,12*64+1)
#    freq = 1/period
#    widths = deg / (2*np.pi*freq)
#    cwtmatr = scipy.signal.cwt( deseason(data), scipy.signal.morlet2, widths=widths, w=deg )
#    psd = np.mean(np.square(np.abs(cwtmatr)),axis=1)
#    return ( period, psd )
longest_period = 20*12
def get_psd_from_wavelet(data):
   deg = 6
   period = np.arange(1, longest_period + 1)
   freq = 1 / period
   widths = deg / (2 * np.pi * freq)
   cwtmatr = scipy.signal.cwt(data, scipy.signal.morlet2, widths=widths, w=deg)
   psd = np.mean( np.square( np.abs(cwtmatr) ), axis=1)
   return (period, psd)
#---------------------------------------------------------------------------------------------------
def get_alpha_sigma2(data_in):
   ## Fit AR1 model to estimate the autocorrelation (alpha) and variance (sigma2)
   mod = ARIMA( data_in, order=(1,0,0) )
   res = mod.fit()
   return (res.params[1],res.params[2])
#---------------------------------------------------------------------------------------------------
def get_tmp_file(case,var):
   return f'data/{tmp_file_head}.{case}.{var}.reg_{region}.lat1_{lat1}.lat2_{lat2}.lon1_{lon1}.lon2_{lon2}.nc'
#---------------------------------------------------------------------------------------------------
for v in range(num_var):
   hc.printline()
   print('  var: '+hc.tcolor.MAGENTA+var[v]+hc.tcolor.ENDC)
   data_list = []
   # glb_avg_list = []
   # lat_list,lon_list = [],[]
   # if 'lev_list' in locals(): lev = lev_list[v]
   period_list, power_list = [],[]
   variance_list = []
   alpha_list = []
   sigma2_list = []
   #----------------------------------------------------------------------------
   scrip_ds = xr.open_dataset('grid_files/ne30pg2_scrip.nc').rename({'grid_size':'ncol'})
   scrip_ds = scrip_ds.rename({'grid_center_lat':'lat','grid_center_lon':'lon'})
   area,lat,lon = scrip_ds['grid_area'],scrip_ds['lat'],scrip_ds['lon']
   mask = xr.full_like(area,True,dtype=bool)
   mask = mask & (lat>=lat1) & (lat<=lat2)
   mask = mask & (lon>=lon1) & (lon<=lon2)
   area = area.where( mask, drop=True)
   #----------------------------------------------------------------------------
   for c in range(num_case):
      print('    case: '+hc.tcolor.GREEN+case[c]+hc.tcolor.ENDC)
      #-------------------------------------------------------------------------
      # tmp_file = f'data/{tmp_file_head}.{case[c]}.{var[v]}.nc'
      tmp_file = get_tmp_file(case[c],var[v])
      if recalculate_timeseries:
         #----------------------------------------------------------------------
         # if case[c]=='ERA5':
         #    obs_root = '/global/cfs/cdirs/m3312/whannah/2023-CPL/ERA5'
         #    obs_file = f'{obs_root}/ERA5_ANN_198501_201412_climo.ne30pg2.nc'
         #    ds = xr.open_dataset(obs_file)
         #    data = ds['ts'] 
         #    print()
         #    print(data)
         #    print()
         #    exit()
         if case[c]=='BEST':
            obs_file = '/global/cfs/cdirs/m3312/whannah/obs_data/BEST/Land_and_Ocean_LatLong1.remap_ne30pg2.nc'
            ds = xr.open_dataset(obs_file)#.rename({'latitude':'lat','longitude':'lon'})
            data = ds['temperature']
            time_list = [None]*len(data['time'])
            # convert anomalies to absolute temperature and build time coord
            for t,tt in enumerate(data['time'].values):
               yr_val = int(np.floor(tt))
               mn_ind = int(np.floor((tt-yr_val)*12))
               time_list[t] = np.datetime64(f'{yr_val}-{(mn_ind+1):02d}')
               data[t,:] = data[t,:] + ds['climatology'].isel(month_number=mn_ind)
            time = xr.DataArray(time_list,dims=('time'))
            time = time.assign_coords(time=time)
            month_length = time.dt.days_in_month
            wgts = month_length.groupby("time.year") / month_length.groupby("time.year").sum()
            data = data.assign_coords(time=time)
            # data = data.resample(time='Y').mean(dim='time')
            # time = time.resample(time='Y').mean(dim='time')
            
            # subset in time
            t_beg_found, t_end_found = False, False
            for t,y in enumerate(data['time.year'].values):
               if not t_beg_found and y==yr1: t_beg=t; t_beg_found = True
               # if not t_end_found and y==yr2: t_end=t; t_end_found = True
               if y==yr2: t_end=t
            data = data.isel(time=slice(t_beg,t_end+1))            

         else:
            file_path = f'{case_dir[c]}/{case[c]}/{case_sub[c]}/*.eam.{htype}.*'
            file_list_all = sorted(glob.glob(file_path))
            file_list1 = [] # data for anomalies
            file_list2 = [] # data for baseline
            for f in range(len(file_list_all)):
               yr = int(file_list_all[f][-7:-7+4])
               if yr>= yr1 and yr<= yr2: file_list1.append(file_list_all[f])
            if file_list1==[]: exit(f'\nERROR: no files found for file_path:\n{file_path}\n')   
            ds_main = xr.open_mfdataset( file_list1 )
            # ds_base = xr.open_mfdataset( file_list2 )
            # data = ds_main[var[v]] - ds_base[var[v]].mean(dim='time')
            data = ds_main[var[v]] - 273.15
         #----------------------------------------------------------------------
         data = data.where( mask, drop=True)
         data_avg = ( (data*area).sum(dim='ncol') / area.sum(dim='ncol') )
         #----------------------------------------------------------------------
         # Write to file 
         if os.path.isfile(tmp_file) : os.remove(tmp_file)
         tmp_ds = xr.Dataset()
         tmp_ds[var[v]] = data
         tmp_ds.to_netcdf(path=tmp_file,mode='w')
      #-------------------------------------------------------------------------
      tmp_ds = xr.open_dataset( tmp_file )
      data = tmp_ds[var[v]]
      #-------------------------------------------------------------------------
      # data = data + 273.15
      #-------------------------------------------------------------------------
      if print_stats: 
         hc.print_stat(data,name=var[v],stat='naxsh',indent=' '*6,compact=True)
      #-------------------------------------------------------------------------
      #-------------------------------------------------------------------------
      if case[c]!='v2.LR.historical':
         ( period, psd ) = get_psd_from_wavelet(data_avg.values)

         variance_list.append(np.var(data_avg.values))

         (alpha,sigma2) = get_alpha_sigma2(data_avg.values)
         alpha_list.append(alpha)
         sigma2_list.append(sigma2)

      #-------------------------------------------------------------------------
      period_list.append( period/12. )
      power_list .append( psd )
   #----------------------------------------------------------------------------
   # # calculate common limits for consistent contour levels
   # data_min = np.min([np.nanmin(d) for d in data_list])
   # data_max = np.max([np.nanmax(d) for d in data_list])

   # if plot_diff:
   #    rmse_list = []
   #    tmp_data_list = []
   #    for c in range(num_case): 
   #       rmse_tmp = np.sqrt( np.mean( np.square( data_list[c] - data_baseline )))
   #       rmse_list.append(rmse_tmp)
   #       tmp_data_list.append( data_list[c] - data_baseline )
   #    diff_data_min = np.min([np.nanmin(d) for d in tmp_data_list])
   #    diff_data_max = np.max([np.nanmax(d) for d in tmp_data_list])
   #----------------------------------------------------------------------------
   # Create plot
   #----------------------------------------------------------------------------
   tres = copy.deepcopy(res)
   tres.trXMinF = np.min([np.min(d) for d in period_list])
   tres.trXMaxF = np.max([np.max(d) for d in period_list])
   tres.trYMinF = 0
   tres.trYMaxF = np.max([np.max(d) for d in power_list])

   # tres.trYMaxF = 15
   tres.trXReverse = True

   if 'ens_psd_max' in locals(): tres.trYMaxF = np.max([tres.trYMaxF,np.max(ens_psd_max)])
   for c in range(num_case):
      
      tres.xyLineColor,tres.xyDashPattern = clr[c],dsh[c]

      tres.xyLineThicknessF = res.xyLineThicknessF
      # if case[c] in ['HadSST']: tres.xyLineThicknessF = res.xyLineThicknessF * 2
      if 'v2.LR.historical_' in case[c]: tres.xyLineThicknessF = res.xyLineThicknessF / 2

      tplot = ngl.xy(wks, period_list[c], power_list[c], tres)

      ip = v
      if c==0:
         plot[ip] = tplot 
      else:
         ngl.overlay(plot[ip],tplot)

      if case[c]=='v2.LR.historical':
         eres = copy.deepcopy(tres)
         ens_spread_data      = np.zeros([2,len(ens_psd_min)])
         ens_spread_data[0,:] = ens_psd_min
         ens_spread_data[1,:] = ens_psd_max
         eres.xyLineColor = [0,0,0,0]
         eres.nglXYAboveFillColors = ens_spread_color
         eres.nglXYBelowFillColors = ens_spread_color
         ngl.overlay(plot[ip], ngl.xy(wks, period_list[c], ens_spread_data, eres) )

      #-------------------------------------------------------------------------
      # Set strings at top of plot
      lat_chk,lon_chk = 'lat1' in locals(), 'lon1' in locals()
      var_str = var[v]
      # if var[v]=="PRECT" : var_str = "Precipitation [mm/day]"
      # if var[v]=="TMQ"   : var_str = "Column Water Vapor [mm]"
      lft_str = ''
      ctr_str = ''
      if not lat_chk and not lon_chk : 
         ctr_str = 'Global'
      else:
         ctr_str += f' {region} '
         if lat_chk:      ctr_str += f' {lat1}:{lat2}N '
         if lon_chk:      ctr_str += f' {lon1}:{lon2}E '

   hs.set_subtitles(wks, plot[ip], '', region, '', font_height=0.015)

   #-------------------------------------------------------------------------
   # # Function to compute the red noise spectrum
   # def get_P(alpha, period, sigma2):
   #    freq = 1/(period_list[c]*12)
   #    alpha_sq = np.square(alpha)
   #    P = ( 1 - alpha_sq ) / ( 1 + alpha_sq - 2*alpha*np.cos(2*np.pi*freq) )
   #    return  P
   #-------------------------------------------------------------------------
   c = 0 

   lres = hs.res_xy()
   lres.xyDashPattern = 1

   alpha_assumed = 0.72
   red_freq = 1/(period_list[c]*12)
   alpha_sq = np.square(alpha_assumed)
   P = ( 1 - alpha_sq ) / ( 1 + alpha_sq - 2*alpha_assumed*np.cos(2*np.pi*red_freq) )
   # q95 = 0.5 * P * 5.99 * sigma2_list[c]
   chi_val = 5.99 
   p = 0.05
   Ws = 0.5 * P * 5.99 * sigma2_list[c]
   # q95 = 2 / ( chi_val*(1-p/2) ) * Ws
   q95 = Ws* 2 / chi_val

   lres.xyLineColor = 'red'
   ngl.overlay(plot[ip], ngl.xy(wks, period_list[c], P, lres) )
   
   lres.xyLineColor = 'blue'
   ngl.overlay(plot[ip], ngl.xy(wks, period_list[c], q95, lres) )

   #-------------------------------------------------------------------------
   # # indicate 95% confidence level
   # if True:

   #    # period = period_list[0]
   #    # red_freq = 1/(period*12)
      
   #    alpha_assumed = 0.72
   #    # alpha_assumed = 0.90

   #    lres = hs.res_xy()
   #    # lres.xyLineColor = 'gray'
   #    lres.xyLineColor = 'red'
   #    lres.xyDashPattern = 1

   #    # plot red noise spectrum
   #    # alpha_sq = np.square(alpha)
   #    # P = ( 1 - alpha_sq ) / ( 1 + alpha_sq - 2*alpha*np.cos(2*np.pi*red_freq) )
   #    # P = get_P(alpha_list[c], red_freq, sigma2)
   #    # ngl.overlay(plot[ip], ngl.xy(wks, period, P, lres) )

   #    # plot 95% confidence level
   #    for c in range(num_case):
   #       lres.xyDashPattern = 2
   #       lres.xyLineColor = clr[c]
   #       # q95 = P*5.99*0.5*variance_list[c]
   #       # Compute the local wavelet power spectrum upper bound eq. 18, 5.99 = 95-th quantile for a Chi² r.v. with 2 DOF
   #       # alpha = alpha_list[c]
   #       # alpha_sq = np.square(alpha)
   #       # P = ( 1 - alpha_sq ) / ( 1 + alpha_sq - 2*alpha*np.cos(2*np.pi*red_freq) )

   #       P = get_P(alpha_assumed, period_list[c], sigma2_list[c])
   #       # P = get_P(alpha_list[c], period_list[c], sigma2_list[c])
   #       lres.xyLineColor = 'red'
   #       ngl.overlay(plot[ip], ngl.xy(wks, period_list[c], P, lres) )
         
   #       q95 = 0.5 * P * 5.99 * sigma2_list[c]
   #       lres.xyLineColor = 'blue'
   #       ngl.overlay(plot[ip], ngl.xy(wks, period_list[c], q95, lres) )

#----------------------------------------------------------------------------
# Add legend
#----------------------------------------------------------------------------
lgres = ngl.Resources()
lgres.vpWidthF, lgres.vpHeightF  = 0.08, 0.1  
lgres.lgLabelFontHeightF = 0.015
lgres.lgLineThicknessF   = 12
lgres.lgMonoLineColor    = False
lgres.lgMonoDashIndex    = True
lgres.lgDashIndex        = 0
lgres.lgLineColors       = clr
lgres.lgLabelJust    = 'CenterLeft'
lbl = [f'  {n}' for n in case_name]
# pid = ngl.legend_ndc(wks, len(case_name), lbl, 0.18, 0.68, lgres)

#---------------------------------------------------------------------------------------------------
# Finalize plot
#---------------------------------------------------------------------------------------------------

# layout = [len(plot),1]
# layout = [num_var,num_case] if var_x_case else [num_case,num_var]

layout = [int(np.ceil(len(plot)/float(num_plot_col))),num_plot_col]

pnl_res = hs.setres_panel()
pnl_res.nglPanelYWhiteSpacePercent = 5

### add panel labels
# pnl_res.nglPanelFigureStrings            = list(string.ascii_lowercase)
# pnl_res.nglPanelFigureStringsJust        = "TopLeft"
# pnl_res.nglPanelFigureStringsFontHeightF = 0.01

# if num_var==1  : layout = [num_case,num_var]
# if num_case==1 : layout = [num_var,num_case]

ngl.panel(wks,plot[0:len(plot)],layout,pnl_res)
ngl.end()



hc.trim_png(fig_file)
#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------

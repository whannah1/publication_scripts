import os, sys, ngl, subprocess as sp, numpy as np, xarray as xr, copy, string, dask, glob
import hapy_common as hc, hapy_E3SM   as he, hapy_setres as hs
import cmocean
import scipy, pywt
from statsmodels.tsa.arima.model import ARIMA
sys.path.append(os.getcwd())
import QBO_diagnostic_methods as QBO_methods
#-------------------------------------------------------------------------------
# based on E3SM diagnostics package:
# https://github.com/E3SM-Project/e3sm_diags/blob/main/e3sm_diags/driver/qbo_driver.py
#-------------------------------------------------------------------------------
case_name,case,case_dir,case_sub,case_grid,clr,dsh,mrk = [],[],[],[],[],[],[],[]
yr1_list = []
yr2_list = []
def add_case(case_in,n=None,p=None,s=None,g=None,d=0,c='black',m=0,yr1=None,yr2=None):
   global name,case,case_dir,case_sub,clr,dsh,mrk
   if n is None:
      tmp_name = ''
   else:
      tmp_name = n
   case.append(case_in); case_name.append(tmp_name)
   case_dir.append(p); case_sub.append(s); case_grid.append(g)
   dsh.append(d) ; clr.append(c) ; mrk.append(m)
   yr1_list.append(yr1)
   yr2_list.append(yr2)
#-------------------------------------------------------------------------------

add_case('ERA5', n='ERA5', c='black')

tmp_path,tmp_sub = '/global/cfs/cdirs/e3smdata/simulations/','archive/atm/hist'
add_case('v2.LR.amip_0101', n='E3SMv2 AMIP 101', d=1, c='red',   p=tmp_path,s=tmp_sub) 
add_case('v2.LR.amip_0201', n='E3SMv2 AMIP 201', d=1, c='green', p=tmp_path,s=tmp_sub)
add_case('v2.LR.amip_0301', n='E3SMv2 AMIP 301', d=1, c='blue',  p=tmp_path,s=tmp_sub)

tmp_path,tmp_sub = '/global/cfs/cdirs/m3312/whannah/e3smv3_amip','archive/atm/hist'
add_case('v3.LR.amip_0101', n='E3SMv3 AMIP 101', d=0, c='red',   p=tmp_path,s=tmp_sub)
add_case('v3.LR.amip_0151', n='E3SMv3 AMIP 151', d=0, c='green', p=tmp_path,s=tmp_sub)
add_case('v3.LR.amip_0201', n='E3SMv3 AMIP 201', d=0, c='blue',  p=tmp_path,s=tmp_sub)

#-------------------------------------------------------------------------------

var = ['U']
pow_spec_lev = 20.

fig_file,fig_type = 'figs/QBO-wavelet-spectra','png'
tmp_file_head     = 'data/QBO.wavelet_spectra'

lat1,lat2 = -5,5

# htype,first_file,num_files = 'h0', (1979-1870)*12, 12*40#65
# htype,yr1,yr2 = 'h0',1995,1999

longest_period = 5*12

print_stats = True
var_x_case = False
use_common_label_bar = True

num_plot_col = 1

recalculate_timeseries = False

# year_start = 1950+first_file/12
# year_end   = 1950+first_file/12+num_files/12+1

#---------------------------------------------------------------------------------------------------
# Set up plot resources
#---------------------------------------------------------------------------------------------------
num_var,num_case = len(var),len(case)

if 'lev' not in vars(): lev = np.array([0])

wkres = ngl.Resources()
npix = 2048; wkres.wkWidth,wkres.wkHeight=npix,npix
wks = ngl.open_wks(fig_type,fig_file,wkres)
plot = [None]*num_var*1

res = hs.res_xy()
res.vpHeightF = 0.3
# res.xyMarkLineMode = "MarkLines"
res.xyMarkerSizeF = 0.008
res.xyMarker = 16
res.xyLineThicknessF = 12
res.tiXAxisString = 'Period [months]'
res.tiYAxisString = 'Variance [m~S~2~N~ s~S~-2~N~]'

# tm_log = np.array([6,8,16,24,32,48,64,96])
tm_log = np.array([8,16,24,28,32,64])
res.tmXBMode      = 'Explicit'
res.tmXBValues    = tm_log
res.tmXBLabels    = tm_log

res.xyXStyle ='Log'

#---------------------------------------------------------------------------------------------------
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
            # Subtract the month's mean over num_years from xraw's data for this month in this year
            # i.e., get the difference between this month's value and it's "usual" value
            x_deseasoned[month_index] = xraw[month_index] - xclim[month]
    return x_deseasoned
#---------------------------------------------------------------------------------------------------
def get_psd_from_wavelet_scipy(data):
   deg = 6
   period = np.arange(1, longest_period + 1)
   freq = 1 / period
   widths = deg / (2 * np.pi * freq)
   cwtmatr = scipy.signal.cwt(data, scipy.signal.morlet2, widths=widths, w=deg)
   psd = np.mean( np.square( np.abs(cwtmatr) ), axis=1)
   return (period, psd)
#---------------------------------------------------------------------------------------------------
def get_psd_from_wavelet_pywt(data):
   deg = 6
   period = np.arange(1, longest_period + 1)
   widths = deg / ( 2 * np.pi / period )
   [cfs, freq] = pywt.cwt(data, scales=widths, wavelet='cmor1.5-1.0')
   psd = np.mean( np.square( np.abs(cfs) ), axis=1)
   period = 1 / freq
   return (period, psd)
#---------------------------------------------------------------------------------------------------
def get_alpha_sigma2(data_in):
   ## Fit AR1 model to estimate the autocorrelation (alpha) and variance (sigma2)
   mod = ARIMA( data_in, order=(1,0,0) )
   res = mod.fit()
   return (res.params[1],res.params[2])
#---------------------------------------------------------------------------------------------------
def get_tmp_file(case,var,yr1,yr2):
   # return f'{tmp_file_head}.{case}.{var}.{yr1}-{yr2}.nc'
   return f'{tmp_file_head}.{case}.{var}.nc'
#---------------------------------------------------------------------------------------------------
def mask_data(ds,data,lat_name='lat',lon_name='lon'):
   global lat1,lat2
   
   if 'ncol' in data.dims:
      tmp_data = np.ones([len(ds['ncol'])],dtype=bool)
      tmp_dims,tmp_coords = ('ncol'),{'ncol':data['ncol']}
   else:
      tmp_data = np.ones([len(ds[lat_name]),len(ds[lon_name])],dtype=bool)
      tmp_dims,tmp_coords = (lat_name,lon_name),{lat_name:ds[lat_name],lon_name:ds[lon_name]}
   mask = xr.DataArray( tmp_data, coords=tmp_coords, dims=tmp_dims )
   mask = mask & (ds[lat_name]>= lat1) & (ds[lat_name]<= lat2)
   data = data.where( mask.compute(), drop=True)
   return data
#---------------------------------------------------------------------------------------------------
for v in range(num_var):
   hc.printline()
   print(f'  var: {hc.tclr.MAGENTA}{var[v]}{hc.tclr.END}')
   tvar = var[v]
   area_name = 'area'
   #----------------------------------------------------------------------------
   # read the data
   data_list,lev_list = [],[]
   wav_power_list, wav_period_list = [],[]
   fft_power_list, fft_period_list = [],[]
   variance_list, alpha_list, sigma2_list = [],[],[]
   for c in range(num_case):
      yr1,yr2 = yr1_list[c],yr2_list[c]
      tmp_file = get_tmp_file(case[c],var[v],yr1,yr2)
      print(f'    case: {hc.tclr.GREEN}{case[c]}{hc.tclr.END}  =>  {tmp_file}')
      if recalculate_timeseries:
         if case[c]=='ERA5':
            lat_name,lon_name = 'lat','lon'
            xy_dims = (lon_name,lat_name)
            obs_root = '/global/cfs/cdirs/e3sm/diagnostics/observations/Atm'
            input_file_name = f'{obs_root}/time-series/ERA5/ua_197901_201912.nc'
            ds = xr.open_dataset( input_file_name )

            # ds = ds.isel(time=slice( (12*(yr1-1979)),(12*(yr2+1-1979)), ))

            # print()
            # print(ds)
            # print()
            # print(ds['time'][0])
            # print(ds['time'][-1])
            # print()
            # exit()

            # ds = ds.where( ds['time.year']>=yr1, drop=True)
            # ds = ds.where( ds['time.year']<=yr2, drop=True)

            # print()
            # print(ds)
            # print()

            # if num_files<=40*12: 
            #    ds = ds.isel(time=slice(0,num_files))
            # else:
            #    exit('ERROR - need to handle ERA5 data for longer time series')
            area = QBO_methods.calculate_area(ds['lon'].values,ds['lat'].values,ds['lon_bnds'].values,ds['lat_bnds'].values)
            area = xr.DataArray( area, coords=[ds['lat'],ds['lon']] )  
            data = ds['ua']
            data = data.rename({'plev':'lev'})
            data['lev'] = data['lev']/1e2
            data = data.sel(lev=pow_spec_lev)
            # tmp_data = np.ones([len(ds[lat_name]),len(ds[lon_name])],dtype=bool)
            # tmp_coords = {lat_name:ds[lat_name],lon_name:ds[lon_name]}
            # mask = xr.DataArray( tmp_data, coords=tmp_coords, dims=(lat_name,lon_name) )
            # mask = mask & (ds[lat_name]>= lat1) & (ds[lat_name]<= lat2)
            # data = data.where( mask, drop=True)
            # area = area.where( mask, drop=True)
            data = mask_data(ds,data)
            area = mask_data(ds,area)
            data_avg = ( (data*area).sum(dim=xy_dims) / area.sum(dim=xy_dims) )
         else:
            data_dir_tmp,data_sub_tmp = None, None
            if case_dir[c] is not None: data_dir_tmp = case_dir[c]
            if case_sub[c] is not None: data_sub_tmp = case_sub[c]
            case_obj = he.Case( name=case[c], atm_comp='eam', 
                                data_dir=data_dir_tmp, data_sub=data_sub_tmp, time_freq=None )
            if 'lat1' in vars() : case_obj.lat1,case_obj.lat2 = lat1,lat2
            if 'lon1' in vars() : case_obj.lon1,case_obj.lon2 = lon1,lon2

            # data = case_obj.load_data(tvar,    component='eam',htype=htype,first_file=first_file,num_files=num_files,lev=pow_spec_lev).isel(lev=0)
            # area = case_obj.load_data(area_name,component='eam',htype=htype,first_file=first_file,num_files=num_files).astype(np.double)

            file_path = f'{case_dir[c]}/{case[c]}/{case_sub[c]}/*.eam.{htype}.*'
            file_list = sorted(glob.glob(file_path))

            # subset files that fall within [yr1:yr2]
            file_list_all = file_list ; file_list = []
            for f in range(len(file_list_all)):
               yr = int(file_list_all[f][-10:-10+4])
               if yr>=yr1 and yr<=yr2: file_list.append(file_list_all[f])

            ds = xr.open_mfdataset( file_list )

            ds = ds.where( ds['time.year']>=yr1, drop=True)
            ds = ds.where( ds['time.year']<=yr2, drop=True)

            # data = ds[tvar]
            area = ds[area_name]

            data = he.interpolate_to_pressure(ds,data_mlev=ds[tvar],lev=np.array([pow_spec_lev]),
                                              ds_ps=ds,ps_var='PS',interp_type=2,extrap_flag=True).isel(lev=0)

            data = mask_data(ds,data)
            area = mask_data(ds,area)

            data_avg = ( (data*area).sum(dim='ncol') / area.sum(dim='ncol') )

         
         data_avg.load()

         ### print stats after time averaging
         # if print_stats: hc.print_stat(data_avg,name=var[v],stat='naxsh',indent='    ',compact=True)

         ds_out = xr.Dataset( coords=data_avg.coords )
         ds_out[var[v]] = data_avg
         ds_out.to_netcdf(path=tmp_file,mode='w')
      #----------------------------------------------------------------------
      else:
      
         tmp_ds = xr.open_dataset( tmp_file, use_cftime=True  )
         data_avg = tmp_ds[var[v]]

         # print(hc.tcolor.RED+'WARNING: artificially halving the data'+hc.tcolor.ENDC)
         # data_avg = data_avg[10*12:]

         # convert to anomalies
         data_avg = data_avg - data_avg.mean(dim='time')

         # detrend in time
         fit = xr.polyval(data_avg['time'], data_avg.polyfit(dim='time', deg=1).polyfit_coefficients)
         data_avg = data_avg - fit

         # hc.print_time_length(data_avg.time,print_span=True, print_length=False,indent=' '*6)
      #-------------------------------------------------------------------------

      # print()
      # print(data_avg)
      # print()

      hc.print_stat(data_avg,name=var[v],stat='naxsh',indent='    ',compact=True)

      # ( period, wavelet_spec ) = get_psd_from_wavelet_scipy(data_avg.values)
      ( period, wavelet_spec ) = get_psd_from_wavelet_pywt(data_avg.values)

      wav_power_list.append( np.sqrt(wavelet_spec) )
      wav_period_list.append( period )

      #-------------------------------------------------------------------------
      # ( period, wavelet_spec ) = get_psd_from_wavelet_pywt(data_avg.values)

      # fft_power_list.append( np.sqrt(wavelet_spec) )
      # fft_period_list.append( period )

      # #-----------------------------------------------------------------------------
      # # wavelet confidence interval
      # variance_list.append(np.var(data_avg.values))
      # (alpha,sigma2) = get_alpha_sigma2(data_avg.values)
      # alpha_list.append(alpha)
      # sigma2_list.append(sigma2)

      # #-----------------------------------------------------------------------------
      # # Calculate FFT for comparison
      # fft_period = np.concatenate( (np.arange(2.0, 33.0), np.arange(34.0, 100.0, 2.0)), axis=0 )
      # fft_psd_spec, fft_pow_spec = QBO_methods.get_psd_from_deseason(data_avg.values, fft_period)
      
      # fft_power_list.append( fft_pow_spec )
      # fft_period_list.append( fft_period)
   
   #----------------------------------------------------------------------------
   # Create plot
   tres = copy.deepcopy(res)
   ip = v
   tres.trXMinF = 8#np.min(tm_log)
   tres.trXMaxF = 64#np.max(tm_log)
   tres.trYMinF = 0
   data_max = np.max([np.max(d) for d in wav_power_list])
   tres.trYMaxF = data_max + 0.02*data_max
   for c in range(num_case):
      tres.xyLineColor,tres.xyDashPattern = clr[c],dsh[c]
      tplot = ngl.xy(wks, wav_period_list[c], wav_power_list[c], tres)
      if c==0: plot[ip] = tplot
      if c!=0: ngl.overlay(plot[ip],tplot)
   hs.set_subtitles(wks, plot[ip], '', '', f'{int(pow_spec_lev)} mb Zonal Wind Wavelet Power Spectrum', font_height=0.015)

   # #----------------------------------------------------------------------------
   # # Create plot
   # tres = copy.deepcopy(res)
   # for i in range(2):
   #    ip = v*2+i
   #    if i==0: tmp_period_list,tmp_power_list = wav_period_list, wav_power_list
   #    if i==1: tmp_period_list,tmp_power_list = fft_period_list, fft_power_list
   #    tres.trXMinF = np.min(tm_log); #np.min([np.min(d) for d in tmp_period_list])
   #    tres.trXMaxF = np.max(tm_log); #np.max([np.max(d) for d in tmp_period_list])
   #    tres.trYMinF = 0
   #    tres.trYMaxF = np.max([np.max(d) for d in tmp_power_list])
   #    for c in range(num_case):
   #       tres.xyLineColor,tres.xyDashPattern = clr[c],dsh[c]
   #       tplot = ngl.xy(wks, tmp_period_list[c], tmp_power_list[c], tres)
   #       if c==0: plot[ip] = tplot
   #       if c!=0: ngl.overlay(plot[ip],tplot)
   #    if i==0: hs.set_subtitles(wks, plot[ip], '', '', 'scipy Wavelet Power Spectrum', font_height=0.015)
   #    # if i==1: hs.set_subtitles(wks, plot[ip], '', '', 'FFT Power Spectrum', font_height=0.015)
   #    if i==1: hs.set_subtitles(wks, plot[ip], '', '', 'pywt Wavelet Power Spectrum', font_height=0.015)
   # #----------------------------------------------------------------------------
   # # Create plot
   # tres = copy.deepcopy(res)
   # for i in range(2):
   #    # ip = v*2+i
   #    ip = 0
   #    if i==0: tmp_period_list,tmp_power_list = wav_period_list, wav_power_list
   #    if i==1: tmp_period_list,tmp_power_list = fft_period_list, fft_power_list
   #    tres.trXMinF = np.min(tm_log); #np.min([np.min(d) for d in tmp_period_list])
   #    tres.trXMaxF = np.max(tm_log); #np.max([np.max(d) for d in tmp_period_list])
   #    tres.trYMinF = 0
   #    tres.trYMaxF = np.max([np.max(d) for d in tmp_power_list])
   #    for c in range(num_case):
   #       # tres.xyLineColor,tres.xyDashPattern = clr[c],dsh[c]
   #       if i==0 : tres.xyLineColor,tres.xyDashPattern = 'blue',0
   #       if i==1 : tres.xyLineColor,tres.xyDashPattern = 'red',1
   #       tplot = ngl.xy(wks, tmp_period_list[c], tmp_power_list[c], tres)
   #       if c==0 and i==0: 
   #          plot[ip] = tplot
   #       else:
   #          ngl.overlay(plot[ip],tplot)
   #    # if i==0: hs.set_subtitles(wks, plot[ip], '', '', 'scipy Wavelet Power Spectrum', font_height=0.015)
   #    # if i==1: hs.set_subtitles(wks, plot[ip], '', '', 'FFT Power Spectrum', font_height=0.015)
   #    if i==1: hs.set_subtitles(wks, plot[ip], '', '', 'Wavelet Power Spectrum', font_height=0.015)
   #-------------------------------------------------------------------------
   # # indicate 95% confidence level
   # if True:
   #    period = period_list[0]
   #    red_freq = 1/(period*12)
   #    alpha = 0.90

   #    lres = hs.res_xy()
   #    lres.xyLineColor = 'gray'
   #    lres.xyDashPattern = 1

   #    # plot 95% confidence level
   #    for c in range(num_case):
   #       lres.xyDashPattern = 2
   #       lres.xyLineColor = clr[c]
   #       # q95 = P*5.99*0.5*variance_list[c]
   #       # Compute the local wavelet power spectrum upper bound eq. 18, 5.99 = 95-th quantile for a Chi² r.v. with 2 DOF
   #       alpha = alpha_list[c]
   #       alpha_sq = np.square(alpha)
   #       P = ( 1 - alpha_sq ) / ( 1 + alpha_sq - 2*alpha*np.cos(2*np.pi*red_freq) )
   #       q95 = 0.5 * P * 5.99 * sigma2_list[c]
   #       ngl.overlay(plot[ip], ngl.xy(wks, period, q95, lres) )
   #-------------------------------------------------------------------------
   # add vertical lines for dominant periods
   lres = hs.res_xy()
   lres.xyLineColor = 'gray'
   lres.xyDashPattern = 2
   
   yy = np.array([-1e3,1e3])
   
   for f in [32,28,24]:
      xx = np.array([1,1]) * f
      ngl.overlay(plot[ip], ngl.xy(wks, xx, yy, lres) )

#----------------------------------------------------------------------------
# Add legend
#----------------------------------------------------------------------------
lgres = ngl.Resources()
lgres.vpWidthF, lgres.vpHeightF  = 0.08, 0.15
lgres.lgLabelFontHeightF = 0.01
lgres.lgLineThicknessF   = 12
lgres.lgMonoLineColor    = False
# lgres.lgMonoDashIndex    = True
# lgres.lgDashIndex        = 0
lgres.lgDashIndexes      = dsh
lgres.lgLineColors       = clr
lgres.lgLabelJust    = 'CenterLeft'
# lbl = [f'  {n}' for n in case_name]
# pid = ngl.legend_ndc(wks, len(case_name), lbl, 0.2, 0.7, lgres)

def add_legend(ind,ypos):
   global clr,dsh,case_name
   lgres.lgDashIndexes      = dsh[ind]
   lgres.lgLineColors       = clr[ind]
   pid = ngl.legend_ndc(wks, 1, [f'  {case_name[ind]}'], 0.2, ypos, lgres)

y_top = 0.75
dy = 0.022
for c in range(num_case):
   add_legend(c, y_top-dy*c-dy*0)
# add_legend(0, y_top-dy*0-dy*0)
# add_legend(1, y_top-dy*1-dy*1)
# add_legend(2, y_top-dy*2-dy*1)
# add_legend(3, y_top-dy*3-dy*1)
# add_legend(4, y_top-dy*4-dy*2)
# add_legend(5, y_top-dy*5-dy*2)
# add_legend(6, y_top-dy*6-dy*2)

#---------------------------------------------------------------------------------------------------
# Finalize plot
#---------------------------------------------------------------------------------------------------
layout = [int(np.ceil(len(plot)/float(num_plot_col))),num_plot_col]

pnl_res = hs.setres_panel()
pnl_res.nglPanelYWhiteSpacePercent = 5

ngl.panel(wks,plot,layout,pnl_res)
ngl.end()

hc.trim_png(fig_file)
#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------

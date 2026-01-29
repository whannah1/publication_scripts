import os, ngl, subprocess as sp, numpy as np, xarray as xr, copy, string, dask
import hapy_common as hc, hapy_E3SM   as he, hapy_setres as hs
import cmocean
import scipy
import warnings
# import statsmodels
from statsmodels.tsa.arima.model import ARIMA
#-------------------------------------------------------------------------------
# based on E3SM diagnostics package:
# https://github.com/E3SM-Project/e3sm_diags/blob/main/e3sm_diags/driver/qbo_driver.py
#-------------------------------------------------------------------------------
case_name,case,case_dir,case_sub,case_grid,clr,dsh,mrk = [],[],[],[],[],[],[],[]
def add_case(case_in,n=None,p=None,s=None,g=None,d=0,c='black',m=0):
   global name,case,case_dir,case_sub,clr,dsh,mrk
   if n is None:
      tmp_name = ''
   else:
      tmp_name = n
   case.append(case_in); case_name.append(tmp_name)
   case_dir.append(p); case_sub.append(s); case_grid.append(g)
   dsh.append(d) ; clr.append(c) ; mrk.append(m)
#-------------------------------------------------------------------------------
scrip_file_path  = os.getenv('HOME')+'/E3SM/data_grid/ne30pg2_scrip.nc'
tmp_path_hst_mmf = '/global/cfs/cdirs/m3312/whannah/2023-CPL'
tmp_path_hst_v2  = '/global/cfs/cdirs/m3312/whannah/e3smv2_historical'
tmp_path_co2_mmf = '/global/cfs/cdirs/m3312/whannah/2023-CPL'
tmp_sub = 'archive/atm/hist'
# add_case('v2.LR.historical_0101',                                  n='E3SMv2',  c='cyan',    p=tmp_path_hst_v2, s='archive/atm/hist')
# add_case('v2.LR.historical_0151',                                  n='E3SMv2',  c='orange', p=tmp_path_hst_v2, s='archive/atm/hist')
# add_case('v2.LR.historical_0201',                                  n='E3SMv2',  c='green',  p=tmp_path_hst_v2, s='archive/atm/hist')
# add_case('v2.LR.historical_0251',                                  n='E3SMv2',  c='purple',   p=tmp_path_hst_v2, s='archive/atm/hist')
# add_case('v2.LR.historical_0301',                                  n='E3SMv2',  c='pink', p=tmp_path_hst_v2, s='archive/atm/hist')
add_case('v2.LR.historical',                                       n='E3SMv2 Ens Mean',  c='red',  p=None, s=None)
add_case('E3SM.INCITE2023-CPL.ne30pg2_EC30to60E2r2.WCYCL20TR-MMF1',n='E3SM-MMF',c='blue', p=tmp_path_hst_mmf,s=tmp_sub)
add_case('HadSST',                                                 n='HadSST',  c='black')
#-------------------------------------------------------------------------------

var = ['TS']

fig_file,fig_type = 'figs/FXX-ENSO-wavelet-spectra','png'
tmp_file_head = 'data/ENSO-power-spectra' # use data from FFT script

# region,lat1,lat2,lon1,lon2 = 'Nino3'  ,-5,5,360-150,360-90
region,lat1,lat2,lon1,lon2 = 'Nino3.4',-5,5,190,240
# region,lat1,lat2,lon1,lon2 = 'Nino4'  ,-5,5,160,360-150


# htype,num_files = 'h0',12*30 ; first_file,first_file_v2 = 12*30,12*80
htype,num_files = 'h0',12*65 ; first_file,first_file_v2 = 12*00,12*50

print_stats = True

var_x_case = False

use_common_label_bar = True

num_plot_col = 1

recalculate_timeseries = False

# year_start = 1950+first_file/12
# year_end   = 1950+first_file/12+num_files/12+1

year_start = 1950
year_end   = 2014

# period = np.concatenate( (np.arange(2.0, 33.0), np.arange(34.0, 120.0, 2.0)), axis=0 )
# period = np.concatenate( (np.arange(6,33), np.arange(34,120,2)), axis=0 )

# period = np.logspace( np.log10(6), np.log10(120), num=40).round(decimals=4)
# period = np.logspace( np.log10(6), np.log10(120), num=20).round(decimals=8)

# print(period/12) ; exit()

#---------------------------------------------------------------------------------------------------
# Set up plot resources
#---------------------------------------------------------------------------------------------------
num_var,num_case = len(var),len(case)

if 'lev' not in vars(): lev = np.array([0])

wkres = ngl.Resources()
npix = 2048; wkres.wkWidth,wkres.wkHeight=npix,npix
wks = ngl.open_wks(fig_type,fig_file,wkres)
plot = [None]*num_var

res = hs.res_xy()
res.vpHeightF = 0.4
# res.xyMarkLineMode = "MarkLines"
res.xyMarkerSizeF = 0.008
res.xyMarker = 16
res.xyLineThicknessF = 12
# res.tmYLLabelFontHeightF         = 0.008
# res.tmXBLabelFontHeightF         = 0.008
# res.tmXBAutoPrecision = False
# res.tmXBPrecision = 2

res.tiXAxisString = 'Period [years]'
res.tiYAxisString = 'Variance [~F34~0~F~C~S~2~N~]'

# res.trXMinF = 0.5
# res.trXMaxF = 10


res.tmXBMode      = 'Explicit'

# tm_log = np.array([1,2,4,8,16,32,64])
# res.tmXBValues    = tm_log
# res.tmXBLabels    = tm_log
# res.xyXStyle ='Log'

tm_lin = np.array([2,4,6,8,10])
res.tmXBValues    = tm_lin
res.tmXBLabels    = tm_lin



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
def get_psd_from_wavelet(data):
   deg = 6
   # period = np.arange(6,12*12+1)
   period = np.arange(6,12*10+1)
   # period = np.arange(6,12*64+1)
   freq = 1/period
   widths = deg / (2*np.pi*freq)
   with warnings.catch_warnings():
      warnings.simplefilter("ignore", category=DeprecationWarning)
      cwtmatr = scipy.signal.cwt( deseason(data), scipy.signal.morlet2, widths=widths, w=deg )
   # dof = len(data) - widths
   # print(dof)
   # exit()
   psd = np.mean(np.square(np.abs(cwtmatr)),axis=1)
   return ( period, psd )

#---------------------------------------------------------------------------------------------------
# def get_alpha_sigma2(data_in):
#    ## Fit AR1 model to estimate the autocorrelation (alpha) and variance (sigma2)
#    mod = ARIMA( data_in, order=(1,0,0) )
#    res = mod.fit()
#    return (res.params[1],res.params[2])
#---------------------------------------------------------------------------------------------------
def get_tmp_file(case,var):
   return f'{tmp_file_head}.{case}.{var}.reg_{region}.lat1_{lat1}.lat2_{lat2}.lon1_{lon1}.lon2_{lon2}.nc'
#---------------------------------------------------------------------------------------------------
for v in range(num_var):
   hc.printline()
   print(f'  var: {hc.tclr.MAGENTA}{var[v]}{hc.tclr.END}')
   tvar = var[v]
   area_name = 'area'
   #----------------------------------------------------------------------------
   # read the data
   #----------------------------------------------------------------------------
   data_list,lev_list = [],[]
   period_list, power_list = [],[]
   variance_list = []
   alpha_list = []
   sigma2_list = []
   for c in range(num_case):

      obs_flag = False
      if case[c] in ['HadSST']: obs_flag = True
      if case[c] in ['BEST']  : obs_flag = True

      # tmp_file = f'{tmp_file_head}.{case[c]}.{var[v]}.reg_{region}.lat1_{lat1}.lat2_{lat2}.lon1_{lon1}.lon2_{lon2}.nc'
      tmp_file = get_tmp_file(case[c],var[v])

      print(f'    case: {hc.tclr.GREEN}{case[c]}{hc.tclr.END}  =>  {tmp_file}')

      # avoid creating large chunks
      with dask.config.set(**{'array.slicing.split_large_chunks': True}):  
         if recalculate_timeseries:
            if obs_flag:
               if case[c]=='HadSST':
                  sst_name, lat_name, lon_name = 'sst','latitude','longitude'
                  file_name = os.getenv('HOME')+'/Data/Obs/HadSST/HadISST_sst.nc'
               if case[c]=='BEST':
                  sst_name, lat_name, lon_name = 'sst','latitude','longitude'
                  file_name = '/global/cfs/cdirs/m3312/whannah/obs_data/BEST/Land_and_Ocean_LatLong1.remap_ne30pg2.nc'
               ds = xr.open_dataset( file_name )
               data = ds[sst_name]
               data = data.where( data['time.year']>=year_start, drop=True)
               data = data.where( data['time.year']<=year_end,   drop=True)
               xy_dims = (lon_name,lat_name)
               xlon, ylat = np.meshgrid(ds[lat_name],ds[lon_name])
               R = hc.earth_radius(ylat)
               dlat = np.deg2rad(np.gradient(ylat, axis=0))
               dlon = np.deg2rad(np.gradient(xlon, axis=1))
               dy,dx = dlat * R , dlon * R * np.cos(np.deg2rad(ylat))
               area = np.absolute(dy*dx) / np.square(R) # calculate area and convert to steridians
               area = xr.DataArray(area,dims=xy_dims).transpose()
               tlon1,tlon2 = lon1,lon2
               if tlon1>180: tlon1 = tlon1 - 360
               if tlon2>180: tlon2 = tlon2 - 360
               tmp_data = np.ones([len(ds[lat_name]),len(ds[lon_name])],dtype=bool)
               tmp_coords = {lat_name:ds[lat_name],lon_name:ds[lon_name]}
               mask = xr.DataArray( tmp_data, coords=tmp_coords, dims=(lat_name,lon_name) )
               mask = mask & (ds[lat_name]>= lat1) & (ds[lat_name]<= lat2)
               mask = mask & (ds[lon_name]>=tlon1) & (ds[lon_name]<=tlon2)
               data = data.where( mask, drop=True)
               area = area.where( mask, drop=True)
               data_avg = ( (data*area).sum(dim=xy_dims) / area.sum(dim=xy_dims) )
            else:
               data_dir_tmp,data_sub_tmp = None, None
               if case_dir[c] is not None: data_dir_tmp = case_dir[c]
               if case_sub[c] is not None: data_sub_tmp = case_sub[c]
               case_obj = he.Case( name=case[c], atm_comp='eam', 
                                   data_dir=data_dir_tmp, data_sub=data_sub_tmp, time_freq=None )
               if 'lat1' in vars() : case_obj.lat1,case_obj.lat2 = lat1,lat2
               if 'lon1' in vars() : case_obj.lon1,case_obj.lon2 = lon1,lon2

               first_file_tmp = first_file_v2 if 'v2.LR.historical' in case[c] else first_file
               data = case_obj.load_data(tvar,    
                                         component='eam',
                                         htype=htype,
                                         first_file=first_file_tmp,
                                         num_files=num_files,
                                         lev=lev)
               area = case_obj.load_data(area_name,
                                         component='eam',
                                         htype=htype,
                                         first_file=first_file_tmp,
                                         num_files=num_files).astype(np.double)
               data_avg = ( (data*area).sum(dim='ncol') / area.sum(dim='ncol') )
            
            data_avg.load()

            ### print stats after time averaging
            # if print_stats: hc.print_stat(data_avg,name=var[v],stat='naxsh',indent='    ',compact=True)

            ds_out = xr.Dataset( coords=data_avg.coords )
            ds_out[var[v]] = data_avg
            ds_out.to_netcdf(path=tmp_file,mode='w')
         #----------------------------------------------------------------------
         else:
            if case[c]=='v2.LR.historical':
               v2_amip_ens_list = []
               v2_amip_ens_list.append('v2.LR.historical_0101')
               v2_amip_ens_list.append('v2.LR.historical_0151')
               v2_amip_ens_list.append('v2.LR.historical_0201')
               v2_amip_ens_list.append('v2.LR.historical_0251')
               v2_amip_ens_list.append('v2.LR.historical_0301')
               ens_cnt = 0
               psd_list = []
               variance_tmp_list = []
               # alpha_tmp_list = []
               # sigma2_tmp_list = []
               for e,ens_member in enumerate(v2_amip_ens_list):
                  tmp_file = get_tmp_file(ens_member,var[v])
                  tmp_ds = xr.open_dataset( tmp_file, use_cftime=True  )
                  ens_member_data = tmp_ds[var[v]]
                  hc.print_time_length(ens_member_data.time,print_span=True, print_length=False,indent=' '*6)
                  
                  ( period, psd_tmp ) = get_psd_from_wavelet(ens_member_data.values)

                  psd_list.append(psd_tmp)

                  variance_tmp_list.append(np.var(ens_member_data.values))

                  # (alpha,sigma2) = get_alpha_sigma2(ens_member_data.values)
                  # alpha_tmp_list.append(alpha)
                  # sigma2_tmp_list.append(sigma2)

                  # if ens_cnt==0: 
                  #    psd = np.zeros_like(psd_tmp)
                  #    ens_psd_min,ens_psd_max = psd_tmp.copy(),psd_tmp.copy()
                  # else:
                  #    for t in range(len(psd)):
                  #       ens_psd_min[t] = np.min([ ens_psd_min[t], psd_tmp[t] ])
                  #       ens_psd_max[t] = np.max([ ens_psd_max[t], psd_tmp[t] ])
                  # psd = ( psd*ens_cnt + psd_tmp ) / (ens_cnt+1)
                  # ens_cnt += 1
               psd_list = np.array(psd_list)
               psd = np.mean(psd_list,axis=0)
               ens_psd_min = np.empty(len(psd))
               ens_psd_max = np.empty(len(psd))
               for i in range(len(psd)):
                  sd = np.std(psd_list[:,i])
                  ens_psd_min[i] = psd[i] - sd
                  ens_psd_max[i] = psd[i] + sd

               variance_list.append(np.mean(variance_tmp_list,axis=0))

               # alpha_list.append(np.mean(alpha_tmp_list,axis=0))
               # sigma2_list.append(np.mean(sigma2_tmp_list,axis=0))
               
               # print(psd)
               # exit()
            else:
               tmp_ds = xr.open_dataset( tmp_file, use_cftime=True  )
               data_avg = tmp_ds[var[v]]

               # convert to anomalies
               data_avg = data_avg - data_avg.mean()

               # detrend in time
               fit = xr.polyval(data_avg['time'], data_avg.polyfit(dim='time', deg=1).polyfit_coefficients)
               data_avg = data_avg - fit

               hc.print_time_length(data_avg.time,print_span=True, print_length=False,indent=' '*6)
      

      # if 'lev' in data_avg.dims: data_avg = data_avg.isel(lev=0)
      #----------------------------------------------------------------------
      if case[c]!='v2.LR.historical':
         ( period, psd ) = get_psd_from_wavelet(data_avg.values)

         variance_list.append(np.var(data_avg.values))

         # (alpha,sigma2) = get_alpha_sigma2(data_avg.values)
         # alpha_list.append(alpha)
         # sigma2_list.append(sigma2)

      #----------------------------------------------------------------------
      period_list.append( period/12. )
      power_list .append( psd )
   #----------------------------------------------------------------------------
   # Create plot
   #----------------------------------------------------------------------------
   tres = copy.deepcopy(res)
   tres.trXMinF = np.min([np.min(d) for d in period_list])
   tres.trXMaxF = np.max([np.max(d) for d in period_list])
   tres.trYMinF = 0
   tres.trYMaxF = np.max([np.max(d) for d in power_list])

   tres.trYMaxF = 16
   tres.trXReverse = False

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
         # fclr = 'orange'
         fclr = 'pink'
         # fclr = 'tan'
         # fclr = ( 1, 0, 0 )
         eres.nglXYAboveFillColors = fclr
         eres.nglXYBelowFillColors = fclr
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

   # hs.set_subtitles(wks, plot[len(plot)-1], case_name[c], '', '', font_height=0.015)
   hs.set_subtitles(wks, plot[len(plot)-1], '', region, '', font_height=0.015)

   #-------------------------------------------------------------------------
   # # Function to compute the red noise spectrum
   # def get_P(alpha, period, sigma2):
   #    freq = 1/(period_list[c]*12)
   #    alpha_sq = np.square(alpha)
   #    P = ( 1 - alpha_sq ) / ( 1 + alpha_sq - 2*alpha*np.cos(2*np.pi*freq) )
   #    return  P
   #-------------------------------------------------------------------------
   for c in range(num_case):

      alpha_assumed = 0.72
      red_freq = 1/(period_list[c]*12)
      alpha_sq = np.square(alpha_assumed)
      P = ( 1 - alpha_sq ) / ( 1 + alpha_sq - 2*alpha_assumed*np.cos(2*np.pi*red_freq) )
      chi_val = 5.99 
      Ws = 0.5 * P * 5.99 * variance_list[c]
      q95 = Ws * 2 / chi_val

      lres = hs.res_xy()
      lres.xyDashPattern = 2
      lres.xyLineColor = clr[c]

      ngl.overlay(plot[ip], ngl.xy(wks, period_list[c], q95, lres) )

   # lres.xyLineColor = 'red'
   # ngl.overlay(plot[ip], ngl.xy(wks, period_list[c], P, lres) )
   
   # lres.xyLineColor = 'blue'
   # ngl.overlay(plot[ip], ngl.xy(wks, period_list[c], q95, lres) )

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
# lgres = ngl.Resources()
# lgres.vpWidthF, lgres.vpHeightF  = 0.08, 0.1  
# lgres.lgLabelFontHeightF = 0.015
# lgres.lgLineThicknessF   = 12
# lgres.lgMonoLineColor    = False
# lgres.lgMonoDashIndex    = True
# lgres.lgDashIndex        = 0
# lgres.lgLineColors       = clr
# lgres.lgLabelJust    = 'CenterLeft'
# lbl = [f'  {n}' for n in case_name]
# # pid = ngl.legend_ndc(wks, len(case_name), lbl, 0.18, 0.68, lgres)

lgres = ngl.Resources()
lgres.vpWidthF           = 0.06
lgres.vpHeightF          = 0.12#*num_case
lgres.lgLabelFontHeightF = 0.014
lgres.lgLabelFont        = "courier"
lgres.lgMonoDashIndex    = False
lgres.lgLineLabelsOn     = False
lgres.lgLineThicknessF   = 40
lgres.lgLabelJust        = 'CenterLeft'
lgres.lgLineColors       = clr
lgres.lgDashIndexes      = dsh

indent = ' '*4
labels = case_name
for i in range(len(labels)): labels[i] = indent+labels[i] 

# if add_obs_TS: labels.insert(obs_pos,indent+'HadCRU')
# if add_obs_TS: labels.insert(obs_pos,indent+'BEST')
if 'v2.LR.historical' in case:
   labels.insert(0,indent+'E3SMv2 Ens Min/Max')
   lgres.lgLineColors.insert(0,'pink')
   lgres.lgDashIndexes.insert(0,0)
   # xyLineOpacities

pid = ngl.legend_ndc(wks, len(labels), labels, 0.58, 0.78, lgres)

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

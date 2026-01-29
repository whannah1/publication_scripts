import os, ngl, subprocess as sp, numpy as np, xarray as xr, copy, string, dask
import hapy_common as hc, hapy_E3SM   as he, hapy_setres as hs
import cmocean
import scipy
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
tmp_path_hst_mmf = '/global/cfs/cdirs/m3312/whannah/2023-CPL'
tmp_path_hst_v2  = '/global/cfs/cdirs/m3312/whannah/e3smv2_historical'
scrip_file_path  = os.getenv('HOME')+'/E3SM/data_grid/ne30pg2_scrip.nc'
add_case('v2.LR.historical_0101',                                  n='E3SMv2',  c='cyan',    p=tmp_path_hst_v2, s='archive/atm/hist')
add_case('v2.LR.historical_0151',                                  n='E3SMv2',  c='orange', p=tmp_path_hst_v2, s='archive/atm/hist')
add_case('v2.LR.historical_0201',                                  n='E3SMv2',  c='green',  p=tmp_path_hst_v2, s='archive/atm/hist')
add_case('v2.LR.historical_0251',                                  n='E3SMv2',  c='purple',   p=tmp_path_hst_v2, s='archive/atm/hist')
add_case('v2.LR.historical_0301',                                  n='E3SMv2',  c='pink', p=tmp_path_hst_v2, s='archive/atm/hist')
# add_case('v2.LR.historical',                                       n='E3SMv2',  c='red',  p=None, s=None)
add_case('HadSST',                                                 n='HadSST',  c='black')
add_case('E3SM.INCITE2023-CPL.ne30pg2_EC30to60E2r2.WCYCL20TR-MMF1',n='E3SM-MMF',c='blue', p=tmp_path_hst_mmf,s='archive/atm/hist')
#-------------------------------------------------------------------------------

var = ['TS']

fig_file,fig_type = 'figs/FXX-ENSO-power-spectra','png'
tmp_file_head = 'data/ENSO-power-spectra'

# region,lat1,lat2,lon1,lon2 = 'Nino3'  ,-5,5,360-150,360-90
# region,lat1,lat2,lon1,lon2 = 'Nino4'  ,-5,5,160,360-150
region,lat1,lat2,lon1,lon2 = 'Nino3.4',-5,5,190,240

# htype,num_files = 'h0',12*30 ; first_file,first_file_v2 = 12*30,12*80
htype,num_files = 'h0',12*65 ; first_file,first_file_v2 = 12*00,12*50

print_stats = True

var_x_case = False

use_common_label_bar = True

num_plot_col = 1

recalculate_timeseries = True

year_start = 1950+first_file/12
year_end   = 1950+first_file/12+num_files/12+1

# period = np.concatenate( (np.arange(2.0, 33.0), np.arange(34.0, 120.0, 2.0)), axis=0 )
# period = np.concatenate( (np.arange(6,33), np.arange(34,120,2)), axis=0 )

period = np.logspace( np.log10(6), np.log10(120), num=40).round(decimals=4)
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
res.vpHeightF = 0.3
# res.xyMarkLineMode = "MarkLines"
res.xyMarkerSizeF = 0.008
res.xyMarker = 16
res.xyLineThicknessF = 8
# res.tmYLLabelFontHeightF         = 0.008
# res.tmXBLabelFontHeightF         = 0.008
# res.tmXBAutoPrecision = False
# res.tmXBPrecision = 2

res.tiXAxisString = 'Period [years]'
res.tiYAxisString = 'Amplitude'

# res.trXMinF = 0.5
res.trXMaxF = 10

tm_log = np.array([1,2,4,8,16])
res.tmXBMode      = 'Explicit'
res.tmXBValues    = tm_log
res.tmXBLabels    = tm_log

res.xyXStyle ='Log'

#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------
def get_comp(case):
   comp = 'eam'
   return comp
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
def ceil_log2(x):
    """
    Given a number, calculate the exponent for the next power of 2.
    Example:
        ceil_log2(16) = 4
        ceil_log2(17) = 5
    """
    return np.ceil(np.log2(x)).astype("int")
#---------------------------------------------------------------------------------------------------
def get_psd_from_deseason(xraw, period_new):
    x_deseasoned = deseason(xraw)

    # Sampling frequency: assumes frequency of sampling = 1 month
    sampling_frequency = 1
    # Calculate the period as a function of frequency
    period0 = 1 / sampling_frequency
    L0 = len(xraw)
    NFFT0 = 2 ** ceil_log2(L0)

    # Apply fft on x_deseasoned with n = NFFT
    x0 = scipy.fftpack.fft(x_deseasoned, n=NFFT0) / L0
    # Frequency (cycles/month). Frequency will be increasing.
    frequency0 = sampling_frequency * np.arange(0, (NFFT0 / 2 + 1)) / NFFT0
    # Period (months/cycle). Calculate as a function of frequency. Thus, period will be decreasing.
    period0 = np.zeros_like(frequency0)
    for f,freq in enumerate(frequency0):
      period0[f] = 0 if freq==0 else 1/freq

    # Calculate amplitude as a function of frequency
    amplitude0 = 2 * abs(x0[0 : int(NFFT0 / 2 + 1)])
    # Calculate power spectral density as a function of frequency
    psd_x0 = amplitude0**2 / L0
    # Total spectral power
    # In the next code block, we will perform an interpolation using the period
    # (interpolating values of amplitude0_flipped and psd_x0_flipped from period0_flipped to period_new).
    # For that interpolation, we want the period to be increasing.
    # Therefore, we will flip the following values:
    period0_flipped = period0[::-1]  # type: ignore
    amplitude0_flipped = amplitude0[::-1]
    psd_x0_flipped = psd_x0[::-1]

    amplitude_new0 = np.interp(period_new, period0_flipped[:-1], amplitude0_flipped[:-1])
    psd_x_new0 = np.interp(period_new, period0_flipped[:-1], psd_x0_flipped[:-1])
    return psd_x_new0, amplitude_new0
#---------------------------------------------------------------------------------------------------
# def get_psd_from_deseason_alt(xraw):
#     x_deseasoned = deseason(xraw)

#     # Sampling frequency: assumes frequency of sampling = 1 month
#     sampling_frequency = 1
#     # Calculate the period as a function of frequency
#     period0 = 1 / sampling_frequency
#     L0 = len(xraw)
#     NFFT0 = 2 ** ceil_log2(L0)

#     # Apply fft on x_deseasoned with n = NFFT
#     x0 = scipy.fftpack.fft(x_deseasoned, n=NFFT0) / L0
#     # Frequency (cycles/month). Frequency will be increasing.
#     frequency0 = sampling_frequency * np.arange(0, (NFFT0 / 2 + 1)) / NFFT0
#     # Period (months/cycle). Calculate as a function of frequency. Thus, period will be decreasing.
    
#     # period0 = 1 / frequency0
#     period0 = frequency0
#     for f,freq in enumerate(frequency0):
#       period0[f] = 0 if freq==0 else 1/freq

#     # Calculate amplitude as a function of frequency
#     amplitude0 = 2 * abs(x0[0 : int(NFFT0 / 2 + 1)])
#     # Calculate power spectral density as a function of frequency
#     # psd_x0 = amplitude0**2 / L0

#     # we want the period to be increasing
#     period_out    = period0[::-1]
#     amplitude_out = amplitude0[::-1]

#     # discard last value
#     period_out    = period_out[:-1]
#     amplitude_out = amplitude_out[:-1]

#     return period_out, amplitude_out
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
   for c in range(num_case):

      obs_flag = False
      if case[c] in ['HadSST']: obs_flag = True

      # tmp_file = f'{tmp_file_head}.{case[c]}.{var[v]}.reg_{region}.lat1_{lat1}.lat2_{lat2}.lon1_{lon1}.lon2_{lon2}.nc'
      tmp_file = get_tmp_file(case[c],var[v])

      print(f'    case: {hc.tclr.GREEN}{case[c]}{hc.tclr.END}  =>  {tmp_file}')

      # avoid creating large chunks
      with dask.config.set(**{'array.slicing.split_large_chunks': True}):  
         if recalculate_timeseries:
            if obs_flag:
               if case[c]=='HadSST': sst_name,lat_name,lon_name,file_name = 'sst','latitude','longitude',os.getenv('HOME')+'/Data/Obs/HadSST/HadISST_sst.nc'
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
               case_obj = he.Case( name=case[c], atm_comp=get_comp(case[c]), 
                                   data_dir=data_dir_tmp, data_sub=data_sub_tmp, time_freq=None )
               if 'lat1' in vars() : case_obj.lat1,case_obj.lat2 = lat1,lat2
               if 'lon1' in vars() : case_obj.lon1,case_obj.lon2 = lon1,lon2

               first_file_tmp = first_file_v2 if 'v2.LR.historical' in case[c] else first_file
               data = case_obj.load_data(tvar,    
                                         component=get_comp(case[c]),
                                         htype=htype,
                                         first_file=first_file_tmp,
                                         num_files=num_files,
                                         lev=lev)
               area = case_obj.load_data(area_name,
                                         component=get_comp(case[c]),
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
               for e,ens_member in enumerate(v2_amip_ens_list):
                  tmp_file = get_tmp_file(ens_member,var[v])
                  tmp_ds = xr.open_dataset( tmp_file, use_cftime=True  )
                  ens_member_data = tmp_ds[var[v]]
                  hc.print_time_length(ens_member_data.time,print_span=True, print_length=False,indent=' '*6)
                  
                  psd_tmp, amp_tmp = get_psd_from_deseason(ens_member_data.values, period)
                  # period, amp_tmp = get_psd_from_deseason_alt(ens_member_data.values)

                  if ens_cnt==0: 
                     amp = np.zeros_like(amp_tmp)
                     ens_amp_min,ens_amp_max = amp_tmp.copy(),amp_tmp.copy()
                  else:
                     for t in range(len(amp)):
                        ens_amp_min[t] = np.min([ ens_amp_min[t], amp_tmp[t] ])
                        ens_amp_max[t] = np.max([ ens_amp_max[t], amp_tmp[t] ])
                  amp = ( amp*ens_cnt + amp_tmp ) / (ens_cnt+1)
                  ens_cnt += 1
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
      # period = np.concatenate( (np.arange(2.0, 33.0), np.arange(34.0, 120.0, 2.0)), axis=0 )
      # psd, amp = get_psd_from_deseason(data_avg.values, period)

      if case[c]!='v2.LR.historical':
         psd, amp = get_psd_from_deseason(data_avg.values, period)
         # period, amp = get_psd_from_deseason_alt(data_avg.values)
         
      #----------------------------------------------------------------------
      period_list.append( period )
      power_list .append( amp )
   #----------------------------------------------------------------------------
   # Create plot
   #----------------------------------------------------------------------------
   tres = copy.deepcopy(res)
   tres.trYMinF = 0
   tres.trYMaxF = np.max([np.nanmax(d) for d in power_list])
   if 'ens_amp_max' in locals(): tres.trYMaxF = np.max([tres.trYMaxF,np.max(ens_amp_max)])
   for c in range(num_case):

      
      tres.xyLineColor,tres.xyDashPattern = clr[c],dsh[c]

      tres.xyLineThicknessF = res.xyLineThicknessF
      if case[c] in ['HadSST']: tres.xyLineThicknessF = res.xyLineThicknessF * 2

      tplot = ngl.xy(wks, period_list[c]/12., power_list[c], tres)

      ip = v
      if c==0:
         plot[ip] = tplot 
      else:
         ngl.overlay(plot[ip],tplot)

      if case[c]=='v2.LR.historical':
         eres = copy.deepcopy(tres)
         ens_spread_data      = np.zeros([2,len(ens_amp_min)])
         ens_spread_data[0,:] = ens_amp_min
         ens_spread_data[1,:] = ens_amp_max
         eres.xyLineColor = [0,0,0,0]
         eres.nglXYAboveFillColors = 'orange'
         eres.nglXYBelowFillColors = 'orange'
         ngl.overlay(plot[ip], ngl.xy(wks, period_list[c]/12., ens_spread_data, eres) )

      #-------------------------------------------------------------------------
      # Set strings at top of plot
      #-------------------------------------------------------------------------
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

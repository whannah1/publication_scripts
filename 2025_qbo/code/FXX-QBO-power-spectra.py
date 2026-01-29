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
#-------------------------------------------------------------------------------
cscratch = '/global/cscratch1/sd/whannah/e3sm_scratch/cori-knl'
gscratch = '/global/cfs/cdirs/m3312/whannah/2022-QBO-TEST'
pscratch = '/pscratch/sd/w/whannah/e3sm_scratch/pm-cpu'

add_case('E3SM.QBO-TEST.F2010.ne30pg2.L72.01',        n='E3SM control',     d=0,c='black', p=gscratch,s='run')
add_case('E3SM.QBO-TEST.F2010.ne30pg2.L72-nsu40.01',  n='E3SM L72 smoothed',d=0,c='red',   p=gscratch,s='run')
# add_case('E3SM.QBO-TEST.F2010.ne30pg2.L80-rsu40.01',  n='E3SM L80 refined', d=0,c='blue' , p=gscratch,s='run')

add_case('E3SM.QBO-TEST-02.F2010.ne30pg2.L72.GW-MOD-0',n='E3SM L72',        d=1,c='black', p=pscratch,s='run')
add_case('E3SM.QBO-TEST-02.F2010.ne30pg2.L72-nsu40',   n='E3SM L72-nsu40',  d=1,c='red',   p=pscratch,s='run')
#add_case('E3SM.QBO-TEST-02.F2010.ne30pg2.L72-rscl',    n='E3SM L72-rscl',   d=1,c='purple',p=pscratch,s='run')
#add_case('E3SM.QBO-TEST-02.F2010.ne30pg2.L72-rlim',    n='E3SM L72-rlim',   d=1,c='pink',  p=pscratch,s='run')

### Runs for QBO paper
# add_case('E3SM.QBO-TEST.F2010.ne30pg2.L72.01',       n='E3SM control')
# add_case('E3SM.QBO-TEST.F2010.ne30pg2.L72-nsu40.01', n='E3SM L72 smoothed')
# add_case('E3SM.QBO-TEST.F2010.ne30pg2.L80-rsu40.01', n='E3SM L80 refined')

### v3 candidate w/ smoothed L72
#add_case('20230127.amip.v3atm.FourthSmoke.chrysalis',                   n='v3 '     ,         c='black',p='/lcrc/group/e3sm/ac.golaz/E3SMv3_dev'  ,s='archive/atm/hist')
#add_case('20230214.amip.v3atm.FourthSmoke.chrysalis.L72-nsu40',         n='v3 L72sm',         c='red',  p='/lcrc/group/e3sm/ac.whannah/E3SMv3_dev',s='archive/atm/hist')
#add_case('20230228.amip.v3atm.FourthSmoke.chrysalis.errgw-02',          n='v3 effgw-02',      c='green',p='/lcrc/group/e3sm/ac.whannah/E3SMv3_dev',s='run')
#add_case('20230228.amip.v3atm.FourthSmoke.chrysalis.L72-nsu40.errgw-02',n='v3 L72sm+effgw-02',c='blue', p='/lcrc/group/e3sm/ac.whannah/E3SMv3_dev',s='run')

# /lcrc/group/e3sm/ac.golaz/E3SMv3_dev/20230127.amip.v3atm.FourthSmoke.chrysalis
# /lcrc/group/e3sm/ac.whannah/E3SMv3_dev/20230214.amip.v3atm.FourthSmoke.chrysalis.L72-nsu40/
# /lcrc/group/e3sm/ac.whannah/E3SMv3_dev/20230228.amip.v3atm.FourthSmoke.chrysalis.errgw-02
# /lcrc/group/e3sm/ac.whannah/E3SMv3_dev/20230228.amip.v3atm.FourthSmoke.chrysalis.L72-nsu40.errgw-02

# ERA5 levels
# lev = np.array([   1.,    2.,    3.,    5.,    7.,   10.,   20.,   30.,   50.,   70.,
#                  100.,  125.,  150.,  175.,  200.,  225.,  250.,  300.,  350.,  400.,
#                  450.,  500.,  550.,  600.,  650.,  700.,  750.,  775.,  800.,  825.,
#                  850.,  875.,  900.,  925.,  950.,  975., 1000.])
# lev = np.array([   1.,    2.,    3.,    5.,    7.,   10.,   20.,   30.,   50.,   70.,  100.,  125.,  150.,])
lev = np.array([20.])
# lev = np.array([40.])

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

var = ['U']

fig_type = 'png'
fig_file = 'figs/FXX-QBO-power-spectra'

lat1,lat2 = -5,5
# lat1,lat2 = -10,10

htype,first_file,num_files = 'h0',0,12*10

print_stats = True

var_x_case = False

use_common_label_bar = True

num_plot_col = 1

recalculate_timeseries = True

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
    period0 = 1 / frequency0

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
#---------------------------------------------------------------------------------------------------
for v in range(num_var):
   hc.printline()
   print('  var: '+hc.tcolor.MAGENTA+var[v]+hc.tcolor.ENDC)
   tvar = var[v]
   area_name = 'area'
   #----------------------------------------------------------------------------
   # read the data
   #----------------------------------------------------------------------------
   data_list,lev_list = [],[]
   period_list, power_list = [],[]
   for c in range(num_case):
      # tmp_file = os.getenv('HOME')+f'/Research/E3SM/data_temp/QBO.power-spectra.v1.{case[c]}.{var[v]}.lat1_{lat1}.lat2_{lat2}.nc'
      tmp_file = os.getenv('HOME')+f'/Research/E3SM/data_temp/QBO.power-spectra.v1.{case[c]}.{var[v]}.lev_{lev[0]}.lat1_{lat1}.lat2_{lat2}.nc'

      print('    case: '+hc.tcolor.GREEN+case[c]+hc.tcolor.ENDC)
      print('    time series file: '+tmp_file)

      data_dir_tmp,data_sub_tmp = None, None
      if case_dir[c] is not None: data_dir_tmp = case_dir[c]
      if case_sub[c] is not None: data_sub_tmp = case_sub[c]

      case_obj = he.Case( name=case[c], atm_comp=get_comp(case[c]), 
                          data_dir=data_dir_tmp, data_sub=data_sub_tmp, time_freq=None )

      if 'lat1' in vars() : case_obj.lat1,case_obj.lat2 = lat1,lat2
      if 'lon1' in vars() : case_obj.lon1,case_obj.lon2 = lon1,lon2

      # avoid creating large chunks
      with dask.config.set(**{'array.slicing.split_large_chunks': True}):  
         if recalculate_timeseries:
            data = case_obj.load_data(tvar,    
                                      component=get_comp(case[c]),
                                      htype=htype,
                                      first_file=first_file,
                                      num_files=num_files,
                                      lev=lev)
            area = case_obj.load_data(area_name,
                                      component=get_comp(case[c]),
                                      htype=htype,
                                      first_file=first_file,
                                      num_files=num_files).astype(np.double)

            # hc.print_time_length(data.time,indent=' '*4)

            data_avg = ( (data*area).sum(dim='ncol') / area.sum(dim='ncol') )
            
            ### print stats after time averaging
            # if print_stats: hc.print_stat(data_avg,name=var[v],stat='naxsh',indent='    ',compact=True)

            ds_out = xr.Dataset( coords=data_avg.coords )
            ds_out[var[v]] = data_avg
            ds_out.to_netcdf(path=tmp_file,mode='w')
         #----------------------------------------------------------------------
         else:
            tmp_ds = xr.open_dataset( tmp_file, use_cftime=True  )
            data_avg = tmp_ds[var[v]]

      if 'lev' in data_avg.dims: data_avg = data_avg.isel(lev=0)
      #----------------------------------------------------------------------
      period = np.concatenate( (np.arange(2.0, 33.0), np.arange(34.0, 100.0, 2.0)), axis=0 )
      psd, amp = get_psd_from_deseason(data_avg.values, period)
      #----------------------------------------------------------------------
      period_list.append( period )
      power_list .append( amp )
   #----------------------------------------------------------------------------
   # Create plot
   #----------------------------------------------------------------------------
   tres = copy.deepcopy(res)
   for c in range(num_case):

      res.tiXAxisString = 'Period (months)'
      res.tiYAxisString = 'Amplitude'

      tres.trXMinF = 0
      tres.trXMaxF = 50

      tres.trYMinF = 0
      tres.trYMaxF = np.max([np.nanmax(d) for d in power_list])

      tres.xyLineColor   = clr[c]
      tres.xyDashPattern = dsh[c]

      tplot = ngl.xy(wks, period_list[c], power_list[c], tres)

      ip = v
      if c==0 :
         plot[ip] = tplot
      else:
         ngl.overlay(plot[ip],tplot)

      #-------------------------------------------------------------------------
      # Set strings at top of plot
      #-------------------------------------------------------------------------
      # hs.set_subtitles(wks, plot[len(plot)-1], case_name[c], '', '', font_height=0.015)

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

import os, glob, ngl, subprocess as sp, numpy as np, xarray as xr, copy, string, dask
import hapy_common as hc, hapy_E3SM   as he, hapy_setres as hs
import cmocean
import scipy, numba
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
var,vlev,vclr,vdsh = [],[],[],[]
def add_var(var_name,lev=-1,c='black',d=0): 
   var.append(var_name); vlev.append(lev)
   vclr.append(c); vdsh.append(d)
#-------------------------------------------------------------------------------
##------------------------------------------------------------------------------
cscratch = '/global/cscratch1/sd/whannah/e3sm_scratch/cori-knl'
gscratch = '/global/cfs/cdirs/m3312/whannah/2022-QBO-TEST'
pscratch = '/pscratch/sd/w/whannah/e3sm_scratch/pm-cpu'
##------------------------------------------------------------------------------
add_case('E3SM.QBO-TEST.F2010.ne30pg2.L72.01',       n='E3SM control',     c='black',p=gscratch,s='run')
add_case('E3SM.QBO-TEST.F2010.ne30pg2.L72-nsu40.01', n='E3SM L72 smoothed',c='red',  p=gscratch,s='run')
add_case('E3SM.QBO-TEST.F2010.ne30pg2.L80-rsu40.01', n='E3SM L80 refined', c='blue', p=gscratch,s='run')

# add_case('E3SM.QBO-TEST-02.F2010.ne30pg2.L72.GW-MOD-0', n='E3SM L72 M0',c='red'  ,d=1,p=pscratch,s='run')
# add_case('E3SM.QBO-TEST-02.F2010.ne30pg2.L72.GW-MOD-1', n='E3SM L72 M1',c='green',d=1,p=pscratch,s='run')
#-------------------------------------------------------------------------------

# ERA5 levels
# lev = np.array([   1.,    2.,    3.,    5.,    7.,   10.,   20.,   30.,   50.,   70.,
#                  100.,  125.,  150.,  175.,  200.,  225.,  250.,  300.,  350.,  400.,
#                  450.,  500.,  550.,  600.,  650.,  700.,  750.,  775.,  800.,  825.,
#                  850.,  875.,  900.,  925.,  950.,  975., 1000.])
lev = np.array([   1.,    2.,    3.,    5.,    7.,   10.,   20.,   30.,   50.,   70.,  100.,  125.,  150.,])


#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

# add_var('u',c='black')
add_var('utendepfd',c='blue')
add_var('BUTGWSPEC',c='magenta')

fig_type = 'png'
fig_file = 'figs/FXX-TEM-amplitude-profile'

tmp_file_head = 'data/TEM-amplitude-profile'

lat1,lat2 = -5,5
# lat1,lat2 = -10,10

# htype,first_file,num_files = 'h0',0,12*2
# htype,first_file,num_files = 'h0',0,12*20
remap_str,search_str = 'remap_90x180','h0.tem.'; first_file,num_files = 12*0,12*20

print_stats = True

var_x_case = False

use_common_label_bar = True

num_plot_col = 2

recalculate_timeseries = True # this is commented out

add_obs = False
obs_clr = 'gray'

pline_lev = 20

#---------------------------------------------------------------------------------------------------
# Set up plot resources
#---------------------------------------------------------------------------------------------------
num_var,num_case = len(var),len(case)

wkres = ngl.Resources()
npix = 2048; wkres.wkWidth,wkres.wkHeight=npix,npix
wks = ngl.open_wks(fig_type,fig_file,wkres)
plot = [None]*num_var

res = hs.res_xy()
res.vpHeightF = 0.5
# res.xyMarkLineMode = "MarkLines"
res.xyMarkerSizeF = 0.008
res.xyMarker = 16
res.xyLineThicknessF = 8
# res.tmYLLabelFontHeightF         = 0.008
# res.tmXBLabelFontHeightF         = 0.008
# res.tmXBAutoPrecision = False
# res.tmXBPrecision = 2
res.trYReverse = True
res.xyYStyle ="Log"

#-------------------------------------------------------------------------------
# plot legend in separate file
#-------------------------------------------------------------------------------
if num_case>1:
   legend_file = fig_file+'.legend'
   wkres = ngl.Resources() #; npix = 1024 ; wkres.wkWidth,wkres.wkHeight=npix,npix
   lgd_wks = ngl.open_wks('png',legend_file,wkres)
   lgres = ngl.Resources()
   lgres.vpWidthF           = 0.05
   lgres.vpHeightF          = 0.03*num_case
   lgres.lgLabelFontHeightF = 0.008
   lgres.lgLabelFont        = "courier"
   lgres.lgMonoDashIndex    = False
   lgres.lgLineLabelsOn     = False
   lgres.lgLineThicknessF   = 16
   lgres.lgLabelJust        = 'CenterLeft'
   labels = case_name
   
   lclr = copy.deepcopy(clr)
   ldsh = copy.deepcopy(dsh)
   if add_obs: 
      labels.insert(0,'ERA-interim')
      lclr.insert(0,obs_clr)
      ldsh.insert(0,0)
      lgres.vpHeightF          = 0.03*(num_case+1)
   
   lgres.lgLineColors       = lclr
   lgres.lgDashIndexes      = ldsh

   for i in range(len(labels)): labels[i] = ' '*4+labels[i] 

   pid = ngl.legend_ndc(lgd_wks, len(labels), labels, 0.5, 0.65, lgres)

   ngl.frame(lgd_wks)
   hc.trim_png(legend_file)

# exit()

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
def get_20to40month_fft_amplitude(qboN, levelN):
    # Calculates the amplitude of wind variations in the 20 - 40 month period
    psd_sumN = np.zeros(levelN.shape,dtype="complex_")
    amplitudeN = np.zeros(levelN.shape)

    for ilev in np.arange(len(levelN)):
        # `qboN[:, ilev]` returns the entire 0th dimension for ilev in the 1st dimension of the array.
        y_input = deseason(np.squeeze(qboN[:, ilev]))
        y = scipy.fftpack.fft(y_input)
        n = len(y)
        frequency = np.arange(n / 2) / n
        period = np.full(len(frequency),np.inf)
        period[1:] = 1 / frequency[1:]
        values = y[0 : int(np.floor(n / 2))]
        fyy = values * np.conj(values)
        # Choose the range 20 - 40 months that captures most QBOs (in nature)
        psd_sumN[ilev] = 2 * np.nansum(fyy[(period <= 40) & (period >= 20)])
        amplitudeN[ilev] = np.real( np.sqrt(2 * psd_sumN[ilev]) * (frequency[1] - frequency[0]) )
    return psd_sumN, amplitudeN
#---------------------------------------------------------------------------------------------------
# @numba.njit()
@numba.jit(nopython=True)
def calculate_obs_area(lon,lat,lon_bnds,lat_bnds):
   re = 6.37122e06  # radius of earth
   nlat,nlon = len(lat),len(lon)
   area = np.empty((nlat,nlon),np.float64)
   for j in range(nlat):
      for i in range(nlon):
         dlon = np.absolute( lon_bnds[j,1] - lon_bnds[j,0] )
         dlat = np.absolute( lat_bnds[j,1] - lat_bnds[j,0] )
         dx = re*dlon*np.pi/180.
         dy = re*dlat*np.pi/180.
         area[j,i] = dx*dy
   return area
#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------
for v in range(num_var):
   lev_list, amp_list = [],[]
   hc.printline()
   print('  var: '+hc.tcolor.MAGENTA+var[v]+hc.tcolor.ENDC)
   tvar = var[v]
   if var[v]=='U': ovar = 'ua'
   #----------------------------------------------------------------------------
   #----------------------------------------------------------------------------
   if var[v]=='U' and add_obs:
      # obs_data_file = '/lcrc/group/e3sm/diagnostics/observations/Atm/time-series/ERA-Interim/ua_197901_201612.nc'
      obs_data_file = '/global/cfs/cdirs/e3sm/diagnostics/observations/Atm/time-series/ERA-Interim/ua_197901_201612.nc'
      ds = xr.open_dataset(obs_data_file)
      area = calculate_obs_area(ds['lon'].values,ds['lat'].values,ds['lon_bnds'].values,ds['lat_bnds'].values)
      area = xr.DataArray( area, coords=[ds['lat'],ds['lon']] )
      data = ds[ovar]
      data = data.sel(lat=slice(lat1,lat2))
      area = area.sel(lat=slice(lat1,lat2))
      data_avg = (data*area).sum(dim=('lon','lat')) / area.sum(dim=('lon','lat'))

      psd, amp = get_20to40month_fft_amplitude(data_avg.values, data_avg['plev'].values)
      lev_list.append( data_avg['plev'].values/1e2 )
      amp_list.append( amp )
   #----------------------------------------------------------------------------
   # read the data
   #----------------------------------------------------------------------------
   for c in range(num_case):
      tmp_file = f'{tmp_file_head}.{case[c]}.{var[v]}'
      if 'lat1' in locals(): tmp_file = tmp_file + f'.lat1_{lat1}.lat2_{lat2}'
      if 'lon1' in locals(): tmp_file = tmp_file + f'.lon1_{lon1}.lon2_{lon2}'
      tmp_file = f'{tmp_file}.nc'

      print('    case: '+hc.tcolor.GREEN+case[c]+hc.tcolor.ENDC)
      print('    time series file: '+tmp_file)

      

      # avoid creating large chunks
      # with dask.config.set(**{'array.slicing.split_large_chunks': True}):  

      if recalculate_timeseries:
         file_path = f'{gscratch}/{case[c]}/data_{remap_str}_tem/*.eam.{search_str}*'
         file_list = sorted(glob.glob(file_path))

         if 'first_file' in locals(): file_list = file_list[first_file:]
         if 'num_files'  in locals(): file_list = file_list[:num_files]

         #----------------------------------------------------------------------
         # manually combine files due to issues with open_mfdataset on Perlmutter
         ds = xr.open_dataset(file_list[0])
         for f in file_list[1:]:
           ds_tmp = xr.open_dataset(f)
           ds = xr.concat([ds, ds_tmp], dim='time')
         #----------------------------------------------------------------------

         # Load data 
         data = ds[var[v]].sel(lat=slice(lat1,lat2),plev=slice(200e2,0))
         data_avg = data.mean(dim=('lat'))
         
         ### print stats after averaging
         if print_stats: hc.print_stat(data_avg,name=var[v],stat='naxsh',indent='    ',compact=True)

      #    ds_out = xr.Dataset( coords=data_avg.coords )
      #    ds_out[var[v]] = data_avg
      #    ds_out.to_netcdf(path=tmp_file,mode='w')
      # #----------------------------------------------------------------------
      # else:
      #    tmp_ds = xr.open_dataset( tmp_file, use_cftime=True  )
      #    data_avg = tmp_ds[var[v]]

      #-------------------------------------------------------------------------
      # Unit conversions

      unit_str = ''
      if var[v]=='dudt'        : data_avg = data_avg*86400.; unit_str = 'm/s/day'
      if var[v]=='RES'         : data_avg = data_avg*86400.; unit_str = 'm/s/day'
      if var[v]=='TEND'        : data_avg = data_avg*86400.; unit_str = 'm/s/day'
      if var[v]=='utendepfd'   : data_avg = data_avg*86400.; unit_str = 'm/s/day'
      if var[v]=='utendepfd_lp': data_avg = data_avg*86400.; unit_str = 'm/s/day'
      if var[v]=='utendvtem'   : data_avg = data_avg*86400.; unit_str = 'm/s/day'
      if var[v]=='utendwtem'   : data_avg = data_avg*86400.; unit_str = 'm/s/day'
      if var[v]=='epfy'        : data_avg = data_avg*86400.; unit_str = 'm3/s/day'
      if var[v]=='epfz'        : data_avg = data_avg*86400.; unit_str = 'm3/s/day'
      if var[v]=='depfydy'     : data_avg = data_avg*86400.; unit_str = ''
      if var[v]=='depfzdp'     : data_avg = data_avg*86400.; unit_str = ''
      if var[v]=='psitem'      : data_avg = data_avg*86400.; unit_str = 'kg/day'
      if var[v]=='vtem'        :                             unit_str = 'm/s'
      if var[v]=='wtem'        :                             unit_str = 'm/s'

      if var[v]=='BUTGWSPEC'   : data_avg = data_avg*86400.; unit_str = 'm/s/day'
      if var[v]=='UTGWSPEC'    : data_avg = data_avg*86400.; unit_str = 'm/s/day'
      if var[v]=='UTGWORO'     : data_avg = data_avg*86400.; unit_str = 'm/s/day'

      #----------------------------------------------------------------------
      psd, amp = get_20to40month_fft_amplitude(data_avg.values, data_avg['plev'].values)


      #----------------------------------------------------------------------
      lev_list.append( data_avg['plev'].values/1e2 )
      amp_list.append( amp )

   #----------------------------------------------------------------------------
   # Create plot
   #----------------------------------------------------------------------------
   tres = copy.deepcopy(res)
   tres.tiYAxisString = 'Pressure (hPa)'
   if var[v]=='u'        : tres.tiXAxisString = 'Amplitude (m/s)'
   if var[v]=='utendepfd': tres.tiXAxisString = 'Amplitude (m/s/day)'
   if var[v]=='BUTGWSPEC': tres.tiXAxisString = 'Amplitude (m/s/day)'
   tres.trXMinF = 0
   # tres.trXMaxF = 25 # 
   tres.trXMaxF = np.max([np.nanmax(d) for d in amp_list])
   tres.trYMinF = 1e0
   tres.trYMaxF = 1e2

   for c in range( (num_case+1) if add_obs else num_case ):

      if add_obs:
         if c==0:
            tres.xyLineColor,tres.xyDashPattern = obs_clr,0
         else:
            tres.xyLineColor,tres.xyDashPattern = clr[c-1],dsh[c-1]
      else:
         tres.xyLineColor,tres.xyDashPattern = clr[c],dsh[c]

      tplot = ngl.xy(wks, amp_list[c], lev_list[c], tres)

      ip = v
      if c==0 :
         plot[ip] = tplot
      else:
         ngl.overlay(plot[ip],tplot)

   #-------------------------------------------------------------------------
   # Set strings at top of plot
   var_str = var[v]
   if var[v]=='u'        : var_str = 'Zonal Wind'
   if var[v]=='utendepfd': var_str = 'EP-Flux Divergence'
   if var[v]=='BUTGWSPEC': var_str = 'Convective GW tendency'
   hs.set_subtitles(wks, plot[ip], '', '', var_str, font_height=0.015)

   #-------------------------------------------------------------------------
   # add refence line 
   tres.xyLineColor,tres.xyDashPattern = 'gray',1
   xx = np.array([tres.trXMinF,tres.trXMaxF])
   yy = np.array([pline_lev,pline_lev])
   ngl.overlay(plot[ip],ngl.xy(wks, xx, yy, tres))
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

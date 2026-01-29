import os, ngl, subprocess as sp, numpy as np, xarray as xr, copy, string, dask
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
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

tmp_path,tmp_sub = '/global/cfs/cdirs/m3312/whannah/e3smv3_amip','archive/atm/hist'
add_case('v3.LR.amip_0101', n='EAMv3', c='blue', p=tmp_path,s=tmp_sub)

tmp_path,tmp_sub = '/global/cfs/cdirs/e3smdata/simulations/','archive/atm/hist'
add_case('v2.LR.amip_0101', n='EAMv2', c='red', p=tmp_path,s=tmp_sub) 

lev = np.array([ 1., 2., 3., 5., 7., 10., 20., 30., 50., 70., 100., 125., 150.,])

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

var = ['U']

fig_file,fig_type = 'figs/QBO-amplitude-profile','png'

lat1,lat2 = -5,5
# lat1,lat2 = -10,10

# htype,first_file,num_files = 'h0',12*1,12*1
# htype,first_file,num_files = 'h0',0,12*10
htype,first_file,num_files = 'h0', (1979-1870)*12, 12*35#65

print_stats = True

var_x_case = False

use_common_label_bar = True

num_plot_col = 1

recalculate_timeseries = False

add_obs = True
obs_clr = 'black'

#---------------------------------------------------------------------------------------------------
# Set up plot resources
#---------------------------------------------------------------------------------------------------
num_var,num_case = len(var),len(case)



# if 'lev' not in vars(): lev = np.array([0])

wkres = ngl.Resources()
# npix = 2048; wkres.wkWidth,wkres.wkHeight=npix,npix
wks = ngl.open_wks(fig_type,fig_file,wkres)
plot = [None]*num_var

if clr is None: clr = np.linspace(2,len( ngl.retrieve_colormap(wks) )-1,num_case,dtype=int).tolist()

res = hs.res_xy()
res.vpHeightF = 0.5
# res.xyMarkLineMode = "MarkLines"
res.xyMarkerSizeF = 0.008
res.xyMarker = 16
res.xyLineThicknessF = 12
# res.tmYLLabelFontHeightF         = 0.008
# res.tmXBLabelFontHeightF         = 0.008
# res.tmXBAutoPrecision = False
# res.tmXBPrecision = 2
res.trYReverse = True
res.xyYStyle ="Log"

#-------------------------------------------------------------------------------
# plot legend in separate file
#-------------------------------------------------------------------------------
# if num_case>1:
#    legend_file = fig_file+'.legend'
#    wkres = ngl.Resources() #; npix = 1024 ; wkres.wkWidth,wkres.wkHeight=npix,npix
#    lgd_wks = ngl.open_wks('png',legend_file,wkres)
#    lgres = ngl.Resources()
#    lgres.vpWidthF           = 0.05
#    lgres.vpHeightF          = 0.03*num_case
#    lgres.lgLabelFontHeightF = 0.008
#    lgres.lgLabelFont        = "courier"
#    lgres.lgMonoDashIndex    = False
#    lgres.lgLineLabelsOn     = False
#    lgres.lgLineThicknessF   = 16
#    lgres.lgLabelJust        = 'CenterLeft'
   
#    labels,lclr,ldsh = case_name[::-1],clr[::-1],dsh[::-1]
   
#    if add_obs: 
#       labels.insert(0,'ERA-interim')
#       lclr.insert(0,obs_clr)
#       ldsh.insert(0,0)
#       lgres.vpHeightF          = 0.03*(num_case+1)
#       lgres.lgLineColors       = lclr
#       lgres.lgDashIndexes      = ldsh

#    for i in range(len(labels)): labels[i] = ' '*4+labels[i] 

#    pid = ngl.legend_ndc(lgd_wks, len(labels), labels, 0.5, 0.65, lgres)

#    ngl.frame(lgd_wks)
#    hc.trim_png(legend_file)

# exit('Stopping after legend creation.')

#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------
def get_comp(case):
   comp = 'eam'
   return comp
#---------------------------------------------------------------------------------------------------
def deseason(xraw):
    # Calculates the deseasonalized data
    months_per_year = 12
    ### Create array to hold climatological values and deseasonalized data
    ### Create months_per_year x 1 array of zeros
    # xclim = np.zeros((months_per_year, 1))
    xclim = np.zeros(months_per_year)
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
      obs_data_file = '/global/cfs/cdirs/e3sm/diagnostics/observations/Atm/time-series/ERA5/ua_197901_201912.nc'
      ds = xr.open_dataset(obs_data_file)
      area = calculate_obs_area(ds['lon'].values,ds['lat'].values,ds['lon_bnds'].values,ds['lat_bnds'].values)
      area = xr.DataArray( area, coords=[ds['lat'],ds['lon']] )
      data = ds[ovar]
      data = data.sel(lat=slice(lat1,lat2))
      area = area.sel(lat=slice(lat1,lat2))
      data_avg = (data*area).sum(dim=('lon','lat')) / area.sum(dim=('lon','lat'))

      hc.print_stat(data_avg,name=var[v]+' (ERAi)',stat='naxsh',indent='    ',compact=True)

      psd, amp = get_20to40month_fft_amplitude(data_avg.values, data_avg['plev'].values)
      lev_list.append( data_avg['plev'].values/1e2 )
      amp_list.append( amp )
   #----------------------------------------------------------------------------
   # read the data
   #----------------------------------------------------------------------------
   for c in range(num_case):

      tmp_file = f'data/QBO.amplitude-profile.v1.{case[c]}.{var[v]}.lev_{lev[0]}.lat1_{lat1}.lat2_{lat2}.nc'

      print('    case: '+hc.tcolor.GREEN+case[c]+hc.tcolor.ENDC)
      print('    time series file: '+tmp_file)

      # avoid creating large chunks
      with dask.config.set(**{'array.slicing.split_large_chunks': True}):  
         if recalculate_timeseries:

            data_dir_tmp,data_sub_tmp = None, None
            if case_dir[c] is not None: data_dir_tmp = case_dir[c]
            if case_sub[c] is not None: data_sub_tmp = case_sub[c]
            case_obj = he.Case( name=case[c], atm_comp=get_comp(case[c]), 
                                data_dir=data_dir_tmp, data_sub=data_sub_tmp, time_freq=None )
            if 'lat1' in vars() : case_obj.lat1,case_obj.lat2 = lat1,lat2
            if 'lon1' in vars() : case_obj.lon1,case_obj.lon2 = lon1,lon2

            data = case_obj.load_data(tvar, htype=htype,
                                      first_file=first_file,
                                      num_files=num_files,
                                      lev=lev)
            # print(); print(data)
            # exit()
            area = case_obj.load_data('area',htype=htype,num_files=1).astype(np.double)

            hc.print_stat(data_avg,name=var[v],stat='naxsh',indent='    ',compact=True)

            # hc.print_time_length(data.time,indent=' '*4)

            data_avg = ( (data*area).sum(dim='ncol') / area.sum(dim='ncol') )
            
            ### print stats after time averaging
            hc.print_stat(data_avg,name=var[v],stat='naxsh',indent='    ',compact=True)

            ds_out = xr.Dataset( coords=data_avg.coords )
            ds_out[var[v]] = data_avg
            ds_out.to_netcdf(path=tmp_file,mode='w')
         #----------------------------------------------------------------------
         else:
            tmp_ds = xr.open_dataset( tmp_file, use_cftime=True  )
            data_avg = tmp_ds[var[v]]

      # if 'lev' in data_avg.dims: data_avg = data_avg.isel(lev=0)
      #----------------------------------------------------------------------
      psd, amp = get_20to40month_fft_amplitude(data_avg.values, data_avg['lev'].values)

      #----------------------------------------------------------------------
      lev_list.append( data_avg['lev'].values )
      amp_list.append( amp )

   #----------------------------------------------------------------------------
   # Create plot
   #----------------------------------------------------------------------------
   tres = copy.deepcopy(res)
   tres.tiYAxisString = 'Pressure [mb]'
   tres.tiXAxisString = 'Amplitude [m/s]'
   tres.trXMinF = 0
   # tres.trXMaxF = 25 # 
   tres.trXMaxF = np.max([np.nanmax(d) for d in amp_list])
   tres.trYMinF = 1e0
   tres.trYMaxF = 1e2

   for c in range( (num_case+1) if add_obs else num_case ):

      if add_obs:
         if c==0:
            tres.xyLineColor,tres.xyDashPattern = obs_clr,0
            tres.xyLineThicknessF = 16
         else:
            tres.xyLineColor,tres.xyDashPattern = clr[c-1],dsh[c-1]
            tres.xyLineThicknessF = 8
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
      #-------------------------------------------------------------------------
      # hs.set_subtitles(wks, plot[len(plot)-1], case_name[c], '', '', font_height=0.015)

#----------------------------------------------------------------------------
# Add legend
#----------------------------------------------------------------------------

if add_obs:
   case_name.insert(0,'ERA5')
   dsh.insert(0,0)
   clr.insert(0,obs_clr)

lgres = ngl.Resources()
lgres.vpWidthF, lgres.vpHeightF  = 0.08, 0.15
lgres.lgLabelFontHeightF = 0.015
lgres.lgLineThicknessF   = 15
lgres.lgMonoLineColor    = False
# lgres.lgMonoDashIndex    = True
# lgres.lgDashIndex        = 0
lgres.lgDashIndexes      = dsh
lgres.lgLineColors       = clr
lgres.lgLabelJust    = 'CenterLeft'
lbl = [f'  {n}' for n in case_name]
pid = ngl.legend_ndc(wks, len(lbl), lbl, 0.8, 0.8, lgres)

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

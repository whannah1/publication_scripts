import os, copy, ngl, xarray as xr, numpy as np
import hapy_common as hc, hapy_E3SM   as he, hapy_setres as hs
#---------------------------------------------------------------------------------------------------
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
#---------------------------------------------------------------------------------------------------
var,mask_flag = [],[]
def add_var(var_name,mask=None): 
   var.append(var_name) ; mask_flag.append(mask)
#---------------------------------------------------------------------------------------------------
tmp_path_co2_mmf = '/global/cfs/cdirs/m3312/whannah/2023-CPL'
tmp_path_hst_mmf = '/global/cfs/cdirs/m3312/whannah/2023-CPL'
tmp_path_hst_v2  = '/global/cfs/cdirs/m3312/whannah/e3smv2_historical'
tmp_sub = 'archive/atm/hist'
add_case('E3SM.2023-CO2-TEST-01.GNUGPU.ne30pg2_EC30to60E2r2.WCYCL1850-MMF1.1xCO2',n='MMF PI 1xCO2',c='blue' ,p=tmp_path_co2_mmf,s=tmp_sub)
add_case('E3SM.2023-CO2-TEST-01.GNUGPU.ne30pg2_EC30to60E2r2.WCYCL1850-MMF1.2xCO2',n='MMF PI 2xCO2',c='green',p=tmp_path_co2_mmf,s=tmp_sub)
add_case('E3SM.2023-CO2-TEST-01.GNUGPU.ne30pg2_EC30to60E2r2.WCYCL1850-MMF1.4xCO2',n='MMF PI 4xCO2',c='red'  ,p=tmp_path_co2_mmf,s=tmp_sub)
# add_case('v2.LR.historical_0101',n='E3SMv2',p=tmp_path_hst_v2,s='archive/atm/hist')
# add_case('E3SM.INCITE2023-CPL.ne30pg2_EC30to60E2r2.WCYCL20TR-MMF1',n='E3SM-MMF',p=tmp_path_hst_mmf,s='archive/atm/hist')
#---------------------------------------------------------------------------------------------------

add_var('TS',mask='ocn')


# htype,years,months,first_file,num_files = 'ha',[],[],50,12
htype,years,months,first_file,num_files = 'ha',[],[],50,70

lat1,lat2 = 0,65
lon1,lon2 = 360-80,360

fig_type,fig_file = 'png','figs/FXX-timeseries-AMO'
tmp_file_head = 'data/timeseries-AMO'

recalculate   = True
print_stats   = True
overlay_cases = True

add_trend = False

num_plot_col  = 2

#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------
num_var = len(var)
num_case = len(case)

if 'clr' not in vars() or clr==[]: clr = ['black']*num_case
if 'dsh' not in vars() or dsh==[]: dsh = np.zeros(num_case)
#-------------------------------------------------------------------------------
# plot legend in separate file
#-------------------------------------------------------------------------------
if num_case>1:
   legend_file = fig_file+'.legend'
   wkres = ngl.Resources() ; npix = 2048 ; wkres.wkWidth,wkres.wkHeight=npix,npix
   lgd_wks = ngl.open_wks('png',legend_file,wkres)
   lgres = ngl.Resources()
   # lgres.vpWidthF           = 0.06
   # lgres.vpHeightF          = 0.04
   lgres.lgLabelFontHeightF = 0.008
   lgres.lgLabelFont        = 'courier'
   lgres.lgMonoDashIndex    = False
   lgres.lgLineLabelsOn     = False
   lgres.lgLineThicknessF   = 4#16
   lgres.lgLabelJust        = 'CenterLeft'
   lgres.lgLineColors       = clr
   lgres.lgDashIndexes      = dsh

   labels = case_name
   for i in range(len(labels)): labels[i] = ' '*4+labels[i] 

   pid = ngl.legend_ndc(lgd_wks, len(labels), labels, 0.5, 0.65, lgres)

   ngl.frame(lgd_wks)
   hc.trim_png(legend_file)
   # exit()

#---------------------------------------------------------------------------------------------------
# Set up plot resources
#---------------------------------------------------------------------------------------------------
if overlay_cases:
   plot = [None]*(num_var)
else:
   plot = [None]*(num_case*num_var)

wkres = ngl.Resources()
npix=1024; wkres.wkWidth,wkres.wkHeight=npix,npix
wks = ngl.open_wks(fig_type,fig_file,wkres)

if 'legend_file' in locals(): lgd_wks = ngl.open_wks('png',legend_file,wkres)


res = hs.res_xy()
# res.vpHeightF = 0.5
res.vpHeightF = 0.2
res.tmYLLabelFontHeightF         = 0.015
res.tmXBLabelFontHeightF         = 0.015
res.tiXAxisFontHeightF           = 0.015
res.tiYAxisFontHeightF           = 0.015
res.xyLineThicknessF = 5

lres = hs.res_xy()
lres.xyDashPattern    = 1
lres.xyLineThicknessF = 1
lres.xyLineColor      = 'black'

#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------
for v in range(num_var):
   
   hc.printline()
   print(' '*2+'var: '+hc.tclr.MAGENTA+var[v]+hc.tclr.END)

   time_list,data_list = [],[]

   for c in range(num_case):

      print(' '*4+'case: '+hc.tclr.CYAN+case[c]+hc.tclr.END)

      tmp_file = f'{tmp_file_head}.{case[c]}.{var[v]}.nc'

      if recalculate:
         #----------------------------------------------------------------------
         data_dir_tmp,data_sub_tmp = None, None
         if case_dir[c] is not None: data_dir_tmp = case_dir[c]
         if case_sub[c] is not None: data_sub_tmp = case_sub[c]

         case_obj = he.Case( name=case[c], data_dir=data_dir_tmp, data_sub=data_sub_tmp  )

         tvar = var[v]

         if 'lat1' in vars() : case_obj.lat1,case_obj.lat2 = lat1,lat2
         if 'lon1' in vars() : case_obj.lon1,case_obj.lon2 = lon1,lon2

         #----------------------------------------------------------------------
         # read the data
         area = case_obj.load_data('area',      htype=htype,num_files=1).astype(np.double)
         time = case_obj.load_data('time_bnds', htype=htype,first_file=first_file,num_files=num_files).isel(nbnd=0)
         data = case_obj.load_data(tvar,        htype=htype,first_file=first_file,num_files=num_files)
         #----------------------------------------------------------------------
         # apply mask
         if mask_flag[v] is not None :
            land_frac = case_obj.load_data('LANDFRAC',htype='h0',first_file=first_file,num_files=1).astype(np.double)
            if 'time' in land_frac.dims: land_frac = land_frac.isel(time=0)
            land_frac = land_frac / land_frac.max()
            mask = xr.DataArray( np.ones(land_frac.shape,dtype=bool), dims=land_frac.dims )
            if mask_flag[v]=='lnd': mask = mask & (land_frac.values>0.5)
            if mask_flag[v]=='ocn': mask = mask & (land_frac.values<0.5)
            data = data.where( mask, drop=True)
         #----------------------------------------------------------------------
         # area weighted spatial mean
         data = ( (data*area).sum(dim='ncol') / area.sum(dim='ncol') )

         # Make time start at zero and convert to years
         data['time'] = ( data['time'] - data['time'][0] ).astype('float') / 86400e9 / 365
      
         # convert to anomalies
         data = data - data.mean()

         # detrend in time
         fit = xr.polyval(data['time'], data.polyfit(dim='time', deg=1).polyfit_coefficients)
         data = data - fit

         # perform rolling average
         data = data.rolling(time=10, center=True).mean()

         #----------------------------------------------------------------------
         # write to file
         tmp_ds = xr.Dataset()
         tmp_ds['data'] = data
         tmp_ds.to_netcdf(path=tmp_file,mode='w')
      else:
         #----------------------------------------------------------------------
         # read from file
         tmp_ds = xr.open_dataset( tmp_file )
         data = tmp_ds['data']
      
      #-------------------------------------------------------------------------
      time_mean = data.mean(dim='time').values
      print(' '*6+f'Area Weighted Time Mean : {hc.tclr.GREEN}{time_mean}{hc.tclr.END}')
      if print_stats: hc.print_stat(data,name='',stat='naxs',indent=' '*6,compact=True)

      data_list.append( np.ma.masked_invalid(data.values))
      time_list.append( data['time'].values )

      #-------------------------------------------------------------------------
      # write to file
      #-------------------------------------------------------------------------
      # if write_file : 
      #    tfile = f'/global/homes/w/whannah/E3SM/scratch/{case[0]}/run/{case[0]}.time_series.{var[v]}.nc'
      #    print('writing to file: '+tfile)
      #    avg_X.name = var[v]
      #    avg_X.to_netcdf(path=tfile,mode='w')
      #    exit()

   #----------------------------------------------------------------------------
   # Create plot
   #----------------------------------------------------------------------------
   tres = copy.deepcopy(res)
   tres.tiXAxisString = 'Time [years]'

   ### Make sure plot bounds are consistent
   # tres.trYMinF = 14.4
   tres.trYMinF = np.min([np.nanmin(d) for d in data_list])
   tres.trYMaxF = np.max([np.nanmax(d) for d in data_list])
   tres.trXMinF = np.min([np.nanmin(d) for d in time_list])
   tres.trXMaxF = np.max([np.nanmax(d) for d in time_list])

   print()
   print(f'tres.trYMinF = {tres.trYMinF}')
   print(f'tres.trYMaxF = {tres.trYMaxF}')
   print(f'tres.trXMinF = {tres.trXMinF}')
   print(f'tres.trXMaxF = {tres.trXMaxF}')
   print()

   if var[v]=='NET_TOA_RAD':
      tres.trYMinF = -20
      tres.trYMaxF =  20

   tmp_num_case = num_case

   for c in range(tmp_num_case):
      ip = c*num_var + v
      if overlay_cases: ip = v
      
      tres.xyLineColor   = clr[c]
      tres.xyDashPattern = dsh[c]
      tplot = ngl.xy(wks, time_list[c], data_list[c], tres)
      
      if overlay_cases: 
         if c==0: 
            plot[ip] = tplot
         else:
            ngl.overlay(plot[ip],tplot)
      else:
         plot[ip] = tplot
   #----------------------------------------------------------------------------
   xx = np.array([-1e8,1e8])
   yy = np.array([0,0])
   ngl.overlay( plot[ip], ngl.xy(wks, xx, yy, lres) )
   #----------------------------------------------------------------------------
   # Set strings at top of plot
   var_str = var[v]
   lft_str = ''
   ctr_str = ''

   lat_chk,lon_chk = 'lat1' in locals(), 'lon1' in locals()

   if not lat_chk and not lon_chk : 
      if mask_flag[v] is None: ctr_str = 'Global'
      if mask_flag[v]=='lnd' : ctr_str = 'Land Only'
      if mask_flag[v]=='ocn' : ctr_str = 'Ocean Only'
      if mask_flag[v]=='lnd/ocn ratio':ctr_str = 'Land / Ocean Ratio'
   else:
      if lat_chk:      ctr_str += f' {lat1}:{lat2}N '
      if lon_chk:      ctr_str += f' {lon1}:{lon2}E '
   
   

   if overlay_cases:
      hs.set_subtitles(wks, plot[ip], var_str, ctr_str, '', font_height=0.012)
   else:
      hs.set_subtitles(wks, plot[ip], case_name[c], ctr_str, var_str, font_height=0.015)

#---------------------------------------------------------------------------------------------------
# Finalize plot
#---------------------------------------------------------------------------------------------------
hc.printline()

if 'num_plot_col' in locals():
   if overlay_cases :
      layout = [int(np.ceil(len(plot)/float(num_plot_col))),num_plot_col]
   else:
      if num_case==1 or num_var==1:
         layout = [int(np.ceil(len(plot)/float(num_plot_col))),num_plot_col]
      else:
         layout = [num_var,num_case]
         # layout = [num_case,num_var]
else:
   layout = [num_var,num_case]


ngl.panel(wks,plot[0:len(plot)],layout,hs.setres_panel())
ngl.end()

hc.trim_png(fig_file)

if 'legend_file' in locals(): hc.trim_png(legend_file)

#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------

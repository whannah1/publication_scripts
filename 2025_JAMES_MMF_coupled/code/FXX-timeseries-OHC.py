import os, copy, glob, ngl, xarray as xr, numpy as np, warnings, numba
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
var,lev_list,mask_flag,var_str = [],[],[],[]
def add_var(var_name,lev=-1,mask=None,vstr=None): 
   var.append(var_name); lev_list.append(lev),mask_flag.append(mask)
   if vstr is None: vstr = var_name
   var_str.append(vstr)
#---------------------------------------------------------------------------------------------------
tmp_path_co2_mmf = '/global/cfs/cdirs/m3312/whannah/2023-CPL'
tmp_path_hst_mmf = '/global/cfs/cdirs/m3312/whannah/2023-CPL'
tmp_path_hst_v2  = '/global/cfs/cdirs/m3312/whannah/e3smv2_historical'
tmp_sub = 'archive/ocn/hist'

add_case('v2.LR.historical_0101',                                  n='E3SMv2',  c='cyan',     p=tmp_path_hst_v2, s=tmp_sub)
add_case('v2.LR.historical_0151',                                  n='E3SMv2',  c='orange',   p=tmp_path_hst_v2, s=tmp_sub)
add_case('v2.LR.historical_0201',                                  n='E3SMv2',  c='green',    p=tmp_path_hst_v2, s=tmp_sub)
add_case('v2.LR.historical_0251',                                  n='E3SMv2',  c='purple',   p=tmp_path_hst_v2, s=tmp_sub)
add_case('v2.LR.historical_0301',                                  n='E3SMv2',  c='pink',     p=tmp_path_hst_v2, s=tmp_sub)
# add_case('v2.LR.historical',                                       n='E3SMv2 Ens Mean',  c='red', p=tmp_path_hst_v2, s=tmp_sub)
add_case('E3SM.INCITE2023-CPL.ne30pg2_EC30to60E2r2.WCYCL20TR-MMF1',n='E3SM-MMF',c='blue',p=tmp_path_hst_mmf,s=tmp_sub)

ens_spread_clr = 'lightpink'
#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------

# add_var('oceanHeatContentSfcTo700m',vstr=f'OHC sfc-700m')
# add_var('oceanHeatContent700mTo2000m',vstr=f'OHC 700-2000m')

add_var('meridionalHeatTransportLat',vstr=f'MHT')


htype,yr1,yr2 = 'ha',1950,2014

fig_file,fig_type = 'figs/FXX-timeseries-OHC','png'
tmp_file_head = 'data/timeseries-OHC'

#---------------------------------------------------------------------------------------------------
write_file    = False
print_stats   = True
overlay_cases = True

recalculate = True

add_trend = False

#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------
num_var = len(var)
num_case = len(case)

obs_pos = num_case

# # overlay_cases with single case causes segfault
# if num_case==1 : overlay_cases = False

if 'clr' not in vars() or clr==[]: clr = ['black']*num_case
if 'dsh' not in vars() or dsh==[]: dsh = np.zeros(num_case)

if 'lev' not in vars(): lev = np.array([0])


#---------------------------------------------------------------------------------------------------
# Set up plot resources
#---------------------------------------------------------------------------------------------------
if overlay_cases:
   plot = [None]*(num_var)
else:
   plot = [None]*(num_case*num_var)

wkres = ngl.Resources()
npix=1024*4; wkres.wkWidth,wkres.wkHeight=npix,npix
wks = ngl.open_wks(fig_type,fig_file,wkres)


res = hs.res_xy()
res.vpHeightF = 0.4
# res.vpHeightF = 0.2
res.tmYLLabelFontHeightF   = 0.015
res.tmXBLabelFontHeightF   = 0.015
res.tiXAxisFontHeightF     = 0.015
res.tiYAxisFontHeightF     = 0.015
res.xyLineThicknessF       = 20
res.tiYAxisString          = 'Temperature Anomaly [K]'
res.tiXAxisString          = 'Time [years]'

lres = hs.res_xy()
lres.xyDashPattern    = 1
lres.xyLineThicknessF = 1
lres.xyLineColor      = "black"

#---------------------------------------------------------------------------------------------------
# @numba.njit()
# def parse_time(xtime):
#    ntime = len(xtime)
#    time_list = [None]*ntime
#    for t in range(len(xtime)): 
#       time_list[t] = xtime[t][:10]
#    return time_list
#---------------------------------------------------------------------------------------------------
for v in range(num_var):
   
   hc.printline()
   print('  var: '+hc.tcolor.MAGENTA+var[v]+hc.tcolor.ENDC)

   if 'lev_list' in locals(): lev = lev_list[v]

   time_list,data_list = [],[]
   for c in range(num_case):

      tmp_file = f'{tmp_file_head}.{case[c]}.{var[v]}.nc'

      print('    case: '+hc.tcolor.CYAN+case[c]+hc.tcolor.ENDC)
      print('    time series file: '+tmp_file)

      if recalculate:

         scrip_file_path = '/global/cfs/cdirs/e3sm/inputdata/ocn/mpas-o/EC30to60E2r2/EC30to60E2r2.scrip.201005.nc'
         scrip_ds = xr.open_mfdataset(scrip_file_path).rename({'grid_size':'ncol'})
         area = scrip_ds['grid_area']
         #----------------------------------------------------------------------
         data_root = f'{case_dir[c]}/{case[c]}/{case_sub[c]}'
         file_path = f'{data_root}/*.mpaso.hist.annual.*'
         file_list = sorted(glob.glob(file_path))
         ds = xr.open_mfdataset( file_list )
         ds = ds.rename({'nCells':'ncol'})
         #----------------------------------------------------------------------
         # data_root = f'{case_dir[c]}/{case[c]}/{case_sub[c]}'
         # if 'HeatTransport' in var[v]: file_path = f'{data_root}/*.mpaso.hist.am.meridionalHeatTransport.*'
         # if 'HeatContent'   in var[v]: file_path = f'{data_root}/*.mpaso.hist.am.oceanHeatContent.*'
         # file_list = sorted(glob.glob(file_path))
         # ds = xr.open_mfdataset( file_list, combine='nested', concat_dim='Time' )
         # ds = ds.rename({'nCells':'ncol'})
         # #----------------------------------------------------------------------
         # # Create time coordinate
         # ntime = len(ds['xtime'])
         # time = [None]*ntime
         # with warnings.catch_warnings():
         #    warnings.simplefilter("ignore", category=UserWarning)
         #    for t in range(len(ds['xtime'])): 
         #       time[t] = np.datetime64(ds['xtime'].values[t][:10])
         # time = xr.DataArray(time,coords={'time':time})
         # ds.rename({'Time':'time'})
         # ds.load()
         # ds['time'] = time
         #----------------------------------------------------------------------
         # ntime = len(ds['xtime'])

         # # print()
         # # print(ds['xtime'].values[:])
         # print()
         # print([ds['xtime'].values[t][:10] for t in range(ntime)])
         # exit()
         # # time_list[t] = ds['xtime'].values[:][:10]
         # time = np.datetime64(ds['xtime'].values[:][:10])
         # # time_list = parse_time(ds['xtime'].values)
         # # time = np.datetime64(time_list)
         # time = xr.DataArray(time,coords={'time':time})
         # ds.rename({'Time':'time'})
         # ds['time'] = time
         #----------------------------------------------------------------------
         # ds = ds.where( ds['time.year']>=yr1, drop=True)
         # ds = ds.where( ds['time.year']<=yr2, drop=True)
         #----------------------------------------------------------------------
         data = ds[var[v]]

         print()
         print(data)
         
         data = ( (data*area).sum(dim='ncol') / area.sum(dim='ncol') )

         print()
         print(data)
         exit()
         #----------------------------------------------------------------------
         tmp_ds = xr.Dataset( coords=data.coords )
         tmp_ds[var[v]] = data
         tmp_ds.to_netcdf(path=tmp_file,mode='w')
      #-------------------------------------------------------------------------
      else:
         # if case[c]=='v2.LR.historical':
         #    v2_amip_ens_list = []
         #    v2_amip_ens_list.append('v2.LR.historical_0101')
         #    v2_amip_ens_list.append('v2.LR.historical_0151')
         #    v2_amip_ens_list.append('v2.LR.historical_0201')
         #    v2_amip_ens_list.append('v2.LR.historical_0251')
         #    v2_amip_ens_list.append('v2.LR.historical_0301')
         #    cnt = 0
         #    for e,ens_member in enumerate(v2_amip_ens_list):
         #       tmp_file = f'{tmp_file_head}.{ens_member}.{var[v]}.nc'
         #       tmp_ds = xr.open_dataset( tmp_file, use_cftime=True  )
         #       ens_member_data = tmp_ds[var[v]]
         #       if cnt==0: 
         #          data = xr.zeros_like(ens_member_data)
         #          ens_min = ens_member_data.copy()
         #          ens_max = ens_member_data.copy()
         #       else:
         #          for t in range(len(ens_member_data)):
         #             ens_min[t] = np.min([ens_min[t].values,ens_member_data[t].values])
         #             ens_max[t] = np.max([ens_max[t].values,ens_member_data[t].values])
         #       data = ( data*cnt + ens_member_data ) / (cnt+1)
         #       cnt += 1
         # else:
         if True:
            tmp_ds = xr.open_dataset( tmp_file, use_cftime=True  )
            data = tmp_ds[var[v]]

      #-------------------------------------------------------------------------
      data = data/1e22
      # data = data - data.mean()
      #-------------------------------------------------------------------------
      # # Make time start at zero
      # print(); print(data['time'])
      # year_start = data['time.year'].values[0]
      # data['time'] = ( data['time'] - data['time'][0] ).astype('float') / 86400e9 / 365 + year_start
      # print(); print(data['time'])
      # exit()
      #-------------------------------------------------------------------------
      time_mean = data.mean(dim='time').values
      print('      Area Weighted Time Mean : '+hc.tcolor.GREEN+f'{time_mean:10.6f}'+hc.tcolor.ENDC)

      if print_stats: hc.print_stat(data,name='',stat='naxs',indent=' '*6,compact=True)

      data_list.append( data.values )
      time_list.append( data['time.year'].values )

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
   # set plot bounds here before loading obs
   res.trXMinF = np.min([np.nanmin(d) for d in time_list])
   res.trXMaxF = np.max([np.nanmax(d) for d in time_list])
   data_min = np.min([np.nanmin(d) for d in data_list])
   data_max = np.max([np.nanmax(d) for d in data_list])
   res.trYMinF = data_min - (data_max-data_min)*0.04
   res.trYMaxF = data_max + (data_max-data_min)*0.04

   #----------------------------------------------------------------------------
   # Create plot
   #----------------------------------------------------------------------------
   tres = copy.deepcopy(res)

   tres.xyCurveDrawOrder = 'PostDraw'

   for c in range(num_case):
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

      #------------------------------------------------
      # add ensemble spread
      #------------------------------------------------
      ens_case = 'v2.LR.historical'
      if case[c]==ens_case:

         n = len(time_list[c])
         ens_spread_data = np.zeros( 2*n+1 )
         ens_spread_data[n*0:n*1] = ens_min[:: 1].values
         ens_spread_data[n*1:n*2] = ens_max[::-1].values
         ens_spread_data[n*2]       = ens_min[0].values
         ens_spread_time = np.zeros( 2*n+1 )
         ens_spread_time[n*0:n*1] = time_list[c][:: 1]
         ens_spread_time[n*1:n*2] = time_list[c][::-1]
         ens_spread_time[n*2]       = time_list[c][0]
         eres = ngl.Resources()
         eres.gsFillColor = ens_spread_clr
         # eres.gsFillOpacityF = 0.1
         # eres.tfPolyDrawOrder = 'Draw'
         dum = ngl.add_polygon(wks, plot[ip], ens_spread_time, ens_spread_data, eres)

      #------------------------------------------------
      # add linear trend
      #------------------------------------------------
      if add_trend:
         px = time_list[c]
         py = data_list[c]
         # simple and fast method for regression coeff and intercept
         a = np.cov( px.flatten(), py.flatten() )[1,0] / np.var( px )
         b = np.mean(py) - a*np.mean(px)

         # print regression info
         # if c==0: print()
         # print(' '*4+f'linear regression a: {a}    b: {b}')
         # if c==(num_case-1): print()

         px_range = np.abs( np.max(px) - np.min(px) )
         lx = np.array([-1e2*px_range,1e2*px_range])

         lres.xyLineColor = clr[c]
         ngl.overlay( plot[ip], ngl.xy(wks, lx, lx*a+b , lres) )

   #----------------------------------------------------------------------------
   # Set strings at top of plot
   #----------------------------------------------------------------------------
   if overlay_cases:
      hs.set_subtitles(wks, plot[ip], '', var_str[v], '', font_height=0.01)
   else:
      hs.set_subtitles(wks, plot[ip], case_name[c], '', var_str[v], font_height=0.01)

   #----------------------------------------------------------------------------
   # Ad legend
   #----------------------------------------------------------------------------
   lgres = ngl.Resources()
   lgres.vpWidthF, lgres.vpHeightF  = 0.08, 0.08  
   lgres.lgLabelFontHeightF = 0.01
   lgres.lgLineThicknessF   = 4
   lgres.lgMonoLineColor    = False
   # lgres.lgMonoDashIndex    = True
   lgres.lgLineColors       = clr
   lgres.lgDashIndexes      = dsh
   lgres.lgLabelJust    = 'CenterLeft'
   # pid = ngl.legend_ndc(wks, len(name), name, 0.5, 0.4, lgres)  # 3x2
   # pid = ngl.legend_ndc(wks, len(name), name, 0.5, 0.1, lgres)  # 3x2
   # pid = ngl.legend_ndc(wks, len(name), name, 0.3, 0.5, lgres)  # 1x2

#-------------------------------------------------------------------------------
# Add legend
#-------------------------------------------------------------------------------

lgres = ngl.Resources()
lgres.vpWidthF           = 0.06
lgres.vpHeightF          = 0.08#*num_case
lgres.lgLabelFontHeightF = 0.008
lgres.lgLabelFont        = "courier"
lgres.lgMonoDashIndex    = False
lgres.lgLineLabelsOn     = False
lgres.lgLineThicknessF   = 30
lgres.lgLabelJust        = 'CenterLeft'
lgres.lgLineColors       = clr
lgres.lgDashIndexes      = dsh

indent = ' '*4
labels = case_name
for i in range(len(labels)): labels[i] = indent+labels[i] 

if 'v2.LR.historical' in case:
   labels.insert(0,indent+'E3SMv2 Ens Min/Max')
   lgres.lgLineColors.insert(0,ens_spread_clr)
   lgres.lgDashIndexes.insert(0,0)
   # xyLineOpacities

pid = ngl.legend_ndc(wks, len(labels), labels, 0.38, 0.63, lgres)

#---------------------------------------------------------------------------------------------------
# Finalize plot
#---------------------------------------------------------------------------------------------------
hc.printline()

layout = [len(plot),1]

ngl.panel(wks,plot[0:len(plot)],layout,hs.setres_panel())
ngl.end()

hc.trim_png(fig_file)

#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------

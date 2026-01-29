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

var,lev_list,var_str = [],[],[]
def add_var(var_name,lev=-1,n=None): 
   var.append(var_name); lev_list.append(lev)
   if n is None: n = var_name
   var_str.append(n)
#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------

### 2024 AQP Cess
# tmp_scratch,tmp_sub = '/gpfs/alpine2/atm146/proj-shared/hannah6/e3sm_scratch/','run'
tmp_scratch,tmp_sub = '/pscratch/sd/w/whannah/2024-AQP-CESS','archive/atm/hist'
add_case('E3SM.2024-AQP-CESS-00.FAQP-MMF1.GNUGPU.ne30pg2_ne30pg2.NN_32.SSTP_0K',    d=1,c='blue',     n='MMF ne30 +0K', p=tmp_scratch,s=tmp_sub)
add_case('E3SM.2024-AQP-CESS-00.FAQP-MMF1.GNUGPU.ne30pg2_ne30pg2.NN_32.SSTP_4K',    d=1,c='red',      n='MMF ne30 +4K', p=tmp_scratch,s=tmp_sub)
add_case('E3SM.2024-AQP-CESS-00.FAQP-MMF1.GNUGPU.ne45pg2_ne45pg2.NN_64.SSTP_0K',    d=1,c='blue',     n='MMF ne45 +0K', p=tmp_scratch,s=tmp_sub)
add_case('E3SM.2024-AQP-CESS-00.FAQP-MMF1.GNUGPU.ne45pg2_ne45pg2.NN_64.SSTP_4K',    d=1,c='red',      n='MMF ne45 +4K', p=tmp_scratch,s=tmp_sub)
add_case('E3SM.2024-AQP-CESS-00.FAQP-MMF1.GNUGPU.ne60pg2_ne60pg2.NN_128.SSTP_0K',   d=2,c='blue',     n='MMF ne60 +0K', p=tmp_scratch,s=tmp_sub)
add_case('E3SM.2024-AQP-CESS-00.FAQP-MMF1.GNUGPU.ne60pg2_ne60pg2.NN_128.SSTP_4K',   d=2,c='red',      n='MMF ne60 +4K', p=tmp_scratch,s=tmp_sub)
add_case('E3SM.2024-AQP-CESS-00.FAQP-MMF1.GNUGPU.ne90pg2_ne90pg2.NN_256.SSTP_0K',   d=2,c='blue',     n='MMF ne90 +0K', p=tmp_scratch,s=tmp_sub)
add_case('E3SM.2024-AQP-CESS-00.FAQP-MMF1.GNUGPU.ne90pg2_ne90pg2.NN_256.SSTP_4K',   d=2,c='red',      n='MMF ne90 +4K', p=tmp_scratch,s=tmp_sub)
add_case('E3SM.2024-AQP-CESS-00.FAQP-MMF1.GNUGPU.ne120pg2_ne120pg2.NN_512.SSTP_0K', d=3,c='blue',     n='MMF ne120 +0K',p=tmp_scratch,s=tmp_sub)
add_case('E3SM.2024-AQP-CESS-00.FAQP-MMF1.GNUGPU.ne120pg2_ne120pg2.NN_512.SSTP_4K', d=3,c='red',      n='MMF ne120 +4K',p=tmp_scratch,s=tmp_sub)

add_case('E3SM.2024-AQP-CESS-00.FAQP.GNUCPU.ne30pg2_ne30pg2.NN_32.SSTP_0K',         d=0,c='cyan',     n='EAM ne30 +0K', p=tmp_scratch,s=tmp_sub)
add_case('E3SM.2024-AQP-CESS-00.FAQP.GNUCPU.ne30pg2_ne30pg2.NN_32.SSTP_4K',         d=0,c='magenta',  n='EAM ne30 +4K', p=tmp_scratch,s=tmp_sub)
add_case('E3SM.2024-AQP-CESS-00.FAQP.GNUCPU.ne45pg2_ne45pg2.NN_64.SSTP_0K',         d=1,c='cyan',     n='EAM ne45 +0K', p=tmp_scratch,s=tmp_sub)
add_case('E3SM.2024-AQP-CESS-00.FAQP.GNUCPU.ne45pg2_ne45pg2.NN_64.SSTP_4K',         d=1,c='magenta',  n='EAM ne45 +4K', p=tmp_scratch,s=tmp_sub)
add_case('E3SM.2024-AQP-CESS-00.FAQP.GNUCPU.ne60pg2_ne60pg2.NN_128.SSTP_0K',        d=2,c='cyan',     n='EAM ne60 +0K', p=tmp_scratch,s=tmp_sub)
add_case('E3SM.2024-AQP-CESS-00.FAQP.GNUCPU.ne60pg2_ne60pg2.NN_128.SSTP_4K',        d=2,c='magenta',  n='EAM ne60 +4K', p=tmp_scratch,s=tmp_sub)
add_case('E3SM.2024-AQP-CESS-00.FAQP.GNUCPU.ne90pg2_ne90pg2.NN_256.SSTP_0K',        d=2,c='cyan',     n='EAM ne90 +0K', p=tmp_scratch,s=tmp_sub)
add_case('E3SM.2024-AQP-CESS-00.FAQP.GNUCPU.ne90pg2_ne90pg2.NN_256.SSTP_4K',        d=2,c='magenta',  n='EAM ne90 +4K', p=tmp_scratch,s=tmp_sub)
add_case('E3SM.2024-AQP-CESS-00.FAQP.GNUCPU.ne120pg2_ne120pg2.NN_512.SSTP_0K',      d=3,c='cyan',     n='EAM ne120 +0K',p=tmp_scratch,s=tmp_sub)
add_case('E3SM.2024-AQP-CESS-00.FAQP.GNUCPU.ne120pg2_ne120pg2.NN_512.SSTP_4K',      d=3,c='magenta',  n='EAM ne120 +4K',p=tmp_scratch,s=tmp_sub)

#---------------------------------------------------------------------------------------------------

# add_var('RESTOM')
add_var('NET_TOM_RAD')

add_var('PRECT')

# add_var('PS')
# add_var('TMQ')
# add_var('TCO')
# add_var('PRECT')
# add_var('SHFLX')
# add_var('LHFLX')
# add_var('P-E')

# add_var('TGCLDLWP')
# add_var('TGCLDIWP')

# add_var('TBOT'); 
# add_var('QBOT'); 
# add_var('SHFLX')
# add_var('LHFLX')

# add_var('UBOT'); 
# add_var('VBOT'); 

#---------------------------------------------------------------------------------------------------

htype,years,months,first_file,num_files = 'h0',[],[],0,0
# htype,years,months,first_file,num_files = 'h1',[],[],0,0

# lat1,lat2 = -30,30

fig_file,fig_type = 'figs/FXX-timeseries-global','png'
tmp_file_head = 'data/timeseries-global'

write_file    = False
print_stats   = True

recalculate = True

convert_to_daily_mean  = False

add_trend = False

num_plot_col  = 2

# year_start = 1950 + first_file/12
# year_start = 1985 + first_file/12
year_start = 00 + first_file/12

#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------
num_var = len(var)
num_case = len(case)

if 'lev' not in vars(): lev = np.array([0])

#---------------------------------------------------------------------------------------------------
# Set up plot resources
#---------------------------------------------------------------------------------------------------

plot = [None]*(num_var)

wkres = ngl.Resources()
npix=1024; wkres.wkWidth,wkres.wkHeight=npix,npix
wks = ngl.open_wks(fig_type,fig_file,wkres)

res = hs.res_xy()
res.vpHeightF = 0.5
# res.vpHeightF = 0.2
res.tmYLLabelFontHeightF         = 0.015
res.tmXBLabelFontHeightF         = 0.015
res.tiXAxisFontHeightF           = 0.015
res.tiYAxisFontHeightF           = 0.015
res.xyLineThicknessF = 5

lres = hs.res_xy()
lres.xyDashPattern    = 1
lres.xyLineThicknessF = 1
lres.xyLineColor      = "black"

#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------
for v in range(num_var):
   
   hc.printline()
   print('  var: '+hc.tcolor.MAGENTA+var[v]+hc.tcolor.ENDC)

   if 'lev_list' in locals(): lev = lev_list[v]

   time_list,data_list = [],[]
   for c in range(num_case):

      tmp_file = f'{tmp_file_head}.{case[c]}.{var[v]}.nc'

      print('    case: '+hc.tcolor.CYAN+case[c]+hc.tcolor.ENDC)


      if recalculate:

         data_dir_tmp,data_sub_tmp = None, None
         if case_dir[c] is not None: data_dir_tmp = case_dir[c]
         if case_sub[c] is not None: data_sub_tmp = case_sub[c]
         case_obj = he.Case( name=case[c], data_dir=data_dir_tmp, data_sub=data_sub_tmp  )
         
         if 'lat1' in vars() : case_obj.lat1,case_obj.lat2 = lat1,lat2
         if 'lon1' in vars() : case_obj.lon1,case_obj.lon2 = lon1,lon2

         #----------------------------------------------------------------------
         # lat = case_obj.load_data('lat',  htype=htype,num_files=1)
         # lon = case_obj.load_data('lon',  htype=htype,num_files=1)
         area = case_obj.load_data('area',htype=htype,num_files=1).astype(np.double)
         #----------------------------------------------------------------------
         tvar = var[v]
         if tvar=='RESTOM':
            FSNT = case_obj.load_data('FSNT', htype=htype, first_file=first_file, num_files=num_files, lev=lev)
            FLNT = case_obj.load_data('FLNT', htype=htype, first_file=first_file, num_files=num_files, lev=lev)
            data = FSNT - FLNT
         else:
            data = case_obj.load_data(tvar, htype=htype, first_file=first_file, num_files=num_files, lev=lev)
         #----------------------------------------------------------------------
         # time = data.time - data.time[0]

         # time_bnds = case_obj.load_data('time_bnds',htype=htype,first_file=first_file,num_files=num_files)
         # time = time_bnds.isel(nbnd=0)

         data = ( (data*area).sum(dim='ncol') / area.sum(dim='ncol') )
         
         # if print_stats: hc.print_stat(data,name='',stat='naxs',indent=' '*6,compact=True)

         # Convert to daily mean
         if htype in ['h1','h2'] and convert_to_daily_mean: 
            data = data.resample(time='D').mean(dim='time')

         # if np.all( lev < 0 ) and 'lev' in data.coords : print(f'    lev value: {data.lev.values}')

         # Get rid of lev dimension
         if 'lev' in data.dims : data = data.isel(lev=0)
         #----------------------------------------------------------------------
         # write to file
         tmp_ds = xr.Dataset( coords=data.coords )
         tmp_ds[var[v]] = data
         tmp_ds.to_netcdf(path=tmp_file,mode='w')
      else:
         tmp_ds = xr.open_dataset( tmp_file, use_cftime=True  )
         data = tmp_ds[var[v]]
      #-------------------------------------------------------------------------
      # Make time start at zero      
      data['time'] = ( data['time'] - data['time'][0] ).astype('float') / 86400e9
      
      #-------------------------------------------------------------------------
      time_mean = data.mean(dim='time').values
      # print('      Area Weighted Time Mean : '+hc.tcolor.GREEN+f'{time_mean:10.6f}'+hc.tcolor.ENDC)
      print('      Area Weighted Time Mean : '+hc.tcolor.GREEN+f'{time_mean}'+hc.tcolor.ENDC)

      if print_stats: hc.print_stat(data,name='',stat='naxs',indent=' '*6,compact=True)

      data_list.append( data.values )
      time_list.append( data['time'].values )

   #----------------------------------------------------------------------------
   # Create plot
   #----------------------------------------------------------------------------
   tres = copy.deepcopy(res)
   # tres.tiXAxisString = 'Time [years]'
   tres.tiXAxisString = 'Time [days]'

   # reset start year for all data
   if 'year_start' in locals():
      for t in range(len(time_list)):
         time_list[t] = time_list[t] + year_start

   ### Make sure plot bounds are consistent
   # tres.trYMinF = 14.4
   tres.trYMinF = np.min([np.nanmin(d) for d in data_list])
   tres.trYMaxF = np.max([np.nanmax(d) for d in data_list])
   tres.trXMinF = np.min([np.nanmin(d) for d in time_list])
   tres.trXMaxF = np.max([np.nanmax(d) for d in time_list])

   # if var[v]=='PRECT': tres.trYMinF,tres.trYMaxF  = 2,5

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

   #----------------------------------------------------------------------------
   # tres.xyExplicitLegendLabels = case_name
   # tres.pmLegendDisplayMode    = 'Always'
   # tres.pmLegendOrthogonalPosF =-1.0 # upper-right
   # tres.pmLegendParallelPosF   = 0.8
   # tres.pmLegendWidthF         = 0.12
   # tres.pmLegendHeightF        = 0.08
   # tres.lgLabelFontHeightF     = 0.02
   # tres.lgLineThicknessF       = 4

   # tres.xyLineColors   = clr
   # tres.xyDashPatterns = dsh

   # ip = v
   # xx = np.array(time_list)
   # yy = np.array(data_list)
   # for c in range(tmp_num_case):

   # plot[ip] = ngl.xy(wks, xx, yy, tres)

   # #----------------------------------------------------------------------------
   # # Add legend
   # lgres = ngl.Resources()
   # lgres.vpWidthF, lgres.vpHeightF  = 0.08, 0.08  
   # lgres.lgLabelFontHeightF = 0.015
   # lgres.lgLineThicknessF   = 4
   # lgres.lgLineColors       = clr
   # lgres.lgDashIndexes      = dsh
   # lgres.lgLabelJust    = 'CenterLeft'
   # pid = ngl.legend_ndc(wks, len(case_name), case_name, 0.6, 0.8, lgres)

   #----------------------------------------------------------------------------
   for c in range(tmp_num_case):
      ip = v
      
      tres.xyLineColor   = clr[c]
      tres.xyDashPattern = dsh[c]
      tplot = ngl.xy(wks, time_list[c], data_list[c], tres)
      
      if c==0: 
         plot[ip] = tplot
      else:
         ngl.overlay(plot[ip],tplot)
      
   #----------------------------------------------------------------------------
   # add linear trend
   if add_trend:
      for c in range(tmp_num_case):
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
   ctr_str = ''
   # lat_chk,lon_chk = 'lat1' in locals(), 'lon1' in locals()
   # if not lat_chk and not lon_chk : 
   #    ctr_str = 'Global Mean'
   # else:
   #    if lat_chk:      ctr_str += f' {lat1}:{lat2}N '
   #    if lon_chk:      ctr_str += f' {lon1}:{lon2}E '
   
   hs.set_subtitles(wks, plot[ip], var_str[v], ctr_str, '', font_height=0.02)
   
   #----------------------------------------------------------------------------
   # Add legend
   #----------------------------------------------------------------------------
   # lgres = ngl.Resources()
   # lgres.vpWidthF, lgres.vpHeightF  = 0.08, 0.04*num_case
   # lgres.lgLabelFontHeightF = 0.015
   # lgres.lgLineThicknessF   = 4
   # lgres.lgLineColors       = clr
   # lgres.lgDashIndexes      = dsh
   # lgres.lgLabelJust    = 'CenterLeft'
   # pid = ngl.legend_ndc(wks, len(case_name), case_name, 0.65, 0.85, lgres)
   # pid = ngl.legend_ndc(wks, len(name), name, 0.5, 0.4, lgres)  # 3x2
   # pid = ngl.legend_ndc(wks, len(name), name, 0.5, 0.1, lgres)  # 3x2
   # pid = ngl.legend_ndc(wks, len(name), name, 0.3, 0.5, lgres)  # 1x2

#---------------------------------------------------------------------------------------------------
# Finalize plot
#---------------------------------------------------------------------------------------------------
hc.printline()

layout = [int(np.ceil(len(plot)/float(num_plot_col))),num_plot_col]

ngl.panel(wks,plot[0:len(plot)],layout,hs.setres_panel())
ngl.end()

hc.trim_png(fig_file)

#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------

import os, copy, string, ngl, xarray as xr, numpy as np, glob, random
import hapy_common as hc, hapy_E3SM   as he, hapy_setres as hs
#---------------------------------------------------------------------------------------------------
case_name,case,case_dir,case_sub,clr,dsh,mrk = [],[],[],[],[],[],[]
def add_case(case_in,n=None,p=None,s=None,d=0,c='black',m=0):
   global name,case,case_dir,case_sub,clr,dsh,mrk
   if n is None:
      tmp_name = ''
   else:
      tmp_name = n
   case.append(case_in); case_name.append(tmp_name)
   case_dir.append(p); case_sub.append(s)
   dsh.append(d) ; clr.append(c) ; mrk.append(m)
#---------------------------------------------------------------------------------------------------
var,lev_list,mask_flag,var_str = [],[],[],[]
def add_var(var_name,lev=-1,mask=None,vstr=None): 
   var.append(var_name); lev_list.append(lev),mask_flag.append(mask)
   if vstr is None: vstr = var_name
   var_str.append(vstr)
#---------------------------------------------------------------------------------------------------
tmp_path_co2_mmf = '/global/cfs/cdirs/m3312/whannah/2023-CPL'
tmp_sub = 'archive/atm/hist'
# add_case('E3SM.2023-CO2-TEST-01.GNUGPU.ne30pg2_EC30to60E2r2.WCYCL1850-MMF1.1xCO2',n='1xCO2',c='blue' ,p=tmp_path_co2_mmf,s=tmp_sub)
# add_case('E3SM.2023-CO2-TEST-01.GNUGPU.ne30pg2_EC30to60E2r2.WCYCL1850-MMF1.2xCO2',n='2xCO2',c='green',p=tmp_path_co2_mmf,s=tmp_sub)
# add_case('E3SM.2023-CO2-TEST-01.GNUGPU.ne30pg2_EC30to60E2r2.WCYCL1850-MMF1.4xCO2',n='4xCO2',c='red'  ,p=tmp_path_co2_mmf,s=tmp_sub)
add_case('E3SM.2023-CO2-TEST-01.GNUGPU.ne30pg2_EC30to60E2r2.WCYCL1850-MMF1.1xCO2',n='1xCO2',c='royalblue',   p=tmp_path_co2_mmf,s=tmp_sub)
add_case('E3SM.2023-CO2-TEST-01.GNUGPU.ne30pg2_EC30to60E2r2.WCYCL1850-MMF1.2xCO2',n='2xCO2',c='springgreen3',p=tmp_path_co2_mmf,s=tmp_sub)
add_case('E3SM.2023-CO2-TEST-01.GNUGPU.ne30pg2_EC30to60E2r2.WCYCL1850-MMF1.4xCO2',n='4xCO2',c='tomato',      p=tmp_path_co2_mmf,s=tmp_sub)

#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------

htype,yr1,yr2 = 'ha', 0, 120
fig_file,fig_type = 'figs/F11-eff-climate-sensitivity','png'
tmp_file_head = 'data/eff-climate-sensitivity'

# htype,yr1,yr2 = 'h0', 0, 120*12
# fig_file,fig_type = 'figs/FXX-eff-climate-sensitivity-h0','png'
# tmp_file_head = 'data/eff-climate-sensitivity-h0'

#---------------------------------------------------------------------------------------------------
print_stats   = True
overlay_cases = True

recalculate = False

num_plot_col  = 2

#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------
num_case = len(case)

#---------------------------------------------------------------------------------------------------
# Set up plot resources
#---------------------------------------------------------------------------------------------------
plot = [None]*2

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

lres = hs.res_xy()
lres.xyDashPattern    = 0
lres.xyLineThicknessF = 2
lres.xyLineColor      = 'black'


reg_label = ['2x / 1x','4x / 2x','4x / 1x']
reg_clr = ['orange','cyan','magenta']

lgres = ngl.Resources()
# lgres.vpWidthF, lgres.vpHeightF  = 0.075, 0.085  
lgres.vpWidthF, lgres.vpHeightF  = 0.04, 0.085  
lgres.lgLabelFontHeightF = 0.01
lgres.lgLineThicknessF   = 30#12
lgres.lgMonoLineColor    = False
lgres.lgMonoDashIndex    = True
lgres.lgDashIndex        = 0
lgres.lgLineColors       = reg_clr
lgres.lgLabelJust        = 'CenterLeft'

#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------
hc.printline()
T_list = []
N_list = []
time_list = []
for c in range(num_case):
   tmp_file = f'{tmp_file_head}.{case[c]}.nc'
   print(' '*4+f'case: {hc.tclr.CYAN}{case[c]}{hc.tclr.END}  =>  {tmp_file}')
   if recalculate:
      #-------------------------------------------------------------------------
      file_path = f'{case_dir[c]}/{case[c]}/{case_sub[c]}/*.eam.{htype}.*'
      file_list = sorted(glob.glob(file_path))

      ds = xr.open_mfdataset( file_list )
      #-------------------------------------------------------------------------
      ds = ds.where( ds['time.year']>=yr1, drop=True)
      ds = ds.where( ds['time.year']<=yr2, drop=True)
      #-------------------------------------------------------------------------
      # read the data
      area = ds['area'].isel(time=0)
      T    = ds['TS']
      FSNT = ds['FSNT']
      FLNT = ds['FLNT']
      N = FSNT - FLNT
      #-------------------------------------------------------------------------
      # area weighted global mean
      T = (T*area).sum(dim='ncol') / area.sum(dim='ncol')
      N = (N*area).sum(dim='ncol') / area.sum(dim='ncol')
      #-------------------------------------------------------------------------
      T.load() ; N.load()
      tmp_ds = xr.Dataset( coords=T.coords )
      tmp_ds['T'] = T
      tmp_ds['N'] = N
      tmp_ds['yr1'] = yr1
      tmp_ds['yr2'] = yr2
      tmp_ds.to_netcdf(path=tmp_file,mode='w')
   #----------------------------------------------------------------------------
   else:
      tmp_ds = xr.open_dataset( tmp_file, use_cftime=True  )
      T = tmp_ds['T']
      N = tmp_ds['N']

   # T = T.isel(time=slice(20,150))
   # N = N.isel(time=slice(20,150))
   #----------------------------------------------------------------------------
   if print_stats: 
      hc.print_stat(T,name='T (global mean)',stat='naxs',indent=' '*6,compact=True)
      hc.print_stat(N,name='N (global mean)',stat='naxs',indent=' '*6,compact=True)
   #-------------------------------------------------------------------------
   # Create time index for time series panels
   # time = ( T['time'] - T['time'][0] ).astype('float') / 86400e9 / 365 
   time = T['time.year']
   #----------------------------------------------------------------------------
   T_list.append(T.values)
   N_list.append(N.values)
   time_list.append(time.values)

#---------------------------------------------------------------------------------------------------
# Calculate differences for ECS via Gregory regression
#---------------------------------------------------------------------------------------------------
dT_list1,dN_list1 = [],[]
# yr1_baseline = 20
for c in range(num_case):
   if c==0: ct,cb = 1,0
   if c==1: ct,cb = 2,1
   if c==2: ct,cb = 2,0
   dT_list1.append( T_list[ct] - T_list[cb] )
   dN_list1.append( N_list[ct] - N_list[cb] )

# if print_stats:
#    for c in range(num_case):
#       print()
#       hc.print_stat(dT_list1[c],name=f'dT {reg_label[c]}',indent=' '*4,compact=True) # stat='naxs'
#       hc.print_stat(dN_list1[c],name=f'dN {reg_label[c]}',indent=' '*4,compact=True) # stat='naxs'
# print()
#----------------------------------------------------------------------------
# Use alternative bootstrap sampling to estimate ECS
#----------------------------------------------------------------------------

# nchoice = 10000

# dT_list2,dN_list2 = [],[]

# for c in range(num_case):
#    if c==0: ct,cb = 1,0
#    if c==1: ct,cb = 2,1
#    if c==2: ct,cb = 2,0
#    i_t = random.choices(np.arange( 0,120,dtype=int), k=nchoice)
#    i_b = random.choices(np.arange(60,120,dtype=int), k=nchoice)
#    dT_list2.append( T_list[ct][i_t] - T_list[cb][i_b] )
#    dN_list2.append( N_list[ct][i_t] - N_list[cb][i_b] )

#----------------------------------------------------------------------------
# Create time series plots
#----------------------------------------------------------------------------
tres = copy.deepcopy(res)
tres.tiXAxisString = 'Time [years]'

tres.xyLineThicknessF = 20
tres.trXMinF = np.min([np.nanmin(d) for d in time_list])
tres.trXMaxF = np.max([np.nanmax(d) for d in time_list])

for n in range(1):
   ip = n
   if n==0: 
      data_list = T_list ; vstr = 'Global Mean Surface Temperature'
      xtick = np.arange(286,296+2,2,dtype=int)
      tres.tmYLMode,tres.tmYLValues,tres.tmYLLabels   = "Explicit",xtick,xtick
      tres.trYMinF = 286
      tres.trYMaxF = 296
      tres.tiYAxisString = f'[K]'
   if n==1: 
      data_list = N_list ; vstr = 'Global Mean TOM Radiative Imbalance'
      xtick = np.arange(-2,6+2,2,dtype=int)
      tres.tmYLMode,tres.tmYLValues,tres.tmYLLabels   = "Explicit",xtick,xtick
      tres.trYMinF = -2
      tres.trYMaxF = 6
      tres.tiYAxisString = f'[W/m2]'
   # tres.trYMinF = np.min([np.nanmin(d) for d in data_list])
   # tres.trYMaxF = np.max([np.nanmax(d) for d in data_list])
   for c in range(num_case):
      tres.xyLineColor   = clr[c]
      tres.xyDashPattern = dsh[c]
      tplot = ngl.xy(wks, time_list[c], data_list[c], tres)
      if c==0:
         plot[ip] = tplot
         ngl.overlay( plot[ip], ngl.xy(wks, np.array([-1e8,1e8]), np.array([0,0]) , lres) )
         hs.set_subtitles(wks, plot[ip], '', vstr, '', font_height=0.015)
      else:
         ngl.overlay(plot[ip],tplot)

#---------------------------------------------------------------------------------------------------
# Create ECS / Gregory plot
#---------------------------------------------------------------------------------------------------

for n in range(1):
   ip = 1+n

   if n==0: dT_list,dN_list = dT_list1,dN_list1
   if n==1: dT_list,dN_list = dT_list2,dN_list2

   # calculate regression coefficients
   a_list,xi_list,yi_list = [],[],[]
   for c in range(num_case):
      px,py = dT_list[c],dN_list[c]
      if c==num_case-1: px,py = px/2,py/2
      # simple and fast method for regression coeff and intercept
      a = np.cov( px.flatten(), py.flatten() )[1,0] / np.var( px )
      yi = np.mean(py) - a*np.mean(px)
      xi = -1*yi/a
      a_list.append(a)
      xi_list.append(xi)
      yi_list.append(yi)

   # use bootstrap sampling to estimate ECS uncertainty
   ngrp = 1000
   nsmp = 60
   print()
   min_xi_list = []
   max_xi_list = []
   for c in range(num_case):
      if c==0: ct,cb = 1,0
      if c==1: ct,cb = 2,1
      if c==2: ct,cb = 2,0
      bs_xi_list = []
      for g in range(ngrp):
         ii = random.choices(np.arange( 0,120,dtype=int), k=nsmp)
         px = T_list[ct][ii] - T_list[cb][ii]
         py = N_list[ct][ii] - N_list[cb][ii]
         if c==2: px,py = px/2,py/2
         a = np.cov( px.flatten(), py.flatten() )[1,0] / np.var( px )
         yi = np.mean(py) - a*np.mean(px)
         xi = -1*yi/a
         bs_xi_list.append(xi)
      bs_xi = np.array(bs_xi_list)
      # min_xi_list.append(np.min(bs_xi))
      # max_xi_list.append(np.max(bs_xi))
      print(f'  c: {c}  {reg_label[c]:10}  min: {np.min(bs_xi):8.4f}  mean: {np.average(bs_xi):8.4f}  max: {np.max(bs_xi):8.4f}  std: {np.std(bs_xi):8.4f}')
      min_xi_list.append(np.mean(bs_xi) - 2*np.std(bs_xi))
      max_xi_list.append(np.mean(bs_xi) + 2*np.std(bs_xi))
   print()
   # exit()

   ### Set consistent plot bounds based on regression lines
   dN_min = 0#np.min([np.min(dN) for dN in yi_list])
   dN_max = np.max([np.max(dN) for dN in yi_list])
   dT_min = 0#np.min([np.min(dT) for dT in xi_list])
   dT_max = np.max([np.max(dT) for dT in xi_list])

   dN_span = dN_max - dN_min
   dT_span = dT_max - dT_min

   tres = copy.deepcopy(res)

   tres.trYMinF = dN_min #- dN_span*0.1
   tres.trYMaxF = 4.2#dN_max + dN_span*0.1
   tres.trXMinF = dT_min #- dT_span*0.1
   tres.trXMaxF = 6.8#dT_max + dT_span*0.1

   delta_symbol = '~F33~D~F21~'
   tres.tiXAxisString = f'{delta_symbol}T~B~2x~N~ [K]'
   tres.tiYAxisString = f'{delta_symbol}N~B~2x~N~ [W/m2]'

   tres.xyMarkLineMode = 'Markers'
   tres.xyMarker = 16
   if n==0: tres.xyMarkerSizeF = 0.006
   if n==1: tres.xyMarkerSizeF = 0.001

   
   print()
   for c in range(num_case):
      # ip = 1

      px,py = dT_list[c],dN_list[c]
      if c==num_case-1: px,py = px/2,py/2

      tres.xyMarkerColor = reg_clr[c]

      tplot = ngl.xy(wks, px, py, tres)
      
      if c==0:
         plot[ip] = tplot

         if n==0: cstr = 'Gregory Regression'
         # if n==0: cstr = 'Traditional Gregory Regression'
         # if n==1: cstr = 'Bootstrapped Gregory Regression'
         hs.set_subtitles(wks, plot[ip], '', cstr, '', font_height=0.015)
      else:
         ngl.overlay(plot[ip],tplot)
      #----------------------------------------------------------------------------
      # add regression line
      a,xi,yi = a_list[c],xi_list[c],yi_list[c]
      
      print(' '*4+f'{reg_label[c]:30}   slope: {a:8.4f}  yi: {yi:8.4f}  xi: {xi:8.4f}')

      px_range = np.abs( np.max(px) - np.min(px) )
      lx = np.array([-1e2*px_range,1e2*px_range])
      
      lres.xyMarkLineMode = 'Lines'
      lres.xyDashPattern = 0
      lres.xyLineThicknessF = 12
      lres.xyLineColor = reg_clr[c]
      # if c==num_case-1: lres.xyDashPattern = 2
      ngl.overlay( plot[ip], ngl.xy(wks, lx, lx*a+yi , lres) )

      #----------------------------------------------------------------------------
      # add brackets to indicate std deviation of x-intercept

      lres.xyLineThicknessF = 20 # 12

      # lres.xyMarkLineMode = 'MarkLines'
      # lres.xyMarker = 16
      # lres.xyMarkerSizeF = 0.002
      # lres.xyMarkerColor = reg_clr[c]
      xx = np.array([min_xi_list[c],max_xi_list[c]])
      yy = np.array([0,0])
      ngl.overlay( plot[ip], ngl.xy(wks, xx, yy, lres) )


      dl = 0.18 # 0.1

      xx = np.array([min_xi_list[c],min_xi_list[c]])
      # yy = np.array([0,dl])
      yy = np.array([-1*dl,dl])
      ngl.overlay( plot[ip], ngl.xy(wks, xx, yy, lres) )

      xx = np.array([max_xi_list[c],max_xi_list[c]])
      # yy = np.array([0,dl])
      ngl.overlay( plot[ip], ngl.xy(wks, xx, yy, lres) )

      #----------------------------------------------------------------------------
      # alt brackets to indicate std deviation of x-intercept using NDC coords



      #----------------------------------------------------------------------------
      # alt brackets to indicate std deviation of x-intercept using NDC coords

      # gsres = ngl.Resources()
      # gsres.gsLineColor = reg_clr[c]
      # gsres.gsLineDashPattern = 0
      # gsres.gsLineThicknessF = 4

      # min_xi_ndcx, min_xi_ndcy = ngl.datatondc(plot[ip], min_xi_list[c], 0 )
      # max_xi_ndcx, max_xi_ndcy = ngl.datatondc(plot[ip], max_xi_list[c], 0 )

      # xx = np.array([min_xi_ndcx,min_xi_ndcx])
      # yy = np.array([min_xi_ndcy-1*dl,min_xi_ndcy+dl])

      # dl = 0.1

      # ngl.polyline_ndc(wks, xx, yy, gsres)

      #----------------------------------------------------------------------------
      if c==0:
         # if 'lgd_label' in locals(): del lgd_label
         lgd_label = copy.deepcopy(reg_label)
         for i in range(len(lgd_label)): 
            lgd_label[i] = f'  {reg_label[i]}    ECS: {xi_list[i]:4.2f} K'
         # if n==0: lgx, lgy = 0.25, 0.45
         # if n==1: lgx, lgy = 0.75,  0.45
         # pid = ngl.legend_ndc(wks, len(lgd_label), lgd_label, lgx, lgy, lgres)
         pid = ngl.legend_ndc(wks, len(lgd_label), lgd_label, 0.78, 0.62, lgres)

#----------------------------------------------------------------------------
# Add legend
#----------------------------------------------------------------------------

# lgd_label = reg_label
# for i in range(len(lgd_label)): lgd_label[i] = f'  {reg_label[i]}    ECS: {xi_list[i]:4.2f} K'
# pid = ngl.legend_ndc(wks, len(lgd_label), lgd_label, 0.76, 0.62, lgres)


# lgres = ngl.Resources()
# lgres.vpWidthF           = 0.03
# lgres.vpHeightF          = 0.06
# lgres.lgLabelFontHeightF = 0.008
# # lgres.lgLabelFont        = "courier"
# lgres.lgMonoDashIndex    = False
# lgres.lgLineLabelsOn     = False
# lgres.lgLineThicknessF   = 20
# lgres.lgLabelJust        = 'CenterLeft'
lgres.lgLineColors       = clr
lgres.lgDashIndexes      = dsh

indent = ' '*2
labels = case_name
for i in range(len(labels)): labels[i] = indent+labels[i] 

pid = ngl.legend_ndc(wks, len(labels), labels, 0.3, 0.58, lgres)

#---------------------------------------------------------------------------------------------------
# Finalize plot
#---------------------------------------------------------------------------------------------------
hc.printline()

if 'num_plot_col' in locals():
   layout = [int(np.ceil(len(plot)/float(num_plot_col))),num_plot_col]
else:
   layout = [1,len(plot)]

pnl_res = hs.setres_panel()
pnl_res.nglPanelFigureStrings            = list(string.ascii_lowercase)
pnl_res.nglPanelFigureStringsJust        = "TopLeft"
pnl_res.nglPanelFigureStringsFontHeightF = 0.01
ngl.panel(wks,plot,layout,pnl_res)
ngl.end()
hc.trim_png(fig_file)

if 'legend_file' in locals(): hc.trim_png(legend_file)

#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------

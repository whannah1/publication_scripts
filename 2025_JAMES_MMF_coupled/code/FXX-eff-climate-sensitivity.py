import os, copy, ngl, xarray as xr, numpy as np
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
add_case('E3SM.2023-CO2-TEST-01.GNUGPU.ne30pg2_EC30to60E2r2.WCYCL1850-MMF1.1xCO2',n='1xCO2',c='blue' ,p=tmp_path_co2_mmf,s=tmp_sub)
add_case('E3SM.2023-CO2-TEST-01.GNUGPU.ne30pg2_EC30to60E2r2.WCYCL1850-MMF1.2xCO2',n='2xCO2',c='green',p=tmp_path_co2_mmf,s=tmp_sub)
add_case('E3SM.2023-CO2-TEST-01.GNUGPU.ne30pg2_EC30to60E2r2.WCYCL1850-MMF1.4xCO2',n='4xCO2',c='red'  ,p=tmp_path_co2_mmf,s=tmp_sub)
#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------

# htype='ha';first_file=0;num_files=120-first_file
# htype,first_file,num_files = 'ha',80,40
htype,first_file,num_files = 'ha',100,10

fig_file,fig_type = 'figs/FXX-eff-climate-sensitivity','png'
tmp_file_head = 'data/eff-climate-sensitivity'

#---------------------------------------------------------------------------------------------------
print_stats   = True
overlay_cases = True

recalculate = True

num_plot_col  = 2

#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------
num_case = len(case)

if 'clr' not in vars() or clr==[]: clr = ['black']*num_case
if 'dsh' not in vars() or dsh==[]: dsh = np.zeros(num_case)

#-------------------------------------------------------------------------------
# plot legend in separate file
#-------------------------------------------------------------------------------
# if num_case>1:
#    legend_file = fig_file+'.legend'
#    wkres = ngl.Resources() #; npix = 1024 ; wkres.wkWidth,wkres.wkHeight=npix,npix
#    lgd_wks = ngl.open_wks('png',legend_file,wkres)
#    lgres = ngl.Resources()
#    lgres.vpWidthF           = 0.06
#    lgres.vpHeightF          = 0.06#*num_case
#    lgres.lgLabelFontHeightF = 0.008
#    lgres.lgLabelFont        = "courier"
#    lgres.lgMonoDashIndex    = False
#    lgres.lgLineLabelsOn     = False
#    lgres.lgLineThicknessF   = 8#16
#    lgres.lgLabelJust        = 'CenterLeft'
#    lgres.lgLineColors       = clr
#    lgres.lgDashIndexes      = dsh

#    indent = ' '*4
#    labels = case_name
#    for i in range(len(labels)): labels[i] = indent+labels[i] 

#    pid = ngl.legend_ndc(lgd_wks, len(labels), labels, 0.5, 0.65, lgres)

#    ngl.frame(lgd_wks)
#    hc.trim_png(legend_file)

#---------------------------------------------------------------------------------------------------
# Set up plot resources
#---------------------------------------------------------------------------------------------------
plot = [None]*3

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
      # set up the case object
      data_dir_tmp,data_sub_tmp = None, None
      if case_dir[c] is not None: data_dir_tmp = case_dir[c]
      if case_sub[c] is not None: data_sub_tmp = case_sub[c]
      case_obj = he.Case( name=case[c], data_dir=data_dir_tmp, data_sub=data_sub_tmp  )
      #-------------------------------------------------------------------------
      # read the data
      area = case_obj.load_data('area', htype=htype, num_files=1).astype(np.double)
      T    = case_obj.load_data('TS',   htype=htype, num_files=num_files, first_file=first_file)
      FSNT = case_obj.load_data('FSNT', htype=htype, num_files=num_files, first_file=first_file)
      FLNT = case_obj.load_data('FLNT', htype=htype, num_files=num_files, first_file=first_file)
      N = FSNT - FLNT
      #-------------------------------------------------------------------------
      # perform area weighted global mean
      area_sum = area.sum(dim='ncol')
      T = (T*area).sum(dim='ncol') / area_sum
      N = (N*area).sum(dim='ncol') / area_sum
      #-------------------------------------------------------------------------
      T.load() ; N.load()
      tmp_ds = xr.Dataset( coords=T.coords )
      tmp_ds['T'] = T
      tmp_ds['N'] = N
      tmp_ds['first_file'] = first_file
      tmp_ds['num_files']  = num_files
      tmp_ds.to_netcdf(path=tmp_file,mode='w')
   #----------------------------------------------------------------------------
   else:
      tmp_ds = xr.open_dataset( tmp_file, use_cftime=True  )
      T = tmp_ds['T']
      N = tmp_ds['N']
   #----------------------------------------------------------------------------
   # if print_stats: 
   #    hc.print_stat(T,name='T (global mean)',stat='naxs',indent=' '*6,compact=True)
   #    hc.print_stat(F,name='F (global mean)',stat='naxs',indent=' '*6,compact=True)
   #-------------------------------------------------------------------------
   # Create time index for time series panels
   time = ( T['time'] - T['time'][0] ).astype('float') / 86400e9 / 365 
   #----------------------------------------------------------------------------
   T_list.append(T.values)
   N_list.append(N.values)
   time_list.append(time.values)

#---------------------------------------------------------------------------------------------------
# Calculate differences for ECS via Gregory regression
#---------------------------------------------------------------------------------------------------
dT_list,dN_list = [],[]
# yr1_baseline = 20
for c in range(num_case):
   if c==0: ct,cb = 1,0
   if c==1: ct,cb = 2,1
   if c==2: ct,cb = 2,0
   dT_list.append( T_list[ct] - T_list[cb] )
   dN_list.append( N_list[ct] - N_list[cb] )

reg_label = ['2x / 1x','4x / 2x','4x / 1x (halved)']

if print_stats:
   for c in range(num_case):
      print()
      hc.print_stat(dT_list[c],name=f'dT {reg_label[c]}',indent=' '*4,compact=True) # stat='naxs'
      hc.print_stat(dN_list[c],name=f'dN {reg_label[c]}',indent=' '*4,compact=True) # stat='naxs'
print()
#----------------------------------------------------------------------------
# Create time series plots
#----------------------------------------------------------------------------
tres = copy.deepcopy(res)
tres.tiXAxisString = 'Time [years]'

tres.xyLineThicknessF = 20
tres.trXMinF = np.min([np.nanmin(d) for d in time_list])
tres.trXMaxF = np.max([np.nanmax(d) for d in time_list])

for n in range(2):
   ip = n
   if n==0: data_list = T_list ; vstr = 'Global Mean Tsfc Anomaly'
   if n==1: data_list = N_list ; vstr = 'Global Mean TOM Radiative Imbalance'
   tres.trYMinF = np.min([np.nanmin(d) for d in data_list])
   tres.trYMaxF = np.max([np.nanmax(d) for d in data_list])
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

### Set consistent plot bounds based on regression lines
dN_min = 0#np.min([np.min(dN) for dN in yi_list])
dN_max = np.max([np.max(dN) for dN in yi_list])
dT_min = 0#np.min([np.min(dT) for dT in xi_list])
dT_max = np.max([np.max(dT) for dT in xi_list])

dN_span = dN_max - dN_min
dT_span = dT_max - dT_min

tres = copy.deepcopy(res)

tres.trYMinF = dN_min #- dN_span*0.1
tres.trYMaxF = dN_max + dN_span*0.1
tres.trXMinF = dT_min #- dT_span*0.1
tres.trXMaxF = dT_max + dT_span*0.1

delta_symbol = '~F33~D~F21~'
tres.tiXAxisString = f'{delta_symbol}T'
tres.tiYAxisString = f'{delta_symbol}N'

tres.xyMarkLineMode = 'Markers'
tres.xyMarkerSizeF = 0.005
tres.xyMarker = 16

reg_clr = ['orange','cyan','magenta']
for c in range(num_case):
   ip = 2

   px,py = dT_list[c],dN_list[c]
   if c==num_case-1: px,py = px/2,py/2

   tres.xyMarkerColor = reg_clr[c]

   tplot = ngl.xy(wks, px, py, tres)
   
   if c==0:
      plot[ip] = tplot
      hs.set_subtitles(wks, plot[ip], '', '', '', font_height=0.015)
   else:
      ngl.overlay(plot[ip],tplot)
   #----------------------------------------------------------------------------
   # add regression line
   a,xi,yi = a_list[c],xi_list[c],yi_list[c]
   
   print(' '*4+f'{reg_label[c]}   slope: {a:8.4f}  yi: {yi:8.4f}  xi: {xi:8.4f}')

   px_range = np.abs( np.max(px) - np.min(px) )
   lx = np.array([-1e2*px_range,1e2*px_range])

   
   lres.xyDashPattern = 0
   lres.xyLineThicknessF = 6
   lres.xyLineColor = reg_clr[c]
   # if c==num_case-1: lres.xyDashPattern = 2
   ngl.overlay( plot[ip], ngl.xy(wks, lx, lx*a+yi , lres) )

#----------------------------------------------------------------------------
# Ad legend
#----------------------------------------------------------------------------
lgres = ngl.Resources()
lgres.vpWidthF, lgres.vpHeightF  = 0.08, 0.08  
lgres.lgLabelFontHeightF = 0.01
lgres.lgLineThicknessF   = 8
lgres.lgMonoLineColor    = False
lgres.lgMonoDashIndex    = True
lgres.lgDashIndex        = 0
lgres.lgLineColors       = reg_clr
lgres.lgLabelJust    = 'CenterLeft'
pid = ngl.legend_ndc(wks, len(reg_label), reg_label, 0.55, 0.45, lgres)

#---------------------------------------------------------------------------------------------------
# Finalize plot
#---------------------------------------------------------------------------------------------------
hc.printline()

if 'num_plot_col' in locals():
   layout = [int(np.ceil(len(plot)/float(num_plot_col))),num_plot_col]
else:
   layout = [1,len(plot)]

pnl_res = hs.setres_panel()
ngl.panel(wks,plot,layout,pnl_res)
ngl.end()
hc.trim_png(fig_file)

if 'legend_file' in locals(): hc.trim_png(legend_file)

#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------

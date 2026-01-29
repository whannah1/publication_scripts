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
add_case('E3SM.2023-CO2-TEST-01.GNUGPU.ne30pg2_EC30to60E2r2.WCYCL1850-MMF1.1xCO2',n='1xCO2',c='blue' ,p=tmp_path_co2_mmf,s=tmp_sub)
add_case('E3SM.2023-CO2-TEST-01.GNUGPU.ne30pg2_EC30to60E2r2.WCYCL1850-MMF1.2xCO2',n='2xCO2',c='green',p=tmp_path_co2_mmf,s=tmp_sub)
add_case('E3SM.2023-CO2-TEST-01.GNUGPU.ne30pg2_EC30to60E2r2.WCYCL1850-MMF1.4xCO2',n='4xCO2',c='red'  ,p=tmp_path_co2_mmf,s=tmp_sub)
#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------

htype,yr1,yr2 = 'ha', 0, 120
fig_file,fig_type = 'figs/FXX-CRE-feedback','png'
# tmp_file_head = 'data/CRE-feedback'
tmp_file_head = 'data/CRE-feedback-alt'

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
plot = [None]*6

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
lgres.vpWidthF, lgres.vpHeightF  = 0.075, 0.085  
lgres.lgLabelFontHeightF = 0.01
lgres.lgLineThicknessF   = 12
lgres.lgMonoLineColor    = False
lgres.lgMonoDashIndex    = True
lgres.lgDashIndex        = 0
lgres.lgLineColors       = reg_clr
lgres.lgLabelJust        = 'CenterLeft'

#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------
hc.printline()
T_list = []
SALL_list = []
LALL_list = []
SCRE_list = []
LCRE_list = []
SCLR_list = []
LCLR_list = []
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
      SCRE = ds['SWCF']
      LCRE = ds['LWCF']
      SALL = ds['FSNT']
      LALL = ds['FLNT']
      SCLR = ds['FSNTC']
      LCLR = ds['FLNTC']
      #-------------------------------------------------------------------------
      # area weighted global mean
      T    = (   T*area).sum(dim='ncol') / area.sum(dim='ncol')
      SCRE = (SCRE*area).sum(dim='ncol') / area.sum(dim='ncol')
      LCRE = (LCRE*area).sum(dim='ncol') / area.sum(dim='ncol')
      SALL = (SALL*area).sum(dim='ncol') / area.sum(dim='ncol')
      LALL = (LALL*area).sum(dim='ncol') / area.sum(dim='ncol')
      SCLR = (SCLR*area).sum(dim='ncol') / area.sum(dim='ncol')
      LCLR = (LCLR*area).sum(dim='ncol') / area.sum(dim='ncol')
      #-------------------------------------------------------------------------
      T.load()
      SCRE.load()
      LCRE.load()
      SALL.load()
      LALL.load()
      SCLR.load()
      LCLR.load()
      tmp_ds = xr.Dataset( coords=T.coords )
      tmp_ds['T'] = T
      tmp_ds['SCRE'] = SCRE
      tmp_ds['LCRE'] = LCRE
      tmp_ds['SALL'] = SALL
      tmp_ds['LALL'] = LALL
      tmp_ds['SCLR'] = SCLR
      tmp_ds['LCLR'] = LCLR
      tmp_ds['yr1'] = yr1
      tmp_ds['yr2'] = yr2
      tmp_ds.to_netcdf(path=tmp_file,mode='w')
   #----------------------------------------------------------------------------
   else:
      tmp_ds = xr.open_dataset( tmp_file, use_cftime=True  )
      T = tmp_ds['T']
      SCRE = tmp_ds['SCRE']
      LCRE = tmp_ds['LCRE']
      SALL = tmp_ds['SALL']
      LALL = tmp_ds['LALL']
      SCLR = tmp_ds['SCLR']
      LCLR = tmp_ds['LCLR']
   #----------------------------------------------------------------------------
   if print_stats: 
      hc.print_stat(T,   name='T    (global mean)',stat='naxs',indent=' '*6,compact=True)
      hc.print_stat(SCRE,name='SCRE (global mean)',stat='naxs',indent=' '*6,compact=True)
      hc.print_stat(LCRE,name='LCRE (global mean)',stat='naxs',indent=' '*6,compact=True)
      hc.print_stat(SALL,name='SALL (global mean)',stat='naxs',indent=' '*6,compact=True)
      hc.print_stat(LALL,name='LALL (global mean)',stat='naxs',indent=' '*6,compact=True)
      hc.print_stat(SCLR,name='SCLR (global mean)',stat='naxs',indent=' '*6,compact=True)
      hc.print_stat(LCLR,name='LCLR (global mean)',stat='naxs',indent=' '*6,compact=True)
   #-------------------------------------------------------------------------
   # Create time index for time series panels
   time = T['time.year']
   #----------------------------------------------------------------------------
   T_list.append(T.values)
   
   SALL_list.append(SALL.values)
   LALL_list.append(LALL.values)
   SCRE_list.append(SCRE.values)
   LCRE_list.append(LCRE.values)
   SCLR_list.append(SCLR.values)
   LCLR_list.append(LCLR.values)

   time_list.append(time.values)

#---------------------------------------------------------------------------------------------------
# Calculate differences for ECS via Gregory regression
#---------------------------------------------------------------------------------------------------
dT_list = []
dS_list = []
dL_list = []

dSALL_list = []
dLALL_list = []
dSCRE_list = []
dLCRE_list = []
dSCLR_list = []
dLCLR_list = []
for c in range(num_case):
   if c==0: ct,cb = 1,0
   if c==1: ct,cb = 2,1
   if c==2: ct,cb = 2,0
   dT_list.append( T_list[ct] - T_list[cb] )
   # dS_list.append( S_list[ct] - S_list[cb] )
   # dL_list.append( L_list[ct] - L_list[cb] )
   dSALL_list.append( SALL_list[ct] - SALL_list[cb] )
   dLALL_list.append( LALL_list[ct] - LALL_list[cb] )
   dSCRE_list.append( SCRE_list[ct] - SCRE_list[cb] )
   dLCRE_list.append( LCRE_list[ct] - LCRE_list[cb] )
   dSCLR_list.append( SCLR_list[ct] - SCLR_list[cb] )
   dLCLR_list.append( LCLR_list[ct] - LCLR_list[cb] )

# if print_stats:
#    for c in range(num_case):
#       print()
#       hc.print_stat(dT_list[c],name=f'dT {reg_label[c]}',indent=' '*4,compact=True) # stat='naxs'
#       hc.print_stat(dS_list[c],name=f'dS {reg_label[c]}',indent=' '*4,compact=True) # stat='naxs'
#       hc.print_stat(dL_list[c],name=f'dL {reg_label[c]}',indent=' '*4,compact=True) # stat='naxs'
# print()

#---------------------------------------------------------------------------------------------------
# Create plot
#---------------------------------------------------------------------------------------------------
tres = copy.deepcopy(res)
# for n in range(2):
for n in range(6):
   ip = n

   # if n==0: dN_list = copy.deepcopy(dS_list)
   # if n==1: dN_list = copy.deepcopy(dL_list)

   if n==0: dN_list = copy.deepcopy(dSALL_list)
   if n==1: dN_list = copy.deepcopy(dLALL_list)
   if n==2: dN_list = copy.deepcopy(dSCRE_list)
   if n==3: dN_list = copy.deepcopy(dLCRE_list)
   if n==4: dN_list = copy.deepcopy(dSCLR_list)
   if n==5: dN_list = copy.deepcopy(dLCLR_list)

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
   dN_min = np.min([np.min(dN) for dN in dN_list])
   dN_max = np.max([np.max(dN) for dN in dN_list])
   dT_min = np.min([np.min(dT) for dT in dT_list])
   dT_max = np.max([np.max(dT) for dT in dT_list])

   dN_span = dN_max - dN_min
   dT_span = dT_max - dT_min

   tres = copy.deepcopy(res)

   tres.trYMinF = dN_min #- dN_span*0.1
   tres.trYMaxF = dN_max #+ dN_span*0.1
   tres.trXMinF = dT_min #- dT_span*0.1
   tres.trXMaxF = dT_max #+ dT_span*0.1

   delta_symbol = '~F33~D~F21~'
   tres.tiXAxisString = f'{delta_symbol}T~B~2x~N~ [K]'
   # tres.tiYAxisString = f'{delta_symbol}N~B~2x~N~ [W/m2]'
   
   if n==0: tres.tiYAxisString = f'{delta_symbol}FSNT~B~2x~N~ [W/m2]'
   if n==1: tres.tiYAxisString = f'{delta_symbol}FLNT~B~2x~N~ [W/m2]'
   if n==2: tres.tiYAxisString = f'{delta_symbol}SWCRE~B~2x~N~ [W/m2]'
   if n==3: tres.tiYAxisString = f'{delta_symbol}LWCRE~B~2x~N~ [W/m2]'
   if n==4: tres.tiYAxisString = f'{delta_symbol}FSNTC~B~2x~N~ [W/m2]'
   if n==5: tres.tiYAxisString = f'{delta_symbol}FLNTC~B~2x~N~ [W/m2]'

   tres.xyMarkLineMode = 'Markers'
   tres.xyMarker = 16
   tres.xyMarkerSizeF = 0.006
   # if n==0: tres.xyMarkerSizeF = 0.006
   # if n==1: tres.xyMarkerSizeF = 0.001

   
   print()
   for c in range(num_case):
      # ip = 1

      px,py = dT_list[c],dN_list[c]
      if c==num_case-1: px,py = px/2,py/2

      tres.xyMarkerColor = reg_clr[c]

      tplot = ngl.xy(wks, px, py, tres)
      
      if c==0:
         plot[ip] = tplot

         cstr = ''
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
      # xx = np.array([min_xi_list[c],max_xi_list[c]])
      # yy = np.array([0,0])
      # ngl.overlay( plot[ip], ngl.xy(wks, xx, yy, lres) )

      # dl = 0.1

      # xx = np.array([min_xi_list[c],min_xi_list[c]])
      # yy = np.array([0,dl])
      # ngl.overlay( plot[ip], ngl.xy(wks, xx, yy, lres) )

      # xx = np.array([max_xi_list[c],max_xi_list[c]])
      # yy = np.array([0,dl])
      # ngl.overlay( plot[ip], ngl.xy(wks, xx, yy, lres) )
      #----------------------------------------------------------------------------
      if c==0:
         # if 'lgd_label' in locals(): del lgd_label
         lgd_label = copy.deepcopy(reg_label)
         for i in range(len(lgd_label)): 
            # lgd_label[i] = f'  {reg_label[i]}    ECS: {xi_list[i]:4.2f} K'
            lgd_label[i] = f'  {reg_label[i]}    FB: {a_list[i]:4.2f} W/m2/K'
         if n==0: lgx, lgy = 0.25, 0.45
         if n==1: lgx, lgy = 0.75,  0.45
         pid = ngl.legend_ndc(wks, len(lgd_label), lgd_label, lgx, lgy, lgres)
         # pid = ngl.legend_ndc(wks, len(lgd_label), lgd_label, 0.76, 0.62, lgres)

#----------------------------------------------------------------------------
# Add legend
#----------------------------------------------------------------------------

# lgd_label = reg_label
# for i in range(len(lgd_label)): lgd_label[i] = f'  {reg_label[i]}    ECS: {xi_list[i]:4.2f} K'
# pid = ngl.legend_ndc(wks, len(lgd_label), lgd_label, 0.76, 0.62, lgres)

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

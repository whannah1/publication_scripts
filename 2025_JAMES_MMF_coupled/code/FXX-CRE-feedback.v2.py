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
# fbk_file_head = 'data/CRE-feedback-fbk'

#---------------------------------------------------------------------------------------------------
print_stats   = True
overlay_cases = True

recalculate = False

num_plot_col  = 1

#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------
num_case = len(case)

#---------------------------------------------------------------------------------------------------
# load CMIP6 data from Zelinka et al. (2020)
f = open('CMIP6_ECS_ERF_fbks.txt', 'r')
CMIP_FB_NET   = []
CMIP_FB_CLD   = []
CMIP_FB_CLDSW = []
CMIP_FB_CLDLW = []
for cnt,line in enumerate(f):
   if cnt>=9:
      if line[0]=='-': break
      line = line.strip() # strip line return character
      line = line.split()
      CMIP_FB_CLD  .append( float(line[11]) )
      CMIP_FB_CLDSW.append( float(line[12]) )
      CMIP_FB_CLDLW.append( float(line[13]) )
      CMIP_FB_NET  .append( float(line[14]) )
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
# Create plot
#---------------------------------------------------------------------------------------------------

num_fb_vars = 4
reg_label = ['2x / 1x','4x / 2x','4x / 1x']
reg_clr = ['orange','cyan','magenta']

wkres = ngl.Resources()
npix=1024*2; wkres.wkWidth,wkres.wkHeight=npix,npix
wks = ngl.open_wks(fig_type,fig_file,wkres)

res = hs.res_xy()
res.vpHeightF = 0.3
# res.vpHeightF = 0.2
res.tmYLLabelFontHeightF   = 0.01
res.tmXBLabelFontHeightF   = 0.01
res.tiXAxisFontHeightF     = 0.01
res.tiYAxisFontHeightF     = 0.01
res.tmXBMode               = 'Explicit'
res.tmXBValues             = np.array([0,1,2,3])+0.5
res.tmXBLabels             = ['Cloud Net FB','Cloud SW FB','Cloud LW FB','Net FB']

res.xyMarkLineMode         = 'Markers'
res.xyMarkerSizeF          = 0.01
res.xyMarkerThicknessF     = 8

res.trYMinF                = -2
res.trYMaxF                = 1.5
res.trXMinF                = 0.25
res.trXMaxF                = num_fb_vars - 0.25

# delta_symbol = '~F33~D~F21~'
# res.tiXAxisString = f'{delta_symbol}T~B~2x~N~ [K]'
res.tiYAxisString = f'[W/m2/K]'


for n in range(num_fb_vars):

   if n==0: CMIP_FB = CMIP_FB_CLD
   if n==1: CMIP_FB = CMIP_FB_CLDSW
   if n==2: CMIP_FB = CMIP_FB_CLDLW
   if n==3: CMIP_FB = CMIP_FB_NET

   

   #----------------------------------------------------------------------------
   # calculate regression coefficients
   a_list,xi_list,yi_list = [],[],[]
   for c in range(num_case):

      if c==0: ct,cb = 1,0
      if c==1: ct,cb = 2,1
      if c==2: ct,cb = 2,0

      px = T_list[ct] - T_list[cb]
      
      if n==0: py = ( SCRE_list[ct] - LCRE_list[ct] )  -  ( SCRE_list[cb] - LCRE_list[cb] )
      if n==1: py =   SCRE_list[ct]                    -    SCRE_list[cb]
      if n==2: py =-1*LCRE_list[ct]                    - -1*LCRE_list[cb]
      if n==3: py = ( SALL_list[ct] - LALL_list[ct] )  -  ( SALL_list[cb] - LALL_list[cb] )

      if c==2: px,py = px/2,py/2

      # simple and fast method for regression coeff and intercept
      a = np.cov( px.flatten(), py.flatten() )[1,0] / np.var( px )
      yi = np.mean(py) - a*np.mean(px)
      xi = -1*yi/a
      a_list.append(a)
      xi_list.append(xi)
      yi_list.append(yi)
   #----------------------------------------------------------------------------

   # tres.xyMarker = 2 # plus
   res.xyMarker = 4 # open circle
   res.xyMarkerColor = 'gray'
   res.xyMarkerThicknessF     = 6
   tplot = ngl.xy(wks, np.zeros(len(CMIP_FB))+n+0.5-0.05, CMIP_FB, res)

   if n==0:
      plot = tplot
   else:
      ngl.overlay(plot,tplot)

   xx = np.array([1,1])
   yy = np.mean(np.array(CMIP_FB))*xx
   res.xyMarkerThicknessF     = 10
   res.xyMarkerColor = 'black'; ngl.overlay(tplot,ngl.xy(wks, (n+0.5-0.05)*xx, yy, res))


   res.xyMarker = 16
   res.xyMarkerColor = reg_clr[0]; ngl.overlay(plot,ngl.xy(wks, (n+0.5+0.05)*xx, a_list[0]*xx, res))
   res.xyMarkerColor = reg_clr[1]; ngl.overlay(plot,ngl.xy(wks, (n+0.5+0.05)*xx, a_list[1]*xx, res))
   res.xyMarkerColor = reg_clr[2]; ngl.overlay(plot,ngl.xy(wks, (n+0.5+0.05)*xx, a_list[2]*xx, res))

   # if n==0:
   #    plot = tplot
   # else:
   #    ngl.overlay(plot,tplot)
#-------------------------------------------------------------------------------
hs.set_subtitles(wks, plot, '', 'Global Mean Feedbacks', '', font_height=0.015)
#-------------------------------------------------------------------------------
# add horizontal line
lres = hs.res_xy()
lres.xyDashPattern    = 0
lres.xyLineThicknessF = 2
lres.xyLineColor      = 'black'

ngl.overlay( plot, ngl.xy(wks, np.array([-1e3,1e3]), np.array([0,0]), lres) )

#-------------------------------------------------------------------------------
# Add legend

lgd_clr = reg_clr
lgd_lbl = reg_label

lgd_clr.insert(0,'gray')
lgd_lbl.insert(0,'CMIP6')

lgres = ngl.Resources()
lgres.vpWidthF, lgres.vpHeightF  = 0.06, 0.12
lgres.lgLabelFontHeightF         = 0.015
lgres.lgLabelJust                = 'CenterLeft'
# lgres.lgLineThicknessF   = 60
# lgres.lgMonoLineColor    = False
# lgres.lgMonoDashIndex    = True
# lgres.lgDashIndex        = 0
# lgres.lgLineColors       = lgd_clr
lgres.lgItemType         = 'Markers'
lgres.lgMarkerIndexes    = [4,16,16,16]
lgres.lgMarkerColors     = lgd_clr
lgres.lgMarkerSizeF      = 0.015
lgres.lgMarkerThicknessF = 15

for i in range(len(lgd_lbl)): lgd_lbl[i] = f'  {lgd_lbl[i]}'

pid = ngl.legend_ndc(wks, len(lgd_lbl), lgd_lbl, 0.78, 0.7, lgres)

#---------------------------------------------------------------------------------------------------
# Finalize plot
#---------------------------------------------------------------------------------------------------
hc.printline()

# if 'num_plot_col' in locals():
#    layout = [int(np.ceil(len(plot)/float(num_plot_col))),num_plot_col]
# else:
#    layout = [1,len(plot)]

# pnl_res = hs.setres_panel()
# pnl_res.nglPanelFigureStrings            = list(string.ascii_lowercase)
# pnl_res.nglPanelFigureStringsJust        = "TopLeft"
# pnl_res.nglPanelFigureStringsFontHeightF = 0.01

# ngl.panel(wks,plot,[1,1],pnl_res)

ngl.draw(plot)
ngl.frame(wks)

ngl.end()
hc.trim_png(fig_file)

if 'legend_file' in locals(): hc.trim_png(legend_file)

#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------

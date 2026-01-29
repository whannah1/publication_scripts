import os, glob, ngl, subprocess as sp, numpy as np, xarray as xr, copy, string, sys, cmocean
import hapy_common as hc, hapy_E3SM   as he, hapy_setres as hs
import logging; logging.basicConfig(level=logging.WARNING)
import wavenumber_frequency_functions as wf
#-------------------------------------------------------------------------------

fig_type,fig_file = 'png',f'figs/FXX-wave-filter-bounds'

use_common_label_bar = True

lat1,lat2 = -15,15
yr1,yr2 = 1995,2004

segsize,noverlap = 96,60

tm_period_values  = np.array([1,2,3,5,10,30,90])

#---------------------------------------------------------------------------------------------------
# Set up plot resources
#---------------------------------------------------------------------------------------------------

subtitle_font_height = 0.01

wkres = ngl.Resources()
npix = 2048; wkres.wkWidth,wkres.wkHeight=npix,npix
wks = ngl.open_wks(fig_type,fig_file,wkres)

plot = [None]*2

res = hs.res_contour_fill()
res.vpHeightF = 0.5
res.lbLabelFontHeightF     = 0.022
res.tiXAxisFontHeightF     = 0.025
res.tiYAxisFontHeightF     = 0.025
res.tmYLLabelFontHeightF   = 0.02
res.tmXBLabelFontHeightF   = 0.02
res.tiXAxisString = 'Zonal Wavenumber'
# res.tiYAxisString = 'Frequency (cpd)'
res.tiYAxisString = 'Period [days]'

if use_common_label_bar:  res.lbLabelBarOn = False

res.trYMinF,res.trYMaxF,res.trXMinF,res.trXMaxF = 0.0,0.5,-25,25
# res.trYMinF,res.trYMaxF,res.trXMinF,res.trXMaxF = 0.0,1./5.,-6,6


# res.cnFillPalette = np.array( cmocean.cm.rain(np.linspace(0,1,256)) )

lres = hs.res_xy()
lres.xyLineThicknessF = 1
lres.xyLineColor      = 'black'
lres.xyDashPattern    = 0

#---------------------------------------------------------------------------------------------------
# routine for adding theoretical dispersion curves and MJO frequency bounds
#---------------------------------------------------------------------------------------------------
# get data for dispersion curves
swfq, swwn = wf.genDispersionCurves(Ahe=[50, 25, 12]) # Ahe = eq depth
# swfq.shape # -->(6, 3, 50)
swf = np.ma.masked_invalid( np.where(swfq==1e20, np.nan, swfq) )
swk = np.ma.masked_invalid( np.where(swwn==1e20, np.nan, swwn) )
def add_dispersion_curves(wks_in,plot_in,type='sym'):
   if type=='asm': wave_idx_list = [0,1,2]
   if type=='sym': wave_idx_list = [3,4,5]
   for ii in wave_idx_list:
      lres.xyDashPattern = 0
      ngl.overlay( plot_in, ngl.xy(wks_in, np.array([0,0]), np.array([-1e3,1e3]), lres) )
      ngl.overlay( plot_in, ngl.xy(wks_in, swk[ii, 0,:], swf[ii,0,:], lres) )
      ngl.overlay( plot_in, ngl.xy(wks_in, swk[ii, 1,:], swf[ii,1,:], lres) )
      ngl.overlay( plot_in, ngl.xy(wks_in, swk[ii, 2,:], swf[ii,2,:], lres) )
   # overlay lines for MJO frequency range
   lres.xyDashPattern = 1
   tfrq=1./30.; ngl.overlay( plot_in, ngl.xy(wks_in, np.array([-1e3,1e3]), np.array([tfrq,tfrq]), lres) )
   tfrq=1./90.; ngl.overlay( plot_in, ngl.xy(wks_in, np.array([-1e3,1e3]), np.array([tfrq,tfrq]), lres) )
   return
#---------------------------------------------------------------------------------------------------
# put white rectangle over satellite aliasing signal
#---------------------------------------------------------------------------------------------------
pgres = ngl.Resources()
pgres.nglDraw,pgres.nglFrame = False,False
pgres.gsLineThicknessF = 1
pgres.gsFillColor = 'white'
pgres.gsLineColor = 'white'
awn1,awn2 = 13,16
afq1,afq2 = 0.095,0.125
def hide_aliasing_artifact(wks_in,plot_in):
   bx = [awn1,awn1,awn2,awn2,awn1]
   by = [afq1,afq2,afq2,afq1,afq1]
   pdum = ngl.add_polygon(wks_in,plot_in,bx,by,pgres)
#---------------------------------------------------------------------------------------------------

tmp_file_head = 'data/wk-wave-spectra'
tcase = 'NOAA'
file_var_str = 'FLUT'

wk_tmp_file = f'{tmp_file_head}.{tcase}.{file_var_str}'
wk_tmp_file+= f'.yr_{yr1}_{yr2}'
wk_tmp_file+= f'.lat_{lat1}_{lat2}'
wk_tmp_file+= f'.seg_{segsize}.novr_{noverlap}'
wk_tmp_file+= f'.daily.nc'

#-------------------------------------------------------------------------------

print(f'      loading pre-calculated spectra... {wk_tmp_file}')
wk_ds = xr.open_mfdataset( wk_tmp_file )
num_days   = wk_ds['num_days'].values
# segsize    = wk_ds['segsize']
# noverlap   = wk_ds['noverlap']
rspec_sym  = wk_ds['rspec_sym']
rspec_asy  = wk_ds['rspec_asy']
bkgd_all   = wk_ds['bkgd_all']
bkgd_sym   = wk_ds['bkgd_sym']
bkgd_asm   = wk_ds['bkgd_asm']

#-------------------------------------------------------------------------------
# normalize by the smoothed "background" spectra
nspec_sym = rspec_sym / bkgd_sym
nspec_asy = rspec_asy / bkgd_asm

spec_data_sym = nspec_sym.transpose().values
spec_data_asm = nspec_asy.transpose().values

freq_data = wk_ds['frequency'].values
wvnm_data = wk_ds['wavenumber'].values

#-------------------------------------------------------------------------------
# calculate limits for common color bar
# data_min = np.min(spec_data)
# data_max = np.max(spec_data)
#-------------------------------------------------------------------------------
# Create plot

tres = copy.deepcopy(res)
# colors and levels from E3SM diags
tres.cnFillPalette = ['white', 'gainsboro','lightgray','gray','paleturquoise',
                      'skyblue','palegreen','mediumseagreen','seagreen','yellow',
                      'orange','red','maroon','pink']
tres.cnLevels = [ 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.7, 1.9, 2.1, 2.4, 2.7, 3.0 ]
tres.cnLevelSelectionMode = 'ExplicitLevels'


tres.sfXArray = wvnm_data
tres.sfYArray = freq_data

tres.tmYLMode     = 'Explicit'
tres.tmYLValues   = 1./tm_period_values
tres.tmYLLabels   = tm_period_values


plot[0] = ngl.contour(wks, spec_data_sym, tres)
plot[1] = ngl.contour(wks, spec_data_asm, tres)


add_dispersion_curves(wks,plot[0],type='sym')
add_dispersion_curves(wks,plot[1],type='asm')
hide_aliasing_artifact(wks,plot[0])
hide_aliasing_artifact(wks,plot[1])

#-------------------------------------------------------------------------------
pdum = []
def draw_box(f1,f2,w1,w2,ip=None,c='black',d=0):
   pgres = ngl.Resources()
   pgres.nglDraw,pgres.nglFrame = False,False
   pgres.gsEdgesOn = True
   pgres.gsEdgeColor       = c
   pgres.gsEdgeDashPattern = d
   pgres.gsEdgeThicknessF  = 4
   # pgres.cnFillPattern     = 4 
   pgres.gsFillIndex       = 4
   bx = np.array([w1,w2,w2,w1,w1])
   by = np.array([f1,f1,f2,f2,f1])
   pdum.append( ngl.add_polygon(wks,plot[ip],bx,by,pgres) )
#-------------------------------------------------------------------------------
# draw boxes around filter regions

# # Wandi's regions
# draw_box(f1=0.0,f2=0.4,w1=  1,w2=20,ip=0,c='blue') # Kelvin
# draw_box(f1=0.0,f2=0.1,w1=-20,w2= 0,ip=0,c='red') # Rossby
# draw_box(f1=0.1,f2=0.5,w1=-20,w2=20,ip=1,c='green') # MRG

# New regions
draw_box(f1=1/ 20,f2=1/ 2.5,w1=  1,w2=20,ip=0,c='blue') # Kelvin
draw_box(f1=0,    f2=1/10,  w1=-20,w2= 0,ip=0,c='red') # Rossby
draw_box(f1=1/ 10,f2=1/ 2,  w1=-20,w2=20,ip=1,c='green') # MRG
draw_box(f1=1/100,f2=1/20,  w1=  0,w2=10,ip=0,c='purple',d=0) # MJO
draw_box(f1=1/100,f2=1/20,  w1=  0,w2=10,ip=1,c='purple',d=0) # MJO

#-------------------------------------------------------------------------------
### wave filter bounds from wavenumer_frequency_functions.py
# if waveName in ['Kelvin','kelvin','KELVIN']:
#   freqInd = (fftData['frequency']<=0.4) & (fftData['frequency']>=0)
#   kInd = (fftData['wavenumber']<=20) & (fftData['wavenumber']>=1)

# elif waveName in ['MRG','IG0','mrg','ig0']:
#   freqInd = (fftData['frequency']<=0.5) & (fftData['frequency']>=0.1)
#   kInd = (fftData['wavenumber']<=0) & (fftData['wavenumber']>=-20)

# elif waveName in ['Rossby','rossby','ROSSBY']:
#   freqInd = (fftData['frequency']<=0.1) & (fftData['frequency']>=0)
#   kInd = (fftData['wavenumber']<=0) & (fftData['wavenumber']>=-20)
#-------------------------------------------------------------------------------

hs.set_subtitles(wks, plot[0], left_string='NOAA', center_string='Symmetric', right_string='OLR', font_height=subtitle_font_height)
hs.set_subtitles(wks, plot[1], left_string='NOAA', center_string='Asymmetric', right_string='OLR', font_height=subtitle_font_height)

#---------------------------------------------------------------------------------------------------
# Finalize plot
#---------------------------------------------------------------------------------------------------

pnl_res = hs.setres_panel()
# pnl_res.nglPanelYWhiteSpacePercent = 5
# pnl_res.nglPanelXWhiteSpacePercent = 5
# pnl_res.nglPanelFigureStrings            = list(string.ascii_lowercase)
# pnl_res.nglPanelFigureStringsJust        = "TopLeft"
# pnl_res.nglPanelFigureStringsFontHeightF = 0.015

if use_common_label_bar: 
   pnl_res.nglPanelLabelBar = True
   pnl_res.nglPanelLabelBarLabelFontHeightF = 0.01
   pnl_res.nglPanelLabelBarWidthF = 0.5
   # pnl_res.lbLeftMarginF      = 0.3
   # pnl_res.lbRightMarginF     = 0.3

ngl.panel(wks,plot,[1,len(plot)],pnl_res)
hc.trim_png(fig_file)

ngl.end()
#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------

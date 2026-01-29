# v4 is similar to v3 except panels are placed 'manually' using ImageMagick
#-------------------------------------------------------------------------------
import os, ngl, copy, string, xarray as xr, numpy as np, glob, dask, numba, cmocean, subprocess as sp
import hapy_common as hc, hapy_E3SM as he, hapy_setres as hs
#-------------------------------------------------------------------------------
'''

'''
#-------------------------------------------------------------------------------
case_name,case, = [],[]
file_sim,file_obs = [],[]
def add_case(case_in,n='',f=None,obs=None):
   global case_name,case,case_file,case_sub,clr,dsh,mrk
   case.append(case_in); case_name.append(n)
   file_sim.append(f); file_obs.append(obs)
#-------------------------------------------------------------------------------
# scrip_file_sim = '/pscratch/sd/w/whannah/e3sm_scratch/files_grid/ne1024pg2_scrip.nc'
# scrip_file_obs = '/pscratch/sd/w/whannah/e3sm_scratch/files_grid/IMERG_1800x3600_scrip.nc'
# scrip_file_era = '~/HICCUP/files_grid/scrip_ERA5_721x1440.nc'
sim_data_root  = '/pscratch/sd/p/paullric/SCREAM_40day'
# obs_data_root  = '/pscratch/sd/w/whannah/Obs'

sim_topo_file = '/global/cfs/cdirs/e3sm/inputdata/atm/cam/topo/USGS-gtopo30_ne1024np4pg2_16xconsistentSGH_20190528_converted.nc'

add_case('DYAMOND2_SCREAMv1' ,n='SCREAMv1 Jan 2020',f=f'{sim_data_root}/dy2/ne1024pg2_ne1024pg2_DY2_AR_climo.nc'     ,obs='data/e5.oper.an.ar_climo.2020012200_2020022823.remap.nc')
add_case('Apr1_2013_SCREAMv1',n='SCREAMv1 Apr 2013',f=f'{sim_data_root}/apr/ne1024pg2_ne1024pg2_Apr2013_AR_climo.nc' ,obs='data/e5.oper.an.ar_climo.2013040300_2013050923.remap.nc')
add_case('DYAMOND1_SCREAMv1' ,n='SCREAMv1 Aug 2016',f=f'{sim_data_root}/dy1/ne1024pg2_ne1024pg2_DY1_AR_climo.nc'     ,obs='data/e5.oper.an.ar_climo.2016080300_2016090923.remap.nc')
add_case('Oct1_2013_SCREAMv1',n='SCREAMv1 Oct 2013',f=f'{sim_data_root}/oct/ne1024pg2_ne1024pg2_40dayrun_AR_climo.nc',obs='data/e5.oper.an.ar_climo.2013100300_2013110923.remap.nc')

#-------------------------------------------------------------------------------
fig_file,fig_type = f'figs/AR-density','png'

lat1,lat2 = -70,70

use_common_labelbar  = True
# subtitle_font_height = 0.008
subtitle_font_height = 0.012
# print_stats          = False

#---------------------------------------------------------------------------------------------------
def run_cmd(cmd,verbose=True,indent='  '):
   cmd_str = '\n' + indent + hc.tcolor.GREEN + cmd + hc.tcolor.ENDC + '\n'
   if verbose: print(cmd_str)
   proc = sp.Popen(cmd, stdout=sp.PIPE, stderr=sp.PIPE, shell=True, universal_newlines=True)
   (msg, err) = proc.communicate()
   if verbose and msg!='': print(f'  msg: {msg}')
   if err!='' and not verbose: print(cmd_str)
   if err!='': print(f'err: {err}'); exit()
   return msg
#---------------------------------------------------------------------------------------------------
num_case = len(case)

res = hs.res_contour_fill_map()
res.vpHeightF = 0.38
# res.vpWidthF  = 0.4
res.tmYLLabelFontHeightF         = 0.01
res.tmXBLabelFontHeightF         = 0.01
# res.tiYAxisString = 'Latitude'
# res.tiXAxisString = 'Longitude'
res.lbLabelFontHeightF    = 0.015
res.lbTitlePosition       = 'Bottom'
res.lbTopMarginF          = -0.2
res.lbBottomMarginF       =  0.2+0.01
# res.lbLeftMarginF     = 0.5
# res.lbRightMarginF    = 0.5
# res.cnFillPalette = 'MPL_viridis'



res.mpLimitMode           = 'LatLon'
res.mpMinLatF             = lat1
res.mpMaxLatF             = lat2
# res.tiYAxisString         = 'Latitude'
# res.tiXAxisString         = 'Longitude'
# res.tmYLLabelFontHeightF  = 0.01
# res.tmXBLabelFontHeightF  = 0.01
# res.lbLabelFontHeightF    = 0.01
# res.lbLeftMarginF         = 0.5
# res.lbRightMarginF        = 0.5
# res.lbTitlePosition       = 'Bottom'
res.tmXBOn                = False
res.tmYLOn                = False
# res.mpGeophysicalLineThicknessF = 2

if use_common_labelbar: res.lbLabelBarOn = False

# cnt_res = hs.res_contour()
# cnt_res.tmXBOn                = False
# cnt_res.tmYLOn                = False
# cnt_res.cnLineThicknessF      = 2.

lres = hs.res_xy()
lres.xyDashPattern    = 0
lres.xyLineThicknessF = 4
lres.xyLineColor      = 'red'

#---------------------------------------------------------------------------------------------------
fig_file_list_str = ''
for n in range(2):
   wkres = ngl.Resources()
   npix = 2048; wkres.wkWidth,wkres.wkHeight=npix,npix
   fig_file_tmp = f'{fig_file}.{n:02d}'
   wks = ngl.open_wks(fig_type,fig_file_tmp,wkres)
   fig_file_list_str += f' {fig_file_tmp}.{fig_type}'
   plot = [None]*num_case
   for c in range(num_case):
      print(' '*4+'case: '+hc.tcolor.CYAN+case[c]+hc.tcolor.ENDC)
      #-------------------------------------------------------------------------
      ds_sim = xr.open_dataset( file_sim[c] )
      data_sim = ds_sim['annualmean_binary_tag'   ].isel(time=0)
      #-------------------------------------------------------------------------
      if n==1:
         ds_obs = xr.open_dataset( file_obs[c] )
         data_obs = ds_obs['annualmean_AR_binary_tag'].isel(time=0)
         data_obs = data_obs.rename({'longitude':'lon'})
         data_obs = data_obs.rename({'latitude':'lat'})
         data_diff = data_sim - data_obs
      #-------------------------------------------------------------------------
      res.sfXArray = ds_sim['lon'].values
      res.sfYArray = ds_sim['lat'].values
      res.cnLevelSelectionMode = 'ExplicitLevels'
      #-------------------------------------------------------------------------
      if n==0:
         tres = copy.deepcopy(res)
         tres.cnFillPalette = np.array( cmocean.cm.rain(np.linspace(0,1,256)) )
         # tres.cnLevels      = np.arange(0.05,0.4+0.05,0.05)
         tres.cnLevels      = np.arange(0.04,0.44+0.08,0.08)
         ip = num_case*0+c
         plot[ip] = ngl.contour_map(wks, data_sim.values, tres)
         hs.set_subtitles(wks, plot[ip], case_name[c], '', '', font_height=subtitle_font_height)
      #-------------------------------------------------------------------------
      if n==1:
         tres = copy.deepcopy(res)
         tres.cnFillPalette = np.array( cmocean.cm.balance(np.linspace(0,1,256)) )
         tres.cnLevels      = np.arange(-0.22,0.22+0.04,0.04)
         ip = num_case*0+c
         plot[ip] = ngl.contour_map(wks, data_diff.values, tres)
         hs.set_subtitles(wks, plot[ip], case_name[c], '', 'Diff from ERA5', font_height=subtitle_font_height)
   #----------------------------------------------------------------------------
   # layout = [1,num_case]
   layout = [num_case,1]
   pnl_res = hs.setres_panel()
   pnl_res.nglPanelFigureStrings = list(string.ascii_lowercase)[num_case*n:num_case*(n+1)]
   pnl_res.nglPanelFigureStringsJust        = 'TopLeft'
   pnl_res.nglPanelFigureStringsFontHeightF = 0.01
   pnl_res.lbTitlePosition                  = 'Bottom'
   pnl_res.lbTitleFontHeightF               = 0.012
   pnl_res.lbTitleString                    = 'Time Fraction of Atmospheric River'
   if use_common_labelbar: 
      pnl_res.nglPanelTop    = 0.9
      pnl_res.nglPanelBottom = 0.1
      pnl_res.nglPanelLabelBar = True
      # pnl_res.nglPanelLabelBarHeightF          = 0.05
      # pnl_res.nglPanelLabelBarWidthF           = 0.4
      # pnl_res.nglPanelLabelBarLabelFontHeightF = 0.01
      pnl_res.nglPanelLabelBarHeightF          = 0.05
      pnl_res.nglPanelLabelBarWidthF           = 0.4
      pnl_res.nglPanelLabelBarLabelFontHeightF = 0.012
   ngl.panel(wks,plot,layout,pnl_res)
   ngl.destroy(wks)
   hc.trim_png(fig_file_tmp)
ngl.end()
#---------------------------------------------------------------------------------------------------
# Combine images
run_cmd(f'montage {fig_file_list_str} -geometry 1024x2048+20 -tile 2x1 {fig_file}.{fig_type}')
# run_cmd(f'montage {fig_file}.{fig_type} figs/AR-histogram.png -geometry 2048x2048+20 -tile 2x1 {fig_file}.{fig_type}')
hc.trim_png(fig_file)
#---------------------------------------------------------------------------------------------------

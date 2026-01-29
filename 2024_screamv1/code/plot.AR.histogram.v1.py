# v4 is similar to v3 except panels are placed 'manually' using ImageMagick
#-------------------------------------------------------------------------------
import os, ngl, copy, string, xarray as xr, numpy as np, glob, dask, numba, cmocean, subprocess as sp
import hapy_common as hc, hapy_E3SM as he, hapy_setres as hs
#-------------------------------------------------------------------------------
'''

'''
#-------------------------------------------------------------------------------
case_name,case,case_dir = [],[],[]
def add_case(case_in,n='',p=None):
   global case_name,case,case_dir
   case.append(case_in); case_name.append(n); case_dir.append(p)
#-------------------------------------------------------------------------------
# scrip_file_sim = '/pscratch/sd/w/whannah/e3sm_scratch/files_grid/ne1024pg2_scrip.nc'
# scrip_file_obs = '/pscratch/sd/w/whannah/e3sm_scratch/files_grid/IMERG_1800x3600_scrip.nc'
# scrip_file_era = '~/HICCUP/files_grid/scrip_ERA5_721x1440.nc'
sim_data_root  = '/pscratch/sd/p/paullric/SCREAM_40day'
# obs_data_root  = '/pscratch/sd/w/whannah/Obs'
# sim_topo_file = '/global/cfs/cdirs/e3sm/inputdata/atm/cam/topo/USGS-gtopo30_ne1024np4pg2_16xconsistentSGH_20190528_converted.nc'

add_case('DYAMOND2_SCREAMv1' ,n='Jan 2020',p=f'{sim_data_root}/dy2')
add_case('Apr1_2013_SCREAMv1',n='Apr 2013',p=f'{sim_data_root}/apr')
add_case('DYAMOND1_SCREAMv1' ,n='Aug 2016',p=f'{sim_data_root}/dy1')
add_case('Oct1_2013_SCREAMv1',n='Oct 2013',p=f'{sim_data_root}/oct')

#-------------------------------------------------------------------------------
fig_file,fig_type = f'figs/AR-histogram-v1','png'

subtitle_font_height = 0.01

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

wkres = ngl.Resources()
npix = 2048; wkres.wkWidth,wkres.wkHeight=npix,npix
wks = ngl.open_wks(fig_type,fig_file,wkres)
plot = [None]*num_case*2

res = hs.res_xy()
res.vpHeightF = 0.4
res.tmYLLabelFontHeightF         = 0.028
res.tmXBLabelFontHeightF         = 0.028
res.tiXAxisFontHeightF           = 0.028
res.tiYAxisFontHeightF           = 0.028
res.trYMinF                      = 0
res.tiYAxisString                = 'Frequency'
res.xyLineThicknessF = 4

# lres = hs.res_xy()
# lres.xyDashPattern    = 0
# lres.xyLineThicknessF = 4
# lres.xyLineColor      = 'red'

#---------------------------------------------------------------------------------------------------
wid_dbin = 0.6; wid_bins=np.arange(wid_dbin/2, 8+wid_dbin/2,wid_dbin) ; wid_nbin=len(wid_bins)
len_dbin = 4  ; len_bins=np.arange(len_dbin/2,40+len_dbin/2,len_dbin) ; len_nbin=len(len_bins)
#---------------------------------------------------------------------------------------------------
wid_sim_list = []
wid_obs_list = []
len_sim_list = []
len_obs_list = []
for c in range(num_case):
   print(' '*4+'case: '+hc.tcolor.CYAN+case[c]+hc.tcolor.ENDC)
   file_sim = glob.glob(f'{case_dir[c]}/ne1024pg2_ne1024pg2_*_stats.txt')[0]
   file_obs = glob.glob(f'{case_dir[c]}/e5.oper.an.ar_stats.*.txt')[0]
   #----------------------------------------------------------------------------
   data_sim = np.loadtxt(file_sim)
   data_obs = np.loadtxt(file_obs)
   ar_wid_sim, ar_len_sim = data_sim[:,6], data_sim[:,7]
   ar_wid_obs, ar_len_obs = data_obs[:,6], data_obs[:,7]
   #----------------------------------------------------------------------------
   wid_cnt_sim = np.zeros(wid_nbin)
   wid_cnt_obs = np.zeros(wid_nbin)
   for b in range( wid_nbin ):
      bin_bot = wid_bins[b] - wid_dbin/2. ; bin_top = bin_bot + wid_dbin
      for n in range( len(ar_wid_sim) ):
         if ( ar_wid_sim[n]>=bin_bot ) & ( ar_wid_sim[n]<bin_top ): wid_cnt_sim[b] = wid_cnt_sim[b]+1
         if ( ar_wid_obs[n]>=bin_bot ) & ( ar_wid_obs[n]<bin_top ): wid_cnt_obs[b] = wid_cnt_obs[b]+1
   wid_cnt_sim = wid_cnt_sim/np.sum(wid_cnt_sim)
   wid_cnt_obs = wid_cnt_obs/np.sum(wid_cnt_obs)
   #----------------------------------------------------------------------------
   len_cnt_sim = np.zeros(len_nbin)
   len_cnt_obs = np.zeros(len_nbin)
   for b in range( len_nbin ):
      bin_bot = len_bins[b] - len_dbin/2. ; bin_top = bin_bot + len_dbin
      for n in range( len(ar_len_sim) ):
         if ( ar_len_sim[n]>=bin_bot ) & ( ar_len_sim[n]<bin_top ): len_cnt_sim[b] = len_cnt_sim[b]+1
         if ( ar_len_obs[n]>=bin_bot ) & ( ar_len_obs[n]<bin_top ): len_cnt_obs[b] = len_cnt_obs[b]+1
   len_cnt_sim = len_cnt_sim/np.sum(len_cnt_sim)
   len_cnt_obs = len_cnt_obs/np.sum(len_cnt_obs)
   #----------------------------------------------------------------------------
   wid_sim_list.append(wid_cnt_sim)
   wid_obs_list.append(wid_cnt_obs)
   len_sim_list.append(len_cnt_sim)
   len_obs_list.append(len_cnt_obs)
#---------------------------------------------------------------------------------------------------
wid_sim_max = np.max([np.nanmax(d) for d in wid_sim_list])
wid_obs_max = np.max([np.nanmax(d) for d in wid_obs_list])
wid_max     = np.max([wid_sim_max,wid_obs_max])
len_sim_max = np.max([np.nanmax(d) for d in len_sim_list])
len_obs_max = np.max([np.nanmax(d) for d in len_obs_list])
len_max     = np.max([len_sim_max,len_obs_max])
for c in range(num_case):
   wid_cnt_sim = wid_sim_list[c]
   wid_cnt_obs = wid_obs_list[c]
   len_cnt_sim = len_sim_list[c]
   len_cnt_obs = len_obs_list[c]
   #----------------------------------------------------------------------------
   # ip = c*2
   ip = num_case*0+c
   tres = copy.deepcopy(res)
   tres.tiXAxisString = 'Width (degrees)'
   tres.trXMinF, tres.trXMaxF = np.min(wid_bins), np.max(wid_bins)
   tres.trYMaxF = wid_max
   tres.xyDashPattern = 0
   plot[ip] = ngl.xy(wks, wid_bins, wid_cnt_sim, tres)
   tres.xyDashPattern = 16
   ngl.overlay( plot[ip], ngl.xy(wks, wid_bins, wid_cnt_obs, tres) )
   hs.set_subtitles(wks, plot[ip], case_name[c], '', '', font_height=subtitle_font_height)
   #----------------------------------------------------------------------------
   # ip = c*2+1
   ip = num_case*1+c
   tres = copy.deepcopy(res)
   tres.trXMinF, tres.trXMaxF = np.min(len_bins), np.max(len_bins)
   tres.trYMaxF = len_max
   tres.tiXAxisString = 'Length (degrees)'
   tres.xyDashPattern = 0
   plot[ip] = ngl.xy(wks, len_bins, len_cnt_sim, tres)
   tres.xyDashPattern = 16
   ngl.overlay( plot[ip], ngl.xy(wks, len_bins, len_cnt_obs, tres) )
   hs.set_subtitles(wks, plot[ip], case_name[c], '', '', font_height=subtitle_font_height)
#---------------------------------------------------------------------------------------------------
# layout = [num_case,2]
layout = [2,num_case]
pres = hs.setres_panel()
pres.nglFrame = False
pres.nglPanelFigureStrings = list(string.ascii_lowercase)
pres.nglPanelFigureStringsJust        = 'TopLeft'
pres.nglPanelFigureStringsFontHeightF = 0.012
ngl.panel(wks,plot,layout,pres)
#-------------------------------------------------------------------------------
# Add legend
lgres = ngl.Resources()
lgres.vpWidthF           = 0.05
lgres.vpHeightF          = 0.03
lgres.lgLabelFontHeightF = 0.008
lgres.lgLabelFont        = "courier"
# lgres.lgMonoDashIndex    = True
lgres.lgMonoLineColor    = True
lgres.lgLineLabelsOn     = False
lgres.lgLineThicknessF   = 3
lgres.lgLabelJust        = 'CenterLeft'
# lgres.lgLineColors       = clr
lgres.lgDashIndexes      = [1,0]

labels = ['ERA5','SCREAM']
max_label_len = max([len(n) for n in labels])+1
for n,nn in enumerate(labels): labels[n] = f' {nn:{max_label_len}}'

ndcx,ndcy = 0.14,0.66

pid = ngl.legend_ndc(wks, len(labels), labels, ndcx, ndcy, lgres)

#-------------------------------------------------------------------------------
ngl.frame(wks)
ngl.end()
hc.trim_png(fig_file)
#---------------------------------------------------------------------------------------------------

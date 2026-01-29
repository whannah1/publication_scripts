import os, glob, ngl, subprocess as sp, numpy as np, xarray as xr, copy, string, sys
import hapy_common as hc, hapy_E3SM   as he, hapy_setres as hs
sys.path.insert(1, os.getenv('HOME')+'/wavenumber_frequency-master')
import wavenumber_frequency_functions as wf
import cmocean
# scratch = '/global/cfs/cdirs/m3312/whannah/2022-QBO-TEST'
#-------------------------------------------------------------------------------
case_dir,case_sub = [],[]
case,name = [],[]
clr,dsh,mrk = [],[],[]
def add_case(case_in,n=None,p=None,s=None,d=0,c='black',m=0):
   global name,case,case_dir,case_sub,clr,dsh,mrk
   case.append(case_in); name.append(n)
   case_dir.append(p); case_sub.append(s)
   dsh.append(d) ; clr.append(c) ; mrk.append(m)
#-------------------------------------------------------------------------------
var,lev_list,vstr,vclr,vdsh = [],[],[],[],[]
unit_list = []
def add_var(var_name,var_str=None,lev=None,c='black',d=0,unit=None): 
   var.append(var_name); lev_list.append(lev)
   if var_str is None: var_str=var_name
   vstr.append(var_str); vclr.append(c); vdsh.append(d)
   unit_list.append(unit)
#-------------------------------------------------------------------------------
tmp_path_co2_mmf = '/global/cfs/cdirs/m3312/whannah/2023-CPL'
tmp_path_hst_mmf = '/global/cfs/cdirs/m3312/whannah/2023-CPL'
tmp_path_hst_v2  = '/global/cfs/cdirs/m3312/whannah/e3smv2_historical'
tmp_sub = 'archive/atm/hist'
# add_case('NOAA',n='NOAA')
add_case('v2.LR.historical_0101',                                  n='E3SMv2',  p=tmp_path_hst_v2, s='archive/atm/hist')
add_case('E3SM.INCITE2023-CPL.ne30pg2_EC30to60E2r2.WCYCL20TR-MMF1',n='E3SM-MMF',p=tmp_path_hst_mmf,s='archive/atm/hist')
# add_case('E3SM.2023-CO2-TEST-01.GNUGPU.ne30pg2_EC30to60E2r2.WCYCL1850-MMF1.1xCO2',n='MMF PI 1xCO2',c='blue' ,p=tmp_path_co2_mmf,s=tmp_sub)
# add_case('E3SM.2023-CO2-TEST-01.GNUGPU.ne30pg2_EC30to60E2r2.WCYCL1850-MMF1.2xCO2',n='MMF PI 2xCO2',c='green',p=tmp_path_co2_mmf,s=tmp_sub)
# add_case('E3SM.2023-CO2-TEST-01.GNUGPU.ne30pg2_EC30to60E2r2.WCYCL1850-MMF1.4xCO2',n='MMF PI 4xCO2',c='red'  ,p=tmp_path_co2_mmf,s=tmp_sub)

#-------------------------------------------------------------------------------

add_var('PRECT','Precipitation'    ,unit='mm/day')
add_var('FLUT', 'OLR'              ,unit='W/m2')
add_var('U850', '850 mb Zonal wind',unit='m/s')

spec_type = 'sym' # sym / asym / tot

fig_type,fig_file = 'png',f'figs/F07-wk-spectra-{spec_type}'
fig_file_diff = f'{fig_file}.diff'
tmp_file_head = 'data/wk-wave-spectra'

lat1,lat2 = -15,15

# num_files = 365*10 # specified from last file
yr1,yr2 = 1980,2014 # 45 years
htype = 'h1'


remap_grid = '90x180'

recalculate_spectra = False

var_x_case  = False

normalize, add_diff   = True, True

# num_plot_col = 2

use_common_label_bar = True

tm_period_values  = np.array([1,2,3,5,10,30,90])

# segsize,noverlap = 180,0
segsize,noverlap = 96,60

#---------------------------------------------------------------------------------------------------
def run_cmd(cmd,verbose=True,indent='  '):
   cmd_str = indent + hc.tcolor.GREEN + cmd + hc.tcolor.ENDC
   if verbose: print('\n'+cmd_str)
   proc = sp.Popen(cmd, stdout=sp.PIPE, stderr=sp.PIPE, shell=True, universal_newlines=True)
   (msg, err) = proc.communicate()
   if verbose and msg!='': print(f'  msg: {msg}')
   if err!='' and not verbose: print(cmd_str)
   if err!='': print(f'err: {err}'); exit()
   return msg
#---------------------------------------------------------------------------------------------------
# Set up plot resources
#---------------------------------------------------------------------------------------------------
num_case,num_var = len(case),len(var)

subtitle_font_height = 0.01

wkres = ngl.Resources()
npix = 2048; wkres.wkWidth,wkres.wkHeight=npix,npix
wks = ngl.open_wks(fig_type,fig_file,wkres)
# if add_diff:
#    plot = [None]*(num_var*(num_case+num_case-1))
# else:
#    plot = [None]*(num_var*num_case)

plot = [None]*(num_var*num_case)
if add_diff: 
   wks_diff = ngl.open_wks(fig_type,fig_file_diff,wkres)
   plot_diff = [None]*(num_var*(num_case-1))

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

res.trYMinF,res.trYMaxF,res.trXMinF,res.trXMaxF = 0.0,0.5,-15,15
# res.trYMinF,res.trYMaxF,res.trXMinF,res.trXMaxF = 0.0,0.5,-6,6


# res.cnFillPalette = np.array( cmocean.cm.rain(np.linspace(0,1,256)) )

lres = hs.res_xy()
lres.xyLineThicknessF = 1
lres.xyLineColor      = 'black'
lres.xyDashPattern    = 0

#---------------------------------------------------------------------------------------------------
# routine for adding theoretical dispersion curves and MJO frequency bounds
#---------------------------------------------------------------------------------------------------
# get data for dispersion curves
swfq, swwn = wf.genDispersionCurves() # swfq.shape # -->(6, 3, 50)
swf = np.ma.masked_invalid( np.where(swfq==1e20, np.nan, swfq) )
swk = np.ma.masked_invalid( np.where(swwn==1e20, np.nan, swwn) )
def add_dispersion_curves(wks_in,plot_in):
   for ii in range(3,6):
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
#---------------------------------------------------------------------------------------------------
spec_raw_list = []
for v in range(num_var):
   print(f'\n  var: '+hc.tcolor.GREEN+f'{var[v]}'+hc.tcolor.ENDC)
   spec_list, freq_list, wvnm_list = [],[],[]

   if lev_list[v]==None:
      # var_str = var[v]
      file_var_str = var[v]
   else:
      # var_str = f'{var[v]} {int(lev_list[v]/1e2)} mb'
      file_var_str = f'{var[v]}{int(lev_list[v]/1e2)}'
   # if var=="PRECT" : var_str = "Precipitation"

   for c in range(num_case):
      print(f'    case: '+hc.tcolor.CYAN+f'{case[c]}'+hc.tcolor.ENDC)

      tvar = var[v]
      # if var[v]=='OLR' and case[c]!='NOAA': tvar = 'FLNT'
      if var[v] in ['OLR','FLNT','FLUT']  and case[c]=='NOAA': tvar = 'olr'

      wk_tmp_file = f'{tmp_file_head}.{case[c]}.{file_var_str}.lat_{lat1}_{lat2}.seg_{segsize}.novr_{noverlap}.daily.nc'

      ##########################################################################
      ##########################################################################
      if recalculate_spectra :
         #----------------------------------------------------------------------
         # identify files
         if case[c]=='NOAA':
            file_list = ['/global/cfs/cdirs/m3312/whannah/obs_data/OLR/olr.day.mean.nc']
         else:
            file_path = f'{case_dir[c]}/{case[c]}/data_remap_{remap_grid}/*.eam.{htype}*'
            file_list = sorted(glob.glob(file_path))
         if file_list==[]:
            print()
            print(f'file_path: {file_path}')
            print(f'file_list: {file_list}')
            exit('ERROR: no files found')
         #----------------------------------------------------------------------
         # read the data
         ds = xr.open_mfdataset(file_list)
         data = ds[tvar]#.sel(plev=lev_list[v])
         #----------------------------------------------------------------------
         # subset the data based on year
         data = data.where( data['time.year']>=yr1, drop=True)
         data = data.where( data['time.year']<=yr2, drop=True)
         #----------------------------------------------------------------------
         # reduce to equatorial region
         mask = xr.DataArray( np.ones([len(data['lat']),len(data['lon'])],dtype=bool), dims=('lat','lon') )
         mask = mask & (data['lat']>=lat1) & (data['lat']<=lat2)
         data = data.where( mask, drop=True)
         #----------------------------------------------------------------------
         # make various other adjustments

         # adjust units
         if var[v]=='PRECT': data = data*86400.*1e3

         # Convert to daily mean
         data = data.resample(time='D').mean(dim='time')

         # normalize by standrad deviation to standardize units for raw spectrum
         # data = data / data.std()

         hc.print_stat(data,name=var[v],compact=True,indent=' '*6)

         data.load()
         #----------------------------------------------------------------------
         #----------------------------------------------------------------------
         # print(); print(data); exit()

         ### check for invalid values
         # num_nan = np.sum(np.isnan(data.values))
         # num_inf = np.sum(np.isinf(data.values))
         # if num_nan>0 or num_inf>0:
         #    print(f'  num_nan: {num_nan}')
         #    print(f'  num_inf: {num_inf}')
         #    exit()
         #----------------------------------------------------------------------
         # Calculate the WK spectra
         #----------------------------------------------------------------------
         print(f'      calculating spectra...')

         spec = wf.spacetime_power( data, segsize=segsize, noverlap=noverlap, spd=1, dosymmetries=True )

         # calculate background from average of symmetric & antisymmetric components
         spec_all = spec.mean(dim='component')
         bkgd_all = wf.smooth_wavefreq(spec_all,    kern=wf.simple_smooth_kernel(), nsmooth=50, freq_name='frequency')
         bkgd_sym = wf.smooth_wavefreq(spec[0,:,:], kern=wf.simple_smooth_kernel(), nsmooth=50, freq_name='frequency')
         bkgd_asm = wf.smooth_wavefreq(spec[1,:,:], kern=wf.simple_smooth_kernel(), nsmooth=50, freq_name='frequency')
         
         # separate components and throw away negative frequencies
         nfreq = len(spec['frequency'])
         rspec_sym  = spec[0,:,int(nfreq/2):]
         rspec_asy  = spec[1,:,int(nfreq/2):]
         bkgd_all = bkgd_all[:,int(nfreq/2):] 
         bkgd_sym = bkgd_sym[:,int(nfreq/2):] 
         bkgd_asm = bkgd_asm[:,int(nfreq/2):] 

         # collect everything into a dataset
         num_days = len(data['time'].values)
         wk_ds = xr.Dataset()
         wk_ds['num_days']  = num_days
         wk_ds['segsize']   = segsize
         wk_ds['noverlap']  = noverlap
         wk_ds['rspec_sym'] = rspec_sym
         wk_ds['rspec_asy'] = rspec_asy
         wk_ds['bkgd_all']  = bkgd_all
         wk_ds['bkgd_sym']  = bkgd_sym
         wk_ds['bkgd_asm']  = bkgd_asm
         #----------------------------------------------------------------------
         # Write to file 
         #----------------------------------------------------------------------
         print(f'      writing to file - {wk_tmp_file}')
         wk_ds.to_netcdf(path=wk_tmp_file,mode='w')
      else:
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

      print(f'      num_days: {num_days}')
      ##########################################################################
      ##########################################################################
      # for difference always use the unnormalized spectra
      if add_diff:
         if spec_type=='sym' : spec_raw_list.append(rspec_sym.transpose().values)
         if spec_type=='asym': spec_raw_list.append(rspec_asy.transpose().values)

      # normalize by the smoothed "background" spectra
      if normalize:
         nspec_sym = rspec_sym / bkgd_sym
         nspec_asy = rspec_asy / bkgd_asm

         if spec_type=='sym' : spec_list.append(nspec_sym.transpose().values)
         if spec_type=='asym': spec_list.append(nspec_asy.transpose().values)
         # if spec_type=='tot' : spec_list.append(nspec_sym.transpose().values)
      else:
         if spec_type=='sym' : spec_list.append(rspec_sym.transpose().values)
         if spec_type=='asym': spec_list.append(rspec_asy.transpose().values)
         # if spec_type=='tot' : spec_list.append(rspec_sym.transpose().values)

      freq_list.append(wk_ds['frequency'].values)
      wvnm_list.append(wk_ds['wavenumber'].values)

   #----------------------------------------------------------------------------
   # calculate limits for common color bar
   data_min = np.min([np.nanmin(d) for d in spec_list])
   data_max = np.max([np.nanmax(d) for d in spec_list])
   #----------------------------------------------------------------------------
   # Create plot
   for c in range(num_case):
      tres = copy.deepcopy(res)

      # tres.cnFillPalette = np.array( cmocean.cm.speed(np.linspace(0,1,256)) )
      # tres.cnFillPalette = np.array( cmocean.cm.rain(np.linspace(0,1,256)) )
      # tres.cnFillPalette = np.array( cmocean.cm.delta(np.linspace(0,1,256)) )
      # tres.cnFillPalette = np.array( cmocean.cm.balance(np.linspace(0,1,256)) )
      # tres.cnFillPalette = np.array( cmocean.cm.thermal(np.linspace(0,1,256)) )

      tres.cnLevelSelectionMode = 'ExplicitLevels'

      if normalize:
         # colors and levels from E3SM diags
         # tres.cnFillPalette = ['white', 'gainsboro','lightgray','gray','paleturquoise',
         #                       'skyblue','palegreen','mediumseagreen','seagreen','yellow',
         #                       'orange','red','maroon','pink','magenta']
         # tres.cnLevels = [ 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.7, 1.9, 2.1, 2.4, 2.7, 3.0, 3.3 ]
         tres.cnFillPalette = ['white', 'gray95', 'gray85','gray75',
                               'paleturquoise','skyblue','palegreen','mediumseagreen','seagreen',
                               'yellow','orange','red','maroon','pink','plum','purple']
         tres.cnLevels = [ 0.9,1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.7, 1.9, 2.1, 2.4, 2.7, 3.0, 3.5, 4.0 ]
         # tres.cnLevels = np.arange(1.0,3.5+0.3,0.3)
         # tres.cnLevels = np.array([5,6,7,8,9,10,12,14,16,18,20,24,28,34,40])/1e1
      else:
         # num_clev,aboutZero = 21,False
         # (cmin,cmax,cint) = ngl.nice_cntr_levels(data_min, data_max, 
         #                                  cint=None, max_steps=num_clev, 
         #                                  returnLevels=False, aboutZero=aboutZero )
         # tres.cnLevels = np.linspace(cmin,cmax,num=num_clev)
         
         # if var[v]=='OLR' : tres.cnLevels = np.linspace(0.01,5.01,30)
         # if spec_type=='sym'  and var[v]=='U': tres.cnLevels = np.linspace(0.1,4.8,10)
         # if spec_type=='asym' and var[v]=='U': tres.cnLevels = np.linspace(0.02,1.2,10)
         if spec_type=='sym'  and var[v]=='U': tres.cnLevels = np.linspace(0.01,4.01,40)
         if spec_type=='asym' and var[v]=='U': tres.cnLevels = np.linspace(0.01,1.01,20)

         # log scales
         # if var[v]=='U'   : tres.cnLevels = np.logspace( -2, 0.8, num=20).round(decimals=4)
         # if var[v]=='OLR ' : tres.cnLevels = np.logspace(  0.0, 2.0, num=30).round(decimals=4)
         # if var[v]=='PRECT': tres.cnLevels = np.logspace( -1.5, 0.0, num=30).round(decimals=4)
         # if var[v]=='U850' : tres.cnLevels = np.logspace( -1.5, 1.0, num=30).round(decimals=4)
         # if var[v]=='U200' : tres.cnLevels = np.logspace( -1.0, 2.0, num=30).round(decimals=4)

      #-------------------------------------------------------------------------
      # zero out values over region with satellite aliasing - helps to set difference colormap
      if case[c]=='NOAA':
         for i in range(len(wvnm_list[c])):
            twn = wvnm_list[c][i]
            if twn>=awn1 and  twn<=awn2:
               for j in range(len(freq_list[c])):
                  tfq = freq_list[c][j]
                  if tfq>=afq1 and  tfq<=afq2:
                     spec_list[c][j,i] = 0
      #-------------------------------------------------------------------------
      #-------------------------------------------------------------------------

      tres.sfXArray = wvnm_list[c]
      tres.sfYArray = freq_list[c]

      tres.tmYLMode     = 'Explicit'
      tres.tmYLValues   = 1./tm_period_values
      tres.tmYLLabels   = tm_period_values

      ip = v*num_case+c if var_x_case else c*num_var+v

      plot[ip] = ngl.contour(wks, spec_list[c], tres)
      add_dispersion_curves(wks,plot[ip])
      if case[c]=='NOAA': hide_aliasing_artifact(wks,plot[ip])
      
      hs.set_subtitles(wks, plot[ip], 
                       left_string=name[c], 
                       center_string='', 
                       right_string=vstr[v], 
                       font_height=subtitle_font_height)

#---------------------------------------------------------------------------------------------------
# Finalize plot
#---------------------------------------------------------------------------------------------------
# if add_diff:
#    # layout = [num_var*2,num_case] if var_x_case else [num_case,num_var*2]
#    layout = [num_var,num_case+(num_case-1)] if var_x_case else [num_case+(num_case-1),num_var]
# else:
#    layout = [num_var,num_case] if var_x_case else [num_case,num_var]

layout = [num_var,num_case] if var_x_case else [num_case,num_var]

pnl_res = hs.setres_panel()
pnl_res.nglPanelYWhiteSpacePercent = 5
pnl_res.nglPanelXWhiteSpacePercent = 5
pnl_res.nglPanelFigureStrings            = list(string.ascii_lowercase)
pnl_res.nglPanelFigureStringsJust        = "TopLeft"
pnl_res.nglPanelFigureStringsFontHeightF = 0.015

if use_common_label_bar: 
   pnl_res.nglPanelLabelBar = True
   pnl_res.nglPanelLabelBarLabelFontHeightF = 0.01
   pnl_res.nglPanelLabelBarWidthF = 0.5
   # pnl_res.lbLeftMarginF      = 0.3
   # pnl_res.lbRightMarginF     = 0.3

ngl.panel(wks,plot,layout,pnl_res)
hc.trim_png(fig_file)

#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------
if add_diff:
   # diff_min = np.min([np.nanmin(d-spec_raw_list[0]) for d in spec_raw_list])
   # diff_max = np.max([np.nanmax(d-spec_raw_list[0]) for d in spec_raw_list])
   # print(); print(f'diff_min: {data_min:.2f}'); print(f'diff_max: {data_max:.2f}');print()
   for v in range(num_var):
      
      ib = v*num_case+0
      ie = v*num_case+num_case
      # print()
      # print(f'ib: {ib}')
      # print(f'ie: {ie}')
      diff_min = np.min([np.nanmin(d-spec_raw_list[ib]) for d in spec_raw_list[ib+1:ie]])
      diff_max = np.max([np.nanmax(d-spec_raw_list[ib]) for d in spec_raw_list[ib+1:ie]])

      # print(); print(f'diff_min: {diff_min:.2f}'); print(f'diff_max: {diff_max:.2f}');print()

      # if var[v]=='PRECT': diff_max = 0.2 ; diff_min = -1*diff_max
      # if var[v]=='FLUT' : diff_max = 12. ; diff_min = -1*diff_max
      # if var[v]=='U850' : diff_max = 1.0 ; diff_min = -1*diff_max

      for c in range(num_case):
         tres = copy.deepcopy(res)
         tres.sfXArray = wvnm_list[c]
         tres.sfYArray = freq_list[c]
         tres.tmYLMode     = 'Explicit'
         tres.tmYLValues   = 1./tm_period_values
         tres.tmYLLabels   = tm_period_values
         tres.lbLabelBarOn = True

         tres.lbLeftMarginF      = -0.2
         tres.lbRightMarginF     = -0.2

         tres.cnFillPalette = np.array( cmocean.cm.balance(np.linspace(0,1,256)) )

         num_clev = 8
         (cmin,cmax,cint) = ngl.nice_cntr_levels(diff_min, diff_max, 
                                    cint=None, max_steps=num_clev, 
                                    returnLevels=False, aboutZero=True )
         tres.cnLevels = np.linspace(cmin,cmax,num=num_clev)
         tres.cnLevelSelectionMode = 'ExplicitLevels'

         if var[v]=='PRECT': dc=0.08;nc=3;lim=nc*dc+dc/2;clev_diff=np.arange(-1*lim,lim+dc,dc)
         if var[v]=='FLUT' : dc=6.0 ;nc=3;lim=nc*dc+dc/2;clev_diff=np.arange(-1*lim,lim+dc,dc)
         if var[v]=='U850' : dc=0.16;nc=3;lim=nc*dc+dc/2;clev_diff=np.arange(-1*lim,lim+dc,dc)

         if 'clev_diff' in locals():
            if c==0: print(f'  clev_diff: {clev_diff}')
            tres.cnLevels = clev_diff
            del clev_diff

         # tres.cnLevels = np.arange(0.22,0.22+0.04,0.04)
         # np.linspace(-12,12,11)
         # np.array([-11,-9,-7,-5,-3,-1,1,3,5,7,9,11])

         # 0.2 
         # 12. 
         # 1.0 


         # ip = v*(num_case+num_case-1)+c if var_x_case else (num_case+c-1)*num_var+v
         ip = v

         if c>0:
            ib = v*num_case+0
            ii = v*num_case+c
            # print(f'  ib: {ib}')
            # print(f'  ii: {ii}')
            nspec_diff = spec_raw_list[ii] - spec_raw_list[0]
            plot_diff[ip] = ngl.contour(wks_diff, nspec_diff, tres)
            add_dispersion_curves(wks_diff,plot_diff[ip])
            hide_aliasing_artifact(wks_diff,plot_diff[ip])
            hs.set_subtitles(wks_diff, plot_diff[ip], 
                             left_string=name[c], 
                             # center_string='(diff)', 
                             center_string='',
                             right_string=vstr[v], 
                             center_sub_string='Un-Normalized Difference',
                             font_height=subtitle_font_height)
   #----------------------------------------------------------------------------
   layout = [num_var,num_case] if var_x_case else [num_case,num_var]

   pnl_res = hs.setres_panel()
   pnl_res.nglPanelYWhiteSpacePercent       = 5
   pnl_res.nglPanelXWhiteSpacePercent       = 5
   pnl_res.nglPanelFigureStrings            = list(string.ascii_lowercase[num_case*num_var:])
   pnl_res.nglPanelFigureStringsJust        = "TopLeft"
   pnl_res.nglPanelFigureStringsFontHeightF = 0.015

   # if use_common_label_bar: 
   #    pnl_res.nglPanelLabelBar = True
   #    pnl_res.lbLeftMarginF      = -0.3
   #    pnl_res.lbRightMarginF     = -0.3

   ngl.panel(wks_diff,plot_diff,layout,pnl_res)
   hc.trim_png(fig_file_diff)
   #----------------------------------------------------------------------------
   for fig_file_tmp in [fig_file,fig_file_diff]:
      cmd = f"convert -trim  {fig_file_tmp}.{fig_type} -format '%wx%h\n' info:"
      (msg, err) = sp.Popen(cmd, stdout=sp.PIPE, stderr=sp.PIPE, shell=True, universal_newlines=True).communicate()
      dim_list = np.array( [ int(d) for d in msg.replace('\n','').split('x') ] )
      # max_dim = np.max(dim_list)
      # nx,ny = dim_list[0],dim_list[1]+100
      nx,ny = dim_list[0],dim_list[1]+50
      run_cmd(f'convert -trim -gravity center -extent {nx}x{ny} {fig_file_tmp}.{fig_type} {fig_file_tmp}.{fig_type}')
   #----------------------------------------------------------------------------
   # run_cmd(f'montage {fig_file}.{fig_type} {fig_file_diff}.{fig_type} -geometry +100 -tile 1x2  -gravity center {fig_file}.{fig_type}')
   run_cmd(f'montage {fig_file}.{fig_type} {fig_file_diff}.{fig_type} -geometry +10 -tile 1x2 {fig_file}.{fig_type}')
   hc.trim_png(fig_file)
#---------------------------------------------------------------------------------------------------
ngl.end()
#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------

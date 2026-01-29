import os, ngl, subprocess as sp, numpy as np, xarray as xr, copy, string, sys
import hapy_common as hc, hapy_E3SM   as he, hapy_setres as hs
sys.path.insert(1, os.getenv('HOME')+'/wavenumber_frequency-master')
import wavenumber_frequency_functions as wf
import cmocean
home = os.getenv('HOME')
#---------------------------------------------------------------------------------------------------
case,name,case_dir,case_sub = [],[],[],[]
clr,dsh,mrk = [],[],[]
def add_case(case_in,n=None,c='black',d=0,m=1,p=None,s=None):
   global name,case,clr,dsh,mrk
   if n is None: n = '' 
   case.append(case_in); name.append(n)
   case_dir.append(p); case_sub.append(s)
   clr.append(c); dsh.append(d); mrk.append(m)
#---------------------------------------------------------------------------------------------------
var,lev_list,var_str = [],[],[]
def add_var(vname,vstr=None,lev=-1): 
   if vstr is None: vstr = vname
   var.append(vname); lev_list.append(lev); var_str.append(vstr)
#-------------------------------------------------------------------------------
# add_case('NOAA',n='NOAA OLR',d=f'{home}/Data/Obs/OLR',s='')
### extended runs
tmp_path,tmp_sub = '/pscratch/sd/w/whannah/e3sm_scratch/pm-gpu','archive/atm/hist'
add_case('E3SM.2023-RAD-SENS-00.GNUGPU.ne30pg2_oECv3.F2010-MMF1.NXY_64x1.RNX_1.MCICA_OFF',   n='nx/rx= 64/ 1 mcica OFF',c='red',      d=1,p=tmp_path,s=tmp_sub)
add_case('E3SM.2023-RAD-SENS-00.GNUGPU.ne30pg2_oECv3.F2010-MMF1.NXY_64x1.RNX_8.MCICA_OFF',   n='nx/rx= 64/ 8 mcica OFF',c='palegreen',d=1,p=tmp_path,s=tmp_sub)
add_case('E3SM.2023-RAD-SENS-00.GNUGPU.ne30pg2_oECv3.F2010-MMF1.NXY_64x1.RNX_64.MCICA_OFF',  n='nx/rx= 64/64 mcica OFF',c='purple',   d=1,p=tmp_path,s=tmp_sub)
# add_case('E3SM.2023-RAD-SENS-00.GNUGPU.ne30pg2_oECv3.F2010-MMF1.NXY_64x1.RNX_1.MCICA_ON',    n='nx/rx= 64/ 1 mcica ON ',c='red',      d=0,p=tmp_path,s=tmp_sub)
# add_case('E3SM.2023-RAD-SENS-00.GNUGPU.ne30pg2_oECv3.F2010-MMF1.NXY_64x1.RNX_8.MCICA_ON',    n='nx/rx= 64/ 8 mcica ON ',c='palegreen',d=0,p=tmp_path,s=tmp_sub)
# add_case('E3SM.2023-RAD-SENS-00.GNUGPU.ne30pg2_oECv3.F2010-MMF1.NXY_64x1.RNX_64.MCICA_ON',   n='nx/rx= 64/64 mcica ON ',c='purple',   d=0,p=tmp_path,s=tmp_sub)
#-------------------------------------------------------------------------------

add_var('FLNT','OLR')
# add_var('PRECT')
# add_var('U850')
# add_var('U200')
# add_var('U',lev=850)

fig_type = 'png'
fig_file = 'figs/FXX-wk-wave-spectra'


lat1,lat2 = -15,15

num_year = 5

htype,first_file,num_files = 'h1',0,365*num_year

use_remap,remap_grid = True,'90x180'
# use_remap,remap_grid = True,'73x144'

recalculate_spectra = True

use_daily   = True
normalize   = True
add_diff    = True
var_x_case  = True

num_plot_col = 2

#---------------------------------------------------------------------------------------------------
# Set up plot resources
#---------------------------------------------------------------------------------------------------
num_case,num_var = len(case),len(var)

wkres = ngl.Resources()
npix = 2048; wkres.wkWidth,wkres.wkHeight=npix,npix
wks = ngl.open_wks(fig_type,fig_file,wkres)
# plot = []
if add_diff:
   # plot = [None]*(num_var*(num_case*2-1))
   plot = [None]*(num_var*(num_case*2))
   # if 'NOAA' in case:
   #    plot = [None]*(num_var*(num_case*2-1-1))
   # else:
   #    plot = [None]*(num_var*(num_case*2-1))
else:
   plot = [None]*(num_var*num_case)

res = hs.res_contour_fill()
# res.vpHeightF = 0.4
res.lbLabelFontHeightF = 0.015
res.tiXAxisFontHeightF = 0.02
res.tiYAxisFontHeightF = 0.02

res.tiXAxisString = 'Zonal Wavenumber'
res.tiYAxisString = 'Frequency (cpd)'

res.tiXAxisString = 'Zonal Wavenumber'
res.tiYAxisString = 'Frequency (cpd)'

# res.cnFillPalette = np.array( cmocean.cm.rain(np.linspace(0,1,256)) )

lres = hs.res_xy()
lres.xyLineThicknessF = 1
lres.xyLineColor      = 'black'
lres.xyDashPattern    = 0

#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------

# get data for dispersion curves
swfq, swwn = wf.genDispersionCurves() # swfq.shape # -->(6, 3, 50)
swf = np.ma.masked_invalid( np.where(swfq==1e20, np.nan, swfq) )
swk = np.ma.masked_invalid( np.where(swwn==1e20, np.nan, swwn) )

#---------------------------------------------------------------------------------------------------
# routine for adding theoretical dispersion curves and MJO frequency bounds
#---------------------------------------------------------------------------------------------------
def add_dispersion_curves(plot):
   # add dispersion curves
   for ii in range(3,6):
      lres.xyDashPattern = 0
      ngl.overlay( plot, ngl.xy(wks, np.array([0,0]), np.array([-1e3,1e3]), lres) )
      ngl.overlay( plot, ngl.xy(wks, swk[ii, 0,:], swf[ii,0,:], lres) )
      ngl.overlay( plot, ngl.xy(wks, swk[ii, 1,:], swf[ii,1,:], lres) )
      ngl.overlay( plot, ngl.xy(wks, swk[ii, 2,:], swf[ii,2,:], lres) )
   # overlay lines for MJO frequency range
   lres.xyDashPattern = 1
   tfrq=1./30.; ngl.overlay( plot, ngl.xy(wks, np.array([-1e3,1e3]), np.array([tfrq,tfrq]), lres) )
   tfrq=1./90.; ngl.overlay( plot, ngl.xy(wks, np.array([-1e3,1e3]), np.array([tfrq,tfrq]), lres) )
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

def hide_aliasing_artifact(plot):
   bx = [awn1,awn1,awn2,awn2,awn1]
   by = [afq1,afq2,afq2,afq1,afq1]
   pdum = ngl.add_polygon(wks,plot,bx,by,pgres)
#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------
for v in range(num_var):
   print(f'\n  var: '+hc.tcolor.GREEN+f'{var[v]}'+hc.tcolor.ENDC)
   spec_list, freq_list, wvnm_list = [],[],[]

   for c in range(num_case):
      print(f'    case: '+hc.tcolor.CYAN+f'{case[c]}'+hc.tcolor.ENDC)

      tvar = var[v]
      if var[v]=='FLNT' and case[c]=='NOAA': tvar = 'OLR'

      use_remap_tmp = use_remap
      if case[c]=='NOAA': use_remap_tmp = False

      data_dir_tmp,data_sub_tmp = None, None
      if use_remap_tmp: 
         data_sub_tmp = f'data_remap_{remap_grid}/'
      else:
         if case_sub[c] is not None: data_sub_tmp = case_sub[c]

      if case_dir[c] is not None: data_dir_tmp = case_dir[c]

      # wk_tmp_file = case_obj.data_dir+'/'+case_obj.name+f'/wk.v1.{var[v]}.{lat1}:{lat2}.nc'
      wk_tmp_file = f'data/wk-wave-spectra.{case[c]}.{var[v]}.lat_{lat1}_{lat2}.nyr_{num_year}.nc'

      if use_daily: wk_tmp_file = wk_tmp_file.replace('.nc','.dailiy.nc')

      ##########################################################################
      ##########################################################################
      if recalculate_spectra :
         case_obj = he.Case( name=case[c], data_dir=data_dir_tmp, data_sub=data_sub_tmp )
         case_obj.set_coord_names(tvar)
         #----------------------------------------------------------------------
         # read the data
         #----------------------------------------------------------------------
         print(f'      loading data...')

         data = case_obj.load_data(tvar,htype='h1',first_file=first_file,num_files=num_files
                                  ,use_remap=use_remap_tmp,remap_str=f'remap_{remap_grid}')

         if case[c]=='NOAA':
            yr1 = 2005 ; yr2 = yr1+num_year-1
            date_beg = f'{yr1}-01-01'
            date_end = f'{yr2}-12-31'
            data = data.sel(time=slice(date_beg,date_end))

         # reduce to equatorial subset
         mask = xr.DataArray( np.ones([len(data['lat']),len(data['lon'])],dtype=bool), dims=('lat','lon') )
         mask = mask & (data['lat']>=lat1) & (data['lat']<=lat2)
         data = data.where( mask, drop=True)

         # Convert to daily mean
         if use_daily: data = data.resample(time='D').mean(dim='time')

         # calculate samples per day for FFT
         dtime = ( data['time'][1]-data['time'][0] ).values / np.timedelta64(1, 's')
         spd = int( 86400./dtime )

         hc.print_stat(data,name=var[v],compact=True,indent=' '*6)

         data.load()

         num_days = len(data['time'].values)
         #----------------------------------------------------------------------
         # Calculate the WK spectra
         #----------------------------------------------------------------------
         print(f'      calculating spectra...')

         # segsize,noverlap = 180*spd,60*spd
         segsize,noverlap = 96*spd,60*spd

         spec = wf.spacetime_power( data, segsize=segsize, noverlap=noverlap, 
                                    spd=spd, dosymmetries=True )

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
         wk_ds = xr.Dataset()
         wk_ds['num_days']  = num_days
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
         rspec_sym  = wk_ds['rspec_sym']
         rspec_asy  = wk_ds['rspec_asy']
         bkgd_all   = wk_ds['bkgd_all']
         bkgd_sym   = wk_ds['bkgd_sym']
         bkgd_asm   = wk_ds['bkgd_asm']

      print(f'      num_days: {num_days}')
      ##########################################################################
      ##########################################################################

      # normalize by the smoothed "background" spectra
      if normalize:
         nspec_sym = rspec_sym / bkgd_sym
         nspec_asy = rspec_asy / bkgd_asm

         spec_list.append(nspec_sym.transpose().values)
      else:
         spec_list.append(rspec_sym.transpose().values)

      freq_list.append(wk_ds['frequency'].values)
      wvnm_list.append(wk_ds['wavenumber'].values)

   #----------------------------------------------------------------------------
   # calculate limits for common color bar
   #----------------------------------------------------------------------------

   data_min = np.min([np.nanmin(d) for d in spec_list])
   data_max = np.max([np.nanmax(d) for d in spec_list])

   print()
   print(f'    data_min: {data_min:.2f}')
   print(f'    data_max: {data_max:.2f}')
   print()
   
   if add_diff:
      # spec_diff_list = np.zeros(len(spec_list))
      # for c in range(num_case):
      #    spec_diff_list[c] = spec_list[c] - spec_list[0]
      diff_min = np.min([np.nanmin(d-spec_list[0]) for d in spec_list])
      diff_max = np.max([np.nanmax(d-spec_list[0]) for d in spec_list])

   #----------------------------------------------------------------------------
   # Create plot
   #----------------------------------------------------------------------------
   for c in range(num_case):
      #-------------------------------------------------------------------------
      #-------------------------------------------------------------------------
      tres = copy.deepcopy(res)

      tres.cnLevelSelectionMode = 'ExplicitLevels'

      if normalize: 
         # tres.cnLevels = np.arange(4,30+2,2)/1e1
         tres.cnLevels = np.array([3,4,5,6,7,8,9,10,11,12,14,17,20,24,28])/1e1
      else:
         # num_clev,aboutZero = 21,False
         # (cmin,cmax,cint) = ngl.nice_cntr_levels(data_min, data_max, 
         #                                  cint=None, max_steps=num_clev, 
         #                                  returnLevels=False, aboutZero=aboutZero )
         # tres.cnLevels = np.linspace(cmin,cmax,num=num_clev)
         tres.cnLevels = np.arange(2,30+2,2)

         # # log scales
         # if var[v]=='OLR ' : tres.cnLevels = np.logspace(  0.0, 2.0, num=30).round(decimals=4)
         # if var[v]=='PRECT': tres.cnLevels = np.logspace( -1.5, 0.0, num=30).round(decimals=4)
         # if var[v]=='U850' : tres.cnLevels = np.logspace( -1.5, 1.0, num=30).round(decimals=4)
         # if var[v]=='U200' : tres.cnLevels = np.logspace( -1.0, 2.0, num=30).round(decimals=4)

      #-------------------------------------------------------------------------
      # zero out values over region with satellite aliasing - helps to set difference colormap
      #-------------------------------------------------------------------------
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
      # tres.cnFillPalette = np.array( cmocean.cm.speed(np.linspace(0,1,256)) )
      # tres.cnFillPalette = np.array( cmocean.cm.rain(np.linspace(0,1,256)) )
      # tres.cnFillPalette = np.array( cmocean.cm.delta(np.linspace(0,1,256)) )
      # tres.cnFillPalette = np.array( cmocean.cm.balance(np.linspace(0,1,256)) )
      # tres.cnFillPalette = np.array( cmocean.cm.thermal(np.linspace(0,1,256)) )

      tres.sfXArray = wvnm_list[c]
      tres.sfYArray = freq_list[c]
      # tres.trYMaxF = 0.25
      # tres.trXMinF,tres.trXMaxF = -15,15
      tres.trYMaxF = 0.4
      tres.trXMinF,tres.trXMaxF = -15,15

      if add_diff:
         ip = v*(num_case*2)+c
      else:
         ip = v*num_case+c if var_x_case else c*num_var+v

      plot[ip] = ngl.contour(wks, spec_list[c], tres)
      add_dispersion_curves(plot[ip])
      if case[c]=='NOAA': hide_aliasing_artifact(plot[ip])
      hs.set_subtitles(wks, plot[ip], name[c], '', var_str[v], font_height=0.01)

      #-------------------------------------------------------------------------
      # Add difference plot
      #-------------------------------------------------------------------------
      if add_diff:
         # num_clev,aboutZero = 21,True
         # (cmin,cmax,cint) = ngl.nice_cntr_levels(diff_min, diff_max, 
         #                            cint=None, max_steps=num_clev, 
         #                            returnLevels=False, aboutZero=aboutZero )
         # tres.cnLevels = np.linspace(cmin,cmax,num=num_clev)

         tres.cnLevels = np.arange(-1,1+0.2,0.2)

         tres.cnFillPalette = np.array( cmocean.cm.balance(np.linspace(0,1,256)) )

         ip = v*(num_case*2)+num_case+c

         if c==0:
            bres = hs.res_default()
            bres.pmTickMarkDisplayMode = 'Never'
            plot[ip] = ngl.blank_plot(wks,bres)
         else:
            nspec_diff = spec_list[c] - spec_list[0]
            plot[ip] = ngl.contour(wks, nspec_diff, tres)
            add_dispersion_curves(plot[ip])
            hide_aliasing_artifact(plot[ip])
            hs.set_subtitles(wks, plot[ip], name[c], '', var_str[v]+' diff', font_height=0.01)

#---------------------------------------------------------------------------------------------------
# Finalize plot
#---------------------------------------------------------------------------------------------------
if add_diff:
   layout = [num_var*2,num_case] if var_x_case else [num_case,num_var*2]
else:
   layout = [num_var,num_case] if var_x_case else [num_case,num_var]

pnl_res = hs.setres_panel()
pnl_res.nglPanelYWhiteSpacePercent = 1
pnl_res.nglPanelXWhiteSpacePercent = 1

# pnl_res.nglPanelFigureStrings            = list(string.ascii_lowercase)
pnl_res.nglPanelFigureStrings            = ['a','b','c','d','','e','f','g']
pnl_res.nglPanelFigureStringsJust        = "TopLeft"
pnl_res.nglPanelFigureStringsFontHeightF = 0.01

ngl.panel(wks,plot,layout,pnl_res)
ngl.end()

hc.trim_png(fig_file)
#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------
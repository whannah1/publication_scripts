import os, glob, ngl, subprocess as sp, numpy as np, xarray as xr, copy, string, sys
import hapy_common as hc, hapy_E3SM   as he, hapy_setres as hs
sys.path.insert(1, os.getenv('HOME')+'/wavenumber_frequency-master')
import wavenumber_frequency_functions as wf
import cmocean
scratch = '/global/cfs/cdirs/m3312/whannah/2022-QBO-TEST'
#-------------------------------------------------------------------------------
case_dir,case_sub = [],[]
case,name = [],[]
clr,dsh,mrk = [],[],[]
def add_case(case_in,n=None,p=None,s=None,d=0,c='black',m=0):
   global name,case,case_dir,case_sub,clr,dsh,mrk
   case.append(case_in); name.append(n)
   case_dir.append(p); case_sub.append(s)
   dsh.append(d) ; clr.append(c) ; mrk.append(m)
var,lev_list,vclr,vdsh = [],[],[],[]
def add_var(var_name,lev=-1,c='black',d=0): 
   var.append(var_name); lev_list.append(lev)
   vclr.append(c); vdsh.append(d)
#-------------------------------------------------------------------------------
add_case('E3SM.QBO-TEST.F2010.ne30pg2.L72.01',       n='L72 control')
add_case('E3SM.QBO-TEST.F2010.ne30pg2.L72-nsu40.01', n='L72 smooth')
add_case('E3SM.QBO-TEST.F2010.ne30pg2.L80-rsu40.01', n='L80')
#-------------------------------------------------------------------------------

add_var('U',lev=400e2)
add_var('U',lev=200e2)
add_var('U',lev=100e2)
add_var('U',lev=70e2)
add_var('U',lev=50e2)
add_var('U',lev=20e2)
add_var('U',lev=10e2)

phase_index_lev = 20e2

# decide what to plot (doesn't affect how things get calculated)
spec_type = 'sym' # sym / asym / tot
plot_phase = 1 # 1=westerly / 2=easterly

fig_type = 'png'
fig_file = f'figs/FXX-wk-spectra-by-phase.{spec_type}.phase_{plot_phase}'

lat1,lat2 = -5,5

first_file,num_files = 0,365*10
# num_files = 180
htype = 'h2'

use_remap,remap_grid = True,'90x180'
# use_remap,remap_grid = True,'73x144'

recalculate_spectra = True

add_diff    = True
var_x_case  = True

# num_plot_col = 2

seg_len     = 90
seg_cnt_min = 4

#---------------------------------------------------------------------------------------------------
# Set up plot resources
#---------------------------------------------------------------------------------------------------
num_case,num_var = len(case),len(var)

subtitle_font_height = 0.006

wkres = ngl.Resources()
# npix = 2048; wkres.wkWidth,wkres.wkHeight=npix,npix
npix = 4096; wkres.wkWidth,wkres.wkHeight=npix,npix
wks = ngl.open_wks(fig_type,fig_file,wkres)
if add_diff:
   plot = [None]*(num_var*(num_case+num_case-1))
else:
   plot = [None]*(num_var*num_case)

res = hs.res_contour_fill()
res.vpHeightF = 0.5
res.lbLabelFontHeightF = 0.015
res.tiXAxisFontHeightF = 0.02
res.tiYAxisFontHeightF = 0.02
res.tiXAxisString = 'Zonal Wavenumber'
res.tiYAxisString = 'Frequency (cpd)'

res.trYMinF,res.trYMaxF = 0.0,0.2
res.trXMinF,res.trXMaxF = -6,6

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
def add_dispersion_curves(plot):
   for ii in range(3,6):
      lres.xyDashPattern = 0
      ngl.overlay( plot, ngl.xy(wks, np.array([0,0]), np.array([-1e3,1e3]), lres) )
      ngl.overlay( plot, ngl.xy(wks, swk[ii, 0,:], swf[ii,0,:], lres) )
      ngl.overlay( plot, ngl.xy(wks, swk[ii, 1,:], swf[ii,1,:], lres) )
      ngl.overlay( plot, ngl.xy(wks, swk[ii, 2,:], swf[ii,2,:], lres) )
   # # overlay lines for MJO frequency range
   # lres.xyDashPattern = 1
   # tfrq=1./30.; ngl.overlay( plot, ngl.xy(wks, np.array([-1e3,1e3]), np.array([tfrq,tfrq]), lres) )
   # tfrq=1./90.; ngl.overlay( plot, ngl.xy(wks, np.array([-1e3,1e3]), np.array([tfrq,tfrq]), lres) )
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
   var_msg = f'  var: '+hc.tcolor.GREEN+f'{var[v]}'+hc.tcolor.ENDC
   if lev_list[v]!=None: var_msg += f'{var_msg} ({int(lev_list[v]/1e2)} mb)'
   print(); print(var_msg)

   spec_list, freq_list, wvnm_list = [],[],[]

   if lev_list[v]==None:
      var_str = var[v]
   else:
      var_str = f'{var[v]} {int(lev_list[v]/1e2)} mb'
      file_var_str = f'{var[v]}{int(lev_list[v]/1e2)}'
   # if var=="PRECT" : var_str = "Precipitation"

   for c in range(num_case):
      print(' '*4+f'case: '+hc.tcolor.CYAN+f'{case[c]}'+hc.tcolor.ENDC)

      tvar = var[v]
      if var[v]=='OLR': tvar = 'FLNT'

      use_remap_tmp = use_remap

      data_dir_tmp,data_sub_tmp = None, None
      if use_remap_tmp: data_sub_tmp = f'data_remap_{remap_grid}/'
      if case_dir[c] is not None: data_dir_tmp = case_dir[c]
      if case_sub[c] is not None: data_sub_tmp = case_sub[c]

      wk_tmp_file = f'data/wk-wave-spectra-by-phase.{case[c]}.{file_var_str}.lat_{lat1}_{lat2}.nc'

      ##########################################################################
      ##########################################################################
      if recalculate_spectra :
         #----------------------------------------------------------------------
         # read the data
         #----------------------------------------------------------------------
         print(' '*6+'loading data...')

         file_path = f'{scratch}/{case[c]}/data_remap_{remap_grid}_prs/*.eam.{htype}*'
         file_list = sorted(glob.glob(file_path))



         if file_list==[]:
            print()
            print(f'file_path: {file_path}')
            print(f'file_list: {file_list}')
            exit('ERROR: no files found')
         else:
            print()
            print(' '*6+f'file_path: {file_path}')

         
         if 'first_file' in locals(): file_list = file_list[first_file:]
         if 'num_files'  in locals(): file_list = file_list[:num_files]
         

         ds = xr.open_mfdataset(file_list)
         area = ds['area'].isel(time=0)

         

         # mask for equatorial subset
         mask = xr.DataArray( np.ones([len(ds['lat']),len(ds['lon'])],dtype=bool), dims=('lat','lon') )
         mask = mask & (ds['lat']>=lat1) & (ds['lat']<=lat2)

         time_dummy = ds['time']

         # data = ds[tvar].sel(plev=lev_list[v]).where( mask, drop=True)
         # data = data.resample(time='D').mean(dim='time') # Convert to daily mean

         # hc.print_stat(data,name=var[v],compact=True,indent=' '*6)

         # calculate mean wind as a baseline to identify QBO phase
         clim = ds[tvar].sel(plev=phase_index_lev).where( mask, drop=True).mean(dim='time')
         # clim = ds[tvar].sel(plev=lev_list[v]).where( mask, drop=True).mean(dim='time')
         clim = ( (clim*area).sum(dim=('lat','lon')) / area.sum(dim=('lat','lon')) ).mean().values

         # also calculate the standard deviation to delineate phases?
         # clim_std = ds[tvar].sel(plev=lev_list[v]).where( mask, drop=True).std(dim='time')

         pthrshld = 0.1

         print(); print(' '*8+f'clim: {clim}')

         #----------------------------------------------------------------------
         # Calculate the WK spectra
         #----------------------------------------------------------------------
         print(); print(f' '*6+'splitting data into segments...')

         nday = len(ds['time'].values)
         # if seg_len>nday: seg_len = nday # reset to 1 segment for debugging
         nseg = int( nday / seg_len)

         print()
         print(' '*8+f'nday_total: {nday}')
         print(' '*8+f'nseg (max): {nseg}')
         print()

         # stitch together separate datasets for each QBO phase
         cnt1,cnt2 = 0,0
         for n in range(nseg):
            t1,t2 = n*seg_len, (n+1)*seg_len
            
            data_tmp = ds[tvar].sel(plev=lev_list[v]).where( mask, drop=True)[t1:t2,:,:]
            data_tmp = data_tmp.resample(time='D').mean(dim='time') # Convert to daily mean

            # data_tmp = data[t1:t2,:,:].copy()
            data_tmp.load()
            glb_mean = ( (data_tmp*area).sum(dim=('lat','lon')) / area.sum(dim=('lat','lon')) ).mean()
            # print(' '*8+f'{n:03}  glb_mean: {glb_mean.values}')
            if glb_mean>(clim+pthrshld): data_p1 = data_tmp.copy() if cnt1==0 else xr.concat([data_p1,data_tmp],'time'); cnt1+=1
            if glb_mean<(clim-pthrshld): data_p2 = data_tmp.copy() if cnt2==0 else xr.concat([data_p2,data_tmp],'time'); cnt2+=1
         

         if cnt1>=seg_cnt_min:
            # print(); print(data_p1)
            # print(); print(data_p1.time.values)
            data_p1['time'] = time_dummy[:len(data_p1.time)]

         if cnt2>=seg_cnt_min:
            # print(); print(data_p2)
            # print(); print(data_p2.time.values)
            data_p2['time'] = time_dummy[:len(data_p2.time)]

         # exit()

         #----------------------------------------------------------------------
         # Calculate the WK spectra
         #----------------------------------------------------------------------
         if cnt1<seg_cnt_min and cnt2<seg_cnt_min: 
            print(); exit(f'ERROR: the count for both phases is less than seg_cnt_min ({seg_cnt_min})')

         if cnt1>=seg_cnt_min:
            print(f' '*6+'calculating phase 1 (westerly) spectra...')
            print(); print(f' '*8+f'PHASE 1 cnt: {cnt1}')
            print(); print(f' '*10+f'data_p1 shape: {data_p1.shape()}'); print()
            spec_p1 = wf.spacetime_power(data_p1,segsize=seg_len,noverlap=1,spd=1,dosymmetries=True)
            # separate components and throw away negative frequencies
            nfreq = len(spec_p1['frequency'])
            spec_p1_sym,spec_p1_asy  = spec_p1[0,:,int(nfreq/2):], spec_p1[1,:,int(nfreq/2):]

         if cnt2>=seg_cnt_min:
            print(f' '*6+'calculating phase 2 (easterly) spectra...')
            print(); print(' '*8+f'PHASE 2 cnt: {cnt2}')
            # print(); print(data_p2)
            spec_p2 = wf.spacetime_power(data_p2,segsize=seg_len,noverlap=1,spd=1,dosymmetries=True)
            # separate components and throw away negative frequencies
            nfreq = len(spec_p2['frequency'])
            spec_p2_sym,spec_p2_asy  = spec_p2[0,:,int(nfreq/2):], spec_p2[1,:,int(nfreq/2):]

         # collect everything into a dataset
         wk_ds = xr.Dataset()
         wk_ds['nday_tot']    = nday
         wk_ds['seg_len']     = seg_len
         wk_ds['seg_cnt1']    = cnt1
         wk_ds['seg_cnt2']    = cnt2
         if cnt1>=seg_cnt_min:
            wk_ds['spec_p1_sym'] = spec_p1_sym
            wk_ds['spec_p1_asy'] = spec_p1_asy
         if cnt2>=seg_cnt_min:
            wk_ds['spec_p2_sym'] = spec_p2_sym
            wk_ds['spec_p2_asy'] = spec_p2_asy
         #----------------------------------------------------------------------
         # Write to file 
         #----------------------------------------------------------------------
         print(f'      writing to file - {wk_tmp_file}')
         wk_ds.to_netcdf(path=wk_tmp_file,mode='w')
      else:
         print(f'      loading pre-calculated spectra... {wk_tmp_file}')
         wk_ds = xr.open_mfdataset( wk_tmp_file )
         nday      = wk_ds['nday_tot']
         seg_len   = wk_ds['seg_len']
         cnt1      = wk_ds['seg_cnt1'].values
         cnt2      = wk_ds['seg_cnt2'].values
         if cnt1>=seg_cnt_min:
            spec_p1_sym = wk_ds['spec_p1_sym']
            spec_p1_asy = wk_ds['spec_p1_asy']
         if cnt2>=seg_cnt_min:
            spec_p2_sym = wk_ds['spec_p2_sym']
            spec_p2_asy = wk_ds['spec_p2_asy']

      print(f'      cnt 1/2: {cnt1}  /  {cnt2}')
      # print(f'      num_days: {num_days}')
      ##########################################################################
      ##########################################################################

      if plot_phase==1 and cnt1<seg_cnt_min: exit(f'ERROR: insufficient count for phase 1 (cnt={cnt1}) - exiting')
      if plot_phase==2 and cnt2<seg_cnt_min: exit(f'ERROR: insufficient count for phase 2 (cnt={cnt2}) - exiting')

      if plot_phase==1 and spec_type=='sym' : spec_list.append(spec_p1_sym.transpose().values)
      if plot_phase==1 and spec_type=='asym': spec_list.append(spec_p1_asy.transpose().values)
      if plot_phase==2 and spec_type=='sym' : spec_list.append(spec_p2_sym.transpose().values)
      if plot_phase==2 and spec_type=='asym': spec_list.append(spec_p2_asy.transpose().values)

      freq_list.append(wk_ds['frequency'].values)
      wvnm_list.append(wk_ds['wavenumber'].values)

   #----------------------------------------------------------------------------
   # calculate limits for common color bar
   #----------------------------------------------------------------------------

   data_min = np.min([np.nanmin(d) for d in spec_list])
   data_max = np.max([np.nanmax(d) for d in spec_list])

   # print()
   # print(f'    data_min: {data_min:.2f}')
   # print(f'    data_max: {data_max:.2f}')
   # print()
   
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

      # num_clev,aboutZero = 21,False
      # (cmin,cmax,cint) = ngl.nice_cntr_levels(data_min, data_max, 
      #                                  cint=None, max_steps=num_clev, 
      #                                  returnLevels=False, aboutZero=aboutZero )
      # tres.cnLevels = np.linspace(cmin,cmax,num=num_clev)
      
      # if var[v]=='OLR' : tres.cnLevels = np.linspace(0.01,5.01,30)
      if spec_type=='sym'  and var[v]=='U': tres.cnLevels = np.linspace(0.1,4.8,10)
      if spec_type=='asym' and var[v]=='U': tres.cnLevels = np.linspace(0.02,1.2,10)

      # log scales
      # if var[v]=='U'   : tres.cnLevels = np.logspace( -2, 0.8, num=20).round(decimals=4)
      # if var[v]=='OLR ' : tres.cnLevels = np.logspace(  0.0, 2.0, num=30).round(decimals=4)
      # if var[v]=='PRECT': tres.cnLevels = np.logspace( -1.5, 0.0, num=30).round(decimals=4)
      # if var[v]=='U850' : tres.cnLevels = np.logspace( -1.5, 1.0, num=30).round(decimals=4)
      # if var[v]=='U200' : tres.cnLevels = np.logspace( -1.0, 2.0, num=30).round(decimals=4)

      #-------------------------------------------------------------------------
      #-------------------------------------------------------------------------
      # tres.cnFillPalette = np.array( cmocean.cm.speed(np.linspace(0,1,256)) )
      # tres.cnFillPalette = np.array( cmocean.cm.rain(np.linspace(0,1,256)) )
      # tres.cnFillPalette = np.array( cmocean.cm.delta(np.linspace(0,1,256)) )
      # tres.cnFillPalette = np.array( cmocean.cm.balance(np.linspace(0,1,256)) )
      # tres.cnFillPalette = np.array( cmocean.cm.thermal(np.linspace(0,1,256)) )

      tres.sfXArray = wvnm_list[c]
      tres.sfYArray = freq_list[c]

      if add_diff:
         # ip = v*(num_case*2)+c
         ip = v*(num_case+num_case-1)+c
      else:
         ip = v*num_case+c if var_x_case else c*num_var+v

      plot[ip] = ngl.contour(wks, spec_list[c], tres)
      add_dispersion_curves(plot[ip])
      
      hs.set_subtitles(wks, plot[ip], 
                       left_string=name[c], 
                       center_string='', 
                       right_string=var_str, 
                       font_height=subtitle_font_height)

      #-------------------------------------------------------------------------
      # Add difference plot
      #-------------------------------------------------------------------------
      if add_diff:
         num_clev,aboutZero = 21,True
         (cmin,cmax,cint) = ngl.nice_cntr_levels(diff_min, diff_max, 
                                    cint=None, max_steps=num_clev, 
                                    returnLevels=False, aboutZero=aboutZero )
         tres.cnLevels = np.linspace(cmin,cmax,num=num_clev)

         # tres.cnLevels = np.arange(-10,10+2,2)

         if spec_type=='sym'  and var[v]=='U':tres.cnLevels=np.linspace(-2,2,11)
         if spec_type=='asym' and var[v]=='U':tres.cnLevels=np.linspace(-1,1,11)

         tres.cnFillPalette = np.array( cmocean.cm.balance(np.linspace(0,1,256)) )

         ip = v*(num_case+num_case-1)+num_case+c-1

         # if c==0:
         #    bres = hs.res_default()
         #    bres.pmTickMarkDisplayMode = 'Never'
         #    plot[ip] = ngl.blank_plot(wks,bres)
         # else:
         if c>0:
            nspec_diff = spec_list[c] - spec_list[0]
            plot[ip] = ngl.contour(wks, nspec_diff, tres)
            add_dispersion_curves(plot[ip])
            hide_aliasing_artifact(plot[ip])
            hs.set_subtitles(wks, plot[ip], 
                             left_string=name[c], 
                             center_string='(diff)', 
                             right_string=var_str, 
                             # center_sub_string='diff',
                             font_height=subtitle_font_height)

#---------------------------------------------------------------------------------------------------
# Finalize plot
#---------------------------------------------------------------------------------------------------
if add_diff:
   # layout = [num_var*2,num_case] if var_x_case else [num_case,num_var*2]
   layout = [num_var,num_case+(num_case-1)] if var_x_case else [num_case+(num_case-1),num_var]
else:
   layout = [num_var,num_case] if var_x_case else [num_case,num_var]

pnl_res = hs.setres_panel()
pnl_res.nglPanelYWhiteSpacePercent = 1
pnl_res.nglPanelXWhiteSpacePercent = 1

# pnl_res.nglPanelFigureStrings            = list(string.ascii_lowercase)
# pnl_res.nglPanelFigureStringsJust        = "TopLeft"
# pnl_res.nglPanelFigureStringsFontHeightF = 0.01

ngl.panel(wks,plot,layout,pnl_res)
ngl.end()

hc.trim_png(fig_file)
#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------

# v2 is modified from v1 to reduce the memory burden by loading each file separately
#-------------------------------------------------------------------------------
import os, ngl, copy, xarray as xr, numpy as np, glob, dask, numba, cmocean, string
import hapy_common as hc, hapy_E3SM as he, hapy_setres as hs
data_dir,data_sub = None,None
#-------------------------------------------------------------------------------
case_name,case,case_dir,case_sub = [],[],[],[]
clr,dsh,mrk = [],[],[]
obs_flag = []
def add_case(case_in,n='',p=None,s='',g=None,d=0,c='black',m=0,r=False,obs=False):
   global case_name,case,case_dir,case_sub,clr,dsh,mrk
   case.append(case_in); case_name.append(n)
   case_dir.append(p); case_sub.append(s)
   dsh.append(d) ; clr.append(c) ; mrk.append(m)
   obs_flag.append(obs)
#-------------------------------------------------------------------------------
var = []
file_type_list = []
obs_var_list = []
obs_file_type_list = []
def add_var(var_name,file_type,obs_var=None,obs_file_type=None): 
   var.append(var_name)
   file_type_list.append(file_type)
   obs_var_list.append(obs_var)
   obs_file_type_list.append(obs_file_type)
#-------------------------------------------------------------------------------
scrip_file = '/pscratch/sd/w/whannah/e3sm_scratch/files_grid/ne1024pg2_scrip.nc'
# sim_data_root = '/global/cfs/cdirs/e3sm/gsharing/EAMxx'
sim_data_root = '/global/cfs/cdirs/e3sm/terai/SCREAM/v1_production/FourSeasons'
obs_data_root = '/pscratch/sd/w/whannah/Obs'
#-------------------------------------------------------------------------------

add_case('ne1024pg2_ne1024pg2_DY2',     n='SCREAMv1 Jan 2020',p=sim_data_root)
add_case('ne1024pg2_ne1024pg2_Apr2013', n='SCREAMv1 Apr 2013',p=sim_data_root)
add_case('ne1024pg2_ne1024pg2_DY1',     n='SCREAMv1 Aug 2016',p=sim_data_root)
add_case('ne1024pg2_ne1024pg2_Oct2013', n='SCREAMv1 Oct 2013',p=sim_data_root)

add_case('ERA5_Jan_2020',n='ERA5 Jan 2020', p=obs_data_root,s='daily',obs=True)
add_case('ERA5_Apr_2013',n='ERA5 Apr 2013', p=obs_data_root,s='daily',obs=True)
add_case('ERA5_Aug_2016',n='ERA5 Aug 2016', p=obs_data_root,s='daily',obs=True)
add_case('ERA5_Oct_2013',n='ERA5 Oct 2013', p=obs_data_root,s='daily',obs=True)

#-------------------------------------------------------------------------------
add_var('horiz_winds','output.scream.HorizWinds.INSTANT','u','ERA5.daily.atm')
#-------------------------------------------------------------------------------

tmp_data_path,tmp_file_head = 'data','zonal_mean.v1.tmp'

fig_file,fig_type = 'figs/zonal-mean.U-wind','png'

# first_day,num_days =  0,40
first_day,num_days =  10,30

# recalculate_bin,create_plot = True,False
# recalculate_bin,create_plot = True,True
recalculate_bin,create_plot = False,True

use_common_label_bar = True
var_x_case = True
num_plot_col = 4

add_case_mean = False
if add_case_mean: num_plot_col+=1

dlat = 0.5 ; lat_bins = np.arange(-88+dlat/2,88-dlat/2,dlat)

#---------------------------------------------------------------------------------------------------
# Set up plot resources
#---------------------------------------------------------------------------------------------------
num_case = len(case)
num_var  = len(var)

wkres = ngl.Resources()
npix = 2048; wkres.wkWidth,wkres.wkHeight=npix,npix
wks = ngl.open_wks(fig_type,fig_file,wkres)

plot = [None]*num_case*num_var
if add_case_mean: plot = [None]*(num_case+2)*num_var
res = hs.res_contour_fill()
res.vpHeightF = 0.5
# res.vpWidthF  = 0.3
res.tiYAxisFontHeightF           = 0.02
res.tiXAxisFontHeightF           = 0.02
res.tmYLLabelFontHeightF         = 0.02
res.tmXBLabelFontHeightF         = 0.02
res.nglYAxisType = 'LinearAxis'
res.trYReverse = True
res.tiXAxisString = 'Latitude'
res.tiYAxisString = 'Pressure [hPa]'

# res.cnFillMode    = 'RasterFill'

res.trYMinF = 50

if use_common_label_bar:
   res.lbLabelBarOn = False
else:
   res.lbLabelFontHeightF    = 0.015
   res.lbTitlePosition       = 'Bottom'
   res.lbTitleFontHeightF    = 0.015
   # res.lbTopMarginF          = -0.2
   # res.lbBottomMarginF       =  0.2+0.01

if 'lev' not in locals(): lev = np.array([0])

#---------------------------------------------------------------------------------------------------
# Create the band-pass filter
#---------------------------------------------------------------------------------------------------
nwgt = 91
fc_lp = 1./( 20.)
fc_hp = 1./(100.)
wgt = hc.filter_wgts_lp_lanczos(nwgt,fc_lp=fc_lp,fc_hp=fc_hp)

@numba.njit()
def filter_numba(data,data_filt,wgt,tvals0,win_width,ntime):   
   for i in range(ntime):
      if i >= win_width and i < (ntime-win_width-1) :
         tvals = tvals0 + i
         data_filt[i] = np.sum( data[tvals] * wgt )
   return data_filt

@numba.njit()
def bin_numba(data, lat, lat_bins, nbin, ntime, ncol, bin_out, cnt):
   for b in range( nbin ):
      bin_bot = lat_bins[b] - dlat/2. 
      bin_top = lat_bins[b] - dlat/2. + dlat
      for n in range( ncol ):
         if ( lat[n] >=bin_bot ) & ( lat[n]  <bin_top ):
            cnt[b] = cnt[b]+1
            for t in range( ntime ):
               hov[t,b] = hov[t,b] + data[t,n]
      hov[:,b] = hov[:,b]/cnt[b]

#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------
for v in range(num_var):
   lev_list = []
   lat_list = []
   var_list = []
   data_list = []
   date_list = []
   for c in range(num_case):
      tmp_case = case[c]
      if obs_flag[c]:
         tfile_type = obs_file_type_list[v]
         tvar = obs_var_list[v]
         tmp_case = 'ERA5'
         if case[c]=='ERA5_Jan_2020': tfile_type = tfile_type+'.2020-0[12]'
         if case[c]=='ERA5_Apr_2013': tfile_type = tfile_type+'.2013-0[45]'
         if case[c]=='ERA5_Aug_2016': tfile_type = tfile_type+'.2016-0[89]'
         if case[c]=='ERA5_Oct_2013': tfile_type = tfile_type+'.2013-1[01]'
      else:
         tfile_type = file_type_list[v]
         tvar = var[v]
      if c==0: print(' '*2+'var: '+hc.tclr.GREEN+tvar+hc.tclr.END)
      print('\n'+' '*4+'case: '+hc.tclr.CYAN+case[c]+hc.tclr.END)
      
      tmp_file = f'{tmp_data_path}/{tmp_file_head}.{tvar}.{case[c]}.nc'
      tmp_file += '.nc'

      group_name = 'Grid' if case[c]=='IMERG' else None
      #-------------------------------------------------------------------------
      # idenfity the files to load
      file_path = f'{case_dir[c]}/{tmp_case}/{case_sub[c]}/{tfile_type}*'
      file_list = sorted(glob.glob(file_path))
      # trim down file list
      f0,nf = first_day,num_days
      # if obs_flag[c]: f0,nf = first_day*48,num_days*48 # IMERG frequency is 30min
      if nf==0: file_list = file_list[f0:     ] # use all files from first
      if nf >0: file_list = file_list[f0:f0+nf] # use initial files
      if nf <0: file_list = file_list[nf:     ] # use latest files
      #-------------------------------------------------------------------------
      # get time information for plot strings
      ds1 = xr.open_dataset( file_list[0] , group=group_name)
      ds2 = xr.open_dataset( file_list[-1], group=group_name)
      yr1=ds1['time.year'][0].values; mn1=ds1['time.month'][0].values; dy1=ds1['time.day'][0].values
      yr2=ds2['time.year'][0].values; mn2=ds2['time.month'][0].values; dy2=ds2['time.day'][0].values
      # date1 = f'{yr1}/{mn1:02d}/{dy1:02d}'
      # date2 = f'{yr2}/{mn2:02d}/{dy2:02d}'
      # date_list.append(f'{date1} - {date2}')
      date_list.append(f'{mn1}/{dy1}-{mn2}/{dy2}  {yr1}')
      #-------------------------------------------------------------------------
      if recalculate_bin :
         print(' '*6+'recalculating...')
         #-------------------------------------------------------------------
         # # idenfity the files to load
         # file_path = f'{case_dir[c]}/{tmp_case}/{case_sub[c]}/{tfile_type}*'
         # file_list = sorted(glob.glob(file_path))
         # # trim down file list
         # f0,nf = first_day,num_days
         # # if obs_flag[c]: f0,nf = first_day*48,num_days*48 # IMERG frequency is 30min
         # if nf==0: file_list = file_list[f0:     ] # use all files from first
         # if nf >0: file_list = file_list[f0:f0+nf] # use initial files
         # if nf <0: file_list = file_list[nf:     ] # use latest files
         #----------------------------------------------------------------------
         print(' '*6+'Loading data ('+obs_var_list[v]+')...')
         cnt = 0
         for f in range(len(file_list)):
            print(f' '*8+f'f: {f:03d}  file: {hc.tclr.YELLOW}{file_list[f]}{hc.tclr.END}')
            #-------------------------------------------------------------------
            data_ds = xr.open_dataset( file_list[f], group=group_name)
            data = data_ds[tvar]
            #-------------------------------------------------------------------
            if 'longitude' in data.dims : data = data.rename({'longitude':'lon'})
            if 'latitude'  in data.dims : data = data.rename({'latitude':'lat'})
            if 'level'     in data.dims : data = data.rename({'level':'lev'})
            if 'plev'      in data.dims : data = data.rename({'plev':'lev'})
            if 'dim2'      in data.dims : data = data.isel(dim2=0)
            if data['lev'].max() > 2e3 : data['lev'] = data['lev']/1e2
            #-------------------------------------------------------------------
            # zonal running average
            data_avg = data.mean(dim='lon')

            if f==0: data_out = 0
            data_out = ( data_out*(cnt) + data_avg.mean(dim='time') ) / (cnt+1)
            cnt+=1
            #-------------------------------------------------------------------
            # sums for running temporal variance of zonal mean
            if f==0: sum1,sum2 = 0,0; hr_cnt = len(data_avg['time'])
            sum1 += data_avg.sum(dim='time')
            sum2 += (data_avg*data_avg).sum(dim='time')

         #----------------------------------------------------------------------
         mean = sum1/(cnt*hr_cnt)
         var_out  = (sum2/(cnt*hr_cnt)) - (mean*mean)
         #----------------------------------------------------------------------
         # print some summary stats after the full hovmoller is created
         hc.print_stat(data_out,name=tvar,compact=True,indent=' '*6)
         hc.print_stat(var_out, name='Variance',compact=True,indent=' '*6)
         #----------------------------------------------------------------------
         # time = time - time[0]
         #----------------------------------------------------------------------
         # Write to file 
         print(' '*6+f'Writing data to file: {tmp_file}')
         data_out.name = tvar
         tmp_ds = xr.Dataset()
         tmp_ds[tvar]        = data_out
         tmp_ds['variance']  = var_out
         tmp_ds['first_day'] = first_day
         tmp_ds['num_days']  = num_days
         tmp_ds.to_netcdf(path=tmp_file,mode='w')
      else:
         print(' '*6+f'Reading pre-calculated data from file: {hc.tclr.MAGENTA}{tmp_file}{hc.tclr.END}')
         tmp_ds = xr.open_dataset( tmp_file )
         lat_bins  = tmp_ds['lat'].values
         data_out  = tmp_ds[tvar]
         var_out   = tmp_ds['variance']
         # first_day = tmp_ds['first_day']
         # num_days  = tmp_ds['num_days']

      print(); print(data_out['lev'].values)
      #-------------------------------------------------------------------------
      #-------------------------------------------------------------------------
      lev_list.append(data_out['lev'].values)
      lat_list.append(data_out['lat'].values)
      var_list.append(  np.ma.masked_invalid(var_out.values) )
      data_list.append( np.ma.masked_invalid(data_out.values) )
   
   if not create_plot: continue
   
   #----------------------------------------------------------------------------
   #----------------------------------------------------------------------------
   def set_clev(res_in, data_min, data_max, var_name, nlev=11, aboutZero=True):
      (cmin,cmax,cint) = ngl.nice_cntr_levels(data_min, data_max, cint=None, max_steps=nlev,returnLevels=False, aboutZero=aboutZero )
      clev = np.linspace(cmin,cmax,num=nlev)
      res_in.cnLevels = clev
      res_in.cnLevelSelectionMode = 'ExplicitLevels'
      return clev

   def neg_dash_contours(res_in, clev):
      dsh = np.zeros((len(clev)),'i')
      for k in range(len(clev)):
         if (clev[k] < 0.): dsh[k] = 6
      res_in.cnLineDashPatterns    = dsh
      res_in.cnMonoLineDashPattern = False
      res_in.cnFillOn  = False
      res_in.cnLinesOn = True
   #-------------------------------------------------------------------------
   # Create plot
   #-------------------------------------------------------------------------
   cres = copy.deepcopy(res)
   cres.tmXBOn               = False
   cres.tmYLOn               = False
   cres.cnFillOn             = False
   cres.cnLinesOn            = True
   cres.cnLineLabelsOn       = False
   cres.cnInfoLabelOn        = False
   cres.cnLineThicknessF     = 1.
   cres.cnLevelSelectionMode = 'ExplicitLevels'
   # cres.cnLevels             = np.arange(8,800+8,8)
   # cres.cnLevels             = np.linspace(-45,45,num=9*2+1)
   # cspc = 8 ; cmax = cspc*6 ; cmax2 = cspc*4
   # cspc = 6 ; cmax = cspc*8 ; cmax2 = cspc*5
   cspc = 5 ; cmax = cspc*10 ; cmax2 = cspc*6
   cres.cnLevels             = np.arange(-1*cmax,cmax+cspc,cspc)

   # clev = set_clev(cres, cdata_min, cdata_max, cvar_name[v])
   neg_dash_contours(cres, cres.cnLevels)

   res.cnLevelSelectionMode  = 'ExplicitLevels'
   # res.cnLevels              = np.linspace(-45,45,num=9*2+1)
   res.cnLevels              = np.arange(6,54+6,6)


   # res.cnFillPalette = 'MPL_viridis'
   # res.cnFillPalette = np.array( cmocean.cm.rain(np.linspace(0,1,256)) )
   res.cnFillPalette = np.array( cmocean.cm.amp(np.linspace(0,1,256)) )
   # res.cnFillPalette = np.array( cmocean.cm.balance(np.linspace(0,1,256)) )
   
   for c in range(num_case):
      
      ip = v*num_case+c if var_x_case else c*num_var+v
      if add_case_mean: 
         # ip = v*(num_case+2)+c if var_x_case else c*num_var+v
         if 'SCREAMv1' in case_name[c]: ip = 0+c
         if 'ERA5'     in case_name[c]: ip = 1+c

      tres = copy.deepcopy(res)
      # tres.tmXBMode      = 'Explicit'
      # tres.tmXBValues    = np.arange( len(lon_bins) )
      # tres.tmXBLabels    = lon_bins

      if not use_common_label_bar:
         if tvar=='u': tres.lbTitleString = '[m/s]'

      # tres.nglXAxisType = 'LinearAxis'
      # tres.tiXAxisString = 'sin( Latitude )'
      # sin_lat_bins = np.sin(lat_list[c]*np.pi/180.)
      # tres.sfYArray = lev_list[c] / 1e2
      # tres.sfXArray = sin_lat_bins
      # lat_tick = np.array([-90,-60,-30,0,30,60,90])
      # tres.tmXBMode = "Explicit"
      # tres.tmXBValues = np.sin( lat_tick*3.14159/180. )
      # tres.tmXBLabels = lat_tick

      tres.sfYArray,tres.sfXArray = lev_list[c],lat_list[c]
      cres.sfYArray,cres.sfXArray = lev_list[c],lat_list[c]

      # plot[ip] = ngl.contour(wks, data_list[c] ,tres) 
      # ngl.overlay( plot[ip], ngl.contour(wks, var_list[c] ,cres) )

      plot[ip] = ngl.contour(wks, var_list[c] ,tres) 
      ngl.overlay( plot[ip], ngl.contour(wks, data_list[c] ,cres) )

      cres2 = copy.deepcopy(cres)
      cres2.cnMonoLineDashPattern = True
      cres2.cnLineThicknessF = 6.
      cres2.cnLevels         = np.arange(cmax2,cmax+cspc,cspc)
      ngl.overlay( plot[ip], ngl.contour(wks, data_list[c] ,cres2) )

      var_str = tvar
      if var_str=='u': var_str = 'Zonal Wind'
      lstr = case_name[c].split(' ')[0]
      # cstr = case_name[c].replace(lstr+' ','')

      # if 'Jan 2020' in name: cstr = ''
      # if 'Apr 2013' in name: cstr = ''
      # if 'Aug 2016' in name: cstr = ''
      # if 'Oct 2013' in name: cstr = ''
      
      # hs.set_subtitles(wks, plot[ip], lstr, cstr, var_str, font_height=0.015)
      hs.set_subtitles(wks, plot[ip], lstr, '', date_list[c], font_height=0.008)

   if add_case_mean:
      for n in range(2):
         cnt = 0
         for c in range(num_case):
            avg_flag = False
            if n==0 and 'SCREAMv1' in case_name[c]: avg_flag = True
            if n==1 and 'ERA5'     in case_name[c]: avg_flag = True
            
            if avg_flag:
               if cnt==0: 
                  case_avg = data_list[c]*0
                  case_var = var_list[c]*0
                  lev_tmp = lev_list[c]
                  lat_tmp = lat_list[c]
               case_avg = ( case_avg*cnt + data_list[c] ) / (cnt+1)
               case_var = ( case_var*cnt +  var_list[c] ) / (cnt+1)
               cnt+=1

         # ip = v*(num_case+2)+num_case+n if var_x_case else (num_case+n)*num_var+v
         ip = (n+1)*5-1
         print(ip)

         tres = copy.deepcopy(res)
         if not use_common_label_bar:
            if tvar=='u': tres.lbTitleString = '[m/s]'
         tres.sfYArray,tres.sfXArray = lev_tmp,lat_tmp
         cres.sfYArray,cres.sfXArray = lev_tmp,lat_tmp
         plot[ip] = ngl.contour(wks, case_var ,tres) 
         ngl.overlay( plot[ip], ngl.contour(wks, case_avg ,cres) )

         cres2 = copy.deepcopy(cres)
         cres2.cnMonoLineDashPattern = True
         cres2.cnLineThicknessF = 6.
         cres2.cnLevels         = np.arange(cmax2,cmax+cspc,cspc)
         ngl.overlay( plot[ip], ngl.contour(wks, case_avg ,cres2) )

         var_str = tvar
         if var_str=='u': var_str = 'Zonal Wind'
         lstr = case_name[n*4].split(' ')[0]
         
         # hs.set_subtitles(wks, plot[ip], lstr, cstr, var_str, font_height=0.015)
         hs.set_subtitles(wks, plot[ip], lstr, '', 'ALL', font_height=0.008)

#---------------------------------------------------------------------------------------------------
# Finalize plot
#---------------------------------------------------------------------------------------------------
if create_plot:
   layout = [int(np.ceil(len(plot)/float(num_plot_col))),num_plot_col]
   # if 'num_plot_col' in locals():
   #    layout = [int(np.ceil(len(plot)/float(num_plot_col))),num_plot_col]
   # else:
   #    layout = [num_var,num_case] if var_x_case else [num_case,num_var]
   pnl_res = hs.setres_panel()
   if use_common_label_bar: 
      pnl_res.nglPanelLabelBar                  = True
      pnl_res.nglPanelLabelBarLabelFontHeightF  = 0.008
      pnl_res.lbTitleFontHeightF                = 0.008
      pnl_res.lbTopMarginF                      =  0.2
      pnl_res.lbBottomMarginF                   = -0.1
      # pnl_res.lbLeftMarginF                     =  1.
      # pnl_res.lbRightMarginF                    =  1.
      pnl_res.lbTitlePosition                   = 'Bottom'
      if tvar=='u': pnl_res.lbTitleString = '[m/s]'
   pnl_res.nglPanelFigureStrings            = list(string.ascii_lowercase)
   pnl_res.nglPanelFigureStringsJust        = "TopLeft"
   pnl_res.nglPanelFigureStringsFontHeightF = 0.005
   ngl.panel(wks,plot,layout,pnl_res)
   ngl.end()

   hc.trim_png(fig_file)
#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------

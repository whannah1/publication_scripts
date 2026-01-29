# v2 is modified from v1 to reduce the memory burden by loading each file separately
#-------------------------------------------------------------------------------
import os, ngl, copy, xarray as xr, numpy as np, glob, dask, numba, cmocean
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
sim_data_root = '/global/cfs/cdirs/e3sm/gsharing/EAMxx'
obs_data_root = '/pscratch/sd/w/whannah/Obs'

# add_case('DYAMOND2_SCREAMv1' ,n='SCREAMv1 Jan 2020',p=data_root)
# add_case('Apr1_2013_SCREAMv1',n='SCREAMv1 Apr 2013',p=data_root)
# add_case('DYAMOND1_SCREAMv1' ,n='SCREAMv1 Aug 2016',p=data_root)
# add_case('Oct1_2013_SCREAMv1',n='SCREAMv1 Oct 2013',p=data_root)

add_case('Apr1_2013_SCREAMv1',n='SCREAMv1 Apr 2013',p=sim_data_root)
add_case('IMERG',             n='IMERG Apr 2013',   p=obs_data_root,s='hourly_nc',obs=True)

#-------------------------------------------------------------------------------

# add_var('LW_flux_up@tom','output.scream.TOMVars.INSTANT')

# add_var('VapWaterPath','output.scream.VertIntegrals.INSTANT')
# add_var('IceWaterPath','output.scream.VertIntegrals.INSTANT')

add_var('precip','output.scream.SurfVars.INSTANT','precipitation','3B-HHR.MS.MRG.3IMERG')

# add_var('precip_liq_surf_mass','output.scream.SurfVars.INSTANT','precipAvg','3B-DAY.MS.MRG.3IMERG.2016')
# add_var('precip_ice_surf_mass','output.scream.SurfVars.INSTANT')
# add_var('surf_evap','output.scream.SurfVars.INSTANT')
# add_var('horiz_winds@bot','output.scream.SurfVars.INSTANT')

#-------------------------------------------------------------------------------
# tmp_data_path = os.getenv('HOME')+'/Research/E3SM/pub_figs/2023_screamv1_4season/data'
tmp_data_path = 'data'

fig_file,fig_type = 'figs/hovmoller.v2.CONUS','png'


first_day,num_days =  0,40

# lat1,lat2 = 30,48 # Carbone & Tuttle 2008
# lat1,lat2 = 32,45 # Surcel et al. 2010
# lat1,lat2 = 38,49 # Clark et al. 2007 (Fig 7)
lat1,lat2 = 35,45
lon1,lon2 = 360-120,360-78
# lon1,lon2 = 360-114,360-78

# dlon = 1; lon_bins = np.arange(50,180,dlon)
# dlon = 1; lon_bins = np.arange(1,360,dlon)
# dlon = 0.5; lon_bins = np.arange(0.5,359.5,dlon)
dlon = 0.1; lon_bins = np.arange(lon1+dlon/2,lon2-dlon/2,dlon)

# recalculate_hov,create_plot = True,False
recalculate_hov,create_plot = False,True
# recalculate_hov,create_plot = True,True

var_x_case = True

#---------------------------------------------------------------------------------------------------
# Set up plot resources
#---------------------------------------------------------------------------------------------------
num_case = len(case)
num_var  = len(var)

wkres = ngl.Resources()
npix = 2048; wkres.wkWidth,wkres.wkHeight=npix,npix
wks = ngl.open_wks(fig_type,fig_file,wkres)

plot = [None]*num_case*num_var
res = hs.res_contour_fill()
# res.vpHeightF = 0.6
res.vpWidthF  = 0.3
# res.tmYLLabelFontHeightF         = 0.008
# res.tmXBLabelFontHeightF         = 0.008

res.tiYAxisString = 'Time [days]'
res.tiXAxisString = 'Longitude'

# res.cnFillMode    = 'RasterFill'

# res.trXMinF,res.trXMaxF = 0, 180
# res.trYMinF = 10
# res.trYMaxF = 20

# res.cnFillPalette = 'MPL_viridis'
res.cnFillPalette = np.array( cmocean.cm.rain(np.linspace(0,1,256)) )
# res.cnFillPalette = np.array( cmocean.cm.amp(np.linspace(0,1,256)) )

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
def hov_numba(data, lon, lon_bins, nbin, ntime, ncol, hov, cnt):
   for b in range( nbin ):
      bin_bot = lon_bins[b] - dlon/2. 
      bin_top = lon_bins[b] - dlon/2. + dlon
      for n in range( ncol ):
         if ( lon[n] >=bin_bot ) & ( lon[n]  <bin_top ):
            cnt[b] = cnt[b]+1
            for t in range( ntime ):
               hov[t,b] = hov[t,b] + data[t,n]
      hov[:,b] = hov[:,b]/cnt[b]

#---------------------------------------------------------------------------------------------------

# if recalculate_hov :
#    grid_ds = xr.open_dataset( scrip_file ).rename({'grid_size':'ncol'})
#    ncol = grid_ds['ncol'].values
#    # define mask for regional subset
#    mask = xr.DataArray( np.ones([len(ncol)],dtype=bool), dims='ncol' )
#    if 'lat1' in globals(): mask = mask & (grid_ds['grid_center_lat']>=lat1)
#    if 'lat2' in globals(): mask = mask & (grid_ds['grid_center_lat']<=lat2)
#    if 'lon1' in globals(): mask = mask & (grid_ds['grid_center_lon']>=lon1)
#    if 'lon2' in globals(): mask = mask & (grid_ds['grid_center_lon']<=lon2)
#    mask = mask.compute()
#    # get longitude values for binning
#    lon = grid_ds['grid_center_lon']
#    lon = lon.where( mask, drop=True).values
scrip_loaded = False
#---------------------------------------------------------------------------------------------------
for v in range(num_var):
   time_list = []
   lon_list = []
   hov_list = []
   for c in range(num_case):
      if obs_flag[c]:
         tfile_type = obs_file_type_list[v]
         tvar = obs_var_list[v]
      else:
         tfile_type = file_type_list[v]
         tvar = var[v]
      if c==0: print(' '*2+'var: '+hc.tcolor.GREEN+tvar+hc.tcolor.ENDC)
      print('\n'+' '*4+'case: '+hc.tcolor.CYAN+case[c]+hc.tcolor.ENDC)
      
      tmp_file = f'{tmp_data_path}/CONUS.hovmoller.horz.v1.tmp.{tvar}.{case[c]}'
      if 'lat1' in globals(): tmp_file += f'.lat1_{lat1}.lat2_{lat2}'
      if 'lon1' in globals(): tmp_file += f'.lon1_{lon1}.lon2_{lon2}'
      tmp_file += '.nc'

      if recalculate_hov :
         print(' '*6+'recalculating...')
         #-------------------------------------------------------------------
         # idenfity the files to load
         file_path = f'{case_dir[c]}/{case[c]}/{case_sub[c]}/{tfile_type}*'
         file_list = sorted(glob.glob(file_path))
         # trim down file list
         f0,nf = first_day,num_days
         if obs_flag[c]: f0,nf = first_day*48,num_days*48 # IMERG frequency is 30min
         if nf==0: file_list = file_list[f0:     ] # use all files from first
         if nf >0: file_list = file_list[f0:f0+nf] # use initial files
         if nf <0: file_list = file_list[nf:     ] # use latest files

         #----------------------------------------------------------------------
         print(' '*6+'Loading data ('+obs_var_list[v]+')...')
         for f in range(len(file_list)):
            print(f' '*8+f'f: {f:03d}  file: {hc.tcolor.YELLOW}{file_list[f]}{hc.tcolor.ENDC}')
            #-------------------------------------------------------------------
            group_name = 'Grid' if case[c]=='IMERG' else None
            data_ds = xr.open_dataset( file_list[f], group=group_name)
            #-------------------------------------------------------------------
            # Load the data
            if tvar=='precip':
               data = data_ds['precip_liq_surf_mass'] + data_ds['precip_ice_surf_mass']
            else:
               data = data_ds[tvar]
            if obs_flag[c]:
               #-------------------------------------------------------------------
               # load data and stack lat/lon dimensions
               data_ds['lon'] = xr.where(data_ds['lon']<0,data_ds['lon']+360,data_ds['lon'])
               data['lon'] = data_ds['lon']
               lat2D =               np.repeat( data_ds['lat'].values[...,None],len(data_ds['lon']),axis=1)
               lon2D = np.transpose( np.repeat( data_ds['lon'].values[...,None],len(data_ds['lat']),axis=1) )
               lat = xr.DataArray(lat2D,dims=('lat','lon'),coords=(data.lat,data.lon)).stack(ncol=('lat','lon'))
               lon = xr.DataArray(lon2D,dims=('lat','lon'),coords=(data.lat,data.lon)).stack(ncol=('lat','lon'))
               data = data.stack(ncol=('lat','lon'))
               #-------------------------------------------------------------------
               # apply mask
               if f==0:
                  obs_mask = xr.ones_like(lat,dtype=bool).load()#.isel(time=0,drop=True).load()
                  if 'lat1' in globals(): obs_mask = np.logical_and( obs_mask, (lat>=lat1))
                  if 'lat2' in globals(): obs_mask = np.logical_and( obs_mask, (lat<=lat2))
                  if 'lon1' in globals(): obs_mask = np.logical_and( obs_mask, (lon>=lon1))
                  if 'lon2' in globals(): obs_mask = np.logical_and( obs_mask, (lon<=lon2))
               data = data.where( obs_mask, drop=True)
               # lat  =  lat.where( obs_mask, drop=True).values
               lon  =  lon.where( obs_mask, drop=True).values
               data_size = len(data.ncol)
            else:
               if not scrip_loaded:
                  grid_ds = xr.open_dataset( scrip_file ).rename({'grid_size':'ncol'})
                  ncol = grid_ds['ncol'].values
                  # define mask for regional subset
                  mask = xr.DataArray( np.ones([len(ncol)],dtype=bool), dims='ncol' )
                  if 'lat1' in globals(): mask = mask & (grid_ds['grid_center_lat']>=lat1)
                  if 'lat2' in globals(): mask = mask & (grid_ds['grid_center_lat']<=lat2)
                  if 'lon1' in globals(): mask = mask & (grid_ds['grid_center_lon']>=lon1)
                  if 'lon2' in globals(): mask = mask & (grid_ds['grid_center_lon']<=lon2)
                  mask = mask.compute()
                  # get longitude values for binning
                  lon = grid_ds['grid_center_lon']
                  lon = lon.where( mask, drop=True).values
                  scrip_loaded = True
               #-------------------------------------------------------------------
               # apply mask
               data = data.where( mask, drop=True)
               if 'dim2' in data.dims : data = data.isel(dim2=0)
               data_size = len(data['ncol'].values)

            #-------------------------------------------------------------------
            # convert units
            if var[v]=='precip':
               if obs_flag[c]:
                  data = data*24. # mm/hr to mm/day
               else:
                  data = (data/100)*86400 # kg/m2 to mm/day using dtime=100 (ne1024)
            #-------------------------------------------------------------------
            # build time coodinate
            time_tmp = data['time.year'].values*365    \
                      +data['time.dayofyear'].values   \
                      +data['time.hour'].values/24       \
                      +data['time.minute'].values/24/60
            #-------------------------------------------------------------------
            # Setup output arrays
            nbins,ntime = len(lon_bins),len(time_tmp)
            cnt = np.zeros(nbins)
            hov_tmp = xr.DataArray( np.zeros( (ntime,nbins) ), \
                                    coords=[('time',time_tmp),('lon',lon_bins)], \
                                    dims=['time','lon'] )
            #-------------------------------------------------------------------
            # bin the data to create the hovmoller
            data.load()
            hov_numba( data.values, lon, lon_bins, nbins, ntime, data_size, hov_tmp.values, cnt )
            #-------------------------------------------------------------------
            if f==0:
               hov  = hov_tmp
               time = time_tmp
            else:
               hov  = xr.concat([hov,hov_tmp], dim='time')
               time = np.append(time,time_tmp)
         #----------------------------------------------------------------------
         # print some summary stats after the full hovmoller is created
         hc.print_stat(hov,name=tvar,compact=True,indent=' '*6)
         #----------------------------------------------------------------------
         # time = time - time[0]
         #----------------------------------------------------------------------
         # Write to file 
         print(' '*6+f'Writing hovmoller data to file: {tmp_file}')
         hov.name = tvar
         tmp_ds = xr.Dataset()
         # tmp_ds['lon_bins'] = lon_bins
         # tmp_ds['time'] = time
         # tmp_ds['cnt']  = cnt
         tmp_ds[tvar] = hov
         tmp_ds.to_netcdf(path=tmp_file,mode='w')
      else:
         print(' '*6+f'Reading pre-calculated hovmoller data from file: {hc.tcolor.MAGENTA}{tmp_file}{hc.tcolor.ENDC}')
         tmp_ds = xr.open_dataset( tmp_file )
         # lon_bins = tmp_ds['lon_bins'].values
         lon_bins = tmp_ds['lon'].values
         time     = tmp_ds['time'].values
         hov      = tmp_ds[tvar]

      #-------------------------------------------------------------------------
      #-------------------------------------------------------------------------
      time = time - time[0]

      time_list.append(time)
      lon_list.append(lon_bins)
      hov_list.append(hov.values)
   
   if not create_plot: continue
  
   #-------------------------------------------------------------------------
   # Create plot
   #-------------------------------------------------------------------------
   ntime = len(time)
   
   # num_clev = 21
   # aboutZero = False
   # # if tvar in ['U','V']: aboutZero = True
   # min_val = np.min([np.min(h) for h in hov_list])
   # max_val = np.max([np.max(h) for h in hov_list])
   # clev_tup = ngl.nice_cntr_levels(min_val, max_val,       \
   #                                 cint=None, max_steps=num_clev, \
   #                                 returnLevels=False, aboutZero=aboutZero )
   
   # # print(f'clev_tup: {clev_tup}')
   # if clev_tup==None: 
   #    print('SWITCHING TO AUTOMATIC CONTOUR LEVELS!')
   #    res.cnLevelSelectionMode = "AutomaticLevels"   
   # else:
   #    cmin,cmax,cint = clev_tup
   #    res.cnLevels = np.linspace(cmin,cmax,num=num_clev)
   #    res.cnLevelSelectionMode = "ExplicitLevels"

   if var[v]=='precip': 
      # res.cnLevels = np.logspace( -1, 1.35, num=30).round(decimals=4)
      res.cnLevels = np.logspace( -1, 1.35, num=10).round(decimals=2)
      res.cnLevelSelectionMode = "ExplicitLevels"

   for c in range(num_case):
      ip = v*num_case+c if var_x_case else c*num_var+v
      tres = copy.deepcopy(res)
      # tres.tmXBMode      = 'Explicit'
      # tres.tmXBValues    = np.arange( len(lon_bins) )
      # tres.tmXBLabels    = lon_bins

      if var[v]=='precip': tres.lbTitleString = 'mm/day'

      tres.sfYArray = time_list[c]
      tres.sfXArray = lon_list[c]
      plot[ip] = ngl.contour(wks, hov_list[c] ,tres) 
      # plot.append( ngl.contour(wks, hov[nwgt:ntime-nwgt-1,:].values ,res) )
      hs.set_subtitles(wks, plot[ip], case_name[c], '', tvar, font_height=0.008)

   # print('\nContour Levels:')
   # print(ngl.get_float_array(plot[len(plot)-1].contour,"cnLevels"))
   # print()

   # else:
      # res.cnFillOn   = False
      # res.cnLinesOn  = True

      # res.cnLineDashPattern = 1

      # for c in range(num_case):
      #    tplot = ngl.contour(wks, hov[nwgt:ntime-nwgt-1,:].values ,res)
      #    ngl.overlay(plot[c],tplot)
      
      #    levels = ngl.get_float_array(tplot.contour,"cnLevels")

      #    # print(levels)
      #    # print(levels[ np.where(levels>=0) ])

      #    res.cnLevelSelectionMode = "ExplicitLevels"
      #    res.cnLevels = levels[ np.where(levels>=0) ]
      #    res.cnLineDashPattern = 0
      #    ngl.overlay(plot[c],ngl.contour(wks, hov[nwgt:ntime-nwgt-1,:].values ,res))

      


#---------------------------------------------------------------------------------------------------
# Finalize plot
#---------------------------------------------------------------------------------------------------
if create_plot:
   # layout = [1,len(plot)]
   layout = [num_var,num_case] if var_x_case else [num_case,num_var]
   ngl.panel(wks,plot[0:len(plot)],layout,hs.setres_panel())
   ngl.end()

   hc.trim_png(fig_file)
#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------

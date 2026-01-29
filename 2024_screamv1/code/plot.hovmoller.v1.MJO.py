import os, ngl, copy, xarray as xr, numpy as np, glob, dask, numba
import hapy_common as hc, hapy_E3SM as he, hapy_setres as hs
data_dir,data_sub = None,None
#-------------------------------------------------------------------------------
case_name,case,case_dir,case_sub = [],[],[],[]
clr,dsh,mrk = [],[],[]
obs_flag = []
def add_case(case_in,n='',p=None,s='',g=None,d=0,c='black',m=0,r=False,obs=False):
   global name,case,case_dir,case_sub,clr,dsh,mrk
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
data_root = '/global/cfs/cdirs/e3sm/gsharing/EAMxx'
obs_data_root = '/pscratch/sd/w/whannah/Obs'

# add_case('DYAMOND2_SCREAMv1' ,n='SCREAMv1 Jan 2020',p=data_root)
# add_case('Apr1_2013_SCREAMv1',n='SCREAMv1 Apr 2013',p=data_root)
# add_case('DYAMOND1_SCREAMv1' ,n='SCREAMv1 Aug 2016',p=data_root)
# add_case('Oct1_2013_SCREAMv1',n='SCREAMv1 Oct 2013',p=data_root)

add_case('IMERG',n='IMERG Jan 2020',p=obs_data_root,s='daily_QC_Jan_2020',obs=True)
add_case('IMERG',n='IMERG Apr 2013',p=obs_data_root,s='daily_QC_Apr_2013',obs=True)
add_case('IMERG',n='IMERG Aug 2016',p=obs_data_root,s='daily_QC_Aug_2016',obs=True)
add_case('IMERG',n='IMERG Oct 2013',p=obs_data_root,s='daily_QC_Oct_2013',obs=True)

#-------------------------------------------------------------------------------

# add_var('LW_flux_up@tom','output.scream.TOMVars.INSTANT')

# add_var('VapWaterPath','output.scream.VertIntegrals.INSTANT')
# add_var('IceWaterPath','output.scream.VertIntegrals.INSTANT')

add_var('precip_liq_surf_mass','output.scream.SurfVars.INSTANT','precipAvg','3B-DAY.MS.MRG.3IMERG.2016')
# add_var('precip_ice_surf_mass','output.scream.SurfVars.INSTANT')
# add_var('surf_evap','output.scream.SurfVars.INSTANT')
# add_var('horiz_winds@bot','output.scream.SurfVars.INSTANT')

#-------------------------------------------------------------------------------
tmp_data_path = os.getenv('HOME')+'/Research/E3SM/pub_figs/2023_screamv1_4season/data'

fig_file,fig_type = 'figs/hovmoller.v1.MJO','png'


first_file,num_files = 0,2
# first_file,num_files = 0,0

lat1,lat2 = -15,15
# lon1,lon2 =  40,180
lon1,lon2 =  0,360

# dlon = 1; lon_bins = np.arange(50,180,dlon)
# dlon = 1; lon_bins = np.arange(1,360,dlon)
# dlon = 0.5; lon_bins = np.arange(0.5,359.5,dlon)
dlon = 0.5; lon_bins = np.arange(lon1+dlon/2,lon2-dlon/2,dlon)


recalculate_hov = True

apply_filter = False
var_x_case = True

#---------------------------------------------------------------------------------------------------
# Set up plot resources
#---------------------------------------------------------------------------------------------------
num_case = len(case)
num_var  = len(var)

wks = ngl.open_wks(fig_type,fig_file)
plot = [None]*num_case*num_var
res = hs.res_contour_fill()
# res.vpHeightF = 0.6
# res.vpWidthF  = 0.3
# res.tmYLLabelFontHeightF         = 0.008
# res.tmXBLabelFontHeightF         = 0.008

res.tiYAxisString = 'Time [days]'
res.tiXAxisString = 'Longitude'

# res.trXMinF,res.trXMaxF = 0, 180

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

#---------------------------------------------------------------------------------------------------
for v in range(num_var):
   if obs_flag[c]:
      tfile_type = obs_file_type_list
      tvar = obs_var_list[v]
   else:
      tfile_type = file_type_list
      tvar = var[v]
   print(' '*2+'var: '+hc.tcolor.GREEN+tvar+hc.tcolor.ENDC)
   time_list = []
   lon_list = []
   hov_list = []
   for c in range(num_case):
      print(' '*4+'case: '+hc.tcolor.CYAN+case[c]+hc.tcolor.ENDC)
      
      tmp_file = f'{tmp_data_path}/MJO.hovmoller.horz.v1.tmp.{tvar}.{case[c]}'
      if 'lat1' in globals(): tmp_file += f'.lat1={lat1}.lat2={lat2}'
      if 'lon1' in globals(): tmp_file += f'.lon1={lon1}.lon2={lon2}'
      tmp_file += '.nc'

      if recalculate_hov :
         print(' '*6+'recalculating...')
         #-------------------------------------------------------------------
         # idenfity the files to load
         file_path = f'{case_dir[c]}/{case[c]}/{case_sub[c]}/{tfile_type}*'
         file_list = sorted(glob.glob(file_path))
         # trim down file list
         f0,nf = first_file,num_files
         if nf==0: file_list = file_list[f0:     ] # use all files from first
         if nf >0: file_list = file_list[f0:f0+nf] # use initial files
         if nf <0: file_list = file_list[nf:     ] # use latest files
         #----------------------------------------------------------------------
         with dask.config.set(**{'array.slicing.split_large_chunks': True}):
            data_ds = xr.open_mfdataset( file_list )
            #-------------------------------------------------------------------
            # Load the data
            print(' '*6+'Loading data ('+obs_var_list[v]+')...')
            data = data_ds[tvar]
            if obs_flag[c]:
               #-------------------------------------------------------------------
               # apply mask
               mask = xr.DataArray( np.ones(data.shape,dtype=bool),coords=data.coords,dims=data.dims )
               if 'lat1' in globals(): mask = mask & (ds['lat']>=lat1)
               if 'lat2' in globals(): mask = mask & (ds['lat']<=lat2)
               if 'lon1' in globals(): mask = mask & (ds['lon']>=lon1)
               if 'lon2' in globals(): mask = mask & (ds['lon']<=lon2)
               data = data.where( mask, drop=True)
               lat  = ds['lat'].where( mask, drop=True)
               lon  = ds['lon'].where( mask, drop=True)
               data_size = len(lat)*len(lon)
               #-------------------------------------------------------------------
               # stack the lat/lon dimensions
               lon2D = np.transpose( np.repeat( lon.values[...,None],len(lat),axis=1) )
               lon = xr.DataArray(lon2D,dims=('lat','lon')).stack(ncol=('lat','lon'))
               lon = lon.values
               data = data.stack(ncol=('lat','lon'))
            else:
               #-------------------------------------------------------------------
               # apply mask
               data = data.where( mask, drop=True)
               if 'dim2' in data.dims : data = data.isel(dim2=0)
               data_size = len(data['ncol'].values)
            #-------------------------------------------------------------------
            # Convert to daily
            data = data.resample(time='D').mean(dim='time')
            #-------------------------------------------------------------------
            # build daily time coodinate
            time = data['time.year'].values*365 + data['time.dayofyear'].values
            time = time - time[0]
            #-------------------------------------------------------------------
            # Setup output arrays
            nbins,ntime = len(lon_bins),len(time)
            cnt = np.zeros(nbins)
            hov = xr.DataArray( np.zeros( (ntime,nbins) ), \
                                coords=[('time',time),('lon',lon_bins)], \
                                dims=['time','lon'] )
            #-------------------------------------------------------------------
            # bin the data to create the hovmoller
            hov_numba( data.values, lon, lon_bins, nbins, ntime, data_size, hov.values, cnt )

            #-------------------------------------------------------------------
            # print some summary stats
            hc.print_stat(hov,name=tvar,compact=True,indent=' '*6)
            #-------------------------------------------------------------------
            # Write to file 
            print(' '*6+f'Writing hovmoller data to file: {tmp_file}')
            hov.name = tvar
            tmp_ds = xr.Dataset()
            tmp_ds['lon_bins'] = lon_bins
            tmp_ds['time'] = time
            tmp_ds['cnt']  = cnt
            tmp_ds[tvar] = hov
            tmp_ds.to_netcdf(path=tmp_file,mode='w')
      else:
         tmp_ds = xr.open_dataset( tmp_file )
         lon_bins = tmp_ds['lon_bins'].values
         time     = tmp_ds['time'].values
         hov      = tmp_ds[tvar]

      #-------------------------------------------------------------------------
      # Filter in time
      if apply_filter:
         ntime = len(hov.time)
         win_width = np.floor(nwgt/2)
         tvals0 = np.arange( 0-win_width, 0+win_width+1, dtype=np.int)
         # hov_mean = hov.mean(dim='time')
         hov = hov - hov.mean(dim='time')
         hov_filt = hov.copy(data=np.full(hov.shape,np.nan))
         for b in range( len(lon_bins) ):
            filter_numba(hov[:,b].values,hov_filt[:,b].values,wgt,tvals0,win_width,ntime)
            # hov_filt[:,b] = ( hov_filt[:,b] + hov_mean[b] ) - hov_mean.mean()
         hov = hov_filt

         # hc.printline()
         # print(hov)
         # hc.print_stat(hov,name=var[v])
         # hc.printline()

         # exit()
      #----------------------------------------------------------------------
      #----------------------------------------------------------------------

      time_list.append(time)
      lon_list.append(lon_bins)
      hov_list.append(hov.values)


   #----------------------------------------------------------------------------
   #----------------------------------------------------------------------------
   # unit_str = ''
   # if var in ['PRECT','PRECC','PRECL']   : unit_str = '[mm/day]'
   # res.tiXAxisString = unit_str
   # var_str = var
   # if var=='PRECT' : var_str = "Precipitation"
   # res.trXMinF = bin_ds['bins'].min().values
   # res.trXMaxF = bin_ds['bins'].max().values
   #-------------------------------------------------------------------------
   # Create plot
   #-------------------------------------------------------------------------
   ntime = len(time)
   
   num_clev = 21
   aboutZero = False
   # if tvar in ['U','V']: aboutZero = True
   if apply_filter: aboutZero = True
   min_val = np.min([np.min(h) for h in hov_list])
   max_val = np.max([np.max(h) for h in hov_list])
   clev_tup = ngl.nice_cntr_levels(min_val, max_val,       \
                                   cint=None, max_steps=num_clev, \
                                   returnLevels=False, aboutZero=aboutZero )
   if clev_tup==None: 
      print('SWITCHING TO AUTOMATIC CONTOUR LEVELS!')
      res.cnLevelSelectionMode = "AutomaticLevels"   
   else:
      cmin,cmax,cint = clev_tup
      res.cnLevels = np.linspace(cmin,cmax,num=num_clev)
      res.cnLevelSelectionMode = "ExplicitLevels"

   for c in range(num_case):
      ip = v*num_case+c if var_x_case else c*num_var+v
      tres = copy.deepcopy(res)
      # tres.tmXBMode      = 'Explicit'
      # tres.tmXBValues    = np.arange( len(lon_bins) )
      # tres.tmXBLabels    = lon_bins

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
# layout = [1,len(plot)]
layout = [num_var,num_case] if var_x_case else [num_case,num_var]
ngl.panel(wks,plot[0:len(plot)],layout,hs.setres_panel())
ngl.end()

hc.trim_png(fig_file)
#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------
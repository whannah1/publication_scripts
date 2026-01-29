#===================================================================================================
# Walter Hannah - Lawrence Livermore National Lab
# 
# Miscellaneous routines
#   get_host               figure out what machine we are using
#   printline              Prints a line
#   trim_png               crops a png image file
# 
# Data query/calculation routines:
#   print_stat             prints simple stats
#   print_time_length      print formatted time information
#   median                 Calculate median
# 
# Spherical Calculation routines:
#   find_nearest_neighbors_on_sphere_numba
#   find_nearest_neighbors_on_sphere
#   calc_great_circle_distance
#   calc_great_circle_bearing
#   calc_sphereical_triangle_area 
# 
# Binning Routines
#   bin_YbyX
#   bin_ZbyYX
# 
# Filtering
#   filter_wgts_lp_lanczos
# 
# Saturation routines
#   calc_Qsat
#   calc_Qsat_Peters
# 
# Vertical Integration Routines
#   calc_dp3d
# 
#===================================================================================================
import os, subprocess as sp, xarray as xr, numpy as np, dask, numba
#---------------------------------------------------------------------------------------------------
# Constants
#---------------------------------------------------------------------------------------------------
pi      = np.pi
cpd     = 1004.         # J / (kg K)       heat capavity at constant pressure
cpv     = 1850.         # J / (kg K)       heat capavity at constant volume
Tf      = 273.15        # K                freezing temperature of water
g       = 9.81          # m / s^2          acceleration due to gravity
Lv      = 2.5104e6      # J / kg           latent heat of vaporization / evaporation
Lf      = 0.3336e6      # J / kg           latent heat of freezing / melting
Ls      = 2.8440e6      # J / kg           latent heat of sublimation
Po      = 100000.       # Pa               reference pressure
Rd      = 287.04        # J / (kg K)       gas constant for dry air
Rv      = 461.5         # J / (kg K)       gas constant for water vapor
esf     = 611.          # Pa               ?
eps     = 0.622         # (epsilon)        ?

days_per_month = [31,28,31,30,31,30,31,31,30,31,30,31]

deg_to_rad = pi/180.
rad_to_deg = 180./pi

class tcolor:
   ENDC      = '\033[0m'
   BLACK     = '\033[30m'
   RED       = '\033[31m'
   GREEN     = '\033[32m'
   YELLOW    = '\033[33m' 
   BLUE      = '\033[34m'
   MAGENTA   = '\033[35m'
   CYAN      = '\033[36m'
   WHITE     = '\033[37m'
   BOLD      = '\033[1m'
   UNDERLINE = '\033[4m'

#---------------------------------------------------------------------------------------------------
# Miscellaneous Routines
#---------------------------------------------------------------------------------------------------
def get_host():
   """
   Get name of current machine
   """
   # First get some info
   try: 
      host = sp.check_output(["dnsdomainname"],universal_newlines=True).strip()
   except:
     host = None
   if host is not None:
      if 'nersc' in host : host = None
   if host is None or host=='' : host = os.getenv('host')
   if host is None or host=='' : host = os.getenv('HOST')
   opsys = os.getenv('os')
   # Make the final setting of host name
   if opsys=='Darwin'      : host = 'mac'
   if 'cori'   in host     : host = 'cori'
   if 'summit' in host     : host = 'olcf'
   if host=='ccs.ornl.gov' : host = 'olcf'   # rhea
   if host=='olcf.ornl.gov': host = 'olcf'   # andes
   return host
#-------------------------------------------------------------------------------   
def printline(n=80):
   """ 
   print a line of n characters 
   """
   print('-'*n)
#-------------------------------------------------------------------------------
def trim_png(fig_file):
   """ crop white space from png file """
   fig_file = fig_file+".png"
   fig_file = fig_file.replace(os.getenv('HOME')+'/Research/E3SM/','')
   if os.path.isfile(fig_file) :
      cmd = "convert -trim +repage "+fig_file+"   "+fig_file
      os.system(cmd)
      print("\n"+fig_file+"\n")
   else:
      raise FileNotFoundError(f'\ntrim_png(): {fig_file} does not exist?!\n')

#---------------------------------------------------------------------------------------------------
# Data Query Routines
#---------------------------------------------------------------------------------------------------
def print_stat(x,name='(no name)',unit='',fmt='f',stat='naxh',indent='',compact=False):
   """ Print min, avg, max, and std deviation of input """
   if fmt=='f' : fmt = '%.4f'
   if fmt=='e' : fmt = '%e'
   if unit!='' : unit = f'[{unit}]'
   name_len = 12 if compact else len(name)
   msg = ''
   line = f'{indent}{name:{name_len}} {unit}'
   # if not compact: print(line)
   if not compact: msg += line+'\n'
   for c in list(stat):
      if not compact: line = indent
      if c=='h' : line += '   shp: '+str(x.shape)
      if c=='a' : line += '   avg: '+fmt%x.mean()
      if c=='n' : line += '   min: '+fmt%x.min()
      if c=='x' : line += '   max: '+fmt%x.max()
      if c=='s' : line += '   std: '+fmt%x.std()
      # if not compact: print(line)
      if not compact: msg += line+'\n'
   # if compact: print(line)
   if compact: msg += line#+'\n'
   print(msg)
   return msg
#-------------------------------------------------------------------------------
def print_time_length(time,indent='    ',print_msg=True):
   """ print length of time """
   if len(time)>1:
      time_span_days   = ( time[-1] - time[0] + (time[1]-time[0]) ).values.astype('timedelta64[D]') 
      time_span_months = time_span_days.astype('timedelta64[M]') + 1
      time_span_years  = time_span_months.astype('timedelta64[Y]') + 1
      msg = indent
      msg = msg+f'Time length: {str(time_span_days)}'
      if time_span_days > np.timedelta64(60, 'D') : 
         msg = msg+f'  /  {time_span_months}'
      if time_span_months > np.timedelta64(12, 'M')  : 
         msg = msg+f'  /  {time_span_years}'
      if print_msg:
         print(msg)
         return
      else:
         return msg
#-------------------------------------------------------------------------------
def median(array, dim, keep_attrs=False, skipna=False, **kwargs):
   """ Runs a median on an dask-backed xarray.

   This function does not scale!
   It will rechunk along the given dimension, so make sure 
   your other chunk sizes are small enough that it 
   will fit into memory.

   :param DataArray array: An xarray.DataArray wrapping a dask array
   :param dim str: The name of the dim in array to calculate the median
   """
   import dask.array
   if type(array) is xr.Dataset: return array.apply(median, dim=dim, keep_attrs=keep_attrs, **kwargs)

   if not hasattr(array.data, 'dask'): return array.median(dim, keep_attrs=keep_attrs, **kwargs)

   array = array.chunk({dim:-1})
   axis = array.dims.index(dim)
   median_func = np.nanmedian if skipna else np.median
   blocks = dask.array.map_blocks(median_func, array.data, dtype=array.dtype, drop_axis=axis, axis=axis, **kwargs)

   new_coords={k: v for k, v in array.coords.items() if k != dim and dim not in v.dims}
   new_dims = tuple(d for d in array.dims if d != dim)
   new_attrs = array.attrs if keep_attrs else None

   return xr.DataArray(blocks, coords=new_coords, dims=new_dims, attrs=new_attrs)
#-------------------------------------------------------------------------------
def calc_t_stat( D0, D1, S0, S1, N0, N1, t_crit=1.96, verbose=True ):
   """ 
   calculate Student's t-statistic 
   t_crit = 1.96   2-tail test w/ inf dof & P=0.05
   t_crit = 2.5    2-tail test w/ 5 dof & P=0.05
   t_crit = 2.2    2-tail test w/ 10 dof & P=0.05
   """

   ### Standard error
   SE = np.sqrt( S0**2/N0 + S1**2/N1 )

   ### t-statistic - aX is the difference now
   t_stat = ( D1 - D0 ) / SE

   ### Degrees of freedom
   # DOF = (S0**2/N0 + S1**2/N1)**2 /( ( (S0**2/N0)**2 / (N0-1) )   \
   #                                  +( (S1**2/N1)**2 / (N1-1) ) )
   # print(f'  DoF min/max: {DOF.min():6.1f} / {DOF.max():6.1f}')

   # hc.print_stat(SE,name='SE',indent='    ')
   # hc.print_stat(t_stat,name='t statistic',indent='    ')

   ### Critical t-statistic
   

   sig_cnt = np.sum( np.absolute(t_stat) > t_crit )
   sig_pct = sig_cnt / t_stat.size *100

   if verbose: print(f'   SIG COUNT: {sig_cnt:8}   ({sig_pct:5.2f}% of {t_stat.size})')
   
   # for i in range(len(lat_bins)):
   #    msg = f'  lat: {lat_bins[i]}   t_stat: {t_stat[i]}   '
   #    if np.absolute(t_stat[i])>t_crit: msg = msg+tcolor.RED+'SIGNIFICANT'+tcolor.ENDC
   #    print(msg)

   return t_stat

#-------------------------------------------------------------------------------

#---------------------------------------------------------------------------------------------------
# Spherical Calculation routines
#---------------------------------------------------------------------------------------------------
# input should be radians
@numba.njit()
def find_nearest_neighbors_on_sphere_numba(lat,lon,num_neighbors,neighbor_id) :
   ncol = len(lat)
   for n in range(0,ncol) :
      cos_dist = np.sin(lat[n])*np.sin(lat[:]) + \
                 np.cos(lat[n])*np.cos(lat[:]) * np.cos( lon[n]-lon[:] )
      for nn in range(0,ncol+1) :
         if cos_dist[nn] >  1.0 : cos_dist[nn] = 1.0
         if cos_dist[nn] < -1.0 : cos_dist[nn] = -1.0
      dist = np.arccos( cos_dist )
      p_vector = np.argsort(dist,kind='mergesort')
      neighbor_id[n,0:num_neighbors] = p_vector[1:num_neighbors+1]
#-------------------------------------------------------------------------------
# input should be degrees
def find_nearest_neighbors_on_sphere(lat,lon,num_neighbors,return_sorted=False) :
   # return_sorted = False
   ncol = len(lat)
   neighbor_id    = np.empty([num_neighbors+1,ncol], dtype=np.int32)
   neighbor_dist  = np.empty([num_neighbors+1,ncol], dtype=np.float32)
   for n in range(0,ncol):
      dist = calc_great_circle_distance(lat[n],lat[:],lon[n],lon[:])
      if return_sorted:
         p_vector = np.argsort(dist,kind='mergesort')
         neighbor_id[:,n]   = p_vector[0:num_neighbors+1]     # include current point
         neighbor_dist[:,n] = dist[p_vector[0:num_neighbors+1]]
      else:
         p_vector = np.argpartition(dist[1:],num_neighbors)+1
         neighbor_id[:,n] = p_vector[0:num_neighbors]
   # return neighbor_id
   neighbor_ds = xr.Dataset()
   neighbor_ds['neighbor_id']   = (('ncol','neighbors'), neighbor_id)
   neighbor_ds['neighbor_dist'] = (('ncol','neighbors'), neighbor_dist)
   return neighbor_ds
#-------------------------------------------------------------------------------
# input should be in degrees
def calc_great_circle_distance(lat1,lat2,lon1,lon2):
   dlon = lon2 - lon1
   cos_dist = np.sin(lat1*deg_to_rad)*np.sin(lat2*deg_to_rad) + \
              np.cos(lat1*deg_to_rad)*np.cos(lat2*deg_to_rad)*np.cos(dlon*deg_to_rad)
   # print( str(cos_dist.min()) +"   "+ str(cos_dist.max()) )
   cos_dist = np.where(cos_dist> 1.0, 1.0,cos_dist)
   cos_dist = np.where(cos_dist<-1.0,-1.0,cos_dist)
   dist = np.arccos( cos_dist )
   return dist
#-------------------------------------------------------------------------------
# @numba.njit
# def calc_great_circle_bearing(lat1_in,lat2_in,lon1_in,lon2_in):
#    lat1 = lat1_in * deg_to_rad
#    lat2 = lat2_in * deg_to_rad
#    lon1 = lon1_in * deg_to_rad
#    lon2 = lon2_in * deg_to_rad

#    dlon = lon1 - lon2

#    atan_tmp1 = sin(lon2-lon1)*cos(lat2)
#    atan_tmp2 = cos(lat1)*sin(lat2)-sin(lat1)*cos(lat2)*cos(lon2-lon1)
#    bearing = atan2( atan_tmp1, atan_tmp2 )

#    bearing = bearing * rad_to_deg
#    return bearing

#-------------------------------------------------------------------------------
# def calc_sphereical_triangle_area (lat1,lat2,lat3,lon1,lon2,lon3):
#     ;;; compute great circle lengths
#     al = calc_great_circle_distance(lat2, lat3, lon2, lon3 ) 
#     bl = calc_great_circle_distance(lat1, lat3, lon1, lon3 ) 
#     cl = calc_great_circle_distance(lat1, lat2, lon1, lon2 ) 

#     ;;; compute angles
#     sina = sin( ( bl + cl - al ) /2. )
#     sinb = sin( ( al + cl - bl ) /2. )
#     sinc = sin( ( al + bl - cl ) /2. ) 
#     sins = sin( ( al + bl + cl ) /2. )

#     a = sqrt( (sinb*sinc) / (sins*sina) ) 
#     b = sqrt( (sina*sinc) / (sins*sinb) ) 
#     c = sqrt( (sina*sinb) / (sins*sinc) ) 

#     a1 = 2.*atan(a)
#     b1 = 2.*atan(b)
#     c1 = 2.*atan(c)

#     if ( a.gt.b .and. a.gt.c ) then
#         ; a1 = -2*atan(1/a)
#         a1 = -2.*atan2(1.,a)
#     else 
#         if (b.gt.c) then
#             ; b1 = -2.*atan(1./b)
#             b1 = -2.*atan2(1.,b)
#         else 
#             ; c1 = -2.*atan(1./c)
#             c1 = -2.*atan2(1.,c)
#         end if
#     end if

#     # apply Girard's theorem
#     area = totype( a1+b1+c1 ,typeof(lat1))
#     return area

#---------------------------------------------------------------------------------------------------
# Binning Routines
#---------------------------------------------------------------------------------------------------
def bin_YbyX (Vy,Vx,bins=[],bin_min=0,bin_max=1,bin_spc=1,bin_spc_log=20,nbin_log=2,
              bin_mode="manual",verbose=False,wgt=[],
              keep_time=False,keep_lev=False):
   """ Average Vy into bins of Vx values according to bins and bin_mode. 
   Manual mode takes an array bin center values and determines thebin spacings. 
   Explicit mode takes a list of bin edge values, which is useful for an 
   irregularly spacing, such as logarithmic.  """
   #----------------------------------------------------------------------------
   # Manual mode - use min, max, and spc (i.e. stride) to define bins
   if bin_mode == "manual":
      # bins = np.arange(bin_min, bin_max+bin_spc, bin_spc)
      # nbin = len(bin_coord)
      nbin    = np.round( ( bin_max - bin_min + bin_spc )/bin_spc ).astype(np.int)
      bins    = np.linspace(bin_min,bin_max,nbin)
      bin_coord = xr.DataArray( bins )
   #----------------------------------------------------------------------------
   # Explicit mode - requires explicitly defined bin edges
   if bin_mode == "explicit":
      bins = xr.DataArray( bins )   # convert the input to a DataArray
      nbin = len(bins)-1
      bin_coord = ( bins[0:nbin-1+1] + bins[1:nbin+1] ) /2.
   #----------------------------------------------------------------------------
   # Log mode - logarithmically spaced bins that increase in width by bin_spc_log [%]
   if bin_mode == "log":
      bin_log_wid = np.zeros(nbin_log)
      bin_log_ctr = np.zeros(nbin_log)
      bin_log_ctr[0] = bin_min
      bin_log_wid[0] = bin_spc
      for b in range(1,nbin_log):
         bin_log_wid[b] = bin_log_wid[b-1] * (1.+bin_spc_log/1e2)  # note - bin_spc_log is in %
         bin_log_ctr[b] = bin_log_ctr[b-1] + bin_log_wid[b-1]/2. + bin_log_wid[b]/2.
      nbin = nbin_log
      bin_coord = xr.DataArray( bin_log_ctr )
      ### Print coords
      # for b in bin_coord.values : print( '{0:8.2f}'.format(b) )
   #----------------------------------------------------------------------------
   # create output data arrays
   nlev  = len(Vy['lev'])  if 'lev'  in Vy.dims else 1
   ntime = len(Vy['time']) if 'time' in Vy.dims else 1

   if ntime==1 and keep_time==True : keep_time = False

   shape,dims,coord = (nbin,),'bin',[('bin', bin_coord)]
   
   if nlev>1 and keep_lev and not keep_time :
      coord = [ ('bin', bin_coord), ('lev', Vy['lev']) ]
      shape = (nbin,nlev)
      dims = ['bin','lev']   
   if nlev==1 and not keep_lev and not keep_time :
      shape,dims,coord = (nbin,),'bin',[('bin', bin_coord)]
   # if nlev>1 and keep_time==True :
   #    coord = [ ('bin', bin_coord), ('time', Vy['time']), ('lev', Vy['lev']) ]
   #    shape = (nbin,ntime,nlev)
   #    dims = ['bin','time','lev']   
   if nlev==1 and keep_time==True :
      coord = [('bin', bin_coord), ('time', Vy['time'])]
      shape = (nbin,ntime)
      dims = ['bin','time']
   
   mval = np.nan
   bin_val = xr.DataArray( np.full(shape,mval,dtype=Vy.dtype), coords=coord, dims=dims )
   bin_std = xr.DataArray( np.full(shape,mval,dtype=Vy.dtype), coords=coord, dims=dims )
   bin_cnt = xr.DataArray( np.zeros(shape,dtype=Vy.dtype), coords=coord, dims=dims )
   #----------------------------------------------------------------------------
   # # set up wgts
   # if len(wgt)==0 : 
   #    wgt = np.ones( Vy.shape )
   # else:
   #    wgt = np.ma.masked_array( wgt, condition)

   lev_chk=False
   if 'lev' in Vy.dims and len(Vy.lev)>1 and keep_lev : lev_chk = True

   # if 'lev' in Vy.dims : 
   #    if len(Vy.lev)>1 : 
   #       lev_chk = True
   #    else:
   #       print("!!!!!!!!!")
   # else:
   #    lev_chk=False

   if lev_chk :
      avg_dims = ['ncol']
      if 'time' in Vy.dims : avg_dims = ['time','ncol']
      avg_dims_wgt = ['ncol']
      # if len(wgt)!=0 : 
         # if 'time' in wgt.dims : 
         #    avg_dims_wgt = ['time','ncol']
         # else:
         #    avg_dims_wgt = ['ncol']

   val_chk = np.isfinite(Vx.values)

   use_masked_arrays = False
   # if dask.is_dask_collection(Vy) and not lev_chk:
   #    # use masked array method to avoid numpy division warnings on dask arrays
   #    use_masked_arrays = True
   #    mask = np.logical_and( np.isfinite(Vy.compute().values), np.isfinite(Vx.compute().values) )
   #    Vy_tmp = Vy.load().data
   #    Vx_tmp = Vx.load().data
   #    Vy_tmp = np.ma.masked_array( np.where(Vy_tmp,Vy_tmp,-999), mask)
   #    Vx_tmp = np.ma.masked_array( np.where(Vx_tmp,Vx_tmp,-999), mask)
   #----------------------------------------------------------------------------
   # Loop through bins
   for b in range(nbin):
      if bin_mode == "manual":
         bin_bot = bin_min - bin_spc/2. + bin_spc*(b  )
         bin_top = bin_min - bin_spc/2. + bin_spc*(b+1)
      if bin_mode == "explicit":
         bin_bot = bins[b]  .values
         bin_top = bins[b+1].values
      if bin_mode == "log":
         bin_bot = bin_log_ctr[b] - bin_log_wid[b]/2.
         bin_top = bin_log_ctr[b] + bin_log_wid[b]/2.

      condition = xr.DataArray( np.full(Vx.shape,False,dtype=bool), coords=Vx.coords )
      condition.values = ( np.where(val_chk,Vx.values,bin_bot-1e3) >=bin_bot ) \
                        &( np.where(val_chk,Vx.values,bin_bot-1e3)  <bin_top )

      if np.sum(condition)>0 :
         if lev_chk :
            ### xarray method
            if len(wgt)==0 : 
               bin_val[b,:] = Vy.where(condition,drop=True).mean( dim=avg_dims, skipna=True )
            else:
               if wgt.dims != Vy.dims : 
                  wgt, *__ = xr.broadcast(wgt, Vy) 
                  if 'time' in Vy.dims :
                     wgt = wgt.transpose('time','lev','ncol')
                  else :
                     wgt = wgt.transpose('lev','ncol')
               if 'time' in Vy.dims : 
                  bin_val[b,:] = ( (Vy*wgt).where(condition,drop=True).sum( dim='ncol', skipna=True ) \
                                      / wgt.where(condition,drop=True).sum( dim='ncol', skipna=True ) ).mean(dim='time', skipna=True )
               else:
                  bin_val[b,:] = ( (Vy*wgt).where(condition,drop=True).sum( dim='ncol', skipna=True ) \
                                      / wgt.where(condition,drop=True).sum( dim='ncol', skipna=True ) )
            bin_std[b,:] = Vy.where(condition,drop=True).std(  dim=avg_dims, skipna=True )
            bin_cnt[b,:] = Vy.where(condition,drop=True).count(dim=avg_dims)

            ### masked array method - Not sure how to make this work with a lev dimension...
            # bin_val[b,:] = np.sum( np.where(condition, Vy_tmp*wgt, 0 ) ) / np.sum( np.where(condition, wgt, 0 ) )
            # bin_val[b,:] = bin_val[b,:] / wgt.where(condition,drop=True).sum( dim=['time','ncol'])
            # bin_std[b,:] = np.where(condition, np.std( Vy_tmp ) )
            # bin_cnt[b,:] = np.sum( condition )
         elif keep_time and 'time' in Vy.dims:
            avg_dims = ['ncol']
            # bin_val[b,:] = Vy.where(condition,drop=True).mean( dim=avg_dims, skipna=True )
            # bin_cnt[b,:] = Vy.where(condition,drop=True).count(dim=avg_dims)
            if wgt.dims != Vy.dims : 
               wgt, *__ = xr.broadcast(wgt, Vy) 
               wgt = wgt.transpose('time','ncol')
            bin_val[b,:] = ( (Vy*wgt).where(condition,drop=True).sum( dim='ncol', skipna=True ) \
                                / wgt.where(condition,drop=True).sum( dim='ncol', skipna=True ) )
         else:
            
            if use_masked_arrays:
               ### numpy masked array method
               with np.errstate(divide='ignore', invalid="ignore"):
                  print("!")
                  # tmp = np.sum( np.where(condition, Vy_tmp*wgt, 0 ) )
                  # print("!")
                  # tmp = np.sum( np.where(condition, wgt, 0 ) )
                  bin_cnt[b] = np.sum( np.where(condition, 1, 0 ) )
                  exit()

               bin_val[b] = np.sum( np.where(condition, Vy_tmp*wgt, 0 ) ) / np.sum( np.where(condition, wgt, 0 ) )
               bin_std[b] = np.where(condition, np.std( Vy_tmp ) )
               bin_cnt[b] = np.sum( condition )
            else:
               ### xarray method - can't avoid divide warnings with dask arrays
               bin_val[b] = Vy.where(condition).mean(skipna=True)
               bin_std[b] = Vy.where(condition).std(skipna=True)
               bin_cnt[b] = np.sum( condition )
            
      # print('b: '+str(b)+'  cnt: '+str(bin_cnt[b])+'  val: '+str(bin_val[b]))
      if not lev_chk and verbose : print('b: {0:8.0f}  cnt: {1:8.0f}  val: {2:12.4e} '.format( b, bin_cnt[b].values, bin_val[b].values ) )
   #----------------------------------------------------------------------------
   ### add mask (doesn't propagate to xarray?)
   # bin_val = np.ma.masked_invalid(bin_val)
   # bin_std = np.ma.masked_invalid(bin_std)
   #----------------------------------------------------------------------------
   # use a dataset to hold all the output
   dims = ('bins',)
   if lev_chk and not keep_time : dims = ('bins','lev')
   if not lev_chk and keep_time : dims = ('bins','time')
   if lev_chk and keep_time     : dims = ('bins','time','lev')
   
   bin_ds = xr.Dataset()
   bin_ds['bin_val'] = (dims, bin_val )
   bin_ds['bin_std'] = (dims, bin_std )
   bin_ds['bin_cnt'] = (dims, bin_cnt )
   bin_ds['bin_pct'] = (dims, bin_cnt/bin_cnt.sum()*1e2 )
   bin_ds.coords['bins'] = ('bins',bin_coord)
   if lev_chk : bin_ds.coords['lev'] = ( 'lev', xr.DataArray(Vy['lev']) )

   #----------------------------------------------------------------------------
   return bin_ds
#-------------------------------------------------------------------------------
def bin_ZbyYX (Vz,Vy,Vx,binsy,binsx,opt):
   """ """
#-------------------------------------------------------------------------------


#---------------------------------------------------------------------------------------------------
# Filtering routines
#---------------------------------------------------------------------------------------------------
def filter_wgts_lp_lanczos(window_len, fc_lp=0, fc_hp=0, sigma_factor=1):
   """
   Calculate weights for a low pass Lanczos filter taken from:
   Duchon C. E. (1979) Lanczos Filtering in One and Two Dimensions. 
      Journal of Applied Meteorology, Vol 18, pp 1016-1022.
   note: 0 < fc_lp < fc_lp < 0.5
   Args:
      window_len: The length of the filter window
      fc_lp: low-pass cutoff frequency in inverse time steps (set to zero to disable)
      fc_hp: high-pass cutoff frequency in inverse time steps (set to zero to disable)
   """
   order = ((window_len - 1) // 2 ) + 1
   nwts = 2 * order + 1
   w = np.zeros([nwts])
   n = nwts // 2
   k = np.arange(1., n)
   # Sigma factor to remove Gibbs oscillation
   sigma = np.sin(np.pi * k / n) * n / (np.pi * k)
   # specify filter weights
   factor_lp = np.sin(2. * np.pi * fc_lp * k) / (np.pi * k)
   factor_hp = np.sin(2. * np.pi * fc_hp * k) / (np.pi * k)
   w[n] = 2 * ( fc_lp - fc_hp )
   w[n-1:0:-1] = ( factor_lp - factor_hp ) * sigma
   w[n+1:-1]   = ( factor_lp - factor_hp ) * sigma
   return w[1:-1]
#---------------------------------------------------------------------------------------------------
# Saturation routines
#---------------------------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# Calculate saturation specific humidity (kg/kg) from temperature (k) and pressure (mb)
# Buck Research Manual (1996)
def calc_Qsat (Ta,Pa) :
   Tc = Ta - 273.
   ew = 6.1121*(1.0007+3.46e-6*Pa)*np.exp((17.502*Tc)/(240.97+Tc))       # in mb
   qs = 0.62197*(ew/(Pa-0.378*ew))                                       # mb -> kg/kg
   return qs
#-------------------------------------------------------------------------------
# Calculate saturation specific humidity (kg/kg) from temperature (k) and pressure (mb)
def calc_Qsat_Peters (Ta,Pa):
    qs = (380./Pa)*np.exp(17.27*(Ta-273.0)/(Ta-36.0))
    return qs

#---------------------------------------------------------------------------------------------------
# Vertical Integration routines
#---------------------------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
def calc_dp3d(ps,lev):
   """
   calculate pressure thickness for vertical integration
   both inputs should be xarray DataArray type
   ps  = surface pressure
   lev = pressure levels of the data (not the model coordinate)
   """

   ps3d,pl3d = xr.broadcast(ps,lev*100.)

   pl3d = pl3d.transpose('time','lev','ncol')
   ps3d = ps3d.transpose('time','lev','ncol')

   nlev = len(lev)
   tvals = slice(0,nlev-2)
   cvals = slice(1,nlev-1)
   bvals = slice(2,nlev-0)

   # Calculate pressure thickness
   dp3d = xr.DataArray( np.full(pl3d.shape,np.nan), coords=pl3d.coords )
   dp3d[:,nlev-1,:] = ps3d[:,nlev-1,:].values - pl3d[:,nlev-1,:].values
   dp3d[:,cvals,:] = pl3d[:,bvals,:].values - pl3d[:,tvals,:].values

   # Deal with cases where levels are below surface pressure
   condition = pl3d[:,cvals,:].values<ps3d[:,cvals,:].values
   new_data  = ps3d[:,bvals,:].values-pl3d[:,tvals,:].values
   dp3d[:,cvals,:] = dp3d[:,cvals,:].where( condition, new_data )

   # Screen out negative dp values
   dp3d = dp3d.where( dp3d>0, np.nan )

   return dp3d

#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------

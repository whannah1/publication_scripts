#===================================================================================================
# Walter Hannah - Lawrence Livermore National Lab
# 
# 
#===================================================================================================
import hapy_common as hc
import os, glob, subprocess as sp, datetime, re, ngl
from dateutil.relativedelta import relativedelta
import xarray as xr, dask, numpy as np, pandas as pd

valid_obs_list = ['TRMM','GPCP','ERAi','ERA5','MAC','GPM','IMERG','NOAA']

default_data_dir = os.getenv('HOME')+'/Data/Obs'
default_data_sub = 'daily'

#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------
def fix_coord_names(ds,case=None,lat_name=None,lon_name=None):
   if case is not None:
      lat_name = case.lat_name
      lon_name = case.lon_name
   else:
      if lat_name is None or lon_name is None: 
         raise ValueError('lat and lon names must be provided if case object is omitted')
   if 'nlat'      in ds.coords : ds = ds.rename({'nlat':lat_name})
   if 'nlon'      in ds.coords : ds = ds.rename({'nlon':lon_name})
   if 'latitude'  in ds.coords : ds = ds.rename({'latitude':lat_name})
   if 'longitude' in ds.coords : ds = ds.rename({'longitude':lon_name})
   if 'level'     in ds.coords : ds = ds.rename({'level':'lev'})
   if case is not None:
      setattr(case, 'lat', ds[lat_name] )
      setattr(case, 'lon', ds[lon_name] )
   return ds
#---------------------------------------------------------------------------------------------------
# Case class
#---------------------------------------------------------------------------------------------------
class Case(object):
   ''' Class for loading Obs data '''
   def __init__(self, name, verbose=False ,lev=-1  \
               ,time_freq=None                \
               ,data_dir=default_data_dir          \
               ,data_sub=default_data_sub          ):
      global default_data_dir, default_data_sub
      ### Initialize case properites
      self.name = name
      self.short_name = self.name
      self.verbose = verbose
      self.data_dir = default_data_dir if data_dir is None else data_dir 
      self.data_sub = default_data_sub if data_sub is None else data_sub 
      self.hist_path = ''#self.data_dir+'/'+self.name+'/'+self.data_sub
      self.lev = lev
      self.obs = True
      self.grid = ''
      self.time_freq = time_freq
      self.lat_name = 'lat'
      self.lon_name = 'lon'
      # self.ncol_name = 'ncol'
      # self.area_name = 'area'
      #----------------------------------------------------
      #----------------------------------------------------
      if name=='GPCP':
         self.data_sub=''
         self.short_name = name
      if name=='TRMM':
         self.data_sub='monthly/'
         self.short_name = name
      if self.name=='NOAA': 
         self.time_freq = 'daily'
      #----------------------------------------------------
      # set the file paths
      #----------------------------------------------------
      if self.time_freq=='monthly': self.file_path_monthly = f'{self.data_dir}/{self.name}/monthly/*'
      if self.time_freq=='daily'  : self.file_path_daily   = f'{self.data_dir}/{self.name}/daily/*'

      if self.name=='GPCP': 
         if self.time_freq=='monthly': self.file_path_monthly = self.data_dir+self.name+'/precip.mon.mean.nc'
         if self.time_freq=='daily'  : self.file_path_daily   = self.data_dir+self.name+'/daily/GPCP.daily.*nc'
      if self.name=='TRMM': 
         if self.time_freq=='monthly': self.file_path_monthly = self.data_dir+self.name+'/monthly/3B43*'
         if self.time_freq=='daily'  : self.file_path_daily   = self.data_dir+self.name+'/daily/3B42_Daily*.nc4'
      if self.name=='ERAi': 
         if self.time_freq=='monthly': self.file_path_monthly = self.data_dir+self.name+'/data/*monthly*'
         if self.time_freq=='daily'  : self.file_path_daily   = self.data_dir+self.name+'/data/*6hour*'
      # if self.name=='ERA5': 
      #    if self.time_freq=='monthly': self.file_path_monthly = self.data_dir+self.name+'/monthly/*monthly*'
      #    if self.time_freq=='daily'  : self.file_path_daily   = self.data_dir+self.name+'/daily/*daily*'
      if self.name=='MAC': 
         if self.time_freq=='monthly': self.file_path_monthly = f'{self.data_dir}/{self.name}/monthly_twp/*'
         if self.time_freq=='daily'  : self.file_path_daily   = f'{self.data_dir}/{self.name}/daily_twp/*'

      self.file_path = self.data_dir+'/'+self.data_sub+'/*'
      # self.file_path = f'{self.data_dir}/{self.name}/{self.data_sub}'

      if self.time_freq is not None:
         if self.time_freq=='monthly': self.file_path = self.file_path_monthly
         if self.time_freq=='daily'  : self.file_path = self.file_path_daily

      # exit(self.file_path)
      
      if self.name=='GPM': self.file_path += '.nc4'
      #----------------------------------------------------
      ### add coordinates as attributes   
      #----------------------------------------------------
      # print()
      # print('  data_dir : '+self.data_dir)
      # print('  data_sub : '+self.data_sub)
      # print('  file_path: '+self.file_path)
      # print(); print(glob.glob(self.file_path))
      # print()
      # exit()
      ds = xr.open_dataset( sorted(glob.glob(self.file_path))[0] )
      ds = fix_coord_names(ds,case=self)
   #----------------------------------------------------------------------------
   def __str__(self):
      # return "%s" % (self.name)
      indent = '    '
      str_out = f'\n{self.name}\n'
      fmt_key_len = 18
      for key in self.__dict__.keys(): 
         attribute = getattr(self,key)
         # if attribute!='' and not isinstance(attribute, xr.Dataset) :
         #    str_out += f'{indent}{key:{fmt_key_len}}:  {attribute}\n'
         str_out += f'{indent}{key:{fmt_key_len}}:  {attribute}\n'
      return str_out
   #----------------------------------------------------------------------------
   def __repr__(self):
      return {'name : ':self.name \
             ,'grid : ':self.grid }
   #----------------------------------------------------------------------------
   def get_hist_file_list(self):
      ''' '''
      files = sp.check_output(f'ls {self.file_path}/*.nc'           \
                              ,shell=True,stderr=sp.STDOUT \
                              ,universal_newlines=True ).split('\n')
      return files
   #----------------------------------------------------------------------------
   def get_hist_freq(self):
      ''' '''
      # ?
   #----------------------------------------------------------------------------
   def get_hist_vars(self):
      ''' '''
      # ?
   #----------------------------------------------------------------------------
   def get_scrip(self):
      ''' return scrip file as Xarray dataset '''
      # sfile = self.grid+'_scrip.nc'
      # sdir = home+'/E3SM/data_grid/'
      # if not os.path.isfile(sdir+sfile) : 
      #    sdir = home+'/Research/E3SM/data_grid/'
      # if not os.path.isfile(sdir+sfile) :
      #    print('\nERROR: get_scrip(): scrip file does not exist! ')
      #    print('case: '+self.name)
      #    print('grid: '+self.grid+'\n')
      #    exit()
      # scripfile = xr.open_dataset(sdir+sfile)
      scripfile = ''
      return scripfile
   #----------------------------------------------------------------------------
   def set_coord_names(self,var):
      # # dynamics grid data
      # if 'DYN_' in var or var in ['VOR','DIV'] : 
      #    physgrid_list = ['pg2','pg3','pg4']
      #    if any(g in self.name for g in physgrid_list) : 
      #       # update grid name
      #       for g in physgrid_list:
      #          if g in self.name : self.grid = self.grid.replace(g,'np4')
      #       # update area and ncol names
      #       self.ncol_name,self.area_name = 'ncol_d','area_d'
      # # regional subset
      # if '_to_' in var:
      #    self.ncol_name = var.replace(var.split('_')[0],'ncol')
      return
   #----------------------------------------------------------------------------
   # def get_mask(self,var):
   #    ''' Create mask that matches the dimensions of input 'var' '''
      # htype = 'h0'
      # file_path = self.hist_path+self.name+'.cam.'+htype+'*'
      # ds = xr.open_mfdataset( file_path )
      
      # mask = xr.DataArray( np.ones(ds[var].shape,dtype=bool), \
      #                      coords=ds[var].coords, dims=ds[var].dims )
      # if hasattr(self,'lat1'): mask = mask & (ds['lat']>=self.lat1)
      # if hasattr(self,'lat2'): mask = mask & (ds['lat']<=self.lat2)
      # if hasattr(self,'lon1'): mask = mask & (ds['lon']>=self.lon1)
      # if hasattr(self,'lon2'): mask = mask & (ds['lon']<=self.lon2)
   #----------------------------------------------------------------------------
   def get_mask(self,ds,var):
      ''' Create mask with only ncol dimension from first available file '''
      mask = xr.DataArray( np.ones(ds[var].shape,   \
                           dtype=bool),               \
                           coords=ds[var].coords,   \
                           dims=ds[var].dims )
      if hasattr(self,'lat1'): mask = mask & (ds['lat']>=self.lat1)
      if hasattr(self,'lat2'): mask = mask & (ds['lat']<=self.lat2)
      if hasattr(self,'lon1'): mask = mask & (ds['lon']>=self.lon1)
      if hasattr(self,'lon2'): mask = mask & (ds['lon']<=self.lon2)

      return mask
   #----------------------------------------------------------------------------
   def extract_time(self,ds):
      if self.name=='TRMM':
         for item in ds.FileHeader.split(';\n'):
            if 'StartGranuleDateTime' in item : time1 = item.split('=')[1]
            if 'StopGranuleDateTime'  in item : time2 = item.split('=')[1]
         fmt = '%Y-%m-%dT%H:%M:%S.%fZ'
         time1 = datetime.datetime.strptime(time1,fmt)
         time2 = datetime.datetime.strptime(time2,fmt)
         avg_time = time1 + (time1-time2)/2.
         ds = ds.assign_coords(time=avg_time)
         # if self.time_freq=='monthly': var_name = 'precipitation'
         # if self.time_freq=='daily'  : var_name = 'precip'
         ds['precipitation'] = ds['precipitation'].expand_dims('time')
         # print(ds)
         # exit()
         return ds
         # return ds.assign(time=avg_time)
   #----------------------------------------------------------------------------
   def extract_time_monthly(self,ds):
      if self.name=='TRMM':
         nmonth = len(ds['month'])
         date = ds.attrs['date']
         fmt = '%Y%m%d'
         time0 = datetime.datetime.strptime(date,fmt)
         new_time = [time0]*nmonth
         for m in range(nmonth) :
            time1 = time0+relativedelta(months=m)
            time2 = time1+relativedelta(months=1) 
            new_time[m] = time1 + (time2-time1)/2.
         del(ds.attrs['date'])
         ds = ds.assign_coords(time=new_time)
         return ds
   #----------------------------------------------------------------------------
   def load_data(self,var,htype='',ps_htype=None,lev=np.array([0]),  \
                 num_years=-1, years=[], months=[],      \
                 num_files=0, first_file=0,              \
                 component='eam',use_remap=False,remap_str='.remap') :
      ''' '''
      import re # why do I need to do this here to avoid errors?
      #----------------------------------------------------
      ### set the file paths
      # if self.time_freq=='monthly': file_path = self.file_path_monthly
      # if self.time_freq=='daily'  : file_path = self.file_path_daily

      file_path = self.file_path

      if self.name=='ERAi': 
         file_path = file_path+f'.{var}.199[0-4]*'
         print('WARNING - only using 1990-1994 for ERAi - WARNING')
      #----------------------------------------------------
      ### adjust the mask
      # if htype=='h0' :
      #    if self.lon1> 180 : self.lon1 = self.lon1-360
      #    if self.lon1<-180 : self.lon1 = self.lon1+360
      #    if self.lon2> 180 : self.lon2 = self.lon2-360
      #    if self.lon2<-180 : self.lon2 = self.lon2+360
      #----------------------------------------------------
      ### open the file
      if len(glob.glob(file_path))>0 : 
         # print(); exit(file_path)
         file_list = sorted(glob.glob(file_path))
         # print(); print(file_list)
         ### Remove "remapped" files unless use_remap is true
         file_list_tmp = file_list.copy()
         for file_name in file_list : 
            if use_remap and remap_str not in file_name: file_list_tmp.remove(file_name)
            if not use_remap and remap_str in file_name: file_list_tmp.remove(file_name)
         file_list = file_list_tmp
         # print(); print(f'use_remap: {use_remap}  remap_str: {remap_str}  file_path: {file_path}')
         # print(); print(file_list); exit()
         # Limit number of files if num_files!=0
         if num_files>0 : file_list = file_list[first_file:first_file+num_files]     # use initial files
         if num_files<0 : file_list = file_list[num_files:]      # use latest files
         if self.name in ['TRMM'] : 
            if htype=='h0' :
               # ds = xr.open_mfdataset( file_list, concat_dim='time', combine='by_coords', preprocess=self.extract_time_monthly )
               ds = xr.open_mfdataset( file_list, concat_dim='time', combine='by_coords', preprocess=self.extract_time )
            else:
               ds = xr.open_mfdataset( file_list, concat_dim='time', combine='by_coords', preprocess=self.extract_time )
         else:
            # print(); print(file_list)
            ds = xr.open_mfdataset( file_list, combine='by_coords'  )
      else :
         print('\nERROR: load_data(): No history files found')
         print('file path: '+file_path+'\n')
         exit()

      ds = fix_coord_names(ds,case=self)
      #----------------------------------------------------
      # Correct timestamp - shift by mean delta t
      if self.name=='TRMM':
         if htype!='h0' :
            dtime = ds['time'].diff(dim='time').values.astype('timedelta64[h]').mean()
            ds['time'] = pd.to_datetime(ds.time.get_index('time')) - dtime*0.5
      #----------------------------------------------------
      ### map variable names to E3SM convention
      tvar = var
      if self.name=='GPCP' :
         if var in ['PRECT','PRECC','PRECL'] : tvar = 'precip'
      if self.name=='TRMM' :
         if 'nlat' in ds.coords : ds = ds.rename({'nlat':self.lat_name,'nlon':self.lon_name})
         if var in ['PRECT','PRECC','PRECL'] : 
            # if self.time_freq=='monthly':  tvar = 'precip'
            if self.time_freq=='monthly':  tvar = 'precipitation'
            if self.time_freq=='daily'  :  tvar = 'precipitation'
      if var=='area'  : tvar = self.lat_name
      
      if self.name=='MAC':
         if var=='TGCLDLWP': tvar = 'cloudlwp'
         # if var=='TGCLDLWP': tvar = 'totallwp'
         if var=='TWP':      tvar = 'totallwp'
         if var=='CWP':      tvar = 'cloudlwp'

      if self.name=='GPM' or self.name=='IMERG':
         if var=='PRECT': tvar = 'precipAvg'

      # if self.name=='ERAi' or self.name=='ERA5' :
      if self.name=='ERAi':
         # check for numeric characters
         if re.search(r'\d', var) is None:
            contains_lev = False
            ivar = var
         else:
            contains_lev = True
            ivar = re.sub('[0-9]', '', var)

         if ivar=='MSE'     : tvar = 'T'
         if ivar=='TS'      : tvar = 'skt'
         if ivar=='PS'      : tvar = 'sp'
         if ivar=='T'       : tvar = 't'
         if ivar=='Q'       : tvar = 'q'
         if ivar=='U'       : tvar = 'u'
         if ivar=='TGCLDLWP': tvar = 'tclw'
         if ivar=='TGCLDIWP': tvar = 'tciw'
         if ivar=='TMQ'     : tvar = 'tcw'
         if ivar=='CWV'     : tvar = 'tcw'
         if ivar=='FLNT'    : tvar = 'ttr'
         if ivar=='FSNT'    : tvar = 'tsr'
         if ivar=='LHFLX'   : tvar = 'mslhf'
         if ivar=='SHFLX'   : tvar = 'msshf'
         # if ivar=='Q850': tvar = ''
         # if ivar=='T850': tvar = ''
         # if ivar=='U850': tvar = ''

         if contains_lev: 
            lev = np.array( int( ''.join(c for c in var if c.isdigit()) ) )
            

      #----------------------------------------------------
      if tvar not in ds :
         print(ds)
         print('\nERROR: load_data: '+var+' ('+tvar+') not found in dataset\n')
         if 'contains_lev' in locals(): 
            print(f'  contains_lev is True - lev = {lev}')
         if 'ivar' in locals():
            print(f'  var:{var}  ivar:{ivar}  tvar:{tvar}')
         else:
            print(f'  var:{var}  tvar:{tvar}')
         exit()
      #----------------------------------------------------
      ### Load and subset the data if lat and lon limits are set
      X = ds[tvar].where(self.get_mask(ds,tvar),drop=True)
      #----------------------------------------------------
      ### derived variable cases
      # if var=='MSE' : 
      #    X = X + ds['Z3'].where( self.get_mask() ,drop=True) * hc.g /hc.cpd \
      #          + ds[ 'Q'].where( self.get_mask() ,drop=True) * hc.Lv/hc.cpd 
      #    X['long_name'] = 'Moist Static Energy'
      #    X['units']     = 'K'
      if var=='area' :
         re   = 6.37122e06                               # radius of earth
         dlon = ( ds[self.lon_name][2] - ds[self.lon_name][1] )*np.pi/180.     # assume dlon is constant
         dlat = ( ds[self.lat_name][2] - ds[self.lat_name][1] )*np.pi/180.     # assume dlat is constant
         dlon = np.absolute(dlon)
         dlat = np.absolute(dlat)
         X = xr.DataArray( np.empty([len(ds[self.lat_name]),len(ds[self.lon_name])]), \
                           coords=[('lat',ds[self.lat_name]),('lon',ds[self.lon_name])], \
                           dims=['lat','lon'] )
         for j in range(len(ds[self.lat_name])):
            X[j,:] = ( re * dlat * dlon * np.cos(ds[self.lat_name][j]*np.pi/180.) ).values
      #----------------------------------------------------
      ### Adjust units   
      if self.name=='TRMM' and tvar=='precipitation' : 
         if self.time_freq=='monthly':
            X = X*24
            X['unts'] = 'mm/day'
      if self.name=='ERA5' and var=='TS' : 
         if self.time_freq=='monthly':
            X = X-273
            X['unts'] = 'DegC'
      #----------------------------------------------------
      ### ignore time dimension for grid variables
      # if var in ['lat','lon','area'] : X = X[0,:]
      #----------------------------------------------------
      ### Subset in time
      if 'time' in X.dims :
         time_mask = xr.DataArray( np.ones([len(X.time)],dtype=bool), coords=[('time', X.time)], dims='time' )
         if num_years>=0 : time_mask = time_mask & ( ( X['time.year']-X['time.year'][0]) <=num_years)
         if years !=[]   : time_mask = time_mask & [ y in years  for y in (X['time.year']-X['time.year'][0]).values ]
         if len(months)>0: time_mask = time_mask & [ m in months for m in X['time.month'].values ]
         X = X.sel( time=X.time.where(time_mask,drop=True) )
      #----------------------------------------------------
      ### adjust dimensions
      if self.name=='TRMM' and var not in ['area','lat','lon'] :
         X = X.transpose('time','lat','lon')
      #----------------------------------------------------
      return X
   #----------------------------------------------------------------------------

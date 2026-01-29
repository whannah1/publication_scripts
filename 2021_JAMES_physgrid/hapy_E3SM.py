#===================================================================================================
# Walter Hannah - Lawrence Livermore National Lab
# 
# 
#===================================================================================================
import os, errno, glob, datetime, ngl, numba
import subprocess as sp
import xarray as xr, dask, numpy as np
import hapy_common as hc, hapy_obs    as ho
home = os.getenv("HOME")
#---------------------------------------------------------------------------------------------------
# Default values
#---------------------------------------------------------------------------------------------------
verbose = False
default_data_dir  = home+"/Data/E3SM/"
default_data_sub  = "atm/"
default_time_freq = 'monthly'
# E3SM_default_lev = np.array([30,50,75,100,125,150,200,250,300,350,400,450,500,550,600,650,700,750,800,825,850,875,900,925,950,975,1000])
E3SM_default_lev = np.array([50,100,150,200,300,400,500,600,700,800,850,925,975,1000])
hist_type_list = ['h0','h1','h2','h3','h4','hgb','hgb0']

#----------------------------------------------------
# machine specific defaults
host = hc.get_host()
if host=='cori':
   default_data_dir='/global/homes/w/whannah/E3SM/scratch/'
   default_data_sub='run/'
if host=='olcf':
   default_data_dir='/ccs/home/hannah6/E3SM/scratch/'
   # default_data_dir='/autofs/nccs-svm1_home1/hannah6/E3SM/scratch/'
   default_data_sub='run/'
if host=='mac':
   default_data_dir=os.getenv('HOME')+'/E3SM/scratch/'
   default_data_sub='run/'

#----------------------------------------------------
# E3SM parameters from cime/src/share/util/shr_const_mod.F90
SHR_CONST_G       = 9.80616
SHR_CONST_AVOGAD  = 6.02214e26
SHR_CONST_BOLTZ   = 1.38065e-23
SHR_CONST_MWDAIR  = 28.966
SHR_CONST_RGAS    = SHR_CONST_AVOGAD*SHR_CONST_BOLTZ
SHR_CONST_RDAIR   = SHR_CONST_RGAS/SHR_CONST_MWDAIR

#---------------------------------------------------------------------------------------------------
# Vertical interpolation
#---------------------------------------------------------------------------------------------------
def interpolate_to_pressure(ds,data_mlev=None,
                           var_name=None,lev=np.array([0]),
                           ds_ps=None,ps_var='PS',
                           interp_type=2,extrap_flag=False,):
   """
   Interpolate given variable to pressure levels specified by lev
   Either set var_name to load from ds or provide a dataArray for data_mlev
   interp_type => 1=LINEAR, 2=LOG, 3=LOG LOG
   """
   if data_mlev is None and var_name is None:
      raise ValueError('var_name and data_mlev inputs cannot both be None')
   if data_mlev is None: data_mlev = ds[var_name]
   if 'lev' not in data_mlev.coords: 
      raise ValueError('lev coordinate is missing from input data')
   if ds_ps is None: ds_ps = ds
   if np.all( lev > 0 ) :
      if type(lev)==type(1) : lev = np.array([float(lev)])       # If lev is single integer then convert to list
      hya, hyb = ds['hyam'], ds['hybm']
      if 'time' in hya.dims : hya = hya.isel(time=0).values
      if 'time' in hyb.dims : hyb = hyb.isel(time=0).values
      #-------------------------------------------------------------------------
      PS_dum = ds_ps[ps_var]
      # PS_dum = PS_dum.where( self.get_mask(ds),drop=True)
      if 'P0' in ds.variables: 
         P0 = ds['P0']/1e2
      else:
         P0 = xr.DataArray(1e3)
      if 'time' in P0.dims : P0 = P0.isel(time=0).values
      #-------------------------------------------------------------------------
      # Create empty array with new lev dim
      data_plev = xr.full_like( data_mlev.isel(lev=0), np.nan ).drop('lev')   
      data_plev = data_plev.expand_dims(dim={'lev':lev}, axis=data_mlev.get_axis_num('lev'))
      data_plev.values = np.full(data_plev.shape,np.nan)
      #-------------------------------------------------------------------------
      # Add dummy dimension if not lat/lon data
      if not ('lat' in data_mlev.dims and 'lon' in data_mlev.dims):
         PS_dum = PS_dum.expand_dims(dim='dummy',axis=len(ds_ps[ps_var].dims))
         data_mlev = data_mlev.expand_dims(dim='dummy',axis=len(data_mlev.dims))
      #-------------------------------------------------------------------------
      # Do the interpolation
      data_tmp = ngl.vinth2p( data_mlev, hya, hyb, lev, PS_dum.values, 
                              interp_type, P0, 1, extrap_flag)
      # Remove the dummy dimension
      if 'dummy' in data_mlev.dims: data_tmp = data_tmp[:,:,:,0]
      # data_plev.data = dask.array.from_array( data_tmp ,chunks='auto')
      # data_plev.data = np.ma.masked_values(data_plev.data,1e30)
      data_plev.values = np.ma.masked_values( data_tmp ,1e30)
      #-------------------------------------------------------------------------
      return data_plev
   else:
      raise ValueError('input levels cannot be negative')
#---------------------------------------------------------------------------------------------------
# Case class initiator
#---------------------------------------------------------------------------------------------------
def Case( name, verbose=False, lev=-1
         ,populate_files=True
         ,time_freq='monthly'
         ,data_dir=None, data_sub=None ,grid=None 
         ,atm_comp='eam'):
   #----------------------------------------------------
   # Divert to Obs module for obs cases
   if name in ['TRMM','GPCP','ERAi','ERA5','MAC'] : 
      return ho.Case( name, verbose, lev, time_freq)
   else:
      if data_dir is None: data_dir = default_data_dir
      if data_sub is None: data_sub = default_data_sub
      # if (lev==-1).all() : lev = E3SM_default_lev
      if np.all(lev==-1) : lev = E3SM_default_lev
      return Case_E3SM( name, verbose, populate_files, lev, time_freq, data_dir, data_sub, atm_comp=atm_comp )

#---------------------------------------------------------------------------------------------------
# Case class
#---------------------------------------------------------------------------------------------------
class Case_E3SM(object):
   """ Class for loading E3SM data """
   def __init__(self, name, verbose=False
               ,populate_files=True
               ,lev=E3SM_default_lev
               ,time_freq=default_time_freq
               ,data_dir=default_data_dir
               ,data_sub=default_data_sub       
               ,grid=None
               ,atm_comp='eam' ):
      #----------------------------------------------------
      # Initialize case properites
      self.name = name
      self.short_name = self.name
      self.verbose = verbose
      self.atm_comp = atm_comp
      self.data_dir = data_dir
      self.data_sub = data_sub
      self.hist_path = self.data_dir+'/'+self.name+'/'+self.data_sub
      self.grid = grid
      # self.lev = lev
      self.time_freq = time_freq
      self.obs = False
      self.ncol_name = 'ncol'
      self.area_name = 'area'
      self.mirror_equator = False  # flag to turn on loading zonal bands on both sides of equator
      #----------------------------------------------------
      if 'CESM' in name: 
         self.data_dir  = self.data_dir.replace('E3SM','CESM')
         self.hist_path = self.hist_path.replace('E3SM','CESM')
         self.atm_comp = 'cam'
      #----------------------------------------------------
      # special cases
      if name=='earlyscience.FC5AV1C-H01A.ne120.sp1_64x1_1000m.20190329':
         if host=='olcf': 
            self.data_dir='/gpfs/alpine/proj-shared/cli115/crjones/e3sm/'
            self.data_sub='run/' 
         self.hist_path = self.data_dir+self.name+'/'+self.data_sub
         self.short_name = 'EarlyScience2018'
      if name=='earlyscience.FC5AV1C-H01A.ne120.sp1_64x1_1000m':
         if host=='cori': 
            self.data_dir='/project/projectdirs/m3312/crjones/e3sm/'
            self.hist_path = self.data_dir+'early_science/monthly_hist/'
         if host=='olcf': 
            self.data_dir='/gpfs/alpine/proj-shared/cli115/crjones/e3sm/'
            self.hist_path = self.data_dir+'earlyscience.FC5AV1C-H01A.ne120.sp1_64x1_1000m*/run/'
         self.data_sub='' 
         self.short_name = 'SP-E3SM ne120np4 (ES)'
      if name=='earlyscience.FC5AV1C-L.ne30.sp1_64x1_1000m.20190415':
         self.data_dir='/project/projectdirs/m3312/crjones/e3sm/'
         self.data_sub='' 
         self.hist_path = self.data_dir+self.name+'/'+self.data_sub
         self.short_name = 'SP-E3SM'
      if name=='earlyscience.FC5AV1C-L.ne30.E3SM.20190519':
         self.data_dir='/project/projectdirs/m3312/whannah/'
         self.data_sub='atm/' 
         self.hist_path = self.data_dir+self.name+'/'+self.data_sub
         self.short_name = 'E3SM'
      if name=='earlyscience.FC5AV1C-H01A.ne120.E3SM.20190329':
         self.data_dir='/project/projectdirs/m3312/crjones/e3sm/'
         self.data_sub='' 
         self.hist_path = self.data_dir+'/early_science_e3sm/'+self.data_sub
         self.short_name = 'E3SM ne120np4 (ES)'
      #----------------------------------------------------
      # use history files to fill in the details
      if populate_files:
         if '*' in self.hist_path :
            path_chk = os.path.exists( glob.glob(self.hist_path)[0] )
         else:
            path_chk = os.path.exists( self.hist_path )
         if path_chk :
            files = self.get_hist_file_list(component=self.atm_comp)
            ds = xr.open_dataset(files[0])
            # Set grid name
            if self.grid is None and 'ne' in ds.attrs : 
               if ds.ne>0:
                  self.grid = 'ne'+str(ds.ne)
                  self.grid = self.grid+'np4'
                  for tnpg in ['pg1','pg2','pg3','pg4'] :
                     if tnpg in self.name : self.grid = self.grid.replace('np4',tnpg)
               else:
                  if 'conusx4v1'    in self.name: self.grid = 'conusx4v1'
                  if 'conusx4v1pg2' in self.name: self.grid = 'conusx4v1pg2'
            # print('\n'+self)
            # print('\n'+'case init - grid: '+self.grid)
            # print('\n'+files[0])
            # print('\n'+ds.attrs)
            # print()
            # exit()
         else:
            err_msg = '\nERROR: Case class init: history path does not exist! '
            err_msg += f'\n  case: {self.name}'
            err_msg += f'\n  path: {self.hist_path}'
            raise ValueError(err_msg)
      #----------------------------------------------------
      # check for history files
      if populate_files:
         component = 'eam'
         data_dir = self.data_dir
         data_sub = self.data_sub
         name = self.name
         for hist in hist_type_list :
            if name=='earlyscience.FC5AV1C-H01A.ne120.sp1_64x1_1000m' :
               if host=='cori':
                  tmp_name = 'early_science'
                  if hist=='h0': setattr(self, hist, f'{data_dir}{tmp_name}/monthly_hist/'   +f'{name}.{component}.{hist}.*' )
                  if hist=='h1': setattr(self, hist, f'{data_dir}{tmp_name}/hourly_2d_hist/' +f'{name}.{component}.{hist}.*' )
                  if hist=='h2': setattr(self, hist, f'{data_dir}{tmp_name}/3hourly_3d_hist/'+f'{name}.{component}.{hist}.*' )
                  if hist=='h3': setattr(self, hist, f'{data_dir}{tmp_name}/hourly_crm_hist/'+f'{name}.{component}.{hist}.*' )
               if host=='olcf': 
                  data_sub = run
                  tmp_name = 'earlyscience.FC5AV1C-H01A.ne120.sp1_64x1_1000m'
                  if hist=='h0': setattr(self, hist, f'{data_dir}{tmp_name}*/{data_sub}/{name}*.{component}.{hist}.*' )
                  if hist=='h1': setattr(self, hist, f'{data_dir}{tmp_name}*/{data_sub}/{name}*.{component}.{hist}.*' )
                  if hist=='h2': setattr(self, hist, f'{data_dir}{tmp_name}*/{data_sub}/{name}*.{component}.{hist}.*' )
                  if hist=='h3': setattr(self, hist, f'{data_dir}{tmp_name}*/{data_sub}/{name}*.{component}.{hist}.*' )
            elif name=='earlyscience.FC5AV1C-L.ne30.sp1_64x1_1000m.20190415' :
               if hist=='h0'  : setattr(self, hist, f'{data_dir}{name}/'                +f'{name}.{component}.{hist}.*' )
               if hist=='h1'  : setattr(self, hist, f'{data_dir}{name}/hourly_2d_hist/' +f'{name}.{component}.{hist}.*' )
               if hist=='h2'  : setattr(self, hist, f'{data_dir}{name}/3hourly_3d_hist/'+f'{name}.{component}.{hist}.*' )
               if hist=='h3'  : setattr(self, hist, f'{data_dir}{name}/hourly_crm_hist/'+f'{name}.{component}.{hist}.*' )
               if hist=='hgb' : setattr(self, hist, f'/project/projectdirs/m3312/whannah/{name}/atm/{name}.{component}.{hist}.*' )
               if hist=='hgb0': setattr(self, hist, f'/project/projectdirs/m3312/whannah/{name}/atm/{name}.{component}.{hist}.*' )
            else:
               setattr(self, hist, f'{data_dir}/{name}/{data_sub}/{name}.{component}.{hist}.*' )

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
   def get_hist_file_list(self,htype=None,component='eam'):
      """ 
      Retreive list of data files. Default looks for all atmos files. 
      """
      # if no htype provided just look for any htype
      if htype is None : htype = 'h'
      # Using a wildcard here helps with the early science case
      msg,err = sp.Popen('ls '+self.hist_path+f'/*.{component}.{htype}*' \
                         ,shell=True,universal_newlines=True \
                         ,stdin=sp.PIPE,stdout=sp.PIPE,stderr=sp.PIPE).communicate()
      files = msg.split('\n')
      # Remove the last item if it's empty
      if files[-1]=='' : files.pop()
      # Add this section for "multi-instance" ensemble cases
      if files==[]:
         msg,err = sp.Popen('ls '+self.hist_path+f'/*.{component}_*.{htype}*' \
                            ,shell=True,universal_newlines=True \
                            ,stdin=sp.PIPE,stdout=sp.PIPE,stderr=sp.PIPE).communicate()
         files = msg.split('\n')
         # Remove the last item if it's empty
         if files[-1]=='' : files.pop()
      # Add this section for CLM spin up case that only have clm history files
      if files==[]: 
         msg,err = sp.Popen('ls '+self.hist_path+f'/*.clm2.{htype}*' \
                            ,shell=True,universal_newlines=True \
                            ,stdin=sp.PIPE,stdout=sp.PIPE,stderr=sp.PIPE).communicate()
         files = msg.split('\n')
         # Remove the last item if it's empty
         if files[-1]=='' : files.pop()
      if files==[]: raise ValueError('No history files found!')
      return files
   #----------------------------------------------------------------------------
   def get_hist_freq(self):
      """
      Determine temporal frequency of history files
      """
      # ?
   #----------------------------------------------------------------------------
   def get_hist_vars(self):
      """ 
      Retreive list of history variables
      """
      # ?
   #----------------------------------------------------------------------------
   def get_scrip(self,scrip_dir=None):
      """ 
      Return scrip file as xarray dataset 
      """
      scrip_file_name = f'{self.grid}_scrip.nc'
      if scrip_dir==None:
         
         scrip_file_path1 = f'{home}/E3SM/data_grid/{scrip_file_name}'
         scrip_file_path2 = f'{home}/Research/E3SM/data_grid/{scrip_file_name}'
      else:
         scrip_file_path1 = f'{scrip_dir}/{scrip_file_name}'
         scrip_file_path2 = None
      file1_missing, file2_missing = False, False
      if scrip_file_path1 is not None:
         if not os.path.isfile(scrip_file_path1) : file1_missing = True
      if scrip_file_path2 is not None:
         if not os.path.isfile(scrip_file_path2) : file2_missing = True
      if file1_missing and file2_missing: 
         print('\nERROR: get_scrip(): scrip file does not exist! ')
         print(f'case: {self.name}')
         print(f'grid: {self.grid}')
         print(f'file attempt 1: {scrip_file_path1}')
         print(f'file attempt 2: {scrip_file_path2}')
         print()
         raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), scrip_file_name)
      if scrip_file_path1 is not None:
         if os.path.isfile(scrip_file_path1): scrip_file_path = scrip_file_path1
      if scrip_file_path2 is not None:
         if os.path.isfile(scrip_file_path2): scrip_file_path = scrip_file_path2
      scripfile = xr.open_dataset(scrip_file_path)
      # print(scrip_file_path)
      return scripfile
   #----------------------------------------------------------------------------
   def set_coord_names(self,var):
      # dynamics grid data
      if 'DYN_' in var or var in ['VOR','DIV'] : 
         physgrid_list = ['pg2','pg3','pg4']
         if any(g in self.name for g in physgrid_list) : 
            # update grid name
            for g in physgrid_list:
               if g in self.name : self.grid = self.grid.replace(g,'np4')
            # update area and ncol names
            self.ncol_name,self.area_name = 'ncol_d','area_d'
      # regional subset
      if '_to_' in var:
         self.ncol_name = var.replace(var.split('_')[0],'ncol')
      return
   #----------------------------------------------------------------------------
   def get_pressure(self,ds, use_interface=False):
      """
      Calculate 3D pressure field from hybrid coordinates
      """
      if use_interface:
        a = ds['hyai']
        b = ds['hybi']
      else:
        a = ds['hyam']
        b = ds['hybm']
      p0 = ds['P0']
      ps = ds['PS']
      pressure = a * p0 + b * ps
      pressure = pressure * 1e-2
      pressure.attrs['units'] = 'hPa'
      pressure.attrs['long_name'] = 'Pressure'
      return pressure
   #----------------------------------------------------------------------------
   def get_mask(self,ds=None,htype=None,lat_name='lat',lon_name='lon'):
      """
      Create mask with only ncol dimension from first available file 
      """
      regional_subset = False
      if ds is None :
         hist_file = self.get_hist_file_list(htype=htype,component=self.atm_comp)[0]
         ds = xr.open_dataset(hist_file)
         ### deal with case when ncol_name not in file?
         # if self.ncol_name not in ds.dims : 
         #    self.get_hist_file_list(htype='h1')[0]
         ### deal with regional output
         for d in ds.dims:
            if 'ncol' in d:
               if '_to_' in d:
                  regional_subset,regional_ncol_name = True,d
      # Separate logic is needed when working with different scenarios 
      # raw vs. remapped data vs. regional subsets
      if regional_subset :
         ncol = ds[regional_ncol_name]
         tmp_data = np.ones([len(ncol)],dtype=bool)
         tmp_coords = [(regional_ncol_name, ncol)]
         mask = xr.DataArray( tmp_data, coords=tmp_coords, dims=self.ncol_name )
      elif 'ncol' in ds.dims:
         if self.ncol_name=='ncol_d' : lat_name,lon_name = 'lat_d','lon_d'
         ncol = ds[self.ncol_name]
         tmp_data = np.ones([len(ncol)],dtype=bool)
         tmp_coords = [(self.ncol_name, ncol)]
         mask = xr.DataArray( tmp_data, coords=tmp_coords, dims=self.ncol_name )
      elif 'ni' in ds.dims: # for CICE data
         lat_name,lon_name = 'TLAT','TLON'
         ncol = ds[self.ncol_name]
         tmp_data = np.ones([len(ncol)],dtype=bool)
         tmp_coords = [(self.ncol_name, ncol)]
         mask = xr.DataArray( tmp_data, coords=tmp_coords, dims=self.ncol_name )
      elif lat_name in ds.dims and lon_name in ds.dims:
         tmp_data = np.ones([len(ds[lat_name]),len(ds[lon_name])],dtype=bool)
         tmp_coords = {lat_name:ds[lat_name],lon_name:ds[lon_name]}
         mask = xr.DataArray( tmp_data, coords=tmp_coords, dims=(lat_name,lon_name) )
      else:
         raise ValueError('Something is wrong with the coordinates of the input data! ')

      if hasattr(self,'lat1'): mask = mask & (ds[lat_name]>=self.lat1)
      if hasattr(self,'lat2'): mask = mask & (ds[lat_name]<=self.lat2)
      if hasattr(self,'lon1'): mask = mask & (ds[lon_name]>=self.lon1)
      if hasattr(self,'lon2'): mask = mask & (ds[lon_name]<=self.lon2)

      # Add equatorially symmetric bands to the mask
      if self.mirror_equator and hasattr(self,'lat1') and hasattr(self,'lat2') :
         lat1_mirror = self.lat1*-1
         lat2_mirror = self.lat2*-1
         mirror_mask = xr.DataArray( np.ones([len(ncol)],dtype=bool), \
                           coords=[(self.ncol_name, ncol)], dims=self.ncol_name )
         mirror_mask = mirror_mask & (ds[lat_name]<=lat1_mirror) & (ds[lat_name]>=lat2_mirror)
         mask = mask | mirror_mask

      return mask
   #----------------------------------------------------------------------------
   def get_var_name(self,ds,var,htype):
      """
      Get initial variable to load (tvar) for requested derived quantity (var)
      """
      tvar = var
      if var=='P-E': 
         if htype=='h0': 
            tvar = 'PRECC'
         else:
            tvar = 'PRECT'
      if var=='MSE'        : tvar = 'T'
      if var=='QT'         : tvar = 'Q'
      if var=='CWV'        : tvar = 'TMQ'
      if var=='RH'         : tvar = 'Q'
      if var=='THETA'      : tvar = 'T'
      if var=='NET_TOA_RAD': tvar = 'FSNT'
      if var=='WSPD'       : tvar = 'U'
      if var=='WSPD850'    : tvar = 'U850'
      if var=='WSPD200'    : tvar = 'U200'
      if var=='WSPD_BOT'   : tvar = 'UBOT'
      if var=='PRECT'and 'PRECT' not in ds : tvar = 'PRECC'
      if var=='area' and 'area'  not in ds : tvar = 'area_p'

      ### SPCAM stuff
      if 'SPCAM' in self.name:
         # if var=='SPTLS' : tvar = 'MMF_TLS'
         # if var=='SPQTLS': tvar = 'MMF_QTLS'
         # if var=='SPDT'  : tvar = 'MMF_DT'
         # if var=='SPDQ'  : tvar = 'MMF_DQ'
         if var=='MMF_TLS' : tvar = 'SPTLS'
         if var=='MMF_QTLS': tvar = 'SPQTLS'
         if var=='MMF_DT'  : tvar = 'SPDT'
         if var=='MMF_DQ'  : tvar = 'SPDQ'
      return tvar
   #----------------------------------------------------------------------------
   def get_var_data(self,ds,var,htype):
      """
      Load requested variable, including derived quantities
      """
      # change variable name for derived variables
      tvar = self.get_var_name(ds,var,htype)
      if tvar not in ds :
         print(ds)
         if 'file_path' in locals(): print('file_path: '+file_path)
         raise ValueError('\nERROR: load_data: '+var+' ('+tvar+') not found in dataset\n')
      #----------------------------------------------------
      # Load the data and subset
      #----------------------------------------------------
      data = ds[tvar]
      if 'ncol' in data.dims or 'ncol_d' in data.dims or 'ni' in data.dims :
         data = data.where( self.get_mask(ds), drop=True)
      #----------------------------------------------------
      # derived variable cases   
      #----------------------------------------------------
      # Total Precipitation
      if var=='PRECT' and 'PRECT' not in ds : 
         data = data + ds['PRECL'].where( self.get_mask(ds), drop=True)
         data['long_name'] = 'Total Precipitation'
      #sfc precipitation minus evaporation
      if var=='P-E' : 
         if htype=='h0': 
            data = (data+ds['PRECL'])*86400.*1e3 - ds['LHFLX'].where( self.get_mask(ds), drop=True) * 86400./hc.Lv
         else:
            data = data*86400.*1e3 - ds['LHFLX'].where( self.get_mask(ds), drop=True) * 86400./hc.Lv
            # num_t = len(data['time'])
            # data = data.isel(time=slice(0,num_t-1))*86400.*1e3 - ds['LHFLX'].isel(time=slice(1,num_t)).where( self.get_mask(ds), drop=True) * 86400./hc.Lv
         data['long_name'] = 'Sfc Precip - Evap'
         data['units']     = 'mm/day'
      # Moist Static Energy
      if var=='MSE' : 
         data = data + ds['Z3'].where( self.get_mask(ds), drop=True) * hc.g /hc.cpd \
               + ds[ 'Q'].where( self.get_mask(ds), drop=True) * hc.Lv/hc.cpd 
         data['long_name'] = 'Moist Static Energy'
         data['units']     = 'K'
      # Total Water
      if var=='QT' : 
         data = data + ds['CLDLIQ'].where( self.get_mask(ds), drop=True) \
               + ds['CLDICE'].where( self.get_mask(ds), drop=True) 
         data['long_name'] = 'Total Water'
         data['units']     = 'kg/kg'
      # Relative Humidity
      if var=='RH' : 
         T = ds['T'].where( self.get_mask(ds), drop=True)
         P = self.get_pressure(ds).transpose('time','lev','ncol').where( self.get_mask(ds), drop=True)
         data_sat = hc.calc_Qsat(T,P)
         data = data / data_sat * 1e2
         data['long_name'] = 'Relative Humidity'
         data['units']     = '%'
      # Potential Temperature
      if var=='THETA' : 
         P = self.get_pressure(ds).transpose('time','lev','ncol').where( self.get_mask(ds), drop=True)
         data = data * ( 1000. / P )**(hc.Rd/hc.cpd)
         data['long_name'] = 'Potential Temperature'
         data['units']     = 'K'
      # TOA net radiative heating
      if var=='NET_TOA_RAD' : 
         # print()
         # print(data)
         # print()
         # print(ds['FLNT'])
         # print()
         data = data - ds['FLNT'].where( self.get_mask(ds), drop=True) 
         data['long_name'] = 'TOA net radiative heating'
         # data['units']     = ''
      if var=='WSPD' : 
         U = data
         V = ds['V'].where( self.get_mask(ds), drop=True)
         data = np.sqrt(U**2+V**2)
         data['long_name'] = 'Wind Speed'
      if var=='WSPD850' : 
         U = data
         V = ds['V850'].where( self.get_mask(ds), drop=True)
         data = np.sqrt(U**2+V**2)
         data['long_name'] = '850mb Wind Speed'
      if var=='WSPD200' : 
         U = data
         V = ds['V200'].where( self.get_mask(ds), drop=True)
         data = np.sqrt(U**2+V**2)
         data['long_name'] = '200mb Wind Speed'
      if var=='WSPD_BOT' : 
         U = data
         V = ds['VBOT'].where( self.get_mask(ds), drop=True)
         data = np.sqrt(U**2+V**2)
         data['long_name'] = 'Wind Speed at lowest model level'
      # Return the data array
      return data
   #----------------------------------------------------------------------------
   def adjust_units(self,var_name=None,data=None):
      """
      Adjust units of input data according to variable name
      """
      if var_name=='PRECT':   data['units'],data.values = 'mm/day',data*86400.*1e3
      if var_name=='PRECC':   data['units'],data.values = 'mm/day',data*86400.*1e3
      if var_name=='PRECL':   data['units'],data.values = 'mm/day',data*86400.*1e3
      if var_name=='MMF_DT':  data['units'],data.values = 'K/day', data*86400.
      if var_name=='MMF_TLS': data['units'],data.values = 'K/day', data*86400.
      if var_name=='MMF_DQ':  data['units'],data.values = 'g/kg/day',data*86400.*1e3
      if var_name=='MMF_QTLS':data['units'],data.values = 'g/kg/day',data*86400.*1e3
      if var_name=='SPDT':    data['units'],data.values = 'K/day', data*86400.
      if var_name=='SPTLS':   data['units'],data.values = 'K/day', data*86400.
      if var_name=='SPDQ':    data['units'],data.values = 'g/kg/day',data*86400.*1e3
      if var_name=='SPQTLS':  data['units'],data.values = 'g/kg/day',data*86400.*1e3
      if var_name=='Q':       data['units'],data.values = 'g/kg',  data*1e3
      if var_name=='CRM_QV':  data['units'],data.values = 'g/kg',  data*1e3
      if var_name=='Q850':    data['units'],data.values = 'g/kg',  data*1e3
      if var_name=='DYN_Q':   data['units'],data.values = 'g/kg',  data*1e3
      if var_name=='CLDLIQ':  data['units'],data.values = 'g/kg',  data*1e3
      if var_name=='CLDICE':  data['units'],data.values = 'g/kg',  data*1e3
      if var_name=='MMF_QI':  data['units'],data.values = 'g/kg',  data*1e3
      if var_name=='TS':      data['units'],data.values = 'C',     data - 273
      if var_name=='QRS':     data['units'],data.values = 'K/day', data*86400.
      if var_name=='QRL':     data['units'],data.values = 'K/day', data*86400.

      # if var_name=='MMF_CVT_T':     data['units'],data.values = '', data.
      if var_name=='MMF_VT_Q':     data['units'],data.values = '', data*1e3
      if var_name=='MMF_VT_TEND_T':data['units'],data.values = '', data*86400.
      if var_name=='MMF_VT_TEND_Q':data['units'],data.values = '', data*86400.*1e3
      if var_name=='MMF_VT_TLS':   data['units'],data.values = '', data*86400.
      if var_name=='MMF_VT_QLS':   data['units'],data.values = '', data*86400.*1e3
      

      if var_name in ['MMF_DU','MMF_DV','ZMMTU','ZMMTV','uten_Cu','vten_Cu']:
         data['units'],data.values = 'm/s/day', data*86400.
      return
   #----------------------------------------------------------------------------
   def load_data(self,var,htype='h0',ps_htype=None,lev=np.array([0]),  \
                 num_years=-1, years=[], months=[],      \
                 num_files=0, first_file=0,              \
                 component='eam',use_remap=False,remap_str='.remap',
                 extrap_flag=False) :
      """ """
      if not isinstance(num_files,int):  raise ValueError('hapy_E3SM: load_data(): num_files must be an integer!')
      if not isinstance(first_file,int): raise ValueError('hapy_E3SM: load_data(): first_file must be an integer!')
      if not isinstance(num_years,int):  raise ValueError('hapy_E3SM: load_data(): num_years must be an integer!')
      if ps_htype is None: ps_htype = htype
      #----------------------------------------------------
      # set the file path
      #----------------------------------------------------
      if htype in hist_type_list and component=='eam' and not use_remap :
         file_path = getattr(self, htype)
      else:
         # data_sub = remap_str if use_remap else self.data_sub
         data_sub = self.data_sub
         file_path = self.data_dir+'/'+self.name+'/'+data_sub+'/' +self.name+'.'+component+'.'+htype+'.*'
      #----------------------------------------------------
      # Subset files by year - this may be a bad idea for files that straddle Dec 31...
      #----------------------------------------------------
      if num_years!=-1 and num_years<10 : 
         file_path = file_path.replace('*','000[1-'+str(num_years)+']*')
      #----------------------------------------------------
      # open the files
      #----------------------------------------------------
      if len(glob.glob(file_path))>0 : 
         file_list = sorted(glob.glob(file_path))
         # Remove "remapped" files unless use_remap is true
         file_list_tmp = file_list.copy()
         for file_name in file_list : 
            if use_remap and remap_str not in file_name: file_list_tmp.remove(file_name)
            if not use_remap and remap_str in file_name: file_list_tmp.remove(file_name)
         file_list = file_list_tmp
         # limit to single file if num_file>0
         if num_files==0 and first_file>0: file_list = file_list[first_file:]
         if num_files>0 : file_list = file_list[first_file:first_file+num_files]     # use initial files
         if num_files<0 : file_list = file_list[num_files:]      # use latest files
         # ds = xr.open_mfdataset( file_list, combine='by_coords', concat_dim='time', data_vars='minimal', coords='minimal' )
         # ds = xr.open_mfdataset( file_list, combine='by_coords', concat_dim='time', compat='override', data_vars=var )
         def drop_rad_band(ds): 
            if 'lwband' in ds.variables: ds = ds.drop(['lwband'])
            if 'swband' in ds.variables: ds = ds.drop(['swband'])
            # hack for time coordinate in CICE data
            if 'TLAT' in ds.coords: ds['time'].values = ds['time_bounds'][:,0].values
            return ds

         ### check for incomplete files
         if len(file_list)>1:
            # sizes = []
            # for file in file_list: 
            #    sizes.append( os.stat(file).st_size )
            #    # msg,err = sp.Popen(['du',file], universal_newlines=True, stdout=sp.PIPE, stderr=sp.PIPE).communicate()
            # if sizes[-1] < sizes(0) - np.stddev(np.array(sizes[0:-2])): file_list.pop()

            if os.stat(file_list[-1]).st_size < os.stat(file_list[0]).st_size: file_list.pop()

         ### truncate file list for debugging or working with partial files
         # file_list = sorted(file_list)[:5]
         # print('WARNING - truncating file list - line 600 in hapy_E3SM.py - WARNING')

         ds = xr.open_mfdataset( file_list, combine='by_coords', concat_dim='time', preprocess=drop_rad_band )
      else :
         print('\nERROR: load_data(): No history files found!')
         print('var  : '+var)
         print('path : '+file_path+'\n')
         raise FileNotFoundError('load_data(): No history files found!')
      #----------------------------------------------------
      # Correct timestamp - shift by mean delta t
      #----------------------------------------------------
      # if htype=='h0':
      #    dtime = ds['time'].diff(dim='time').values.astype('timedelta64[h]').mean()
      #    # print(ds['time'].values)
      #    ds['time'] = ds.time.get_index('time') \
      #                -datetime.timedelta(hours=dtime.astype(np.int)*0.5)
      #    # print(ds['time'].values)
      #----------------------------------------------------
      # Load and subset the data - including derived quantities
      #----------------------------------------------------
      data = self.get_var_data(ds,var,htype)
      #----------------------------------------------------
      # Adjust units - if not loading with precalculated data
      #----------------------------------------------------
      if 'hgb' not in htype: self.adjust_units(var_name=var,data=data)
      #----------------------------------------------------
      # Get dataset for surface pressure if interpolating
      #----------------------------------------------------
      if ps_htype!=htype and np.all(lev>0) and 'lev' in data.coords :
         ps_file_list = [ f.replace(f'.{htype}.',f'.{ps_htype}.') for f in file_list]
         ds_ps = xr.open_mfdataset( ps_file_list, combine='by_coords', concat_dim='time', preprocess=drop_rad_band )

         # Extract time frequency in hours
         dt = (ds['time'][1] - ds['time'][0]).values
         dt_hr_freq = f'{dt/86400e9*24.}H'

         # resample PS to match data
         ds_ps = ds_ps.resample(time=dt_hr_freq).mean(dim='time')
      else:
         ds_ps = None
      #----------------------------------------------------
      # vertical subset or interpolation
      #----------------------------------------------------
      # interpolate to pressure if lev is positive
      if np.all( lev > 0 ) and 'lev' in data.coords :
         data = interpolate_to_pressure(ds,data_mlev=data,lev=lev,ds_ps=ds_ps,ps_var='PS'
                                       ,interp_type=2,extrap_flag=extrap_flag)
      
      # if lev<0 use model levels without interpolation
      if np.all(lev<0) and 'lev' in data.coords :
         if  'lev' in data.coords: lev_coord = 'lev'
         if 'ilev' in data.coords: lev_coord = 'ilev'
         data = data.isel({lev_coord:np.absolute(lev)})
      #----------------------------------------------------
      # ignore time dimension for grid variables
      #----------------------------------------------------
      if var in ['lat','lon','lat_d','lon_d','area','area_p','area_d'] : 
         if 'time' in data.dims: data = data.isel(time=0)
      #----------------------------------------------------
      # Subset in time
      #----------------------------------------------------
      if 'time' in data.dims and num_files==0 :
         time_mask = xr.DataArray( np.ones([len(data.time)],dtype=bool), coords=[('time', data.time)], dims='time' )
         if num_years>=0 : time_mask = time_mask & ( ( data['time.year']-data['time.year'][0]) <=num_years)
         if years !=[]   : time_mask = time_mask & [ y in years  for y in (data['time.year']-data['time.year'][0]).values ]
         if len(months)>0: time_mask = time_mask & [ m in months for m in data['time.month'].values ]
         data = data.sel( time=data.time.where(time_mask,drop=True) )
      #----------------------------------------------------
      #----------------------------------------------------
      return data
   #----------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------
# Misc Routines
#---------------------------------------------------------------------------------------------------
def get_GLL_area_bins(case_name):
   if "ne4"   in case_name : area_bins = [ 0.0000, 0.006000, 0.015000, 0.030000 ]
   if "ne30"  in case_name : area_bins = [ 0.0000, 0.000100, 0.000300, 0.000600 ]
   if "ne120" in case_name : area_bins = [ 0.0000, 0.000007, 0.000015, 0.000050 ]
   return area_bins
#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------
def calc_PSL(pmid, phis, ps, T) :
   """ Calculate sea level pressure using method from E3SM (components/eam/src/physics/cam/cpslec.F90) 
   Method: CCM2 hybrid coord version using ECMWF formulation Algorithm: 
   See section 3.1.b in NCAR NT-396 \"Vertical Interpolation and Truncation of Model-Coordinate Data\" """
   gravit = SHR_CONST_G
   rair = SHR_CONST_RDAIR
   xlapse = 6.5e-3   # Temperature lapse rate (K/m)
   time = len( T['time'] )
   ncol = len( T['ncol'] )
   pver = len( T['lev'] )

   alpha = rair*xlapse/gravit

   psl = xr.full_like(ps,np.nan)

   for t in range(0,time) :
      for i in range(0,ncol) :
         if ( np.abs(phis.isel(ncol=i)/gravit) < 1.e-4 ) :
            psl[t,i] = ps[t,i]
         else :
            Tstar = T.isel(time=t,ncol=i,lev=pver-1) * (1.0+alpha*( ps[t,i] /pmid.isel(time=t,ncol=i,lev=pver-1)-1.0) )  # pg 7 eq 5

            TT0 = Tstar + xlapse*phis.isel(ncol=i)/gravit                  # pg 8 eq 13

            if ( Tstar<=290.5 and TT0>290.5 )  :                           # pg 8 eq 14.1
               alph = rair / phis.isel(ncol=i) * (290.5-Tstar)  
            elif (Tstar>290.5 and TT0>290.5) :                            # pg 8 eq 14.2
               alph = 0.0
               Tstar= 0.5 * (290.5 + Tstar)  
            else :
               alph = alpha  
               if (Tstar<255.) : Tstar = 0.5 * (255. + Tstar)              # pg 8 eq 14.3

            beta = phis.isel(ncol=i)/(rair*Tstar)
            psl[t,i] = ps[t,i] * np.exp( beta*( 1.0 - alph*beta/2.0 + (alph*beta)**2 / 3.0 ) )

   return psl

@numba.njit()
def calc_PSL_numba( ntime, nlev, ncol, pmid, phis, ps, T, psl) :
   gravit = SHR_CONST_G
   rair = SHR_CONST_RDAIR
   xlapse = 6.5e-3   # Temperature lapse rate (K/m)
   alpha = rair*xlapse/gravit
   for t in range(0,ntime) :
      for i in range(0,ncol) :
         if ( np.abs(phis[i]/gravit) < 1.e-4 ) :
            psl[t,i] = ps[t,i]
         else :
            Tstar = T[t,i] * ( 1.0+alpha*( ps[t,i]/pmid[t,i]-1.0 ) )       # pg 7 eq 5
            TT0 = Tstar + xlapse*phis[i]/gravit                  # pg 8 eq 13
            if ( Tstar<=290.5 and TT0>290.5 )  :                           # pg 8 eq 14.1
               alph = rair / phis[i] * (290.5-Tstar)  
            elif (Tstar>290.5 and TT0>290.5) :                            # pg 8 eq 14.2
               alph = 0.0
               Tstar= 0.5 * (290.5 + Tstar)  
            else :
               alph = alpha  
               if (Tstar<255.) : Tstar = 0.5 * (255. + Tstar)              # pg 8 eq 14.3
            beta = phis[i]/(rair*Tstar)
            psl[t,i] = ps[t,i] * np.exp( beta*( 1.0 - alph*beta/2.0 + (alph*beta)**2 / 3.0 ) )

#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------
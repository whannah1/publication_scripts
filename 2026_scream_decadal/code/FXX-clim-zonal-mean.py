import os, glob, ngl, subprocess as sp, numpy as np, xarray as xr, copy, string, sys, cmocean
# import hapy_common as hc, hapy_E3SM   as he, hapy_setres as hs
from hapy_setres import set_subtitles
import hapy
#-------------------------------------------------------------------------------
# case_dir,case_sub = [],[]
# case,name = [],[]
# htype_list = []
# comp_list = []
# clr,dsh,mrk = [],[],[]
# def add_case(case_in,n=None,comp=None,htype=None,p=None,s=None,d=0,c='black',m=0):
#    global name,case,case_dir,case_sub,clr,dsh,mrk
#    # if comp  is None: raise ValueError(f'ERROR - add_case: comp argument cannot be None')
#    # if htype is None: raise ValueError(f'ERROR - add_case: htype argument cannot be None')
#    # if p     is None: raise ValueError(f'ERROR - add_case: p argument cannot be None')
#    # if s     is None: raise ValueError(f'ERROR - add_case: s argument cannot be None')
#    case.append(case_in); name.append(n); 
#    comp_list.append(comp); htype_list.append(htype)
#    case_dir.append(p); case_sub.append(s)
#    dsh.append(d) ; clr.append(c) ; mrk.append(m)
#-------------------------------------------------------------------------------
case,opts_list = [],[]
def add_case(case_in,**kwargs):
   case.append(case_in)
   case_opts = {}
   for k, val in kwargs.items(): case_opts[k] = val
   opts_list.append(case_opts)
#-------------------------------------------------------------------------------
var_list = []
eam_var_list = []
eamxx_var_list = []
lev_list = []
unit_list = []
def add_var(var_name,eam_var=None,eamxx_var=None,var_str=None,lev=None,unit=None): 
   var_list.append(var_name)
   eam_var_list.append(eam_var)
   eamxx_var_list.append(eamxx_var)
   lev_list.append(lev)
   unit_list.append(unit)
#-------------------------------------------------------------------------------
tmp_path_ne1024 = '/global/cfs/cdirs/e3smdata/simulations/scream-decadal'
tmp_path_ne256  = '/global/cfs/cdirs/e3smdata/simulations/scream-decadal-ne256'
tmp_path_hst_v3 = '/global/cfs/cdirs/m3312/whannah/e3smv3_amip'
tmp_path_qbo_bm = '/pscratch/sd/w/whannah/e3sm_scratch/pm-cpu/'

# add_case('ERA5',n='ERA5')
# add_case('v3.LR.amip_0101.QBObenchmark.20241008',                                      n='EAMv3 AMIP',          comp='eam',  htype='eam.h0', p=tmp_path_qbo_bm, s='data_remap_90x180')
# add_case('v3.LR.amip_0201',                                                            n='v3.LR.amip_0201',     comp='eam',  htype='eam.h2', p=tmp_path_hst_v3, s='data_remap_90x180')
# add_case('decadal-production-run6',                                                      n='SCREAM ne1024',       comp='eamxx',htype='_ANN_199501_200412_climo.nc',p=tmp_path_ne1024,s='climo_180x360')
add_case('ne256pg2_ne256pg2.F20TR-SCREAMv1.May-12.with.rain.frac.n0128',                 n='SCREAM 13-km control',comp='eamxx',htype='_ANN_199501_200412_climo.nc',p=tmp_path_ne256, s='climo_180x360')
add_case('ne256pg2_ne256pg2.F20TR-SCREAMv1.July-1.spanc800.2xauto.acc150.n0032.test2.1', n='SCREAM 13-km tuned',  comp='eamxx',htype='_ANN_199501_200412_climo.nc',p=tmp_path_ne256, s='climo_180x360')

'''
ERA5 hourly data:
/global/cfs/projectdirs/m3522/cmip6/ERA5/e5.oper.an.pl/201111/e5.oper.an.pl.128_130_t.ll025sc.2011111900_2011111923.nc
/global/cfs/cdirs/m3312/whannah/ERA5/daily
'''
#-------------------------------------------------------------------------------

add_var('Q',eam_var='Q', eamxx_var='qv',     unit='kg/kg')
# add_var('U',eam_var='U', eamxx_var='U',     unit='m/s')
# add_var('T',eam_var='T', eamxx_var='T_mid', unit='K')

num_plot_col = 1

#-------------------------------------------------------------------------------
# ERA5 levels
# lev = np.array([   1.,    2.,    3.,    5.,    7.,   10.,   20.,   30.,   50.,   70.,
#                  100.,  125.,  150.,  175.,  200.,  225.,  250.,  300.,  350.,  400.,
#                  450.,  500.,  550.,  600.,  650.,  700.,  750.,  775.,  800.,  825.,
#                  850.,  875.,  900.,  925.,  950.,  975., 1000.])
lev = np.array([ 50.,   70.,   100.,  125.,  150.,  175.,  200.,  225.,  250.,  300.,  350.,  400.,
                 450.,  500.,  550.,  600.,  650.,  700.,  750.,  775.,  800.,  825.,
                 850.,  875.,  900.,  925.,  950.,  975., 1000.])
lev = lev*1e2
#-------------------------------------------------------------------------------

fig_type,fig_file = 'png',f'figs/FXX-clim-zonal-mean'
# fig_file_diff = f'{fig_file}.diff'

lat1,lat2 = -15,15

# yr1,yr2 = 1975,2020
yr1,yr2 = 1995,2004
# yr1,yr2 = 1995,1995

recalculate = True

var_x_case  = False

# num_plot_col = 2

use_common_label_bar = True

#---------------------------------------------------------------------------------------------------
def run_cmd(cmd,verbose=True,indent='  '):
   cmd_str = indent + hapy.tclr.GREEN + cmd + hapy.tclr.END
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
num_case,num_var = len(case),len(var_list)

subtitle_font_height = 0.01

wkres = ngl.Resources()
npix = 2048; wkres.wkWidth,wkres.wkHeight=npix,npix
wks = ngl.open_wks(fig_type,fig_file,wkres)

plot = [None]*(num_var*num_case)

res = ngl.Resources()
res.nglDraw                      = False
res.nglFrame                     = False
res.tmXTOn                       = False
res.tmXBMajorOutwardLengthF      = 0.
res.tmXBMinorOutwardLengthF      = 0.
res.tmYLMajorOutwardLengthF      = 0.
res.tmYLMinorOutwardLengthF      = 0.
res.tmYLLabelFontHeightF         = 0.015
res.tmXBLabelFontHeightF         = 0.015
res.tiXAxisFontHeightF           = 0.015
res.tiYAxisFontHeightF           = 0.015
res.tmXBMinorOn                  = False
res.tmYLMinorOn                  = False
res.cnFillOn                     = False
res.cnLinesOn                    = True
res.cnLineLabelsOn               = False
res.cnInfoLabelOn                = False
res.cnLineThicknessF             = 3.
res.cnFillPalette = "ncl_default"
res.cnFillOn                     = True
res.cnLinesOn                    = False
res.cnLineLabelsOn               = False
res.cnInfoLabelOn                = False
res.lbOrientation                = "Horizontal"
res.lbLabelFontHeightF           = 0.008

res.vpHeightF = 0.3
res.trYReverse = True
res.tiXAxisString = 'Latitude'
res.tiYAxisString = 'Pressure [hPa]'

# res.tmYLLabelFontHeightF    = 0.008
# res.tmXBLabelFontHeightF    = 0.008
# res.lbLabelFontHeightF      = 0.01

res.tmYLLabelFontHeightF   = 0.01
res.tmXBLabelFontHeightF   = 0.01
res.lbLabelFontHeightF     = 0.015


lres = ngl.Resources()
lres.nglDraw                      = False
lres.nglFrame                     = False
lres.tmXTOn                       = False
lres.tmXBMajorOutwardLengthF      = 0.
lres.tmXBMinorOutwardLengthF      = 0.
lres.tmYLMajorOutwardLengthF      = 0.
lres.tmYLMinorOutwardLengthF      = 0.
lres.tmYLLabelFontHeightF         = 0.015
lres.tmXBLabelFontHeightF         = 0.015
lres.tiXAxisFontHeightF           = 0.015
lres.tiYAxisFontHeightF           = 0.015
lres.tmXBMinorOn                  = False
lres.tmYLMinorOn                  = False
lres.xyLineThicknessF             = 6.
lres.xyLineThicknessF = 1
lres.xyLineColor      = 'black'
lres.xyDashPattern    = 0

#---------------------------------------------------------------------------------------------------
def calculate_area(lon,lat,lon_bnds,lat_bnds):
  re = 6.37122e06  # radius of earth
  nlat,nlon = len(lat),len(lon)
  area = np.empty((nlat,nlon),np.float64)
  for j in range(nlat):
    for i in range(nlon):
      dlon = np.absolute( lon_bnds[j,1] - lon_bnds[j,0] )
      dlat = np.absolute( lat_bnds[j,1] - lat_bnds[j,0] )
      dx = re*dlon*np.pi/180.
      dy = re*dlat*np.pi/180.
      area[j,i] = dx*dy
  return area
#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------
for v in range(num_var):
   print()
   print(f'  var: {hapy.tclr.MAGENTA}{var_list[v]}{hapy.tclr.END}')
   data_list = []
   lat_list = []
   lev_list = []
   for c in range(num_case):
      case_opts = opts_list[c]
      print('    case: '+hapy.tclr.GREEN+case[c]+hapy.tclr.END)

      if recalculate:

         if case[c]=='ERA5':
            #-------------------------------------------------------------------
            # lat_name,lon_name = 'lat','lon'
            # xy_dims = (lon_name,lat_name)
            #-------------------------------------------------------------------
            # obs_root = '/global/cfs/cdirs/e3sm/diagnostics/observations/Atm/time-series/ERA5'
            # if var_list[v]=='U' : obs_var,input_file_name =  'ua', f'{obs_root}/ua_197901_201912.nc'
            # if var_list[v]=='V' : obs_var,input_file_name =  'va', f'{obs_root}/va_197901_201912.nc'
            # if var_list[v]=='T' : obs_var,input_file_name =  'ta', f'{obs_root}/ta_197901_201912.nc'
            # if var_list[v]=='Q' : obs_var,input_file_name = 'hus', f'{obs_root}/hus_197901_201912.nc'
            # if var_list[v]=='O3': obs_var,input_file_name = 'tro3',f'{obs_root}/tro3_197901_201912.nc'
            #-------------------------------------------------------------------
            obs_root = '/pscratch/sd/w/whannah/Obs/ERA5/e3sm_diags_remap'
            if var_list[v]=='U' : obs_var,input_file_name =  'ua', f'{obs_root}/ua_197901_201912.remap_180x360.nc'
            if var_list[v]=='V' : obs_var,input_file_name =  'va', f'{obs_root}/va_197901_201912.remap_180x360.nc'
            if var_list[v]=='T' : obs_var,input_file_name =  'ta', f'{obs_root}/ta_197901_201912.remap_180x360.nc'
            if var_list[v]=='Q' : obs_var,input_file_name = 'hus', f'{obs_root}/hus_197901_201912.remap_180x360.nc'
            if var_list[v]=='O3': obs_var,input_file_name = 'tro3',f'{obs_root}/tro3_197901_201912.remap_180x360.nc'
            #-------------------------------------------------------------------
            # obs_root = '/global/cfs/cdirs/e3sm/diagnostics/observations/Atm/climatology/ERA5'
            # input_file_name = f'{obs_root}/ERA5_ANN_197901_201912_climo.nc'
            #-------------------------------------------------------------------
            ds = xr.open_dataset( input_file_name )
            ds = ds.where( ds['time.year']>=yr1, drop=True)
            ds = ds.where( ds['time.year']<=yr2, drop=True)
            #-------------------------------------------------------------------
            # lat  = ds['lat']
            # area = calculate_area(ds['lon'].values,ds['lat'].values,ds['lon_bnds'].values,ds['lat_bnds'].values)
            # area = xr.DataArray( area, coords=[ds['lat'],ds['lon']] )
            data = ds[obs_var].mean(dim='time')
            # data = data.rename({'plev':'lev'})
            # data['lev'] = data['lev']/1e2
            data = data.sel(plev=lev)
            #-------------------------------------------------------------------
            # data = mask_data(ds,data)
            # area = mask_data(ds,area)
            # data_avg = ( (data*area).sum(dim=xy_dims) / area.sum(dim=xy_dims) )
            #-------------------------------------------------------------------
         else:
            case_dir = case_opts['p']
            case_sub = case_opts['s']
            htype    = case_opts['htype']
            #-------------------------------------------------------------------
            file_path = f'{case_dir}/{case[c]}/{case_sub}/*{htype}*'
            # file_path = f'{case_dir[c]}/{case[c]}/{tmp_sub}/{case[c]}.eam.{htype}.*'
            file_list = sorted(glob.glob(file_path))
            #-------------------------------------------------------------------
            # print(f'    case_dir  : {case_dir}')
            # print(f'    case_sub  : {case_sub}')
            # print(f'    htype     : {htype}')
            # print(f'    file_path : {file_path}')
            # print();print(file_list)
            # exit()
            #-------------------------------------------------------------------
            # # subset files that fall within [yr1:yr2]
            # file_list_all = file_list ; file_list = []
            # for f in range(len(file_list_all)):
            #    yr = int(file_list_all[f][-10:-10+4])
            #    if yr>=yr1 and yr<=yr2: file_list.append(file_list_all[f])
            #-------------------------------------------------------------------
            # if 'first_file' in locals(): file_list = file_list[first_file:]
            # if 'num_files'  in locals(): file_list = file_list[:num_files]
            #-------------------------------------------------------------------
            # for f in file_list: print(f)
            # exit()
            #-------------------------------------------------------------------
            ds = xr.open_mfdataset( file_list )
            if 'time' in ds.dims:
               if len(ds.time)==1: ds = ds.isel(time=0,drop=True)
               # ds = ds.where( ds['time.year']>=yr1, drop=True)
               # ds = ds.where( ds['time.year']<=yr2, drop=True)
            #-------------------------------------------------------------------
            tvar = None
            if case_opts['comp']=='eam'  : tvar = eam_var_list[v]
            if case_opts['comp']=='eamxx': tvar = eamxx_var_list[v]
            data = ds[tvar]
            # area = ds['area']
            # lat  = ds['lat']
            #-------------------------------------------------------------------
            # data = mask_data(ds,data)
            # area = mask_data(ds,area)
            # data_avg = ( (data*area).sum(dim='ncol') / area.sum(dim='ncol') )
            #-------------------------------------------------------------------
            # adjust units
            # if var_list[v]=='utendepfd': data = data*86400. # m/s/s => m/s/day
            #-------------------------------------------------------------------
            data = hapy.vinth2p_simple(data, ds['hyam'], ds['hybm'], lev, ds['ps'], \
                                       interp_type='linear',extrapolate=False)
            if len(data.plev)==1: data = data.isel(plev=0)
            #-------------------------------------------------------------------
            # # adjust time to represent the middle of the month instead of the end
            # time = data_avg.time.values
            # time_orig = copy.deepcopy(time)
            # for i,t in enumerate(time):
            #    if i==0: dt = pd.Timedelta('15 days')
            #    if i>=1: dt = ( time_orig[i] - time_orig[i-1] ) / 2
            #    time[i] = time_orig[i] - dt
            # data['time'] = time
         #-------------------------------------------------------------------------
         data = data.mean(dim='lon')
         data_list.append( data.values )
         lat_list.append( data['lat'].values/1e2 )
         lev_list.append( data['plev'].values )
         # #-------------------------------------------------------------------------
         # # Calculate time and zonal mean
         # if 'ncol' in ds.dims:
         #    with warnings.catch_warnings():
         #       warnings.simplefilter("ignore", category=RuntimeWarning)
         #       bin_ds = hapy.bin_YbyX( data, lat, bin_min=lat1, bin_max=lat2, bin_spc=dlat, \
         #                             wgt=area, keep_lev=True )
         #    lat_bins = bin_ds['bins'].values
         #    # data_binned = np.ma.masked_invalid( bin_ds['bin_val'].transpose().values )
         #    data_binned = bin_ds['bin_val'].transpose()
      #-------------------------------------------------------------------------
      #    tmp_file = get_tmp_file_name(case[c],var_list[v])
      #    print(' '*4+'writing to file: '+tmp_file)
      #    ds_tmp = xr.Dataset() # coords=data_binned.coords )
      #    ds_tmp[var_list[v]] = data_binned
      #    ds_tmp['lat']  = lat_bins
      #    ds_tmp['lev']  = data['lev'].values
      #    ds_tmp.to_netcdf(path=tmp_file,mode='w')
      # else:
      #    tmp_file = get_tmp_file_name(case[c],var_list[v])
      #    print(' '*4+'reading from file: '+tmp_file)
      #    ds_tmp = xr.open_dataset( tmp_file )
      #    data_binned = ds_tmp[var_list[v]]

      # data_list.append( data_binned.values )
      # lat_list.append( ds_tmp['lat'].values )
      # lev_list.append( ds_tmp['lev'].values )


   #----------------------------------------------------------------------------
   # calculate limits for common color bar
   data_min = np.min([np.nanmin(d) for d in data_list])
   data_max = np.max([np.nanmax(d) for d in data_list])
   if data_min==data_max: raise ValueError(hapy.tclr.RED+'WARNING: Difference is zero!'+hapy.tclr.END)
   #----------------------------------------------------------------------------
   tres = copy.deepcopy(res)

   # tres.cnLevels = np.arange(-30,30+3,3)/1e3
   tres.cnLevelSelectionMode = "ExplicitLevels"
   if not hasattr(tres,'cnLevels') : 
      cmin,cmax,cint,clev = ngl.nice_cntr_levels(data_min, data_max, cint=None, max_steps=21, \
                                                 returnLevels=True, aboutZero=False )
      tres.cnLevels = np.linspace(cmin,cmax,num=21)
   #----------------------------------------------------------------------------
   # Create plot
   for c in range(num_case):
      case_opts = opts_list[c]

      tres.cnFillMode = "AreaFill"
      # tres.cnFillMode = "RasterFill"
      tres.sfXArray = lat_list[c]
      tres.sfYArray = lev_list[c]

      ip = v*num_case+c if var_x_case else c*num_var+v

      plot[ip] = ngl.contour(wks,np.ma.masked_invalid(data_list[c]),tres)
   
      set_subtitles(wks, plot[ip], 
                    left_string=case_opts['n'], 
                    center_string='', 
                    right_string=f'{var_list[v]}',
                    font_height=subtitle_font_height)

# #---------------------------------------------------------------------------------------------------
# # Finalize plot
# #---------------------------------------------------------------------------------------------------

layout = [num_var,num_case] if var_x_case else [num_case,num_var]
# layout = [int(np.ceil(len(plot)/float(num_plot_col))),num_plot_col]

pnl_res = ngl.Resources()
pnl_res.nglPanelYWhiteSpacePercent       = 5
pnl_res.nglPanelXWhiteSpacePercent       = 5
pnl_res.nglPanelFigureStrings            = list(string.ascii_lowercase)
pnl_res.nglPanelFigureStringsJust        = "TopLeft"
pnl_res.nglPanelFigureStringsFontHeightF = 0.015

ngl.panel(wks,plot,layout,pnl_res)
hapy.trim_png(fig_file)

ngl.end()
#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------

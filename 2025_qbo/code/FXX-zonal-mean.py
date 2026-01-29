import os, ngl, xarray as xr, numpy as np
import hapy_common as hc, hapy_E3SM as he, hapy_setres as hs
import copy, cftime, warnings
import cmocean
#-------------------------------------------------------------------------------
case_dir,case_sub = [],[]
case,name = [],[]
clr,dsh,mrk = [],[],[]
def add_case(case_in,n=None,p=None,s=None,d=0,c='black',m=0):
   global name,case,case_dir,case_sub,clr,dsh,mrk
   case.append(case_in); name.append(n)
   case_dir.append(p); case_sub.append(s)
   dsh.append(d) ; clr.append(c) ; mrk.append(m)
var,lev_list = [],[]
def add_var(var_name,lev=-1): var.append(var_name); lev_list.append(lev)
#-------------------------------------------------------------------------------
cscratch = '/global/cscratch1/sd/whannah/e3sm_scratch/cori-knl'
gscratch = '/global/cfs/cdirs/m3312/whannah/2022-QBO-TEST'
pscratch = '/pscratch/sd/w/whannah/e3sm_scratch/pm-cpu'

# add_case('MAC-PG',    n='MAC')
# add_case('GPM-PG',    n='IMERG')

add_case('E3SM.QBO-TEST.F2010.ne30pg2.L72.01',       n='control',     c='red'  ,p=gscratch,s='run')
add_case('E3SM.QBO-TEST.F2010.ne30pg2.L72-nsu40.01', n='L72 smoothed',c='green',p=gscratch,s='run')
# add_case('E3SM.QBO-TEST.F2010.ne30pg2.L80-rsu40.01', n='L80 refined', c='blue' ,p=gscratch,s='run')

add_case('E3SM.QBO-TEST-02.F2010.ne30pg2.L72.GW-MOD-0', n='E3SM L72',      d=1,c='black', p=pscratch,s='run')
add_case('E3SM.QBO-TEST-02.F2010.ne30pg2.L72-nsu40',    n='E3SM L72-nsu40',d=1,c='red',   p=pscratch,s='run')
# add_case('E3SM.QBO-TEST-02.F2010.ne30pg2.L72-rscl',     n='E3SM L72-rscl', d=1,c='purple',p=pscratch,s='run')
# add_case('E3SM.QBO-TEST-02.F2010.ne30pg2.L72-rlim',     n='E3SM L72-rlim', d=1,c='pink',  p=pscratch,s='run')
#-------------------------------------------------------------------------------

lev = np.array([1,2,3,5,7,10,20,35,50,75,100,125,150,200,250,300,350,400,450,500,
               550,600,650,700,750,800,825,850,875,900,925,950,975])

# lev = np.array([100,200,300,400,500,600,700,800,900,950])

# add_var('U')
# add_var('OMEGA')
add_var('T')
add_var('Q')
# add_var('CLDLIQ')
# add_var('CLDICE')
# add_var('CLOUD')
# add_var('QRS')
# add_var('QRL')
# add_var('O3')


plot_diff = False

var_x_case = True

htype,years,months,first_file,num_files = 'h0',[],[],0,10*12
# htype,years,months,first_file,num_files = 'h0',[],[],0,20*12


fig_file = 'figs/FXX-zonal-mean'

lat1, lat2, dlat = -88., 88., 2

recalculate = True

tmp_file_head = 'data/zonal-mean'

#---------------------------------------------------------------------------------------------------
# Set up plot resources
#---------------------------------------------------------------------------------------------------
num_var = len(var)
num_case = len(case)

wkres = ngl.Resources()
npix = 2048; wkres.wkWidth,wkres.wkHeight=npix,npix

wks = ngl.open_wks('png',fig_file,wkres)
plot = [None]*(num_var*num_case)
res = hs.res_contour_fill()
res.vpHeightF = 0.3
res.trYReverse = True
# res.tiXAxisString = 'Latitude'
res.tiXAxisString = 'sin( Latitude )'
res.tiYAxisString = 'Pressure [hPa]'

res.tmYLLabelFontHeightF   = 0.01
res.tmXBLabelFontHeightF   = 0.01
res.lbLabelFontHeightF     = 0.015


#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------
# def get_data_dir(c):
#    global case_sub
#    data_dir_tmp = None
#    # if use_remap: data_sub_tmp = f'data_remap_{remap_grid}/'
#    if case_dir[c] is not None: data_dir_tmp = case_dir[c]
#    return data_dir_tmp

# def get_data_sub(c):
#    global case_dir
#    data_sub_tmp = None
#    if case_sub[c] is not None: data_sub_tmp = case_sub[c]
#    return data_sub_tmp

#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------
for v in range(num_var):
   print('  var: '+hc.tcolor.CYAN+var[v]+hc.tcolor.ENDC)
   data_list = []
   bins_list = []
   levs_list = []
   for c in range(num_case):
      print('    case: '+hc.tcolor.GREEN+case[c]+hc.tcolor.ENDC)

      data_dir_tmp,data_sub_tmp = None, None
      # if use_remap: data_sub_tmp = f'data_remap_{remap_grid}/'
      if case_dir[c] is not None: data_dir_tmp = case_dir[c]
      if case_sub[c] is not None: data_sub_tmp = case_sub[c]

      case_obj = he.Case( name=case[c], data_dir=data_dir_tmp, data_sub=data_sub_tmp  )

      # case_obj = he.Case( name=case[c])

      case_obj.set_coord_names(var[v])

      ip = v*num_case+c if var_x_case else c*num_var+v

      #-------------------------------------------------------------------------
      # read the data
      #-------------------------------------------------------------------------
      tmp_file = f'{tmp_file_head}.{case[c]}.{var[v]}.nc'
      print(f'      tmp_file: {tmp_file}')
      if recalculate :
         if 'lon1' in vars() : case_obj.lon1 = lon1
         if 'lon2' in vars() : case_obj.lon2 = lon2

         lat  = case_obj.load_data('lat',  htype=htype)
         area = case_obj.load_data('area', htype=htype).astype(np.double)
         data = case_obj.load_data(var[v],htype=htype, years=years, months=months, first_file=first_file, num_files=num_files, lev=lev )

         # hc.print_stat(data)

         if 'time' in data.dims : 
            time = data.time
            hc.print_time_length(data.time,indent=' '*6)

         if 'ilev' in data.dims : data = data.rename({'ilev':'lev'})

         #----------------------------------------------------------------------
         # Calculate time and zonal mean
         #----------------------------------------------------------------------
         with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=RuntimeWarning)

            bin_ds = hc.bin_YbyX( data.mean(dim='time', skipna=True), lat, \
                                  bin_min=lat1, bin_max=lat2, \
                                  bin_spc=dlat, wgt=area, keep_lev=True )

         lat_bins = bin_ds['bins'].values
         sin_lat_bins = np.sin(lat_bins*np.pi/180.)

         tmp_lev = bin_ds['lev']
         tmp_data = bin_ds['bin_val']
         #----------------------------------------------------------------------
         # write to temporary file
         #----------------------------------------------------------------------
         tmp_ds = xr.Dataset()
         # tmp_ds['bins']         = ( ('bins'), sin_lat_bins )
         # tmp_ds['lev']          = ( ('lev'), tmp_lev )
         # tmp_ds['data_binned']  = ( ('bins','lev'), tmp_data )
         # if 'time'in locals(): tmp_ds['time'] = time
         tmp_ds['bins']         = sin_lat_bins
         tmp_ds['lev']          = tmp_lev
         tmp_ds['data_binned']  = tmp_data
         if 'time'in locals(): tmp_ds['time'] = time
         print(f'      writing to file: {tmp_file}')
         tmp_ds.to_netcdf(path=tmp_file,mode='w')

      else:
         tmp_ds = xr.open_dataset( tmp_file )
         tmp_data     = tmp_ds['data_binned']
         sin_lat_bins = tmp_ds['bins'].values
         tmp_lev      = tmp_ds['lev']

      
      data_list.append( np.ma.masked_invalid( tmp_data.transpose().values ) )
      bins_list.append( sin_lat_bins )
      levs_list.append( tmp_lev.values )

      if 'time'in locals(): del time
   #------------------------------------------------------------------------------------------------
   # Plot zonally averaged data
   #------------------------------------------------------------------------------------------------
   if plot_diff:
      # Find min and max difference values for setting color bar
      tmp_data = data_list - data_list[0]
      for c in range(num_case): tmp_data[c] = data_list[c] - data_list[0]
      diff_data_min = np.min([np.min(d) for d in tmp_data])
      diff_data_max = np.max([np.max(d) for d in tmp_data])
      # diff_data_min = -2. * np.std(tmp_data)
      # diff_data_max =  2. * np.std(tmp_data)

      # print(f'diff_data_min: {diff_data_min}')
      # print(f'diff_data_max: {diff_data_max}')

   for c in range(num_case):
      data_dir_tmp,data_sub_tmp = None, None
      # if use_remap: data_sub_tmp = f'data_remap_{remap_grid}/'
      if case_dir[c] is not None: data_dir_tmp = case_dir[c]
      if case_sub[c] is not None: data_sub_tmp = case_sub[c]

      case_obj = he.Case( name=case[c], data_dir=data_dir_tmp, data_sub=data_sub_tmp  )

      ip = v*num_case+c if var_x_case else c*num_var+v

      #-------------------------------------------------------------------------
      # Set colors, contour levels, and plot strings
      #-------------------------------------------------------------------------
      tres = copy.deepcopy(res)
      # if var[v] in ['OMEGA']  : tres.cnFillPalette = "BlueWhiteOrangeRed"
      # if var[v]=="OMEGA"      : tres.cnLevels = np.linspace(-1,1,21)*0.1
      # if var[v] in ['Q','DYN_Q']: tres.cnLevels = np.linspace(1,16,16)

      # if plot_diff and c>0 : 
      #    if var[v]=="U" : tres.cnLevels = np.linspace(-10,10,11)

      if hasattr(tres,'cnLevels') and not (plot_diff and c>0) : 
         tres.cnLevelSelectionMode = "ExplicitLevels"
      else:
         if (var[v] in ['U','V']) or (plot_diff and c>0) : 
            aboutZero = True
         else:
            aboutZero = False
         
         if plot_diff and c>0 : 
            aboutZero = True
            data_min = diff_data_min
            data_max = diff_data_max
         else:
            data_min = np.min([np.min(d) for d in data_list])
            data_max = np.max([np.max(d) for d in data_list])
         
         cmin,cmax,cint,clev = ngl.nice_cntr_levels(data_min, data_max, cint=None, max_steps=21, \
                                                    returnLevels=True,aboutZero=aboutZero )
         tres.cnLevels = np.linspace(cmin,cmax,num=21)
         tres.cnLevelSelectionMode = "ExplicitLevels"

      # if plot_diff and c>0 : 
         # if var[v]=="U" : tres.cnLevels = np.linspace(-10,10,11)

      var_str = var[v]
      # if var[v]=="PRECT" : var_str = "Precipitation [mm/day]"

      #-------------------------------------------------------------------------
      # Create plot
      #-------------------------------------------------------------------------
      # tres.sfXCStartV = min( lat_bins )
      # tres.sfXCEndV   = max( lat_bins )
      tres.sfXArray = bins_list[c]
      # tres.sfXCStartV = -1. #np.min( sin_lat_bins )
      # tres.sfXCEndV   =  1. #np.max( sin_lat_bins )

      lat_tick = np.array([-90,-60,-30,0,30,60,90])
      res.tmXBMode = "Explicit"
      res.tmXBValues = np.sin( lat_tick*3.14159/180. )
      res.tmXBLabels = lat_tick

      # if 'lev' in vars() : tres.sfYCStartV = min( bin_ds['lev'].values )
      # if 'lev' in vars() : tres.sfYCEndV   = max( bin_ds['lev'].values )
      tres.sfYCStartV = min( levs_list[c] )
      tres.sfYCEndV   = max( levs_list[c] )

      # hc.print_stat(data_list[c],compact=True)
      # hc.print_stat(np.sin( lat_tick*3.14159/180. ),compact=True)
      # hc.print_stat(bin_ds['lev'].values)

      if plot_diff and c >0 :
         plot[ip] = ngl.contour(wks, data_list[c] - data_list[0], tres)
      else:
         plot[ip] = ngl.contour(wks, data_list[c], tres)

      var_str = var[v]
      if plot_diff and c>0 : var_str = var_str+' (diff)'

      # name_str = case_obj.short_name
      name_str = name[c]
      if plot_diff and c==0 : name0 = name_str
      # if plot_diff and c >0 : name_str = name_str+' - '+name0

      hs.set_subtitles(wks, plot[ip], name_str, '', var_str, font_height=0.01)

#---------------------------------------------------------------------------------------------------
# Finalize plot
#---------------------------------------------------------------------------------------------------

layout = [num_var,num_case] if var_x_case else [num_case,num_var]

# layout = [int(np.ceil(len(plot)/float(num_plot_col))),num_plot_col]

# ngl.panel(wks,plot,layout,hs.setres_panel())
ngl.panel(wks,plot,layout,ngl.Resources())
ngl.end()

hc.trim_png(fig_file)
#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------

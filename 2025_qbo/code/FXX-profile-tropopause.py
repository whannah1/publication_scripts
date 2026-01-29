import os, ngl, copy, xarray as xr, numpy as np
import hapy_common as hc, hapy_E3SM as he, hapy_setres as hs
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
##------------------------------------------------------------------------------
cscratch = '/global/cscratch1/sd/whannah/e3sm_scratch/cori-knl'
gscratch = '/global/cfs/cdirs/m3312/whannah/2022-QBO-TEST'
pscratch = '/pscratch/sd/w/whannah/e3sm_scratch/pm-cpu'

scratch = '/global/cfs/cdirs/m4310/whannah/E3SM'
# add_case('E3SM.2023-SCIDAC-v2-AMIP.ne30pg2_EC30to60E2r2.L72' ,n='E3SM L72', c='red',p=scratch,s='archive/atm/hist')
add_case('E3SM.2023-SCIDAC-v2-AMIP.ne30pg2_EC30to60E2r2.L80' ,n='E3SM L80', c='red',p=scratch,s='archive/atm/hist')
# add_case('E3SM.2023-SCIDAC-v2-AMIP.ne30pg2_EC30to60E2r2.L128',n='E3SM L128',p=scratch,s='run')

scratch = '/global/cfs/cdirs/m4310/wandiyu/E3SM/'
# add_case('E3SM.2023-SCIDAC.ne30pg2_EC30to60E2r2.F20TR-MMF1.L64',n='MMF L64', c='',  p=scratch,s='archive/atm/hist')
add_case('E3SM.2023-SCIDAC.ne30pg2_EC30to60E2r2.F20TR-MMF1.L72',n='MMF L72', c='blue',  p=scratch,s='archive/atm/hist')

add_obs = True
obs_case = 'ERA5'
#-------------------------------------------------------------------------------
# lev = np.array([5,8,10,15,20,25,30,40,50,100,150,200])
# plev = np.array([0.1,0.5,1,2,3,4,5,6,7,8,9,10,15,20,25,30,40,50,100,150,200,250,300,400,500,600,700,800,850,900])
# plev = np.array([0.1,0.5,1,2,3,4,5,6,7,8,9,10,15,20,25,30,40,50,100,150,200])
plev = np.array([20,30,40,50,60,70,80,90,100,110,120,130,140,150])

# var = []
add_var('T')
# add_var('RH')
# add_var('OMEGA')
# add_var('THETA')
# add_var('CLOUD')
# add_var('CLDLIQ')
# add_var('CLDICE')
# add_var('QRL')
# add_var('QRS')
# add_var('O3')


num_plot_col = len(var)


fig_type = "png"
fig_file = 'figs/FXX-profile-tropopause'

lat1,lat2 = -5,5

htype,first_file,num_files = 'h0',0,5*12

# Y-axis options
use_height_coord = False
print_stats      = False
print_profile    = False

tmp_file_head = 'data/profile'

#---------------------------------------------------------------------------------------------------
# Set up plot resources
#---------------------------------------------------------------------------------------------------
num_var,num_case = len(var),len(case)

if add_obs:
   clr.insert(0,'black')
   dsh.insert(0,0)

if 'lev' not in vars(): lev = np.array([0])

wks = ngl.open_wks(fig_type,fig_file)
plot = [None]*(num_var)
res = hs.res_xy()
# res.vpWidthF = 0.4
res.xyMarkerSizeF = 0.008
res.xyMarker = 16
res.xyLineThicknessF = 8
# res.tmYLLabelFontHeightF         = 0.008
# res.tmXBLabelFontHeightF         = 0.008

res.tmXBAutoPrecision = False
res.tmXBPrecision = 2

if use_height_coord: 
   res.tiYAxisString = 'Height [km]'
else:
   res.tiYAxisString = 'Pressure [hPa]'
   res.trYReverse = True
   # res.xyYStyle = 'Log'

# print(hc.tcolor.RED+'WARNING - limiting plot to 100mb and above'+hc.tcolor.ENDC)
# res.trYMaxF = 100.

def get_comp(case):
   comp = 'eam'
   if 'CESM' in case: comp = 'cam'
   return comp
#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------
data_list_list,lev_list_list = [],[]
for v in range(num_var):
   hc.printline()
   print(hc.tcolor.GREEN+'  var: '+var[v]+hc.tcolor.ENDC)
   data_list,lev_list = [],[]
   #-------------------------------------------------------------------------------
   # Load obs data
   if add_obs:
      obs_root = '/global/cfs/cdirs/e3sm/diagnostics/observations/Atm/climatology_1985-2014'
      if obs_case=='ERA5': 
         obs_grid_file = os.getenv('HOME')+f'/E3SM/data_grid/721x1440_ERA5_scrip.nc'
         obs_data_file = f'{obs_root}/ERA5/ERA5_ANN_198501_201412_climo.nc'
         if var[v]=='T': obs_var='ta'
      grd_ds = xr.open_mfdataset( obs_grid_file )
      obs_ds = xr.open_mfdataset( obs_data_file )
      lat  = grd_ds['grid_center_lat'].rename({'grid_size':'ncol'})
      lon  = grd_ds['grid_center_lon'].rename({'grid_size':'ncol'})
      area = grd_ds['grid_area'].rename({'grid_size':'ncol'})
      data = obs_ds[obs_var].isel(time=0).stack(ncol=('lat','lon'))
      
      mask = xr.DataArray( np.ones(area.shape,dtype=bool), dims=('ncol') )
      mask = mask & (lat>=lat1) & (lat<=lat2)
      area = area.where( mask, drop=True)
      data = data.where( mask, drop=True)
      data = ( (data*area).sum(dim=['ncol']) / area.sum(dim=['ncol']) )#.values 

      data.load()

      data['plev'] = data['plev']/1e2
      plev_min = np.min(plev)
      plev_max = np.max(plev)
      i_min = None
      i_max = None
      for i,p in enumerate(data['plev'].values):
         if p>=plev_min: i_max = i
         if p<=plev_max and i_min is None: i_min = i
      data = data.isel(plev=slice(i_min,i_max+1))
      # print(i_min)
      # print(i_max)
      # print(); print(data)
      # exit()

      hc.print_stat(data,        name=var[v],compact=True,indent=' '*6)
      hc.print_stat(data['plev'],name=var[v],compact=True,indent=' '*6)

      data_list.append( data.values )
      lev_list.append( data['plev'].values )
      # if use_height_coord:
      #    lev_list.append( Z.values )
      # else:
      #    lev_list.append( data['lev'].values )
   #-------------------------------------------------------------------------------
   # Load simulation data
   for c in range(num_case):
      print(hc.tcolor.CYAN+'    case: '+case[c]+hc.tcolor.ENDC)

      data_dir_tmp,data_sub_tmp = None, None
      # if use_remap: data_sub_tmp = f'data_remap_{remap_grid}/'
      if case_dir[c] is not None: data_dir_tmp = case_dir[c]
      if case_sub[c] is not None: data_sub_tmp = case_sub[c]
      
      #-------------------------------------------------------------------------
      #-------------------------------------------------------------------------
      # tmp_file = f'{tmp_file_head}.{case[c]}.{var[v]}'
      # if 'lat1' in vars(): tmp_file = tmp_file + f'.lat1_{lat1}.lat2_{lat2}'
      # if 'lon1' in vars(): tmp_file = tmp_file + f'.lon1_{lon1}.lon2_{lon2}'
      # tmp_file = f'{tmp_file}.nc'
      # print(f'      tmp_file: {tmp_file}')

      # if recalculate :
      if True:

         case_obj = he.Case( name=case[c], atm_comp=get_comp(case[c]), data_dir=data_dir_tmp, data_sub=data_sub_tmp  )

         #-------------------------------------------------------------------------
         # read the data   
         if 'lat1' in vars(): case_obj.lat1 = lat1; case_obj.lat2 = lat2
         if 'lon1' in vars(): case_obj.lon1 = lon1; case_obj.lon2 = lon2
         
         data = case_obj.load_data(var[v],htype=htype,first_file=first_file,num_files=num_files,lev=plev)
         area = case_obj.load_data('area',htype=htype,first_file=first_file,num_files=1,).astype(np.double)   
         Z = case_obj.load_data('Z3',htype=htype,first_file=first_file,num_files=num_files)

         #-------------------------------------------------------------------------
         if 'time' in data.dims : 
            time = data.time
            hc.print_time_length(data.time,indent=' '*6,print_span=True, print_length=False)

         data = data.mean(dim='time')
         Z = Z.mean(dim='time')
         #-------------------------------------------------------------------------
         # area weighted spatial average 
         Z    = ( (Z   *area).sum(dim='ncol') / area.sum(dim='ncol') )#.values 
         data = ( (data*area).sum(dim='ncol') / area.sum(dim='ncol') )#.values 

      #    #----------------------------------------------------------------------
      #    # write to temporary file
      #    #----------------------------------------------------------------------
      #    tmp_ds = xr.Dataset()
      #    tmp_ds['data'] = data
      #    tmp_ds['time'] = time
      #    if 'lat1' in vars(): tmp_ds['lat1']=lat1; tmp_ds['lat2']=lat2
      #    if 'lon1' in vars(): tmp_ds['lon1']=lon1; tmp_ds['lon2']=lon2
      #    tmp_ds['Z']    = Z
      #    print(f'      writing to file: {tmp_file}')
      #    tmp_ds.to_netcdf(path=tmp_file,mode='w')

      # else:
      #    tmp_ds = xr.open_dataset( tmp_file )
      #    data = tmp_ds['data']
      #    time = tmp_ds['time']
      #    Z    = tmp_ds['Z']
      
      hc.print_stat(data,       name=var[v],compact=True,indent=' '*6)
      hc.print_stat(data['lev'],name=var[v],compact=True,indent=' '*6)

      #-------------------------------------------------------------------------
      if print_profile:
         print()
         for xx in data.values: 
            print(f'    {xx}')
         print()
      #-------------------------------------------------------------------------
      # append final data to list
      data_list.append( data.values )
      if use_height_coord:
         lev_list.append( Z.values )
      else:
         lev_list.append( data['lev'].values )

   data_list_list.append(data_list)
   lev_list_list.append(lev_list)
   
#-------------------------------------------------------------------------------
# Create plot - overlay all cases for each var
#-------------------------------------------------------------------------------
for v in range(num_var):
   data_list = data_list_list[v]
   lev_list = lev_list_list[v]
   
   tres = copy.deepcopy(res)
   
   data_min = np.min([np.nanmin(d) for d in data_list])
   data_max = np.max([np.nanmax(d) for d in data_list])
   tres.trXMinF = data_min
   tres.trXMaxF = data_max

   ip = v


   num_plot_case = num_case+1 if add_obs else num_case
   for c in range(num_plot_case):
      tres.xyLineColor   = clr[c]
      tres.xyMarkerColor = clr[c]
      tres.xyDashPattern = dsh[c]

      tplot = ngl.xy(wks, data_list[c], lev_list[c], tres)

      if c==0:
         plot[ip] = tplot
      else:
         ngl.overlay(plot[ip],tplot)

   ### add vertical line
   lres = hs.res_xy()
   lres.xyLineThicknessF = 1
   lres.xyDashPattern = 0
   lres.xyLineColor = 'black'
   ngl.overlay(plot[ip],ngl.xy(wks, np.array([0,0]), np.array([-1e3,1e8]), lres))

   reg_str = ''
   var_str = var[v]


   if 'lat1' in locals(): 
      lat1_str = f'{lat1}N' if lat1>=0 else f'{(lat1*-1)}S'
      lat2_str = f'{lat2}N' if lat2>=0 else f'{(lat2*-1)}S'
      reg_str += f' {lat1_str}:{lat2_str} '
   if 'lon1' in locals(): 
      lon1_str = f'{lon1}E' #if lon1>=0 and lon1<=360 else f'{(lon1*-1)}S'
      lon2_str = f'{lon2}E' #if lon2>=0 and lon2<=360 else f'{(lon2*-1)}S'
      reg_str += f' {lon1_str}:{lon2_str} '


   hs.set_subtitles(wks, plot[ip], var_str, '', reg_str, font_height=0.01)


#-------------------------------------------------------------------------------
# Add legend
#-------------------------------------------------------------------------------
if num_case>1:
   lgres = ngl.Resources()
   lgres.vpWidthF           = 0.05
   lgres.vpHeightF          = 0.08
   lgres.lgLabelFontHeightF = 0.012
   lgres.lgMonoDashIndex    = True
   lgres.lgLineLabelsOn     = False
   lgres.lgLineThicknessF   = 8
   lgres.lgLabelJust        = 'CenterLeft'
   lgres.lgLineColors       = clr
   lgres.lgDashIndexes      = dsh

   lx,ly = 0.5,0.45
   if num_var==2: lx,ly = 0.3,0.45
   if num_var==4: lx,ly = 0.05,0.5

   # pid = ngl.legend_ndc(wks, len(case_name), case_name, lx, ly, lgres)

#---------------------------------------------------------------------------------------------------
# Finalize plot
#---------------------------------------------------------------------------------------------------

layout = [int(np.ceil(len(plot)/float(num_plot_col))),num_plot_col]

#-- draw a common title string on top of the panel
textres               =  ngl.Resources()
# textres.txFontHeightF =  0.01                  #-- title string size
# ngl.text_ndc(wks,f'time step = {ss_t}',0.5,.97,textres)  #-- add title to plot
textres.txFontHeightF =  0.02                  #-- title string size
if layout[0]==1: y_pos = 0.7
if layout[0]>=2: y_pos = 0.9
# ngl.text_ndc(wks,f'time step = {ss_t}',0.5,y_pos,textres)  #-- add title to plot

pres = hs.setres_panel()
pres.nglPanelTop      =  0.93

ngl.panel(wks,plot,layout,pres)
ngl.end()

hc.trim_png(fig_file)
#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------

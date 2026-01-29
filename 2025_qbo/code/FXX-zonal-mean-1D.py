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
add_case('E3SM.QBO-TEST.F2010.ne30pg2.L80-rsu40.01', n='L80 refined', c='blue' ,p=gscratch,s='run')

# add_case('E3SM.QBO-TEST-02.F2010.ne30pg2.L72.GW-MOD-0', n='E3SM L72 M0',c='red'  ,d=1,p=pscratch,s='run')
# add_case('E3SM.QBO-TEST-02.F2010.ne30pg2.L72.GW-MOD-1', n='E3SM L72 M1',c='green',d=1,p=pscratch,s='run')
#-------------------------------------------------------------------------------

# add_var('PRECT')
# add_var('TMQ')
# add_var('LHFLX')
# add_var('SHFLX')
# add_var('P-E')
# add_var('TS')
# add_var('TGCLDLWP')
# add_var('TGCLDIWP')
# add_var('SWCF')
# add_var('LWCF')


add_var('FSNT')
add_var('FLNT')
# add_var('FSNTOA')  # use with CERES!
# add_var('FLUT')    # use with CERES!
# add_var('NET_TOA_RAD')

# add_var('FSNS');add_var('FLNS')
# add_var('LWCF');add_var('SWCF')
# add_var('TAUX')
# add_var('TAUY')
# add_var('WSPD_BOT')
# add_var('CLDLOW')
# add_var('CLDHGH')
# add_var('CLDTOT')

### use for h0
# add_var('U',lev=850)
# add_var('U',lev=200)

num_plot_col = 2


htype,years,months,first_file,num_files = 'h0',[],[],0,20*12

fig_file = 'figs/FXX-zonal-mean-1D'

bin_dlat = 2

plot_diff = False

print_rmse  = False
print_stats = True

add_legend = True

#---------------------------------------------------------------------------------------------------
# Set up plot resources
#---------------------------------------------------------------------------------------------------
num_var = len(var)
num_case = len(case)

wkres = ngl.Resources()
npix=4096; wkres.wkWidth,wkres.wkHeight=npix,npix
wks = ngl.open_wks('png',fig_file,wkres)
plot = [None]*num_var
res = hs.res_xy()
res.vpHeightF = 0.3
res.xyLineThicknessF = 10

if 'clr' not in locals(): 
   if num_case>1 : clr = np.linspace(2,len( ngl.retrieve_colormap(wks) )-1,num_case,dtype=int)
   else : clr = ['black']

# if num_case>1 and 'dsh' not in vars(): dsh = np.arange(0,num_case,1)
if 'dsh' not in locals(): 
   if num_case>1 : dsh = np.zeros(num_case)
   else : dsh = [0]
# res.xyLineColors   = clr
# res.xyDashPatterns = dsh

# res.tiXAxisString = 'Latitude'
res.tiXAxisString = 'sin( Latitude )'

#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------
msg_list = []
for v in range(num_var):
   print('  var: '+hc.tcolor.MAGENTA+var[v]+hc.tcolor.ENDC)
   data_list,std_list,cnt_list = [],[],[]
   bin_list = []
   glb_avg_list = []
   if 'lev_list' in locals(): lev = lev_list[v]
   for c in range(num_case):
      print('    case: '+hc.tcolor.GREEN+case[c]+hc.tcolor.ENDC)
      if case[c]!='CERES-EBAF':
         # case_obj = he.Case( name=case[c] )
         data_dir_tmp,data_sub_tmp = None, None
         # if use_remap: data_sub_tmp = f'data_remap_{remap_grid}/'
         if case_dir[c] is not None: data_dir_tmp = case_dir[c]
         if case_sub[c] is not None: data_sub_tmp = case_sub[c]

         case_obj = he.Case( name=case[c], data_dir=data_dir_tmp, data_sub=data_sub_tmp  )
      #-------------------------------------------------------------------------
      # read the data
      #-------------------------------------------------------------------------
      if case[c]=='CERES-EBAF':
         if use_obs_remap:
            file_list = sorted(glob.glob(f'{obs_data_path}/{obs_remap_sub}/*'))
         else:
            file_list = sorted(glob.glob(f'{obs_data_path}/clim_180x360/*'))
         file_list = file_list[:num_files]
         ds = xr.open_mfdataset(file_list,combine='by_coords',concat_dim='time').isel(nv=0)

         tvar = var[v]
         if tvar=='NET_TOA_RAD': tvar = 'RESTOA'
         if var[v]=='NET_CF': tvar = 'SWCF'
         
         data = ds[tvar]
         
         if var[v]=='NET_CF': data = data + ds['LWCF']

         if use_obs_remap:
            lat,lon = ds['lat'],ds['lon']
            # scrip_ds = xr.open_dataset(os.getenv('HOME')+f'/E3SM/data_grid/{CERES_remap_scrip}.nc')
            area = scrip_ds['grid_area'].rename({'grid_size':'ncol'})
         else:
            data = data.stack(ncol=('lat','lon'))
            lon2D = np.transpose( np.repeat( ds['lon'].values[...,None],len(ds['lat']),axis=1) )
            lat2D =               np.repeat( ds['lat'].values[...,None],len(ds['lon']),axis=1)
            lon = xr.DataArray(lon2D,dims=('lat','lon')).stack(ncol=('lat','lon'))
            lat = xr.DataArray(lat2D,dims=('lat','lon')).stack(ncol=('lat','lon'))
            fv_scrip_ds = xr.open_dataset(os.getenv('HOME')+'/E3SM/data_grid/cmip6_180x360_scrip.20181001.nc')
            area = fv_scrip_ds['grid_area'].rename({'grid_size':'ncol'})
            # convert ncol to DataArray - multiindex causes problems
            data['ncol'] = np.arange(data.shape[1])
            area['ncol'] = np.arange(data.shape[1])
            lat['ncol']  = np.arange(data.shape[1])
            lon['ncol']  = np.arange(data.shape[1])


         # xlon, ylat = np.meshgrid(lat1D,lon1D)
         # R = hc.earth_radius(ylat)
         # dlat = np.deg2rad(np.gradient(ylat, axis=0))
         # dlon = np.deg2rad(np.gradient(xlon, axis=1))
         # dy,dx = dlat * R , dlon * R * np.cos(np.deg2rad(ylat))
         # area = np.absolute(dy*dx) / np.square(R) # calculate area and convert to steridians
         # area = xr.DataArray(area,dims=('lat','lon')).stack(ncol=('lat','lon'))
         
         # area_alt = area; area_alt['ncol'] = np.arange(data.shape[1])


         # scrip_ds = xr.open_dataset(os.getenv('HOME')+'/E3SM/data_grid/cmip6_180x360_scrip.20181001.nc')
         area = scrip_ds['grid_area'].rename({'grid_size':'ncol'})

         # avoid dask array
         data.load()

         # convert ncol to DataArray - multiindex causes problems
         data['ncol'] = np.arange(data.shape[1])
         area['ncol'] = np.arange(data.shape[1])
         lat['ncol']  = np.arange(data.shape[1])

         
         # print(); print(area_alt)
         # print(); print(area)
         # print()
         # print(); hc.print_stat(area_alt)
         # print(); hc.print_stat(area)
         # print()
         # print(np.sum(area_alt.values))
         # print(np.sum(area.values))
         # exit()

         
      else:
         if 'lon1' in locals() : case_obj.lon1 = lon1
         if 'lon2' in locals() : case_obj.lon2 = lon2
         if 'lev'  in locals() : case_obj.lev  = lev

         lat  = case_obj.load_data('lat',  htype=htype,num_files=1)
         area = case_obj.load_data('area', htype=htype,num_files=1).astype(np.double)
         data = case_obj.load_data(var[v],htype=htype,years=years,months=months,
                                   first_file=first_file,num_files=num_files)

      hc.print_time_length(data.time,indent=(' '*6))

      
      if 'area' in locals() :
         gbl_mean = ( (data.mean(dim='time')*area).sum() / area.sum() ).values 
         ### this also works
         # dummy,area = xr.broadcast(data,area); del dummy
         # gbl_mean = ( (data*area).sum() / area.sum() ).values 
         ### but why doesn't this work...?
         # weights = np.cos(np.deg2rad(lat.values))
         # gbl_mean_alt = np.average(np.multiply(data.values,weights))
         # print(f'\n  global mean :  '+hc.tcolor.RED+f'{gbl_mean    }'+hc.tcolor.ENDC
         #                   +'  /  '+hc.tcolor.RED+f'{gbl_mean_alt}'+hc.tcolor.ENDC+'\n')
         glb_avg_list.append(gbl_mean)
      else:
         lat1D,lon1D = ds['lat'],ds['lon']
         lat2D = np.repeat( lat1D.values[...,None],len(lon1D),axis=1)
         weights = np.cos( np.deg2rad( lat2D ) )
         gbl_mean = np.average( np.multiply( data.mean(dim='time').values, weights ) )
         glb_avg_list.append(gbl_mean)

      if print_stats:
         msg = hc.print_stat(data,name=var[v],stat='naxsh',indent=(' '*6),compact=True)
         msg_list.append('  case: '+case[c]+'\n'+msg)
         if 'area' in locals() :
            print(f'      Area Weighted Global Mean : '+hc.tcolor.CYAN+f'{gbl_mean:6.4}'+hc.tcolor.ENDC)
            

      if print_rmse:
         if c==0:baseline = data
         if c>0:
            rmse = np.sqrt( np.mean( np.square( data.to_masked_array() - baseline.to_masked_array() )))
            print(f'      Root Mean Square Error    : {rmse:6.4}')
            # exit()
      #-------------------------------------------------------------------------
      # Calculate time and zonal mean
      #-------------------------------------------------------------------------

      bin_ds = hc.bin_YbyX( data.mean(dim='time'), lat, 
                           bin_min=-90, bin_max=90, 
                           bin_spc=bin_dlat, wgt=area )

      data_list.append( bin_ds['bin_val'].values )
      
      # std_list.append( bin_ds['bin_std'].values )
      # cnt_list.append( bin_ds['bin_cnt'].values )

      lat_bins = bin_ds['bins'].values


      sin_lat_bins = np.sin(lat_bins*np.pi/180.)

      bin_list.append(sin_lat_bins)

      if 'area' in locals(): del area

   #----------------------------------------------------------------------------
   # Take difference from first case
   #----------------------------------------------------------------------------
   if plot_diff :
      data_tmp = data_list
      data_baseline = data_list[0]
      for c in range(num_case): data_list[c] = data_list[c] - data_baseline
   #----------------------------------------------------------------------------
   # Create plot
   #----------------------------------------------------------------------------
   unit_str = ''
   if var[v] in ['PRECT','PRECC','PRECL']   : unit_str = '[mm/day]'
   if var[v] in ['LHFLX','SHFLX']           : unit_str = '[W/m2]'
   res.tiYAxisString = unit_str

   res.trXMinF = -1. #np.min( sin_lat_bins )
   res.trXMaxF =  1. #np.max( sin_lat_bins )

   res.trYMinF = np.min([np.nanmin(d) for d in data_list])
   res.trYMaxF = np.max([np.nanmax(d) for d in data_list])

   lat_tick = np.array([-90,-60,-30,0,30,60,90])
   res.tmXBMode = "Explicit"
   res.tmXBValues = np.sin( lat_tick*3.14159/180. )
   res.tmXBLabels = lat_tick

   # plot[v] = ngl.xy(wks, np.stack(bin_list), np.ma.masked_invalid(  np.stack(data_list) ), res)

   for c in range(num_case):
      res.xyLineColor   = clr[c]
      res.xyDashPattern = dsh[c]
      tplot = ngl.xy(wks, bin_list[c], np.ma.masked_invalid( data_list[c] ), res)
      if c==0:
         plot[v] = tplot
      else:
         ngl.overlay( plot[v], tplot )

   var_str = var[v]
   if var[v]=="PRECT" : var_str = "Precipitation"
   hs.set_subtitles(wks, plot[v], "", "", var_str, font_height=0.005)


#---------------------------------------------------------------------------------------------------
# Finalize plot
#---------------------------------------------------------------------------------------------------
# layout = [num_var,1]
layout = [int(np.ceil(len(plot)/float(num_plot_col))),num_plot_col]
# layout = [1,num_var]
# if num_var==4 : layout = [2,2]
# if num_var==6 : layout = [3,2]
pres = hs.setres_panel()
pres.nglFrame = False
pres.nglPanelRight = 0.5
ngl.panel(wks,plot,layout,pres)
#-------------------------------------------------------------------------------
# Add legend
#-------------------------------------------------------------------------------
if num_case>1 and add_legend:
   lgres = ngl.Resources()
   lgres.vpWidthF           = 0.05
   lgres.vpHeightF          = 0.01*num_case
   lgres.lgLabelFontHeightF = 0.005
   lgres.lgLabelFont        = "courier"
   lgres.lgMonoDashIndex    = False
   lgres.lgLineLabelsOn     = False
   lgres.lgLineThicknessF   = 16
   lgres.lgLabelJust        = 'CenterLeft'
   lgres.lgLineColors       = clr
   lgres.lgDashIndexes      = dsh

   labels = name
   max_label_len = max([len(n) for n in name])+2
   for n,nn in enumerate(name):
      labels[n] = f'  {nn:{max_label_len}}'
      if num_var==1: labels[n] += f'  {glb_avg_list[n]:6.1f}'

   # ndc_T,ndc_B,ndc_L,ndc_R = ngl.get_bounding_box(plot[0])
   # ndcx = ndc_R + 0.02
   # ndcy = np.average(np.array([ndc_T,ndc_B]))

   ndcx,ndcy = 0.5,0.55

   pid = ngl.legend_ndc(wks, len(labels), labels, ndcx, ndcy, lgres)

   # legend_id_list = [None]*num_var
   # for v in range(num_var):
   #    # ndcx, ndcy = ngl.datatondc(plot[v], 1., 0. )
   #    ndc_T,ndc_B,ndc_L,ndc_R = ngl.get_bounding_box(plot[v])
   #    ndcx = ndc_R + 0.02
   #    ndcy = np.average(np.array([ndc_T,ndc_B]))
   #    print(f'v: {v}   ndcx / ndcy: {ndcx} / {ndcy}')
   #    legend_id_list[v] = ngl.legend_ndc(wks, len(labels), labels, ndcx, ndcy, lgres)
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
ngl.frame(wks)
ngl.end()

# print()
# for msg in msg_list: print(msg)

hc.trim_png(fig_file)
# print(f'\n{fig_file}.png\n')
#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------
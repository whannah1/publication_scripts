import os, ngl, subprocess as sp, numpy as np, xarray as xr, copy, string
import hapy_common as hc, hapy_E3SM   as he, hapy_setres as hs
host = hc.get_host()
data_dir,data_sub = None,None
case,name,clr,dsh = [],[],[],[]
var,lev_list = [],[]
def add_case(case_in,n='',c='black',d=0):
   global name,case
   case.append(case_in); name.append(n); clr.append(c); dsh.append(d)

def add_var(var_name,lev=-1): var.append(var_name); lev_list.append(lev)
#-------------------------------------------------------------------------------
add_case('E3SM.CHX.ne30pg2_ne30pg2.F-MMF1.00',            n='MMF')
# add_case('E3SM.CHX.ne30pg2_ne30pg2.FC5AV1C-L.00',         n='E3SM')
# add_case('E3SM.CHX.ne30pg2_ne30pg2.FC5AV1C-L.no_dcape.00',n='E3SM (no dCAPE)')
#-------------------------------------------------------------------------------

chx_var = 'TGCLDLWP'

add_var('FLNT')
add_var('FSNT')
add_var('NET_TOA_RAD')

fig_type = 'png'
fig_file = 'figs/FXX-chx-composite.zonal-mean'



subset_min_length = 4
scratch = '/global/cscratch1/sd/whannah/e3sm_scratch/cori-knl/chx-detection-data'
tmp_file=f'{scratch}/occurence-partial-chx-over-time'
# tmp_file = 'data/chx-occurrence'

# lat1,lat2 = -40,40 ; lon1,lon2 = 90,260
# lat1,lat2 = 0,30 ; lon1,lon2 = 140,180+40

htype,first_file,num_files = 'h1',0,1*365
# htype,first_file,num_files = 'h1',0,30

print_stats = True
var_x_case = True

scripfile = xr.open_dataset('scrip_files/ne30pg2_scrip.nc')

#---------------------------------------------------------------------------------------------------
# Set up plot resources
#---------------------------------------------------------------------------------------------------
num_var,num_case = len(var),len(case)

wkres = ngl.Resources()
npix = 2048; wkres.wkWidth,wkres.wkHeight=npix,npix

wks = ngl.open_wks(fig_type,fig_file,wkres)
plot = [None]*(num_var)

#---------------------------------------------------------------------------------------------------
# mask for regional subset
#---------------------------------------------------------------------------------------------------
if 'lat1' in locals():
   ncol = scripfile['grid_size'].values
   mask = xr.DataArray( np.ones(len(ncol),dtype=bool), dims=('grid_size') )
   mask = mask & (scripfile['grid_center_lat']>=lat1) & (scripfile['grid_center_lat']<=lat2)
   mask = mask & (scripfile['grid_center_lon']>=lon1) & (scripfile['grid_center_lon']<=lon2)

   scripfile = scripfile.where(mask,drop=True)

   mask = mask.rename({'grid_size':'ncol'})

#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------
for v in range(num_var):
   hc.printline()
   print('  var: '+hc.tcolor.MAGENTA+var[v]+hc.tcolor.ENDC)
   data_list,pchx_list,nchx_list = [],[],[]
   for c in range(num_case):
      print('    case: '+hc.tcolor.GREEN+case[c]+hc.tcolor.ENDC)

      case_obj = he.Case( case[c], time_freq='daily' )
      case_obj.set_coord_names(var[v])
      #-------------------------------------------------------------------------
      # load the data
      #-------------------------------------------------------------------------   
      lat  = case_obj.load_data('lat', htype=htype,num_files=1)
      area = case_obj.load_data('area',htype=htype,num_files=1).astype(np.double)
      data = case_obj.load_data(var[v],htype=htype,first_file=first_file,num_files=num_files)

      ### Convert to daily mean
      data = data.resample(time='D').mean(dim='time')
      
      case_chx_file = f'{tmp_file}.daily.{case[c]}.{chx_var}.sml_{subset_min_length}.nc'
      cnt_ds = xr.open_dataset( case_chx_file )
      pchx_data = cnt_ds['cnt'][1,:,:]

      # if len(data['time']) < len(pchx_data['time']):
      #    ntime = len(data['time'])
      #    pchx_data = pchx_data[:ntime,:]

      if 'lat1' in locals():
         data = data.where(mask,drop=True)
         pchx_data = pchx_data.where(mask,drop=True)
      #-------------------------------------------------------------------------
      # average the data based on amount of checkerboard
      #-------------------------------------------------------------------------
      data.load()

      data_avg = data.mean(dim='time')
      pchx_avg = data.where(pchx_data==1).mean(dim='time')
      nchx_avg = data.where(pchx_data==0).mean(dim='time')

      data_avg.compute()
      pchx_avg.compute()
      nchx_avg.compute()

      #-------------------------------------------------------------------------
      # calculate zonal mean
      #-------------------------------------------------------------------------
      with np.errstate(divide='ignore', invalid="ignore"):
         data_bin_ds = hc.bin_YbyX( data_avg, lat, bin_min=-60, bin_max=60, bin_spc=2, wgt=area )
         nchx_bin_ds = hc.bin_YbyX( nchx_avg, lat, bin_min=-60, bin_max=60, bin_spc=2, wgt=area )
         pchx_bin_ds = hc.bin_YbyX( pchx_avg, lat, bin_min=-60, bin_max=60, bin_spc=2, wgt=area )

      data_avg = data_bin_ds['bin_val']
      nchx_avg = nchx_bin_ds['bin_val']
      pchx_avg = pchx_bin_ds['bin_val']

      if c==0: 
         lat_bins = data_bin_ds['bins'].values
         sin_lat_bins = np.sin(lat_bins*np.pi/180.)
      #-------------------------------------------------------------------------
      # print diagnostic stats
      #-------------------------------------------------------------------------
      if print_stats: 
         hc.print_stat(data_avg,name=f'{var[v]} - data_avg',stat='nax',indent=' '*6,compact=True)
         hc.print_stat(nchx_avg,name=f'{var[v]} - nchx_avg',stat='nax',indent=' '*6,compact=True)
         hc.print_stat(pchx_avg,name=f'{var[v]} - pchx_avg',stat='nax',indent=' '*6,compact=True)

      #-------------------------------------------------------------------------
      # Create plot
      #-------------------------------------------------------------------------
      res = hs.res_xy()
      res.vpHeightF = 0.3

      data_min = np.nanmin([ data_avg.values, nchx_avg.values, pchx_avg.values ])
      data_max = np.nanmax([ data_avg.values, nchx_avg.values, pchx_avg.values ])

      lat_tick = np.array([-90,-60,-30,0,30,60,90])
      res.tmXBMode = "Explicit"
      res.tmXBValues = np.sin( lat_tick*3.14159/180. )
      res.tmXBLabels = lat_tick
      res.trXMinF,res.trXMaxF = -1.,1.
      res.trYMinF,res.trYMaxF = data_min,data_max

      res.xyLineThicknessF = 8
      res.xyDashPatterns   = 0

      ip = v

      res.xyLineColor = 'black'
      plot[ip] = ngl.xy(wks, sin_lat_bins, np.ma.masked_invalid(data_avg.values), res) 
      
      res.xyLineColor = 'blue'
      ngl.overlay(plot[ip], ngl.xy(wks, sin_lat_bins, np.ma.masked_invalid(nchx_avg.values), res) )

      res.xyLineColor = 'red'
      ngl.overlay(plot[ip], ngl.xy(wks, sin_lat_bins, np.ma.masked_invalid(pchx_avg.values), res) )

      #-------------------------------------------------------------------------
      #-------------------------------------------------------------------------
      var_str = var[v]
      if var[v]=='PRECT':     var_str = 'Precipitation'
      if var[v]=='TGCLDLWP':  var_str = 'Liquid Water Path'
      if var[v]=='TGCLDIWP':  var_str = 'Ice Water Path'
      
      hs.set_subtitles(wks, plot[ip], name[c], '', var_str, font_height=0.015)

#---------------------------------------------------------------------------------------------------
# Finalize plot
#---------------------------------------------------------------------------------------------------

# layout = [num_var,num_case] if var_x_case else [num_case,num_var]
layout = [num_var,1]
   
pnl_res = hs.setres_panel()

### add panel labels
pnl_res.nglPanelFigureStrings            = list(string.ascii_lowercase)
pnl_res.nglPanelFigureStringsJust        = "TopLeft"
pnl_res.nglPanelFigureStringsFontHeightF = 0.01

pnl_res.nglPanelYWhiteSpacePercent = 5

ngl.panel(wks,plot,layout,pnl_res)
ngl.end()

hc.trim_png(fig_file)
#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------

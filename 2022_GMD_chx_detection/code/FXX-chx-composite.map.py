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
fig_file = 'figs/FXX-chx-composite.map'



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
plot = [None]*(num_var*4)

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
      data_avg = data.mean(dim='time')
      pchx_avg = data.where(pchx_data==1).mean(dim='time')
      nchx_avg = data.where(pchx_data==0).mean(dim='time')

      if print_stats: 
         hc.print_stat(data_avg,name=f'{var[v]} - pchx_avg',stat='naxh',indent=' '*6,compact=True)
         hc.print_stat(pchx_avg,name=f'{var[v]} - pchx_avg',stat='naxh',indent=' '*6,compact=True)
         hc.print_stat(nchx_avg,name=f'{var[v]} - nchx_avg',stat='naxh',indent=' '*6,compact=True)

      #-------------------------------------------------------------------------
      # append to data lists
      #-------------------------------------------------------------------------
      data_list.append( data_avg.values )
      pchx_list.append( pchx_avg.values )
      nchx_list.append( nchx_avg.values )

   #------------------------------------------------------------------------------------------------
   # determine color levels
   #------------------------------------------------------------------------------------------------
   pchx_min = np.min([np.nanmin(d) for d in pchx_list])
   pchx_max = np.max([np.nanmax(d) for d in pchx_list])

   nchx_min = np.min([np.nanmin(d) for d in nchx_list])
   nchx_max = np.max([np.nanmax(d) for d in nchx_list])

   data_min = np.min([pchx_min,nchx_min])
   data_max = np.max([pchx_max,nchx_max])

   aboutZero = False
   (cmin,cmax,cint) = ngl.nice_cntr_levels(data_min, data_max, cint=None, max_steps=21, 
                                          returnLevels=False, aboutZero=False )
   
   #------------------------------------------------------------------------------------------------
   # Define plot resources
   #------------------------------------------------------------------------------------------------
   res = hs.res_contour_fill_map()
   res.tmYLLabelFontHeightF         = 0.008
   res.tmXBLabelFontHeightF         = 0.008
   res.lbLabelFontHeightF           = 0.018
   res.tmXBOn                       = False
   res.tmYLOn                       = False
   if 'lat1' in vars() : res.mpMinLatF = lat1; res.mpMaxLatF = lat2
   if 'lon1' in vars() : res.mpMinLonF = lon1; res.mpMaxLonF = lon2

   res.cnLevels = np.linspace(cmin,cmax,num=21)
   res.cnLevelSelectionMode = 'ExplicitLevels'
   res.cnFillPalette = "MPL_viridis"

   res.cnFillMode    = "CellFill"
   res.sfXArray      = scripfile['grid_center_lon'].values
   res.sfYArray      = scripfile['grid_center_lat'].values
   res.sfXCellBounds = scripfile['grid_corner_lon'].values
   res.sfYCellBounds = scripfile['grid_corner_lat'].values

   for c in range(num_case):
      #-------------------------------------------------------------------------
      # Create plot
      #-------------------------------------------------------------------------
      # ip = v*num_case+c if var_x_case else c*num_var+v
      ip = v*4
      res.cnLevelSelectionMode = 'ExplicitLevels'
      plot[ip+0] = ngl.contour_map(wks,np.ma.masked_invalid(data_list[c]),res) 
      plot[ip+1] = ngl.contour_map(wks,np.ma.masked_invalid(nchx_list[c]),res) 
      plot[ip+2] = ngl.contour_map(wks,np.ma.masked_invalid(pchx_list[c]),res) 

      diff = pchx_list[c] - data_list[c]
      res.cnLevelSelectionMode = 'AutomaticLevels'
      plot[ip+3] = ngl.contour_map(wks,np.ma.masked_invalid(diff),res) 

      #-------------------------------------------------------------------------
      #-------------------------------------------------------------------------
      var_str = var[v]
      if var[v]=='PRECT':     var_str = 'Precipitation'
      if var[v]=='TMQ':       var_str = 'CWV'
      if var[v]=='TGCLDLWP':  var_str = 'Liquid Water Path'
      if var[v]=='TGCLDIWP':  var_str = 'Ice Water Path'
      
      hs.set_subtitles(wks, plot[ip], name[c], '', var_str, font_height=0.015)

#---------------------------------------------------------------------------------------------------
# Finalize plot
#---------------------------------------------------------------------------------------------------

# layout = [num_var,num_case] if var_x_case else [num_case,num_var]
layout = [num_var,4]
   
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

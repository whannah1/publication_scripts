# plot the fractional occurence of all possible patterns
import os, ngl, sys, numba, copy, xarray as xr, numpy as np, re, string
import hapy_common as hc, hapy_E3SM as he, hapy_setres as hs
import pg_checkerboard_utilities as pg
import cmocean

case,name,clr,dsh,pat = [],[],[],[],[]
var,var_name,lev_list = [],[],[]
def add_case(case_in,n='',c='black',d=0,p=0):
   global name,case
   case.append(case_in); name.append(n); clr.append(c); dsh.append(d); pat.append(p)

def add_var(v,n=None,lev=-1): var.append(v); var_name.append(n); lev_list.append(lev)
#-------------------------------------------------------------------------------
# add_case('MAC-FV',  n='MAC 1x1 deg'                                          ,c='black',p=0)
# add_case('MAC-PG',  n='MAC ne30pg2'                                          ,c='gray',p=0)
# add_case('ERA5-FV', n='ERA5 1x1 degree'                                      ,c='gray', p=0)
# add_case('ERA5-PG', n='ERA5 ne30pg2'                                         ,c='gray', p=0)
# add_case('E3SM.CHX.ne30pg2_ne30pg2.FC5AV1C-L.no_dcape.00',n='E3SM (no DCAPE)',c='green',p=0)
add_case('E3SM.CHX.ne30pg2_ne30pg2.F-MMF1.00',            n='E3SM-MMF'       ,c='blue' ,p=0)
# add_case('E3SM.CHX.ne30pg2_ne30pg2.FC5AV1C-L.00',         n='E3SM'           ,c='red'  ,p=0)
#-------------------------------------------------------------------------------

var = []
# add_var('PRECT')
add_var('TGCLDLWP',n='Liq Water Path')
# add_var('TGCLDIWP',n='Ice Water Path')
# add_var('TMQ')
# add_var('LHFLX')
# add_var('FLNT')
# add_var('FSNT')
# add_var('Q850')
# add_var('T850')
# add_var('U850')


# lat1,lat2 = -60,60
# lat1,lat2,lon1,lon2 = -20,40,80,180+30


fig_type = 'png'
fig_file = 'figs/F07-chx-time-series'

scratch = '/global/cscratch1/sd/whannah/e3sm_scratch/cori-knl/chx-detection-data'
tmp_file=f'{scratch}/occurence-partial-chx-over-time'


ocean_only = False

var_x_case           = True    # controls plot panel arrangement
print_stats          = True    #
verbose              = True    #


ilat,ilon,idx = -20,360-115,10
lat1,lat2 = ilat-idx,ilat+idx
lon1,lon2 = ilon-idx, ilon+idx

#---------------------------------------------------------------------------------------------------
# generate sets of possible neighbor states
#---------------------------------------------------------------------------------------------------

### only pure checkboard
# tmp_file='data/occurrence-chx-only'; rotate_sets=False; sets=pg.chx_only_sets

### find index of pure checkerboard
# chx_idx = []
# for s,tset in enumerate(sets):
#    for chx in chx_sets:
#       if all(tset==chx): chx_idx.append(s)

# (num_set,set_len) = sets.shape
# set_coord,nn_coord = np.arange(num_set),np.arange(set_len)
# sets.assign_coords(coords={'set':set_coord,'neighbors':nn_coord})

# set_labels = pg.get_set_labels(sets)

set_labels = ['no checkerboard','partial checkerboard']

#---------------------------------------------------------------------------------------------------
# Set up plot resources
#---------------------------------------------------------------------------------------------------
num_case,num_var = len(case),len(var)

wkres = ngl.Resources()
# npix=2**13; wkres.wkWidth,wkres.wkHeight=npix,npix # use this for plotting all patterns w/ rotation
wks = ngl.open_wks(fig_type,fig_file,wkres)
plot = [None]*(num_var)

res = hs.res_xy()
# res.tmXBOn = False
# res.tmYLOn = False

#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------
cnt_ds_list = []
for c in range(num_case):
   print('\n  case: '+hc.tcolor.CYAN+case[c]+hc.tcolor.ENDC)

   tname = case[c]
   if 'MAC'  in case[c]: tname = 'MAC'
   if 'ERA5' in case[c]: tname = 'ERA5'
   case_obj = he.Case( name=tname, time_freq='daily' )
   if 'lev' not in vars() : lev = np.array([-1])

   comp = 'eam'
   if case[c]=='EAR5': comp = None
   if case[c]=='E3SM.RGMA.ne30pg2_r05_oECv3.FC5AV1C-L.00': comp = 'cam'

   use_remap = False
   remap_str=f'remap_ne30pg2'
   if case[c]=='ERA5-FV': use_remap = False
   if case[c]=='ERA5-PG': use_remap = True
   if case[c]=='MAC-FV' : use_remap = False
   if case[c]=='MAC-PG' : use_remap = True

   #----------------------------------------------------------------------------
   #----------------------------------------------------------------------------
   for v in range(num_var) :
      case_tmp_file = f'{tmp_file}.daily.{case[c]}.{var[v]}.nc'
      
      print('    var: '+hc.tcolor.GREEN+f'{var[v]:10}'+hc.tcolor.ENDC+'    '+case_tmp_file)
      
      cnt_ds = xr.open_dataset( case_tmp_file )

      case_obj = he.Case( name=case[c] )
      lat = case_obj.load_data('lat',htype='h1',num_files=1)
      lon = case_obj.load_data('lon',htype='h1',num_files=1)

      mask = xr.DataArray( np.ones(len(cnt_ds['ncol']),dtype=bool), coords=[cnt_ds['ncol']] )
      mask = mask & (lat>=lat1) & (lat<=lat2)
      mask = mask & (lon>=lon1) & (lon<=lon2)
      cnt_ds = cnt_ds['cnt'].where(mask,drop=True).to_dataset()

      #-------------------------------------------------------------------------
      # print the time length
      #-------------------------------------------------------------------------
      # msg = '    num_time: '
      # msg += str(cnt_ds['num_time'].values)
      # msg += ' days'
      # print(msg)
      #-------------------------------------------------------------------------
      #-------------------------------------------------------------------------
      
      # cnt_ds['cnt'] = cnt_ds['cnt'] / cnt_ds['num_valid'] 

      # print(); print(cnt_ds); print()
      ### replace NaNs with -1
      # cnt_ds['cnt'] = np.where( cnt_ds['cnt'].values==np.nan, np.zeros(cnt_ds['cnt'].shape), cnt_ds['cnt'].values )
      # cnt_ds = cnt_ds.where( cnt_ds==np.nan, 0. )

      # if print_stats: hc.print_stat(cnt_ds['cnt'],name='final count dataset',stat='naxs',indent='    ')

      cnt_ds_list.append(cnt_ds)

#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------
list_cnt = 0
# plot = [None]*(num_var*num_case*num_set)
plot = [None]*(num_var*num_case)

for c in range(num_case):
   
   # tname = case[c]
   # if 'MAC'  in case[c]: tname = 'MAC'
   # if 'ERA5' in case[c]: tname = 'ERA5'
   # case_obj = he.Case( name=tname, time_freq='daily' )

   # scrip_file_path = 'scrip_files/ne30pg2_scrip.nc'
   # if case[c]=='MAC-FV' or case[c]=='ERA5': scrip_file_path = 'scrip_files/cmip6_180x360_scrip.20210722.nc'

   # scripfile = xr.open_dataset(scrip_file_path)

   # if ocean_only: 
   #    if case[c]=='MAC-FV' : 
   #       print('ocean_only not implemented for MAC-FV!')
   #    else:
   #       landfrac_ds = xr.open_dataset('data/E3SM_landfrac_ne30pg2.nc').isel(time=0)
   #       mask = xr.DataArray( np.ones(len(cnt_ds['ncol']),dtype=bool), coords=[cnt_ds['ncol']] )
   #       mask = mask & (landfrac_ds['LANDFRAC'].values<0.5)
   #       scripfile = scripfile.where(mask.rename({'ncol':'grid_size'}),drop=True)
   

   for v in range(num_var) :

      tres = copy.deepcopy(res)
   
      cnt_ds = cnt_ds_list[list_cnt]
      list_cnt += 1

      # for s in range(num_set):

      if num_var>1:
         ip = s*num_var+v if var_x_case else v*num_set+s
      else:
         # ip = s*num_case+c if var_x_case else c*num_set+s
         ip = c

      tmp_data = cnt_ds['cnt']#.isel(set=s)

      # exit(tmp_data)

      # if ocean_only: tmp_data = tmp_data.where(mask,drop=True)

      ### replace NaNs with zero
      # # tmp_data = np.where( tmp_data.values==np.nan, 0., tmp_data.values )
      # tmp_data = np.where( np.isfinite(tmp_data.values), tmp_data.values, 0. )

      # exit(tmp_data)

      tmp_data = tmp_data.isel(set=1).sum(dim='ncol')
      # tmp_data = tmp_data.isel(set=1).transpose()  # plot count of partial checkerboard each column 

      # Make time start at zero
      tmp_time = ( cnt_ds['time'] - cnt_ds['time'][0] ).astype('float') / 86400e9

      plot[ip] = ngl.xy(wks, tmp_time.values, tmp_data.values, tres) 


      # set_str = set_labels[s]
      hs.set_subtitles(wks, plot[ip], name[c], '', var[v], font_height=0.015)


#---------------------------------------------------------------------------------------------------
# Finalize plot
#---------------------------------------------------------------------------------------------------
# hs.set_plot_labels(wks, plot, font_height=0.01, justify='left')

layout = [len(plot),1]

# if num_var>1:
#    layout = [num_set,num_var] if var_x_case else [num_var,num_set]
# else:
#    layout = [num_set,num_case] if var_x_case else [num_case,num_set]

# if num_case==1 and num_set>3: 
#    num_plot_col = 6
#    layout = [int(np.ceil(len(plot)/float(num_plot_col))),num_plot_col]

pnl_res = hs.setres_panel()
# pnl_res.nglPanelYWhiteSpacePercent = 10
# pnl_res.nglPanelFigureStrings            = list(string.ascii_lowercase)
# pnl_res.nglPanelFigureStringsJust        = "TopLeft"
# pnl_res.nglPanelFigureStringsFontHeightF = 0.015   

ngl.panel(wks,plot,layout,pnl_res)
ngl.end()

hc.trim_png(fig_file)
#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------
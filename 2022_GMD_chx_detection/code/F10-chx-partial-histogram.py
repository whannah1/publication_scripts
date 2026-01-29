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
add_case('OBS-PG',    n='OBS ne30pg2'                                        ,c='black',p=0)
add_case('E3SM.CHX.ne30pg2_ne30pg2.FC5AV1C-L.no_dcape.00',n='E3SM (no DCAPE)',c='green',p=0)
add_case('E3SM.CHX.ne30pg2_ne30pg2.FC5AV1C-L.00',         n='E3SM'           ,c='red'  ,p=0)
add_case('E3SM.CHX.ne30pg2_ne30pg2.F-MMF1.00',            n='E3SM-MMF'       ,c='blue' ,p=0)
#-------------------------------------------------------------------------------

var = []
add_var('TGCLDLWP',n='Liq Water Path')
add_var('PRECT',n='Precipitation')


fig_type = 'png'
fig_file = 'figs/F10-chx-partial-histogram'

scratch = '/global/cscratch1/sd/whannah/e3sm_scratch/cori-knl/chx-detection-data'
tmp_file=f'{scratch}/occurence-partial-chx-over-time'

event_file_head = 'data/event_duration.v1'


ocean_only = True

var_x_case           = True    # controls plot panel arrangement
print_stats          = True    #
verbose              = True    #

recalculate = False

lat1,lat2 = -60,60

# lat1,lat2,lon1,lon2 = 0,30,140,180+40

# ilat,ilon,idx = -20,360-115,20
# lat1,lat2 = ilat-idx,ilat+idx
# lon1,lon2 = ilon-idx, ilon+idx

subset_min_length = 4       # threshold for partial checkerboard when combine_mode=1

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

res.tiXAxisString = 'Event Length [days]'
res.tiYAxisString = 'Number of Partial Checkerboard Events'

res.xyXStyle = 'Log'
res.xyYStyle = 'Log'


scrip_file_path = 'scrip_files/ne30pg2_scrip.nc'
scripfile = xr.open_dataset(scrip_file_path)
lat = scripfile['grid_center_lat'].rename({'grid_size':'ncol'})
lon = scripfile['grid_center_lon'].rename({'grid_size':'ncol'})

#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------
@numba.njit()
def count_event_duration(chx_flag,ncol,ntime):
   cnt = 0
   event = False
   event_len = [2]
   event_cnt = [0]
   for n in range(ncol):
      for t in range(1,ntime):
         beg,end = False,False
         if chx_flag[t,n]==1 and chx_flag[t-1,n]==0: beg,event = True,True
         if chx_flag[t,n]==0 and chx_flag[t-1,n]==1: end,event = True,False
         if end:
            # if cnt>1:
               if cnt not in event_len: 
                  event_len.append(cnt)
                  event_cnt.append(0)
               pos = event_len.index(cnt)
               event_cnt[pos] = event_cnt[pos]+1
         else:
            if beg: cnt = 0
            if event: cnt = cnt+1
   return (event_cnt,event_len)
#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------
def get_alt_name(case,var):
   tname = case
   if 'OBS' in case and var=='TGCLDLWP': tname = 'MAC'
   if 'OBS' in case and var=='PRECT'   : tname = 'GPM'
   alt_name = case
   if 'OBS' in case: alt_name = alt_name.replace('OBS',tname)
   return alt_name

def get_data_file_name(case,var):
   case_tmp = get_alt_name(case,var)
   data_file = f'{tmp_file}.daily.{case_tmp}.{var}.sml_{subset_min_length}.nc'
   return data_file

def get_event_file_name(case,var):
   case_tmp = get_alt_name(case,var)
   event_file = f'{event_file_head}.{case_tmp}.{var}'
   if 'lat1' in globals(): event_file+= f'.lat1_{lat1}.lat2_{lat2}'
   if 'lon1' in globals(): event_file+= f'.lon1_{lat1}.lon2_{lat2}'
   event_file+= f'.sml_{subset_min_length}.nc'
   return event_file

def print_file_status(tfile):
   msg = tfile
   msg = msg.replace('/global/cscratch1/sd/whannah/e3sm_scratch/cori-knl/chx-detection-data','<...>')
   msg = f'  {msg:140}  '
   if os.path.exists(tfile): 
      print(msg+hc.tcolor.GREEN+'OK'+hc.tcolor.ENDC)
   else:
      print(msg+hc.tcolor.RED+'MISSING!'+hc.tcolor.ENDC)
#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------
# check if any files are missing
for i in range(2):
   print()
   for v in range(num_var):
      for c in range(num_case):
         if i==0: print_file_status( get_data_file_name(case[c],var[v]) )
         if i==1: print_file_status( get_event_file_name(case[c],var[v]) )
print()

#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------
list_cnt = 0
# plot = [None]*(num_var*num_case*num_set)
plot = [None]*(num_var)
# plot = [None]

for v in range(num_var):
   print('\n    var: '+hc.tcolor.GREEN+f'{var[v]:10}'+hc.tcolor.ENDC)

   # if v==0 : subset_min_length = 4
   # if v==1 : subset_min_length = 6

   event_cnt_list = []
   event_len_list = []

   for c in range(num_case):

      chx_data_file = get_data_file_name(case[c],var[v])
      event_file    = get_event_file_name(case[c],var[v])
      
      # print('\n  case: '+hc.tcolor.CYAN+case[c]+hc.tcolor.ENDC+'    '+chx_data_file)
      print('\n      case: '+hc.tcolor.CYAN+case[c]+hc.tcolor.ENDC)

      #-------------------------------------------------------------------------
      # Calculate event duration
      #------------------------------------------------------------------------- 
      # event_file = f'{event_file_head}.{case_tmp}.{var[v]}.lat1_{lat1}.lat2_{lat2}.sml_{subset_min_length}.nc'
      
      if recalculate:
         #----------------------------------------------------------------------
         # Load the teporally resolved chx identification data
         #----------------------------------------------------------------------
         cnt_ds = xr.open_dataset( chx_data_file )

         #----------------------------------------------------------------------
         # Ocean mask
         #----------------------------------------------------------------------
         if ocean_only:    
            landfrac_ds = xr.open_dataset('data/E3SM_landfrac_ne30pg2.nc').isel(time=0)
            landfrac = landfrac_ds['LANDFRAC']
            ocn_mask = xr.DataArray( np.ones(len(cnt_ds['ncol']),dtype=bool), coords=[cnt_ds['ncol']] )
            ocn_mask = ocn_mask & (landfrac.values<0.5)
            cnt_ds['cnt'] = cnt_ds['cnt'].where(ocn_mask,drop=True)
         #----------------------------------------------------------------------
         # Apply regional mask
         #----------------------------------------------------------------------
         mask = xr.DataArray( np.ones(len(cnt_ds['ncol']),dtype=bool), coords=[cnt_ds['ncol']] )
         if 'lat1' in locals(): mask = mask & (lat>=lat1) & (lat<=lat2)
         if 'lon1' in locals(): mask = mask & (lon>=lon1) & (lon<=lon2)
         cnt_ds = cnt_ds['cnt'].where(mask,drop=True).to_dataset()
         #----------------------------------------------------------------------
         # misc data prep
         #----------------------------------------------------------------------
         # cnt_ds = cnt_ds_list[list_cnt]
         list_cnt += 1

         ### just use the partial chx flag
         tmp_data = cnt_ds['cnt'][1,:,:]

         ### subset the time dimension for debugging
         # tmp_data.isel(time=slice(0,60))

         ### Make time start at zero
         tmp_time = ( cnt_ds['time'] - cnt_ds['time'][0] ).astype('float') / 86400e9
         
         #----------------------------------------------------------------------
         # partial checkerboard event detection
         #----------------------------------------------------------------------
         (event_cnt,event_len) = count_event_duration(tmp_data[:,:].values,len(cnt_ds['ncol']),len(tmp_time))
         
         # for c,e in enumerate(event_cnt): print(f'  {event_len[c]}   {event_cnt[c]}')

         event_cnt = np.array(event_cnt)
         event_len = np.array(event_len)

         #----------------------------------------------------------------------
         # sort results of event detection
         #----------------------------------------------------------------------
         sort_ind = np.argsort(event_len)
         event_cnt = event_cnt[sort_ind]
         event_len = event_len[sort_ind]

         event_ds = xr.Dataset()
         event_ds['event_cnt'] = ( ('event_len'), event_cnt )
         event_ds['event_len'] = ( ('event_len'), event_len )

         #----------------------------------------------------------------------
         # write to file
         #----------------------------------------------------------------------
         print(f'      writing to file: {event_file}')
         event_ds.to_netcdf(path=event_file,mode='w')

      else:
         event_ds = xr.open_dataset( event_file )

         event_cnt = event_ds['event_cnt'].values
         event_len = event_ds['event_len'].values

         # for e in range(len(event_cnt)): 
         #    print(' '*8+f'{event_len[e]}  {event_cnt[e]}')
         # exit()

      #-------------------------------------------------------------------------
      # expand the data to add zeros
      #-------------------------------------------------------------------------
      new_event_len = []
      new_event_cnt = []
      for c in range(1,max(event_len)+1):
         idx = -1 ; msg = ' '*8+f'{c:5d}  '
         
         if c in event_len:
            idx = np.argwhere(event_len==c)
            if np.all(idx==[None]): idx = -1
            if idx!=-1:
               if len(idx.shape)>0: idx = idx[0]
               new_cnt = event_cnt[idx][0]
         if idx==-1: new_cnt = 0
         # if idx==-1: new_cnt = np.nan
         # print(msg+f'{new_cnt}')
         new_event_len.append(c)
         new_event_cnt.append(new_cnt)
      
      event_len = np.array(new_event_len)
      event_cnt = np.array(new_event_cnt)

      # print()
      # for e in range(len(event_cnt)):
      # # for e in range(min(30,len(event_cnt))): 
      #    print(' '*8+f'{event_len[e]}  {event_cnt[e]}')
      # exit()
      #-------------------------------------------------------------------------
      #-------------------------------------------------------------------------
      event_cnt_list.append(event_cnt)
      event_len_list.append(event_len)
   #----------------------------------------------------------------------------
   #----------------------------------------------------------------------------
   cnt_min,cnt_max = 1e10,0
   len_min,len_max = 1e10,0
   for c in range(num_case):
      for event_cnt in event_cnt_list:
         cnt_min = np.min([ cnt_min, np.min(event_cnt) ])
         cnt_max = np.max([ cnt_max, np.max(event_cnt) ])
         len_min = np.min([ len_min, np.min(event_len) ])
         len_max = np.max([ len_max, np.max(event_len) ])

   # span = np.absolute( data_max - data_min )
   # res.trXMinF,res.trXMaxF = len_min, len_max
   res.trYMinF,res.trYMaxF = cnt_min-cnt_min*0.05, cnt_max

   if res.trYMinF==0: res.trYMinF = 1

   res.trXMinF = 1
   res.trXMaxF = 100
   # if len_max>10 and len_max<100: res.trXMaxF = 100


   #----------------------------------------------------------------------------
   #----------------------------------------------------------------------------
   ip = v

   for c in range(num_case):

      tres = copy.deepcopy(res)
      tres.xyLineColor = clr[c]
      # tres.xyDashPattern = v

      tmp_cnt = event_cnt_list[c]
      tmp_cnt = np.ma.masked_array(tmp_cnt, mask=tmp_cnt==0)

      tplot = ngl.xy(wks, event_len_list[c], tmp_cnt, tres) 
      if c==0:
         plot[ip] = tplot
      else:
         ngl.overlay( plot[ip], tplot )

   var_name = var[v]
   if var_name=='TGCLDLWP': var_name = 'Liq Water Path'
   if var_name=='PRECT':    var_name = 'Precipitation'
   hs.set_subtitles(wks, plot[ip], var_name, '', '', font_height=0.015)

#---------------------------------------------------------------------------------------------------
# Add legend
#---------------------------------------------------------------------------------------------------
lgres = ngl.Resources()
lgres.vpWidthF, lgres.vpHeightF  = 0.09, 0.11
lgres.lgLabelFontHeightF = 0.01
lgres.lgLineThicknessF   = 20
lgres.lgMonoLineColor,lgres.lgLineColors  = False, clr
lgres.lgMonoDashIndex,lgres.lgDashIndexes = True, 0
lgres.lgLabelJust    = 'CenterLeft'

lname = [f' {n}' for n in name]
for (n,l) in enumerate(lname):
   if 'OBS' in lname[n]: 
      lname[n] = lname[n].replace('OBS','MAC')
xpos = 0.25
ypos = 0.7
pid = ngl.legend_ndc(wks, len(name), lname, xpos, ypos, lgres)

lname = [f' {n}' for n in name]
for (n,l) in enumerate(lname):
   if 'OBS' in lname[n]: 
      lname[n] = lname[n].replace('OBS','IMERG')
xpos = 0.75
ypos = 0.7
pid = ngl.legend_ndc(wks, len(name), lname, xpos, ypos, lgres)
#---------------------------------------------------------------------------------------------------
# Finalize plot
#---------------------------------------------------------------------------------------------------
layout = [1,num_var]
# layout = [1,len(plot)]

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
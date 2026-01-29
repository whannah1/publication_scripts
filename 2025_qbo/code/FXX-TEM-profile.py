import os, ngl, copy, xarray as xr, numpy as np, glob
import hapy_common as hc, hapy_E3SM as he, hapy_setres as hs
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
var,lev_list,vclr,vdsh = [],[],[],[]
def add_var(var_name,lev=-1,c='black',d=0): 
   var.append(var_name); lev_list.append(lev)
   vclr.append(c); vdsh.append(d)
##------------------------------------------------------------------------------
cscratch = '/global/cscratch1/sd/whannah/e3sm_scratch/cori-knl'
gscratch = '/global/cfs/cdirs/m3312/whannah/2022-QBO-TEST'
pscratch = '/pscratch/sd/w/whannah/e3sm_scratch/pm-cpu'
##------------------------------------------------------------------------------
add_case('E3SM.QBO-TEST.F2010.ne30pg2.L72.01',       n='E3SM control',     c='red'  ,p=gscratch,s='run')
add_case('E3SM.QBO-TEST.F2010.ne30pg2.L72-nsu40.01', n='E3SM L72 smoothed',c='green',p=gscratch,s='run')
add_case('E3SM.QBO-TEST.F2010.ne30pg2.L80-rsu40.01', n='E3SM L80 refined', c='blue' ,p=gscratch,s='run')

#add_case('E3SM.QBO-TEST-02.F2010.ne30pg2.L72.GW-MOD-0', n='E3SM L72',      d=1,c='black', p=pscratch,s='run')
#add_case('E3SM.QBO-TEST-02.F2010.ne30pg2.L72-nsu40',    n='E3SM L72-nsu40',d=1,c='red',   p=pscratch,s='run')
#add_case('E3SM.QBO-TEST-02.F2010.ne30pg2.L72-rscl',     n='E3SM L72-rscl', d=1,c='purple',p=pscratch,s='run')
#add_case('E3SM.QBO-TEST-02.F2010.ne30pg2.L72-rlim',     n='E3SM L72-rlim', d=1,c='pink',  p=pscratch,s='run')
#-------------------------------------------------------------------------------

lev = np.array([1,2,3,4,5,6,7,8,9,10,15,20,25,30,40,50,100,150,200])

# add_var('u',        c='black')
add_var('utendepfd',c='blue') # u-wind tendency due to TEM EP flux divergence
add_var('BUTGWSPEC',c='magenta')


num_plot_col = 2


fig_file,fig_type = 'figs/FXX-TEM-profile','png'

lat1,lat2 = -15,15

# htype,first_file,num_files = 'h2',0,10#365*1
remap_str,search_str = 'remap_90x180','h0.tem.'; first_file,num_files = 0,12*20
# remap_str,search_str = 'remap_90x180','h2.tem.'; first_file,num_files = 0,365*20

recalculate = True

plot_diff = False

omit_bot,bot_k     = False,-2
omit_top,top_k     = False,30

print_stats        = False
print_profile      = False

tmp_file_head = 'data/TEM.profile'


p_min = 10
p_max = 100

#---------------------------------------------------------------------------------------------------
# Set up plot resources
#---------------------------------------------------------------------------------------------------
num_var,num_case = len(var),len(case)

if 'lev' not in vars(): lev = np.array([0])

wks = ngl.open_wks(fig_type,fig_file)
plot = [None]*(num_case*num_var)

wks = ngl.open_wks(fig_type,fig_file)
plot = [None]*(num_var)
res = hs.res_xy()
# res.vpWidthF = 0.4
# res.xyMarkLineMode = "MarkLines"
res.xyMarkerSizeF = 0.008
res.xyMarker = 16
res.xyLineThicknessF = 8
# res.tmYLLabelFontHeightF         = 0.008
# res.tmXBLabelFontHeightF         = 0.008

res.tmXBAutoPrecision = False
res.tmXBPrecision = 2

res.trYReverse = True
res.xyYStyle = 'Log'

#---------------------------------------------------------------------------------------------------
def get_comp(case):
   comp = 'eam'
   return comp
#---------------------------------------------------------------------------------------------------
# @numba.jit(nopython=True)
# def calculate_equiangular_area(lon,lat,lon_bnds,lat_bnds):
#    re = 6.37122e06  # radius of earth
#    nlat,nlon = len(lat),len(lon)
#    area = np.empty((nlat,nlon),np.float64)
#    for j in range(nlat):
#       for i in range(nlon):
#          dlon = np.absolute( lon_bnds[j,1] - lon_bnds[j,0] )
#          dlat = np.absolute( lat_bnds[j,1] - lat_bnds[j,0] )
#          dx = re*dlon*np.pi/180.
#          dy = re*dlat*np.pi/180.
#          area[j,i] = dx*dy
#    return area
#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------
data_list_list,lev_list_list = [],[]
for v in range(num_var):
   hc.printline()
   print(hc.tcolor.GREEN+'  var: '+var[v]+hc.tcolor.ENDC)
   data_list,lev_list = [],[]
   for c in range(num_case):
      print(hc.tcolor.CYAN+'    case: '+case[c]+hc.tcolor.ENDC)
      #-------------------------------------------------------------------------
      #-------------------------------------------------------------------------
      tmp_file = f'{tmp_file_head}.{case[c]}.{var[v]}'
      if 'lat1' in locals(): tmp_file = tmp_file + f'.lat1_{lat1}.lat2_{lat2}'
      if 'lon1' in locals(): tmp_file = tmp_file + f'.lon1_{lon1}.lon2_{lon2}'
      tmp_file = f'{tmp_file}.nc'
      print(f'      tmp_file: {tmp_file}')

      if recalculate :
         file_path = f'{gscratch}/{case[c]}/data_{remap_str}_tem/*.eam.{search_str}*'
         file_list = sorted(glob.glob(file_path))

         if 'first_file' in locals(): file_list = file_list[first_file:]
         if 'num_files'  in locals(): file_list = file_list[:num_files]

         # for ff in file_list: print(ff)

         # ds = xr.open_mfdataset(file_list)
         #----------------------------------------------------------------------
         # manually combine files due to issues with open_mfdataset on Perlmutter
         ds = xr.open_dataset(file_list[0])
         for f in file_list[1:]:
           ds_tmp = xr.open_dataset(f)
           ds = xr.concat([ds, ds_tmp], dim='time')
         #----------------------------------------------------------------------
         time = ds['time']

         data = ds[var[v]].sel(lat=slice(lat1,lat2),plev=slice(200e2,0))

         hc.print_stat(data,name=f'{var[v]}         ',indent=' '*6,compact=True)
         data = data.mean(dim='lat')
         data = data.mean(dim='time')
         hc.print_stat(data,name=f'{var[v]} variance',indent=' '*6,compact=True)

         #-------------------------------------------------------------------------
         # options for omitting top or bottom levels
         if omit_bot: data = data[:bot_k]
         if omit_top: data = data[top_k:]

         #----------------------------------------------------------------------
         # write to temporary file
         #----------------------------------------------------------------------
         tmp_ds = xr.Dataset()
         tmp_ds['data'] = data
         tmp_ds['time'] = time
         if 'lat1' in vars(): tmp_ds['lat1']=lat1; tmp_ds['lat2']=lat2
         if 'lon1' in vars(): tmp_ds['lon1']=lon1; tmp_ds['lon2']=lon2
         print(f'      writing to file: {tmp_file}')
         tmp_ds.to_netcdf(path=tmp_file,mode='w')

      else:
         tmp_ds = xr.open_dataset( tmp_file )
         data = tmp_ds['data']
         time = tmp_ds['time']
      
      #-------------------------------------------------------------------------
      # append final data to list
      data_list.append( data.values )
      lev_list.append( data['plev'].values/1e2 )

   data_list_list.append(data_list)
   lev_list_list.append(lev_list)

#-------------------------------------------------------------------------------
# Create plot - overlay all cases for each var
#-------------------------------------------------------------------------------
for v in range(num_var):
   data_list = data_list_list[v]
   lev_list = lev_list_list[v]
   
   tres = copy.deepcopy(res)

   # ip = v*num_case+c
   ip = c*num_var+v
   
   baseline = data_list[0]
   if plot_diff:
      for c in range(num_case): 
         data_list[c] = data_list[c] - baseline

   for c in range(num_case): 
      for i in range(len(data_list[c])):
         if lev_list[c][i]>p_max: data_list[c][i] = np.nan
         if lev_list[c][i]<p_min: data_list[c][i] = np.nan
   tres.trYMinF = p_min
   tres.trYMaxF = p_max
   
   data_min = np.min([np.nanmin(d) for d in data_list])
   data_max = np.max([np.nanmax(d) for d in data_list])
   tres.trXMinF = data_min
   tres.trXMaxF = data_max
   ip = v


   for c in range(num_case):
      tres.xyLineColor   = clr[c]
      tres.xyMarkerColor = clr[c]
      tres.xyDashPattern = dsh[c]

      tplot = ngl.xy(wks, np.ma.masked_invalid(data_list[c]), lev_list[c], tres)

      
      if (c==1 and plot_diff) or (c==0 and not plot_diff) :
         plot[ip] = tplot
      elif (plot_diff and c>0) or not plot_diff:
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

   if plot_diff: var_str += ' (diff)'

   hs.set_subtitles(wks, plot[ip], var_str, '', reg_str, font_height=0.01)

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

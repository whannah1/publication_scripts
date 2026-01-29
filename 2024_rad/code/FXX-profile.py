import os, ngl, copy, xarray as xr, numpy as np
import hapy_common as hc, hapy_E3SM as he, hapy_setres as hs
#-------------------------------------------------------------------------------
case      = []
case_name = []
case_dir  = []
case_sub  = []
case_clr  = []
case_dsh  = []
case_mrk  = []
def add_case(case_in,n=None,p=None,s=None,d=0,c='black',m=0):
   global case,case_name,case_dir,case_sub,case_clr,case_dsh,case_mrk
   case     .append(case_in)
   case_name.append(n)
   case_dir .append(p)
   case_sub .append(s)
   case_dsh .append(d)
   case_clr .append(c)
   case_mrk .append(m)
#-------------------------------------------------------------------------------
group_name      = []
group_case      = []
group_case_name = []
group_case_dir  = []
group_case_sub  = []
group_case_clr  = []
group_case_dsh  = []
group_case_mrk  = []
# def add_group(name_in,case,case_name,case_dir,case_sub,case_clr,case_dsh,case_mrk):
def add_group(name_in):
   global case,case_name,case_dir,case_sub,case_clr,case_dsh,case_mrk
   group_name     .append(name_in)
   group_case     .append(case)
   group_case_name.append(case_name)
   group_case_dir .append(case_dir)
   group_case_sub .append(case_sub)
   group_case_clr .append(case_clr)
   group_case_dsh .append(case_dsh)
   group_case_mrk .append(case_mrk)
   case      = []
   case_name = []
   case_dir  = []
   case_sub  = []
   case_clr  = []
   case_dsh  = []
   case_mrk  = []
#-------------------------------------------------------------------------------
var,lev_list = [],[]
def add_var(var_name,lev=-1): var.append(var_name); lev_list.append(lev)
##------------------------------------------------------------------------------
cscratch = '/global/cscratch1/sd/whannah/e3sm_scratch/cori-knl'
gscratch = '/global/cfs/cdirs/m3312/whannah/2022-QBO-TEST'
pscratch = '/pscratch/sd/w/whannah/e3sm_scratch/pm-gpu'

tmp_path,tmp_sub = pscratch,'archive/atm/hist'

add_case('E3SM.2023-RAD-SENS-00.GNUGPU.ne30pg2_oECv3.F2010-MMF1.NXY_64x1.RNX_1.MCICA_OFF',   n='nx/rx= 64/ 1 mcica OFF',  c='red',      d=0,p=tmp_path,s=tmp_sub)
add_case('E3SM.2023-RAD-SENS-00.GNUGPU.ne30pg2_oECv3.F2010-MMF1.NXY_64x1.RNX_2.MCICA_OFF',   n='nx/rx= 64/ 2 mcica OFF',  c='orange',   d=0,p=tmp_path,s=tmp_sub)
add_case('E3SM.2023-RAD-SENS-00.GNUGPU.ne30pg2_oECv3.F2010-MMF1.NXY_64x1.RNX_4.MCICA_OFF',   n='nx/rx= 64/ 4 mcica OFF',  c='gold',     d=0,p=tmp_path,s=tmp_sub)
add_case('E3SM.2023-RAD-SENS-00.GNUGPU.ne30pg2_oECv3.F2010-MMF1.NXY_64x1.RNX_8.MCICA_OFF',   n='nx/rx= 64/ 8 mcica OFF',  c='palegreen',d=0,p=tmp_path,s=tmp_sub)
add_case('E3SM.2023-RAD-SENS-00.GNUGPU.ne30pg2_oECv3.F2010-MMF1.NXY_64x1.RNX_16.MCICA_OFF',  n='nx/rx= 64/16 mcica OFF',  c='green',    d=0,p=tmp_path,s=tmp_sub)
add_case('E3SM.2023-RAD-SENS-00.GNUGPU.ne30pg2_oECv3.F2010-MMF1.NXY_64x1.RNX_32.MCICA_OFF',  n='nx/rx= 64/32 mcica OFF',  c='blue',     d=0,p=tmp_path,s=tmp_sub)
add_case('E3SM.2023-RAD-SENS-00.GNUGPU.ne30pg2_oECv3.F2010-MMF1.NXY_64x1.RNX_64.MCICA_OFF',  n='nx/rx= 64/64 mcica OFF',  c='purple',   d=0,p=tmp_path,s=tmp_sub)
add_group('nx=64 MCICA-OFF')

add_case('E3SM.2023-RAD-SENS-00.GNUGPU.ne30pg2_oECv3.F2010-MMF1.NXY_64x1.RNX_1.MCICA_ON',    n='nx/rx= 64/ 1 mcica ON ',  c='red',      d=0,p=tmp_path,s=tmp_sub)
add_case('E3SM.2023-RAD-SENS-00.GNUGPU.ne30pg2_oECv3.F2010-MMF1.NXY_64x1.RNX_2.MCICA_ON',    n='nx/rx= 64/ 2 mcica ON ',  c='orange',   d=0,p=tmp_path,s=tmp_sub)
add_case('E3SM.2023-RAD-SENS-00.GNUGPU.ne30pg2_oECv3.F2010-MMF1.NXY_64x1.RNX_4.MCICA_ON',    n='nx/rx= 64/ 4 mcica ON ',  c='gold',     d=0,p=tmp_path,s=tmp_sub)
add_case('E3SM.2023-RAD-SENS-00.GNUGPU.ne30pg2_oECv3.F2010-MMF1.NXY_64x1.RNX_8.MCICA_ON',    n='nx/rx= 64/ 8 mcica ON ',  c='palegreen',d=0,p=tmp_path,s=tmp_sub)
add_case('E3SM.2023-RAD-SENS-00.GNUGPU.ne30pg2_oECv3.F2010-MMF1.NXY_64x1.RNX_16.MCICA_ON',   n='nx/rx= 64/16 mcica ON ',  c='green',    d=0,p=tmp_path,s=tmp_sub)
add_case('E3SM.2023-RAD-SENS-00.GNUGPU.ne30pg2_oECv3.F2010-MMF1.NXY_64x1.RNX_32.MCICA_ON',   n='nx/rx= 64/32 mcica ON ',  c='blue',     d=0,p=tmp_path,s=tmp_sub)
add_case('E3SM.2023-RAD-SENS-00.GNUGPU.ne30pg2_oECv3.F2010-MMF1.NXY_64x1.RNX_64.MCICA_ON',   n='nx/rx= 64/64 mcica ON ',  c='purple',   d=0,p=tmp_path,s=tmp_sub)
add_group('nx=64 MCICA-ON')

add_case('E3SM.2023-RAD-SENS-00.GNUGPU.ne30pg2_oECv3.F2010-MMF1.NXY_128x1.RNX_1.MCICA_OFF',  n='nx/rx= 128/  1 mcica OFF',c='red',      d=0,p=tmp_path,s=tmp_sub)
add_case('E3SM.2023-RAD-SENS-00.GNUGPU.ne30pg2_oECv3.F2010-MMF1.NXY_128x1.RNX_2.MCICA_OFF',  n='nx/rx= 128/  2 mcica OFF',c='orange',   d=0,p=tmp_path,s=tmp_sub)
add_case('E3SM.2023-RAD-SENS-00.GNUGPU.ne30pg2_oECv3.F2010-MMF1.NXY_128x1.RNX_4.MCICA_OFF',  n='nx/rx= 128/  4 mcica OFF',c='gold',     d=0,p=tmp_path,s=tmp_sub)
add_case('E3SM.2023-RAD-SENS-00.GNUGPU.ne30pg2_oECv3.F2010-MMF1.NXY_128x1.RNX_8.MCICA_OFF',  n='nx/rx= 128/  8 mcica OFF',c='palegreen',d=0,p=tmp_path,s=tmp_sub)
add_case('E3SM.2023-RAD-SENS-00.GNUGPU.ne30pg2_oECv3.F2010-MMF1.NXY_128x1.RNX_16.MCICA_OFF', n='nx/rx= 128/ 16 mcica OFF',c='green',    d=0,p=tmp_path,s=tmp_sub)
add_case('E3SM.2023-RAD-SENS-00.GNUGPU.ne30pg2_oECv3.F2010-MMF1.NXY_128x1.RNX_32.MCICA_OFF', n='nx/rx= 128/ 32 mcica OFF',c='blue',     d=0,p=tmp_path,s=tmp_sub)
add_case('E3SM.2023-RAD-SENS-00.GNUGPU.ne30pg2_oECv3.F2010-MMF1.NXY_128x1.RNX_64.MCICA_OFF', n='nx/rx= 128/ 64 mcica OFF',c='purple',   d=0,p=tmp_path,s=tmp_sub)
add_case('E3SM.2023-RAD-SENS-00.GNUGPU.ne30pg2_oECv3.F2010-MMF1.NXY_128x1.RNX_128.MCICA_OFF',n='nx/rx= 128/128 mcica OFF',c='pink',     d=0,p=tmp_path,s=tmp_sub)
add_group('nx=128 MCICA-OFF')

add_case('E3SM.2023-RAD-SENS-00.GNUGPU.ne30pg2_oECv3.F2010-MMF1.NXY_128x1.RNX_1.MCICA_ON',   n='nx/rx= 128/  1 mcica ON ',c='red',      d=0,p=tmp_path,s=tmp_sub)
add_case('E3SM.2023-RAD-SENS-00.GNUGPU.ne30pg2_oECv3.F2010-MMF1.NXY_128x1.RNX_2.MCICA_ON',   n='nx/rx= 128/  2 mcica ON ',c='orange',   d=0,p=tmp_path,s=tmp_sub)
add_case('E3SM.2023-RAD-SENS-00.GNUGPU.ne30pg2_oECv3.F2010-MMF1.NXY_128x1.RNX_4.MCICA_ON',   n='nx/rx= 128/  4 mcica ON ',c='gold',     d=0,p=tmp_path,s=tmp_sub)
add_case('E3SM.2023-RAD-SENS-00.GNUGPU.ne30pg2_oECv3.F2010-MMF1.NXY_128x1.RNX_8.MCICA_ON',   n='nx/rx= 128/  8 mcica ON ',c='palegreen',d=0,p=tmp_path,s=tmp_sub)
add_case('E3SM.2023-RAD-SENS-00.GNUGPU.ne30pg2_oECv3.F2010-MMF1.NXY_128x1.RNX_16.MCICA_ON',  n='nx/rx= 128/ 16 mcica ON ',c='green',    d=0,p=tmp_path,s=tmp_sub)
add_case('E3SM.2023-RAD-SENS-00.GNUGPU.ne30pg2_oECv3.F2010-MMF1.NXY_128x1.RNX_32.MCICA_ON',  n='nx/rx= 128/ 32 mcica ON ',c='blue',     d=0,p=tmp_path,s=tmp_sub)
add_case('E3SM.2023-RAD-SENS-00.GNUGPU.ne30pg2_oECv3.F2010-MMF1.NXY_128x1.RNX_64.MCICA_ON',  n='nx/rx= 128/ 64 mcica ON ',c='purple',   d=0,p=tmp_path,s=tmp_sub)
add_case('E3SM.2023-RAD-SENS-00.GNUGPU.ne30pg2_oECv3.F2010-MMF1.NXY_128x1.RNX_128.MCICA_ON', n='nx/rx= 128/128 mcica ON ',c='pink',     d=0,p=tmp_path,s=tmp_sub)
add_group('nx=128 MCICA-ON')

#-------------------------------------------------------------------------------

# lev = np.array([30,50,75,100,125,150,200,250,300,350,400,450,500,550,600,650,700,750,800,825,850,875,900,925,950,975,1000])
# lev = np.array([50,100,150,200,300,400,500,600,700,750,800,850,875,900,925,975])
# lev = np.array([5,8,10,15,20,25,30,40,50,100,150,200])
# plev = np.array([0.1,0.5,1,2,3,4,5,6,7,8,9,10,15,20,25,30,40,50,100,150,200,250,300,400,500,600,700,800,850,900])
# plev = np.array([0.1,0.5,1,2,3,4,5,6,7,8,9,10,15,20,25,30,40,50,100,150,200])
# plev = np.array([50,60,70,80,90,100,110,120,130,140,150,160,170,180,190,200,220,240,260,280,300])
plev = np.arange(50,400+5,5)

# var = []
add_var('T')
add_var('Q')
# add_var('RH')
# add_var('U')
# add_var('V')
# add_var('OMEGA')
# add_var('THETA')
# add_var('CLOUD')
# add_var('CLDLIQ')
# add_var('CLDICE')
# add_var('QRL')
# add_var('QRS')

num_plot_col = len(var)
# num_plot_col = len(group_name)
# num_plot_col = 2


fig_file,fig_type = 'figs/FXX-profile','png'


# lat1,lat2 = -10,10
lat1,lat2 = -30,30
# lat1,lat2 = -40,40

htype,first_file,num_files = 'h0',0,12

recalculate = False

plot_diff = False

# Y-axis options
use_height_coord = False
omit_bot,bot_k   = False,-2
omit_top,top_k   = False,30
print_stats      = False
print_profile    = False

tmp_file_head = 'data/profile'

#---------------------------------------------------------------------------------------------------
# Set up plot resources
#---------------------------------------------------------------------------------------------------
num_var = len(var)
num_group = len(group_name)

if 'diff_mode' not in locals(): diff_mode = 0

if 'lev' not in vars(): lev = np.array([0])

wks = ngl.open_wks(fig_type,fig_file)
plot = [None]*(num_var*num_group)
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
   return comp
#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------
data_list_list_list = []
lev_list_list_list = []
for g in range(num_group):
   case      = group_case[g]
   case_name = group_case_name[g]
   case_dir  = group_case_dir[g]
   case_sub  = group_case_sub[g]
   case_clr  = group_case_clr[g]
   case_dsh  = group_case_dsh[g]
   case_mrk  = group_case_mrk[g]
   num_case = len(case)
   #------------------------------------------------------------------------------------------------
   data_list_list,lev_list_list = [],[]
   for v in range(num_var):
      hc.printline()
      print(hc.tcolor.GREEN+'  var: '+var[v]+hc.tcolor.ENDC)
      data_list,lev_list = [],[]
      for c in range(num_case):
         print(hc.tcolor.CYAN+'    case: '+case[c]+hc.tcolor.ENDC)

         data_dir_tmp,data_sub_tmp = None, None
         # if use_remap: data_sub_tmp = f'data_remap_{remap_grid}/'
         if case_dir[c] is not None: data_dir_tmp = case_dir[c]
         if case_sub[c] is not None: data_sub_tmp = case_sub[c]
         
         #-------------------------------------------------------------------------
         #-------------------------------------------------------------------------
         tmp_file = f'{tmp_file_head}.{case[c]}.{var[v]}'
         if 'lat1' in vars(): tmp_file = tmp_file + f'.lat1_{lat1}.lat2_{lat2}'
         if 'lon1' in vars(): tmp_file = tmp_file + f'.lon1_{lon1}.lon2_{lon2}'
         tmp_file = f'{tmp_file}.nc'
         print(f'      tmp_file: {tmp_file}')

         if recalculate :

            case_obj = he.Case( name=case[c], atm_comp=get_comp(case[c]), data_dir=data_dir_tmp, data_sub=data_sub_tmp  )

            #-------------------------------------------------------------------------
            # read the data   
            if 'lat1' in vars(): case_obj.lat1 = lat1; case_obj.lat2 = lat2
            if 'lon1' in vars(): case_obj.lon1 = lon1; case_obj.lon2 = lon2
            
            data = case_obj.load_data(var[v],htype=htype,first_file=first_file,num_files=num_files,lev=plev)
            area = case_obj.load_data('area',htype=htype,first_file=first_file,num_files=1,).astype(np.double)   
            Z = case_obj.load_data('Z3',htype=htype,first_file=first_file,num_files=num_files)


            if 'ilev' in data.dims : data = data.rename({'ilev':'lev'})

            #-------------------------------------------------------------------------
            if 'time' in data.dims : 
               time = data.time
               hc.print_time_length(data.time,indent=' '*6,print_span=True, print_length=False)

            data = data.mean(dim='time')
            Z = Z.mean(dim='time')

            #-------------------------------------------------------------------------
            # area weighted spatial average 
            if case[c]=='ERA5':
               data = ( (data*area).sum(dim=['lat','lon']) / area.sum(dim=['lat','lon']) )#.values 
            else:
               Z    = ( (Z   *area).sum(dim='ncol') / area.sum(dim='ncol') )#.values 
               data = ( (data*area).sum(dim='ncol') / area.sum(dim='ncol') )#.values 

            #-------------------------------------------------------------------------
            # options for omitting top or bottom levels
            if omit_bot:
               data = data[:bot_k]
               Z = Z[:bot_k]

            if omit_top:
               data = data[top_k:]
               Z = Z[top_k:]

            #----------------------------------------------------------------------
            # write to temporary file
            #----------------------------------------------------------------------
            tmp_ds = xr.Dataset()
            tmp_ds['data'] = data
            tmp_ds['time'] = time
            if 'lat1' in vars(): tmp_ds['lat1']=lat1; tmp_ds['lat2']=lat2
            if 'lon1' in vars(): tmp_ds['lon1']=lon1; tmp_ds['lon2']=lon2
            tmp_ds['Z']    = Z
            print(f'      writing to file: {tmp_file}')
            tmp_ds.to_netcdf(path=tmp_file,mode='w')

         else:
            tmp_ds = xr.open_dataset( tmp_file )
            data = tmp_ds['data']
            time = tmp_ds['time']
            Z    = tmp_ds['Z']
         
         #-------------------------------------------------------------------------
         # print stuff
         if print_stats:
            hc.print_stat(data,name=f'{var[v]} after averaging',indent='    ',compact=True)
            # hc.print_stat(Z,name=f'Z after averaging',indent='    ')

         if print_profile:
            print()
            for xx in data.values: print(f'    {xx}')
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


   data_list_list_list.append(data_list_list)
   lev_list_list_list .append(lev_list_list)
#---------------------------------------------------------------------------------------------------
# determine data range for all groups and variables
data_min_list = [None]*num_var
data_max_list = [None]*num_var
for g in range(num_group):
   case      = group_case[g]
   case_name = group_case_name[g]
   case_dir  = group_case_dir[g]
   case_sub  = group_case_sub[g]
   case_clr  = group_case_clr[g]
   case_dsh  = group_case_dsh[g]
   case_mrk  = group_case_mrk[g]
   num_case = len(case)
   #------------------------------------------------------------------------------------------------
   for v in range(num_var):
      if plot_diff: 
         baseline = copy.deepcopy(data_list_list_list[g][v][0])
         # baseline = data_list_list_list[g][v][-1]
         for c in range(num_case): 
            data_list_list_list[g][v][c] -= baseline
      data_min_tmp = np.min( [ np.nanmin(d) for d in data_list_list_list[g][v] ] )
      data_max_tmp = np.max( [ np.nanmax(d) for d in data_list_list_list[g][v] ] )
      if g==0:
         data_min_list[v] = data_min_tmp
         data_max_list[v] = data_max_tmp
      else:
         data_min_list[v] = np.min( [ data_min_list[v], data_min_tmp ] )
         data_max_list[v] = np.max( [ data_max_list[v], data_max_tmp ] )
#---------------------------------------------------------------------------------------------------
for g in range(num_group):
   case      = group_case[g]
   case_name = group_case_name[g]
   case_dir  = group_case_dir[g]
   case_sub  = group_case_sub[g]
   case_clr  = group_case_clr[g]
   case_dsh  = group_case_dsh[g]
   case_mrk  = group_case_mrk[g]
   num_case = len(case)
   #------------------------------------------------------------------------------------------------
   data_list_list = data_list_list_list[g]
   lev_list_list = lev_list_list_list[g]
   #-------------------------------------------------------------------------------
   # Create plot - overlay all cases for each var
   #-------------------------------------------------------------------------------
   for v in range(num_var):
      data_list = data_list_list[v]
      lev_list = lev_list_list[v]
      
      tres = copy.deepcopy(res)
      
      # baseline = data_list[-1]
      # if plot_diff: data_list[c] = data_list[c] - baseline
      # if diff_mode==1 :
      #    baseline1 = data_list[0]
      #    baseline2 = data_list[1]
      # if diff_mode==2 :
      #    baseline1 = data_list[0]
      #    baseline2 = data_list[2]
      # if plot_diff:
      #    for c in range(num_case): 
      #       if diff_mode==1 :
      #          if c==0 or c==2 : baseline = baseline1
      #          if c==1 or c==3 : baseline = baseline2
      #       if diff_mode==2 :
      #          if c==0 or c==1 : baseline = baseline1
      #          if c==2 or c==3 : baseline = baseline2
      #       data_list[c] = data_list[c] - baseline
      
      data_min = data_min_list[v]
      data_max = data_max_list[v]

      # data_min = np.min([np.nanmin(d) for d in data_list])
      # data_max = np.max([np.nanmax(d) for d in data_list])
      tres.trXMinF = data_min
      tres.trXMaxF = data_max

      if var[v]=='BUTGWSPEC': tres.tiXAxisString = 'm/s/day'

      # ip = v
      # ip = v*num_group+g
      ip = g*num_var+v

      for c in range(num_case):
         tres.xyLineColor   = case_clr[c]
         tres.xyMarkerColor = case_clr[c]
         tres.xyDashPattern = case_dsh[c]

         tplot = ngl.xy(wks, data_list[c], lev_list[c], tres)

         if diff_mode==0 :
            if (c==1 and plot_diff) or (c==0 and not plot_diff) :
               plot[ip] = tplot
            elif (plot_diff and c>0) or not plot_diff:
               ngl.overlay(plot[ip],tplot)

         # if diff_mode==1 :
         #    if (c==2 and plot_diff) or (c==0 and not plot_diff) :
         #       plot[ip] = tplot
         #    elif (plot_diff and c!=0 and c!=1) or not plot_diff:
         #       ngl.overlay(plot[ip],tplot)

         # if diff_mode==2 :
         #    if (c==1 and plot_diff) or (c==0 and not plot_diff) :
         #       plot[ip] = tplot
         #    elif (plot_diff and c!=0 and c!=2) or not plot_diff:
         #       ngl.overlay(plot[ip],tplot)

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


      hs.set_subtitles(wks, plot[ip], group_name[g], reg_str, var_str, font_height=0.005)


#-------------------------------------------------------------------------------
# Add legend
#-------------------------------------------------------------------------------
# if num_case>1:
#    lgres = ngl.Resources()
#    lgres.vpWidthF           = 0.05
#    lgres.vpHeightF          = 0.08
#    lgres.lgLabelFontHeightF = 0.012
#    lgres.lgMonoDashIndex    = True
#    lgres.lgLineLabelsOn     = False
#    lgres.lgLineThicknessF   = 8
#    lgres.lgLabelJust        = 'CenterLeft'
#    lgres.lgLineColors       = case_clr
#    lgres.lgDashIndexes      = case_dsh

#    lx,ly = 0.5,0.45
#    if num_var==2: lx,ly = 0.3,0.45
#    if num_var==4: lx,ly = 0.05,0.5

#    # pid = ngl.legend_ndc(wks, len(case_name), case_name, lx, ly, lgres)

#---------------------------------------------------------------------------------------------------
# Finalize plot
#---------------------------------------------------------------------------------------------------

layout = [int(np.ceil(len(plot)/float(num_plot_col))),num_plot_col]

textres               =  ngl.Resources()
# textres.txFontHeightF =  0.01
# ngl.text_ndc(wks,f'time step = {ss_t}',0.5,.97,textres)
textres.txFontHeightF =  0.02
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

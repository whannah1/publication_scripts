#---------------------------------------------------------------------------------------------------
import os, ngl, xarray as xr, numpy as np, glob, copy
import hapy_common as hc, hapy_E3SM as he, hapy_setres as hs
np.seterr(divide='ignore', invalid='ignore')
np.errstate(divide='ignore', invalid="ignore")
#---------------------------------------------------------------------------------------------------
case1_list, case2_list, name_list = [],[],[]
case_root, case_sub = [],[]
grp_id_list,grp_str_list = [],[]
res_str_list = []
def add_case(case1,case2,n=None,p=None,s=None,grp_id=None,grp_str='',res_str=''):
   case1_list.append(case1); case2_list.append(case2); name_list.append(n)
   case_root.append(p); case_sub.append(s)
   grp_id_list.append(grp_id), grp_str_list.append(grp_str)
   res_str_list.append(res_str)
   # if res_str== 'ne30': dx = 360*111/( 30*4*2) ; res_str_list.append(f'dx~{int(dx)}km')
   # if res_str== 'ne45': dx = 360*111/( 45*4*2) ; res_str_list.append(f'dx~{int(dx)}km')
   # if res_str== 'ne60': dx = 360*111/( 60*4*2) ; res_str_list.append(f'dx~{int(dx)}km')
   # if res_str== 'ne90': dx = 360*111/( 90*4*2) ; res_str_list.append(f'dx~{int(dx)}km')
   # if res_str=='ne120': dx = 360*111/(120*4*2) ; res_str_list.append(f'dx~{int(dx)}km')
#---------------------------------------------------------------------------------------------------
var_list,var_name,unit_str = [],[],[]
def add_var(var,n=None,u=None):
   var_list.append(var); var_name.append(n if n is not None else var)
   unit_str.append(u)
#---------------------------------------------------------------------------------------------------

# tmp_scratch,tmp_sub = '/gpfs/alpine2/atm146/proj-shared/hannah6/e3sm_scratch/','run'
tmp_scratch,tmp_sub = '/pscratch/sd/w/whannah/2024-AQP-CESS','archive/atm/hist'

# add_var('ECS')
# add_var('RESTOM', n='RESTOM all-sky')
# add_var('RESTOMC',n='RESTOM clear-sky')
# add_var('PRECT',   n='Precipitation')
# add_var('TMQ',     n='Column Water Vapor')
# add_var('FSNT',    n='TOA Net Shortwave Flux')
# add_var('FLNT',    n='TOA Net Longwave Flux')
# add_var('TGCLDLWP',n='Liquid Water Path')
# add_var('TGCLDIWP',n='Ice Water Path')
# # add_var('CLDTOT')
# add_var('CLDLOW',  n='Low Cloud Fraction')
# # add_var('CLDMED',  n='Mid Cloud Fraction')
# add_var('CLDHGH',  n='High Cloud Fraction')
# add_var('LWCF')
# add_var('SWCF')

# add_var('PRECT',   n='Precipitation',           u='mm/day')
# add_var('FSNT',    n='TOA Net Shortwave Flux',  u='W/m2')
# add_var('TGCLDLWP',n='Liquid Water Path',       u='g/m2')
add_var('CLDLOW',  n='Low Cloud Fraction',      u='%')

# add_var('TMQ',     n='Column Water Vapor',      u='kg/m2')
# add_var('FLNT',    n='TOA Net Longwave Flux',   u='W/m2')
# add_var('TGCLDIWP',n='Ice Water Path',          u='g/m2')
add_var('CLDHGH',  n='High Cloud Fraction',     u='%')



# num_plot_col = len(var_list)
num_plot_col = 4

recalculate = False

fig_file,fig_type = 'figs/FXX-glb-mean-vs-res','png'

# htype,first_file,last_file = 'h0',int(100/10),int(500/10)
# tmp_file_head     = 'data/glb-mean-vs-res'

htype,first_file,last_file = 'h0',int(100/10),int(800/10)
tmp_file_head     = 'data/sensitivity-vs-res-alt'

#---------------------------------------------------------------------------------------------------

case1 = 'E3SM.2024-AQP-CESS-00.FAQP-MMF1.GNUGPU.ne30pg2_ne30pg2.NN_32.SSTP_0K'
case2 = 'E3SM.2024-AQP-CESS-00.FAQP-MMF1.GNUGPU.ne30pg2_ne30pg2.NN_32.SSTP_4K'
add_case(case1,case2, n='MMF ne30', p=tmp_scratch,s=tmp_sub,grp_id=0,grp_str='MMF',res_str='ne30')
case1 = 'E3SM.2024-AQP-CESS-00.FAQP-MMF1.GNUGPU.ne45pg2_ne45pg2.NN_64.SSTP_0K'
case2 = 'E3SM.2024-AQP-CESS-00.FAQP-MMF1.GNUGPU.ne45pg2_ne45pg2.NN_64.SSTP_4K'
add_case(case1,case2, n='MMF ne45', p=tmp_scratch,s=tmp_sub,grp_id=0,grp_str='MMF',res_str='ne45')
case1 = 'E3SM.2024-AQP-CESS-00.FAQP-MMF1.GNUGPU.ne60pg2_ne60pg2.NN_128.SSTP_0K'
case2 = 'E3SM.2024-AQP-CESS-00.FAQP-MMF1.GNUGPU.ne60pg2_ne60pg2.NN_128.SSTP_4K'
add_case(case1,case2, n='MMF ne60', p=tmp_scratch,s=tmp_sub,grp_id=0,grp_str='MMF',res_str='ne60')
case1 = 'E3SM.2024-AQP-CESS-00.FAQP-MMF1.GNUGPU.ne90pg2_ne90pg2.NN_256.SSTP_0K'
case2 = 'E3SM.2024-AQP-CESS-00.FAQP-MMF1.GNUGPU.ne90pg2_ne90pg2.NN_256.SSTP_4K'
add_case(case1,case2, n='MMF ne90', p=tmp_scratch,s=tmp_sub,grp_id=0,grp_str='MMF',res_str='ne90')
# case1 = 'E3SM.2024-AQP-CESS-00.FAQP-MMF1.GNUGPU.ne120pg2_ne120pg2.NN_512.SSTP_0K'
# case2 = 'E3SM.2024-AQP-CESS-00.FAQP-MMF1.GNUGPU.ne120pg2_ne120pg2.NN_512.SSTP_4K'
# add_case(case1,case2, n='MMF ne120',p=tmp_scratch,s=tmp_sub,grp_id=0,grp_str='MMF',res_str='ne120')



case1 = 'E3SM.2024-AQP-CESS-00.FAQP.GNUCPU.ne30pg2_ne30pg2.NN_32.SSTP_0K'
case2 = 'E3SM.2024-AQP-CESS-00.FAQP.GNUCPU.ne30pg2_ne30pg2.NN_32.SSTP_4K'
add_case(case1,case2, n='EAM ne30', p=tmp_scratch,s=tmp_sub,grp_id=1,grp_str='EAM',res_str='ne30')
case1 = 'E3SM.2024-AQP-CESS-00.FAQP.GNUCPU.ne45pg2_ne45pg2.NN_64.SSTP_0K'
case2 = 'E3SM.2024-AQP-CESS-00.FAQP.GNUCPU.ne45pg2_ne45pg2.NN_64.SSTP_4K'
add_case(case1,case2, n='EAM ne45', p=tmp_scratch,s=tmp_sub,grp_id=1,grp_str='EAM',res_str='ne45')
case1 = 'E3SM.2024-AQP-CESS-00.FAQP.GNUCPU.ne60pg2_ne60pg2.NN_128.SSTP_0K'
case2 = 'E3SM.2024-AQP-CESS-00.FAQP.GNUCPU.ne60pg2_ne60pg2.NN_128.SSTP_4K'
add_case(case1,case2, n='EAM ne60', p=tmp_scratch,s=tmp_sub,grp_id=1,grp_str='EAM',res_str='ne60')
case1 = 'E3SM.2024-AQP-CESS-00.FAQP.GNUCPU.ne90pg2_ne90pg2.NN_256.SSTP_0K'
case2 = 'E3SM.2024-AQP-CESS-00.FAQP.GNUCPU.ne90pg2_ne90pg2.NN_256.SSTP_4K'
add_case(case1,case2, n='EAM ne90', p=tmp_scratch,s=tmp_sub,grp_id=1,grp_str='EAM',res_str='ne90')
# case1 = 'E3SM.2024-AQP-CESS-00.FAQP.GNUCPU.ne120pg2_ne120pg2.NN_512.SSTP_0K'
# case2 = 'E3SM.2024-AQP-CESS-00.FAQP.GNUCPU.ne120pg2_ne120pg2.NN_512.SSTP_4K'
# add_case(case1,case2, n='EAM ne120',p=tmp_scratch,s=tmp_sub,grp_id=1,grp_str='EAM',res_str='ne120')


case1 = 'E3SM.2024-AQP-CESS-00.FAQP.GNUCPU.ne30pg2_ne30pg2.NN_32.SSTP_0K.ALT-NCPL_72'
case2 = 'E3SM.2024-AQP-CESS-00.FAQP.GNUCPU.ne30pg2_ne30pg2.NN_32.SSTP_4K.ALT-NCPL_72'
add_case(case1,case2, n='EAM ne30', p=tmp_scratch,s=tmp_sub,grp_id=2,grp_str='EAM (alt)',res_str='ne30')
case1 = 'E3SM.2024-AQP-CESS-00.FAQP.GNUCPU.ne45pg2_ne45pg2.NN_64.SSTP_0K.ALT-NCPL_72'
case2 = 'E3SM.2024-AQP-CESS-00.FAQP.GNUCPU.ne45pg2_ne45pg2.NN_64.SSTP_4K.ALT-NCPL_72'
add_case(case1,case2, n='EAM ne45', p=tmp_scratch,s=tmp_sub,grp_id=2,grp_str='EAM (alt)',res_str='ne45')
case1 = 'E3SM.2024-AQP-CESS-00.FAQP.GNUCPU.ne60pg2_ne60pg2.NN_128.SSTP_0K.ALT-NCPL_72'
case2 = 'E3SM.2024-AQP-CESS-00.FAQP.GNUCPU.ne60pg2_ne60pg2.NN_128.SSTP_4K.ALT-NCPL_72'
add_case(case1,case2, n='EAM ne60', p=tmp_scratch,s=tmp_sub,grp_id=2,grp_str='EAM (alt)',res_str='ne60')
case1 = 'E3SM.2024-AQP-CESS-00.FAQP.GNUCPU.ne90pg2_ne90pg2.NN_256.SSTP_0K.ALT-NCPL_72'
case2 = 'E3SM.2024-AQP-CESS-00.FAQP.GNUCPU.ne90pg2_ne90pg2.NN_256.SSTP_4K.ALT-NCPL_72'
add_case(case1,case2, n='EAM ne90', p=tmp_scratch,s=tmp_sub,grp_id=2,grp_str='EAM (alt)',res_str='ne90')
# case1 = 'E3SM.2024-AQP-CESS-00.FAQP.GNUCPU.ne120pg2_ne120pg2.NN_512.SSTP_0K.ALT-NCPL_72'
# case2 = 'E3SM.2024-AQP-CESS-00.FAQP.GNUCPU.ne120pg2_ne120pg2.NN_512.SSTP_4K.ALT-NCPL_72'
# add_case(case1,case2, n='EAM ne120',p=tmp_scratch,s=tmp_sub,grp_id=2,grp_str='EAM (alt)',res_str='ne120')

#---------------------------------------------------------------------------------------------------
num_case = len(case1_list)
num_var = len(var_list)

wkres = ngl.Resources()
npix=2048; wkres.wkWidth,wkres.wkHeight=npix,npix
wks = ngl.open_wks(fig_type,fig_file,wkres)

plot = [None]*num_var
#---------------------------------------------------------------------------------------------------
def load_time(case,case_root,case_sub):
   global htype,first_file,last_file
   file_path = f'{case_root}/{case}/{case_sub}/*eam.{htype}*'
   file_list = sorted(glob.glob(file_path))
   file_list = file_list[first_file:last_file]
   if len(file_list)==0:
      print( '\nERROR: load_time(): No files found!')
      print(f'  file_path: {file_path}')
      exit()
   ds = xr.open_mfdataset(file_list)
   return ds['time']
#---------------------------------------------------------------------------------------------------
def load_data(var,case,case_root,case_sub):
   global htype,first_file,last_file
   file_list = sorted(glob.glob(f'{case_root}/{case}/{case_sub}/*eam.{htype}*'))
   file_list = file_list[first_file:last_file]
   ds = xr.open_mfdataset(file_list)
   ds = ds.mean(dim='time')
   if var=='PRECT':
      da = ds['PRECC'] + ds['PRECL']
   else:
      da = ds[var]
   gbl_mean = ( (da*ds['area']).sum(dim='ncol') / ds['area'].sum(dim='ncol') )
   # gbl_mean = gbl_mean.values 
   if var=='PRECC': gbl_mean = gbl_mean*86400*1e3
   if var=='PRECT': gbl_mean = gbl_mean*86400*1e3
   return gbl_mean
#---------------------------------------------------------------------------------------------------
print()
for v in range(num_var):
   print(f'  var: {var_list[v]}')
   data_list = []
   for c in range(num_case):

      tmp_file = f'{tmp_file_head}.{case1_list[c]}.{var_list[v]}.nc'

      if recalculate:

         # if var_list[v]=='ECS':

         #    time  = load_time(case1_list[c],case_root[c],case_sub[c])
         #    dt = (time[1] - time[0])/86400e9
         #    num_day = float( (len(time)*dt).values )

         #    T1    = load_data('TS'  ,case1_list[c],case_root[c],case_sub[c])
         #    FSNT1 = load_data('FSNT',case1_list[c],case_root[c],case_sub[c])
         #    FLNT1 = load_data('FLNT',case1_list[c],case_root[c],case_sub[c])
         #    T2    = load_data('TS'  ,case2_list[c],case_root[c],case_sub[c])
         #    FSNT2 = load_data('FSNT',case2_list[c],case_root[c],case_sub[c])
         #    FLNT2 = load_data('FLNT',case2_list[c],case_root[c],case_sub[c])

         #    N1 = FSNT1 - FLNT1
         #    N2 = FSNT2 - FLNT2

         #    ECS = ( N2 - N1 ) / ( T2 - T1 )
         #    # ECS = ( T2 - T1 ) / ( N2 - N1 )

         #    msg = f'  {name_list[c]:10}'
         #    msg+= f'  nday : {num_day:6.1f}  '
         #    msg+= f'  ECS  : {ECS.values:8.4f}  '
         #    # msg+= f'  N2   : {N2:8.2f}  N1: {N1:8.2f}  dN: {(N2-N1):8.4f}'
         #    # msg+= f'  T2   : {T2:8.2f}  T1: {T1:8.2f}  dT: {(T2-T1):8.4f}'
         #    print(msg)

         #    data = ECS

         if var_list[v]=='RESTOM':
            FSNT1 = load_data('FSNT',case1_list[c],case_root[c],case_sub[c])
            FLNT1 = load_data('FLNT',case1_list[c],case_root[c],case_sub[c])
            # FSNT2 = load_data('FSNT',case2_list[c],case_root[c],case_sub[c])
            # FLNT2 = load_data('FLNT',case2_list[c],case_root[c],case_sub[c])
            # data = ( FSNT2 - FLNT2 ) - ( FSNT1 - FLNT1 )
            data = FSNT1 - FLNT1
         elif var_list[v]=='RESTOMC':
            FSNT1 = load_data('FSNTC',case1_list[c],case_root[c],case_sub[c])
            FLNT1 = load_data('FLNTC',case1_list[c],case_root[c],case_sub[c])
            # FSNT2 = load_data('FSNTC',case2_list[c],case_root[c],case_sub[c])
            # FLNT2 = load_data('FLNTC',case2_list[c],case_root[c],case_sub[c])
            # data = ( FSNT2 - FLNT2 ) - ( FSNT1 - FLNT1 )
            data = FSNT1 - FLNT1
         else:
            data = load_data(var_list[v] ,case1_list[c],case_root[c],case_sub[c])
            # data1 = load_data(var_list[v] ,case1_list[c],case_root[c],case_sub[c])
            # data2 = load_data(var_list[v] ,case2_list[c],case_root[c],case_sub[c])
            # data = data2 - data1

         #----------------------------------------------------------------------
         # write to file
         tmp_ds = xr.Dataset( coords=data.coords )
         tmp_ds[var_list[v]] = data
         tmp_ds.to_netcdf(path=tmp_file,mode='w')
      else:
         tmp_ds = xr.open_dataset( tmp_file, use_cftime=True  )
         data = tmp_ds[var_list[v]]

      
      # if var_list[v]=='TMQ'     : data = data*1e3 # kg/m2 => g/m2
      if var_list[v]=='TGCLDLWP': data = data*1e3 # kg/m2 => g/m2
      if var_list[v]=='TGCLDIWP': data = data*1e3 # kg/m2 => g/m2
      if var_list[v]=='CLDLOW'  : data = data*1e2 # frac => %
      if var_list[v]=='CLDHGH'  : data = data*1e2 # frac => %

      data_list.append(data.values)

   # print()

   
   #------------------------------------------------------------------------------------------------
   # regroup
   grp_id_set = list(set(grp_id_list))
   num_group = len(grp_id_set)
   data_array = np.full( [num_group,int(num_case/num_group)], np.nan )
   pos_array = np.full( [num_group,int(num_case/num_group)], np.nan )
   lbl_list = []
   lgd_list = []
   for g in range(num_group):
      cnt = 0
      for c in range(num_case):
         if grp_id_list[c]==grp_id_set[g]:
            data_array[g,cnt] = data_list[c]
            pos_array[g,cnt] = cnt
            if g==0: lbl_list.append(res_str_list[c])
            if cnt==0: lgd_list.append(grp_str_list[c])
            cnt += 1

   for g in range(num_group):
      data_array[g,:] = data_array[g,:] - data_array[g,0]

   # print(); print(pos_array)
   # print(); print(lgd_list)
   # exit()
   #------------------------------------------------------------------------------------------------
   # Create plot

   # ngl.define_colormap(wks,'BlAqGrYeOrReVi200')
   # clr = np.linspace(2,len( ngl.retrieve_colormap(wks) )-1,num_group,dtype=int)

   clr = ['red','blue','green']

   res = hs.res_xy()
   res.vpHeightF = 0.5
   res.tmYLLabelFontHeightF   = 0.02
   res.tmXBLabelFontHeightF   = 0.02
   res.tiXAxisFontHeightF     = 0.025
   res.tiYAxisFontHeightF     = 0.025
   res.xyLineThicknessF       = 10
   res.xyMarkLineMode         = 'MarkLines'
   res.xyMarkerSizeF          = 0.02
   res.xyMarker               = 16

   res.tiXAxisString          =  'Resolution'
   res.tiYAxisString          = f'Anomaly from ne30 [{unit_str[v]}]'
   
   # if var_list[v] in ['ECS']:
   #    res.tiYAxisString          = var_name[v]
   # else:
   #    # delta_symbol = '~F8~D~F21~'
   #    res.tiYAxisString          = var_name[v]

   res.tmXBMode               = 'Explicit'
   res.tmXBValues             = pos_array[0,:]
   res.tmXBLabels             = lbl_list

   res.xyLineColors           = clr
   res.xyMarkerColors         = clr

   x_span = np.max(pos_array)  - np.min(pos_array)
   y_span = np.max(data_array) - np.min(data_array)
   res.trXMinF = np.min(pos_array)  - 0.08*x_span
   res.trXMaxF = np.max(pos_array)  + 0.08*x_span
   res.trYMinF = np.min(data_array) - 0.08*y_span
   res.trYMaxF = np.max(data_array) + 0.08*y_span

   if v==0:
      res.xyExplicitLegendLabels = lgd_list
      res.pmLegendDisplayMode    = 'Always'
      # res.pmLegendOrthogonalPosF =-1.0 # upper-left
      # res.pmLegendParallelPosF   = 0.2
      # res.pmLegendOrthogonalPosF =-1.0 # upper-right
      # res.pmLegendParallelPosF   = 0.8
      # res.pmLegendOrthogonalPosF =-0.7 # middle-right
      # res.pmLegendParallelPosF   = 0.85
      res.pmLegendOrthogonalPosF =-0.55 # lower-right
      res.pmLegendParallelPosF   = 0.8
      res.pmLegendWidthF         = 0.22
      res.pmLegendHeightF        = 0.18

   lres = hs.res_xy()
   lres.xyDashPattern         = 1
   lres.xyLineThicknessF      = 1
   lres.xyLineColor           = "black"
   #----------------------------------------------------------------------------

   ip = v

   plot[ip] = ngl.xy(wks, pos_array, data_array, res)


   # lx = np.array([-1e5,1e5])
   # ngl.overlay( plot[ip], ngl.xy(wks, lx, lx*0 , lres) )

   hs.set_subtitles(wks, plot[ip], '', '', var_name[v], font_height=0.008)

#-------------------------------------------------------------------------------

layout = [int(np.ceil(len(plot)/float(num_plot_col))),num_plot_col]

pres = hs.setres_panel()
pres.nglPanelYWhiteSpacePercent = 10
pres.nglPanelXWhiteSpacePercent = 5
ngl.panel(wks,plot,layout,pres)

# ngl.draw(plot)
# ngl.frame(wks)
ngl.end()

hc.trim_png(fig_file)

#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------
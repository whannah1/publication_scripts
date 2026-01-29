import os, ngl, copy, xarray as xr, numpy as np, glob
# import hapy_common as hc, hapy_E3SM as he, hapy_setres as hs
#---------------------------------------------------------------------------------------------------
class tclr: END, GREEN, MAGENTA, CYAN = '\033[0m','\033[32m','\033[35m','\033[36m'
#---------------------------------------------------------------------------------------------------
case_name,case,case_dir,case_sub,case_grid,clr,dsh,mrk = [],[],[],[],[],[],[],[]
remap_flag = []
def add_case(case_in,n=None,p=None,s=None,g=None,d=0,c='black',m=0,r=False):
   global name,case,case_dir,case_sub,clr,dsh,mrk
   if n is None:
      tmp_name = ''
   else:
      tmp_name = n
   case.append(case_in); case_name.append(tmp_name)
   case_dir.append(p); case_sub.append(s); case_grid.append(g)
   dsh.append(d) ; clr.append(c) ; mrk.append(m)
   remap_flag.append(r)
#---------------------------------------------------------------------------------------------------
var,vclr,vdsh,var_str = [],[],[],[]
var_unit = []
def add_var(var_name,s=None,u='',clr='black',dsh=0): 
   var.append(var_name); vclr.append(clr); vdsh.append(dsh)
   var_str.append( var_name if s is None else s )
   var_unit.append(u)
#---------------------------------------------------------------------------------------------------
lev = np.array([30,50,75,100,125,150,200,250,300,350,400,450,500,550,600,650,700,750,800,825,850,875,900,925,950,975,1000])
# lev = np.array([1,2,3,4,5,6,7,8,9,10,15,20,25,30,40,50,100,150,200])
# lev = np.array([50,100,150,200,300,400,500,600,700,750,800,850,875,900,925,975])
# lev = np.array([0.1,0.5,1,2,3,4,5,6,7,8,9,10,15,20,25,30,40,50,100,150,200])
# lev = np.array([1,2,3,4,5,6,7,8,9,10,15,20,25,30,40,50,100,150,200,250,300,350,400,450,500,550,600,650,700,750,800,825,850,875,900,925,950,975,1000])
#---------------------------------------------------------------------------------------------------

pscratch,psub = '/pscratch/sd/w/whannah/e3sm_scratch/pm-gpu','run'
add_case('E3SM.2024-RCEMIP-DOMAIN-TEST-01.FRCE-MMF1_300.NX_256x1.DX_4000',n='256x1 4km', d=0, c='purple', p=pscratch,s=psub)
add_case('E3SM.2024-RCEMIP-DOMAIN-TEST-01.FRCE-MMF1_300.NX_128x1.DX_4000',n='128x1 4km', d=0, c='blue',   p=pscratch,s=psub)
add_case('E3SM.2024-RCEMIP-DOMAIN-TEST-01.FRCE-MMF1_300.NX_64x1.DX_4000', n=' 64x1 4km', d=0, c='green',  p=pscratch,s=psub)
add_case('E3SM.2024-RCEMIP-DOMAIN-TEST-01.FRCE-MMF1_300.NX_32x1.DX_4000', n=' 32x1 4km', d=0, c='red',    p=pscratch,s=psub)

#---------------------------------------------------------------------------------------------------

add_var('RELHUM',   s='RH', u='[%]')
add_var('MMF_MCUP' ,s='Upward Cloud Mass Flux (2 m/s threshold)',   u='[kg/m2/s]')
add_var('MMF_MCUP2',s='Upward Cloud Mass Flux (0.1 m/s threshold)', u='[kg/m2/s]')

# add_var('T')
# add_var('Q')
# add_var('CLDLIQ')
# add_var('CLDICE')
# add_var('U')
# add_var('V')
# add_var('OMEGA')
# add_var('CLDLIQ')
# add_var('CLDICE')

#---------------------------------------------------------------------------------------------------
# lat1,lat2 = -10,10
# lat1,lat2 = -30,30
#---------------------------------------------------------------------------------------------------

num_plot_col = len(var)

fig_file,fig_type = 'figs/profile.v1','png'
# fig_file,fig_type = 'figs/profile.v1','pdf'

htype,first_file,num_files = 'h0',0,1


plot_diff = False

use_height_coord   = False

#---------------------------------------------------------------------------------------------------
# Set up plot resources
#---------------------------------------------------------------------------------------------------
num_case,num_var = len(case),len(var)

wkres = ngl.Resources() ; npix = 2048*2 ; wkres.wkWidth,wkres.wkHeight=npix,npix
wks = ngl.open_wks(fig_type,fig_file,wkres)
plot = [None]*(num_var)
res = ngl.Resources()
res.nglDraw                      = False
res.nglFrame                     = False
res.tmXTOn                       = False
res.tmXBMajorOutwardLengthF      = 0.
res.tmXBMinorOutwardLengthF      = 0.
res.tmYLMajorOutwardLengthF      = 0.
res.tmYLMinorOutwardLengthF      = 0.
res.tmYLLabelFontHeightF         = 0.015
res.tmXBLabelFontHeightF         = 0.015
res.tiXAxisFontHeightF           = 0.015
res.tiYAxisFontHeightF           = 0.015
res.tmXBMinorOn                  = False
res.tmYLMinorOn                  = False
res.xyLineThicknessF             = 6.
# res.vpWidthF = 0.4
res.xyLineThicknessF = 15
# res.tmYLLabelFontHeightF = 0.008
# res.tmXBLabelFontHeightF = 0.008

res.tmXBAutoPrecision = False
res.tmXBPrecision = 2

if use_height_coord: 
   res.tiYAxisString = 'Height [km]'
else:
   res.tiYAxisString = 'Pressure [hPa]'
   res.trYReverse = True


ngl.define_colormap(wks,'MPL_coolwarm')
clr_tmp = np.linspace(2,len( ngl.retrieve_colormap(wks) )-1,6,dtype=int)
clr = clr_tmp[1:4+1]
clr = clr[::-1]

#---------------------------------------------------------------------------------------------------
def print_stat(x,name='(no name)',unit='',fmt='f',stat='naxh',indent='',compact=False):
   if fmt=='f' : fmt = '%.4f'
   if fmt=='e' : fmt = '%e'
   if unit!='' : unit = f'[{unit}]'
   name_len = 12 if compact else len(name)
   msg = ''
   line = f'{indent}{name:{name_len}} {unit}'
   if not compact: msg += line+'\n'
   for c in list(stat):
      if not compact: line = indent
      if c=='h' : line += '   shp: '+str(x.shape)
      if c=='a' : line += '   avg: '+fmt%x.mean()
      if c=='n' : line += '   min: '+fmt%x.min()
      if c=='x' : line += '   max: '+fmt%x.max()
      if c=='s' : line += '   std: '+fmt%x.std()
      if not compact: msg += line+'\n'
   if compact: msg += line
   print(msg)
   return msg
#---------------------------------------------------------------------------------------------------
def set_subtitles(wks, plot, left_string='', center_string='', right_string='', font_height=0.01):
   ttres         = ngl.Resources()
   ttres.nglDraw = False
   ttres.txFontHeightF = font_height
   # Use plot extent to call ngl.text(), otherwise you will see this error:
   # GKS ERROR NUMBER   51 ISSUED FROM SUBROUTINE GSVP  : --RECTANGLE DEFINITION IS INVALID
   strx = ngl.get_float(plot,"trXMinF")
   stry = ngl.get_float(plot,"trYMinF")
   # Set annotation resources to describe how close text is to be attached to plot
   amres = ngl.Resources()
   if not hasattr(ttres,"amOrthogonalPosF"):
      amres.amOrthogonalPosF = -0.50-0.02   # Top of plot plus a little extra to stay off the border
   else:
      amres.amOrthogonalPosF = ttres.amOrthogonalPosF
   # Add strings to the top of the plot
   if left_string != '':
      tx_id_l = ngl.text(wks, plot, left_string, strx, stry, ttres)
      amres.amJust,amres.amParallelPosF = "BottomLeft", -0.5   # Left-justified
      anno_id_l = ngl.add_annotation(plot, tx_id_l, amres)
   if center_string != '':
      tx_id_c = ngl.text(wks, plot, center_string, strx, stry, ttres)
      amres.amJust,amres.amParallelPosF = "BottomCenter", 0.0   # Centered
      anno_id_c = ngl.add_annotation(plot, tx_id_c, amres)
   if right_string != '':
      tx_id_r = ngl.text(wks, plot, right_string, strx, stry, ttres)
      amres.amJust,amres.amParallelPosF = "BottomRight", 0.5   # Right-justifed
      anno_id_r = ngl.add_annotation(plot, tx_id_r, amres)
   return
#---------------------------------------------------------------------------------------------------
def trim_png(fig_file,verbose=True):
   """ crop white space from png file """
   fig_file = fig_file+".png"
   fig_file = fig_file.replace(os.getenv('HOME')+'/Research/E3SM/','')
   if os.path.isfile(fig_file) :
      cmd = "convert -trim +repage "+fig_file+"   "+fig_file
      os.system(cmd)
      if verbose: print("\n"+fig_file+"\n")
   else:
      raise FileNotFoundError(f'\ntrim_png(): {fig_file} does not exist?!\n')
#---------------------------------------------------------------------------------------------------
data_list_list,lev_list_list = [],[]
for v in range(num_var):
   print('-'*80)
   print(' '*4+f'var: {tclr.GREEN}{var[v]}{tclr.END}')
   data_list,lev_list = [],[]
   for c in range(num_case):
      print(' '*4+f'case: {tclr.CYAN}{case[c]}{tclr.END}')
      #-------------------------------------------------------------------------
      file_path = f'{case_dir[c]}/{case[c]}/{case_sub[c]}/*eam.{htype}.*'
      file_list = sorted(glob.glob(file_path))
      if 'first_file' in globals(): file_list = file_list[first_file:]
      if 'num_files' in globals(): file_list = file_list[:num_files]
      #-------------------------------------------------------------------------
      ds = xr.open_mfdataset( file_list )
      #-------------------------------------------------------------------------
      ncol = ds['ncol'].values
      mask = xr.DataArray( np.ones([len(ncol)],dtype=bool), coords=[('ncol', ncol)], dims='ncol' )
      if 'lat1' in locals(): mask = mask & (ds[lat_name]>=lat1)
      if 'lat2' in locals(): mask = mask & (ds[lat_name]<=lat2)
      if 'lon1' in locals(): mask = mask & (ds[lon_name]>=lon1)
      if 'lon2' in locals(): mask = mask & (ds[lon_name]<=lon2)
      mask.compute()
      #-------------------------------------------------------------------------
      area = ds['area'].isel(time=0,missing_dims='ignore')
      data = ds[var[v]]
      if use_height_coord:
         Z = ds['Z3']
      #-------------------------------------------------------------------------
      area = area.where( mask, drop=True)
      data = data.where( mask, drop=True)
      if use_height_coord:
         Z = Z.where( mask, drop=True)
      #-------------------------------------------------------------------------
      data = data.mean(dim='time')
      if use_height_coord: Z = Z.mean(dim='time')
      #-------------------------------------------------------------------------
      if use_height_coord: Z = ( (Z*area).sum(dim='ncol') / area.sum(dim='ncol') )
      data = ( (data*area).sum(dim='ncol') / area.sum(dim='ncol') )
      #-------------------------------------------------------------------------
      print_stat(data,name=f'{var[v]} after averaging',indent=' '*4,compact=True)
      #-------------------------------------------------------------------------
      data_list.append( data.values )
      if use_height_coord:
         lev_list.append( Z.values )
      else:
         lev_list.append( data['lev'].values )
      #-------------------------------------------------------------------------
   data_list_list.append(data_list)
   lev_list_list.append(lev_list)

#-------------------------------------------------------------------------------
# Create plot
#-------------------------------------------------------------------------------
for v in range(num_var):
   data_list = data_list_list[v]
   lev_list = lev_list_list[v]
   
   tres = copy.deepcopy(res)

   tres.tiXAxisString = f'{var_unit[v]}'

   # ip = v*num_case+c
   ip = c*num_var+v
   
   baseline = data_list[0]
   if plot_diff:
      for c in range(num_case): 
         data_list[c] = data_list[c] - baseline

   data_min = np.min([np.nanmin(d) for d in data_list])
   data_max = np.max([np.nanmax(d) for d in data_list])
   tres.trXMinF = data_min
   tres.trXMaxF = data_max
   ip = v

   if 'MMF_MCUP' in var[v]:
      tres.trXMinF = 0
      tres.trXMaxF = 0.016

   for c in range(num_case):
      tres.xyLineColor   = clr[c]
      tres.xyMarkerColor = clr[c]
      tres.xyDashPattern = dsh[c]

      tplot = ngl.xy(wks, data_list[c], lev_list[c], tres)

      if (c==1 and plot_diff) or (c==0 and not plot_diff) :
         plot[ip] = tplot
      elif (plot_diff and c>0) or not plot_diff:
         ngl.overlay(plot[ip],tplot)

   ### add vertical line
   lres = ngl.Resources()
   lres.nglDraw                      = False
   lres.nglFrame                     = False
   lres.tmXTOn                       = False
   lres.tmXBMajorOutwardLengthF      = 0.
   lres.tmXBMinorOutwardLengthF      = 0.
   lres.tmYLMajorOutwardLengthF      = 0.
   lres.tmYLMinorOutwardLengthF      = 0.
   lres.tmYLLabelFontHeightF         = 0.015
   lres.tmXBLabelFontHeightF         = 0.015
   lres.tiXAxisFontHeightF           = 0.015
   lres.tiYAxisFontHeightF           = 0.015
   lres.tmXBMinorOn                  = False
   lres.tmYLMinorOn                  = False
   lres.xyLineThicknessF             = 6.
   lres.xyLineThicknessF = 1
   lres.xyDashPattern = 0
   lres.xyLineColor = 'black'
   ngl.overlay(plot[ip],ngl.xy(wks, np.array([0,0]), np.array([-1e3,1e8]), lres))

   cstr = ''
   if 'lat1' in locals(): 
      lat1_str = f'{lat1}N' if lat1>=0 else f'{(lat1*-1)}S'
      lat2_str = f'{lat2}N' if lat2>=0 else f'{(lat2*-1)}S'
      cstr += f' {lat1_str}:{lat2_str} '
   if 'lon1' in locals(): 
      lon1_str = f'{lon1}E' #if lon1>=0 and lon1<=360 else f'{(lon1*-1)}S'
      lon2_str = f'{lon2}E' #if lon2>=0 and lon2<=360 else f'{(lon2*-1)}S'
      cstr += f' {lon1_str}:{lon2_str} '

   lstr = var_str[v]
   if plot_diff: lstr += ' (diff)'

   set_subtitles(wks, plot[ip], '', cstr, lstr, font_height=0.01)

#-------------------------------------------------------------------------------
# Add legend
#-------------------------------------------------------------------------------
if num_case>1:
   lgres = ngl.Resources()
   lgres.vpWidthF           = 0.05
   lgres.vpHeightF          = 0.08
   lgres.lgLabelFontHeightF = 0.01
   lgres.lgMonoDashIndex    = True
   lgres.lgLineLabelsOn     = False
   lgres.lgLineThicknessF   = 15
   lgres.lgLabelJust        = 'CenterLeft'
   lgres.lgLineColors       = clr
   lgres.lgDashIndexes      = dsh

   lx,ly = 0.2,0.5

   pid = ngl.legend_ndc(wks, len(case_name), case_name, lx, ly, lgres)

#---------------------------------------------------------------------------------------------------
# Finalize plot
#---------------------------------------------------------------------------------------------------
layout = [int(np.ceil(len(plot)/float(num_plot_col))),num_plot_col]

pres = ngl.Resources()
pres.nglPanelYWhiteSpacePercent = 5
pres.nglPanelXWhiteSpacePercent = 5
pres.nglPanelTop      =  0.93

ngl.panel(wks,plot,layout,pres)
ngl.end()

trim_png(fig_file)
#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------

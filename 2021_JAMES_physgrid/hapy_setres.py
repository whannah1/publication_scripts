#===================================================================================================
# Walter Hannah - Lawrence Livermore National Lab
# 
# basic plot defaults
#  res_default()
#  res_xy()
# 
# contour plot defaults
#  res_contour()
#  res_contour_fill()
#  res_contour_fill_map()
# 
# panel plot defaults
#  setres_panel()
#  setres_panel_labeled()
# 
# plot string routines
#  set_subtitles(wks, plot, left_string, center_string, right_string, tres)
# 
#===================================================================================================
import xarray as xr
import numpy as np
import ngl
#---------------------------------------------------------------------------------------------------
# basic plot defaults
#---------------------------------------------------------------------------------------------------
def res_default():
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
   # res.gsnLeftStringFontHeightF    = 0.015
   # res.gsnRightStringFontHeightF   = 0.015
   # res.gsnLeftString               = ""
   # res.gsnrightString              = ""
   return res
#-------------------------------------------------------------------------------
def res_xy():
   res = res_default()
   res.xyLineThicknessF             = 6.
   return res
#---------------------------------------------------------------------------------------------------
# contour plot defaults
#---------------------------------------------------------------------------------------------------
def res_contour():
   res = res_default()
   res.cnFillOn                     = False
   res.cnLinesOn                    = True
   res.cnLineLabelsOn               = False
   res.cnInfoLabelOn                = False
   res.cnLineThicknessF             = 3.
   # res.gsnContourNegLineDashPattern = 1
   # res.gsnContourZeroLineThicknessF = 6
   # res.cnLevelSelectionMode         = "ExplicitLevels"
   return res
#-------------------------------------------------------------------------------
def res_contour_fill():
   res = res_contour()
   res.cnFillPalette = "ncl_default"
   res.cnFillOn                     = True
   res.cnLinesOn                    = False
   res.cnLineLabelsOn               = False
   res.cnInfoLabelOn                = False
   res.lbOrientation                = "Horizontal"
   res.lbLabelFontHeightF           = 0.008
   return res
#-------------------------------------------------------------------------------
def res_contour_fill_map():
   res = res_contour_fill()
   res.mpGridAndLimbOn = False
   res.mpCenterLonF     = 180
   res.mpLimitMode = "LatLon" 
   return res
#---------------------------------------------------------------------------------------------------
# panel plot defaults
#---------------------------------------------------------------------------------------------------
def setres_panel():
   res = ngl.Resources()
   res.nglPanelYWhiteSpacePercent = 5
   res.nglPanelXWhiteSpacePercent = 5
   return res
#-------------------------------------------------------------------------------
def setres_panel_labeled():
   res = setres_panel()
   res.nglPanelFigureStrings               = ['a','b','c','d','e','f','g','h','i'\
                                             ,'j','k','l','m','n','o','p','q','r'\
                                             ,'s','t','u','v','w','x','y','z']
   res.amJust                              = "TopLeft"
   res.nglPanelFigureStringsFontHeightF    = 0.015
   return res
#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------
def set_cell_fill(res,case_obj=None,lat=None,lon=None,htype=None,scrip_dir=None):
   if hasattr(res,'cnFillMode')    : del res.cnFillMode
   if hasattr(res,'sfXCStartV')    : del res.sfXCStartV
   if hasattr(res,'sfXCEndV')      : del res.sfXCEndV
   if hasattr(res,'sfYCStartV')    : del res.sfYCStartV
   if hasattr(res,'sfYCEndV')      : del res.sfYCEndV
   if hasattr(res,'sfXArray')      : del res.sfXArray
   if hasattr(res,'sfYArray')      : del res.sfYArray
   if hasattr(res,'sfXCellBounds') : del res.sfXCellBounds
   if hasattr(res,'sfYCellBounds') : del res.sfYCellBounds
   # if case_obj.obs or use_remap :
   if lat is not None and lon is not None :
      # res.cnFillMode    = "AreaFill"
      res.cnFillMode    = "RasterFill"
      if lat is None: lat = case_obj.lat
      if lon is None: lon = case_obj.lon
      if isinstance(lat, xr.DataArray):
         lat2D = np.repeat( lat.values[...,None],len(lon),axis=1)
         lon2D = np.repeat( lon.values[...,None],len(lat),axis=1)
      elif isinstance(lat,np.ndarray):
         lat2D = np.repeat( lat[...,None],len(lon),axis=1)
         lon2D = np.repeat( lon[...,None],len(lat),axis=1)
      else:
         raise AssertionError('input lat and lon data must be numpy array or xarray DataArray!')
      lon2D = np.transpose( lon2D )
      if case_obj.name=='TRMM' and not hasattr(case_obj,'lon1'): 
         lat2D,lon2D = ngl.add_cyclic(lat2D),ngl.add_cyclic(lon2D)
      res.sfXArray      = lon2D
      res.sfYArray      = lat2D
      res.sfXCStartV    = np.min(lon2D[0,:])
      res.sfXCEndV      = np.max(lon2D[0,:])
      res.sfYCStartV    = np.min(lat2D[:,0])
      res.sfYCEndV      = np.max(lat2D[:,0])
   else:
      res.cnFillMode    = "CellFill"
      scripfile = case_obj.get_scrip(scrip_dir)
      mask = case_obj.get_mask(htype=htype)
      
      #-------------------------------------------------------------------------
      # Deal with case of regional subset data
      if len(mask.values) < len(scripfile['grid_size'].values):

         ### get lat/lon bounds from coordinates
         ds = xr.open_dataset(case_obj.get_hist_file_list(htype=htype,component=case_obj.atm_comp)[0])
         lat = ds[ case_obj.ncol_name.replace('ncol','lat') ].values
         lon = ds[ case_obj.ncol_name.replace('ncol','lon') ].values
         lat1,lat2 = np.min(lat),np.max(lat)
         lon1,lon2 = np.min(lon),np.max(lon)

         # try comparing scrip coords to bounds of data coords
         lat_cond1 = scripfile['grid_center_lat'].values>=lat1
         lat_cond2 = scripfile['grid_center_lat'].values<=lat2
         lon_cond1 = scripfile['grid_center_lon'].values>=lon1
         lon_cond2 = scripfile['grid_center_lon'].values<=lon2
         lat_mask = np.logical_and(lat_cond1,lat_cond2)
         lon_mask = np.logical_and(lon_cond1,lon_cond2)
         scrip_mask = xr.DataArray( np.logical_and(lat_mask,lon_mask), \
                              coords=[scripfile['grid_size']], dims='grid_size' )

         if np.sum(scrip_mask) != np.sum(mask.values):
            # since the above approach didn't work,
            # try using ceil/floor values for bounds
            lat_cond1 = scripfile['grid_center_lat'].values>=np.floor(lat1)
            lat_cond2 = scripfile['grid_center_lat'].values<=np.ceil( lat2)
            lon_cond1 = scripfile['grid_center_lon'].values>=np.floor(lon1)
            lon_cond2 = scripfile['grid_center_lon'].values<=np.ceil( lon2)
            lat_mask = np.logical_and(lat_cond1,lat_cond2)
            lon_mask = np.logical_and(lon_cond1,lon_cond2)
            scrip_mask = xr.DataArray( np.logical_and(lat_mask,lon_mask), \
                                 coords=[scripfile['grid_size']], dims='grid_size' )

         if np.sum(scrip_mask) != np.sum(mask.values):
            # give up and throw an error
            raise AssertionError('Cannot match mask to scrip file coordinates!')
         else:
            mask = scrip_mask
      #-------------------------------------------------------------------------
      # apply the mask to the scrip coordinate data
      if 'ncol'   in mask.dims: mask = mask.rename({'ncol'  :'grid_size'})
      if 'ncol_d' in mask.dims: mask = mask.rename({'ncol_d':'grid_size'})
      res.sfXArray      = scripfile['grid_center_lon'].where(mask,drop=True).values
      res.sfYArray      = scripfile['grid_center_lat'].where(mask,drop=True).values
      res.sfXCellBounds = scripfile['grid_corner_lon'].where(mask,drop=True).values
      res.sfYCellBounds = scripfile['grid_corner_lat'].where(mask,drop=True).values
#---------------------------------------------------------------------------------------------------
# plot string routines
#---------------------------------------------------------------------------------------------------
# add three subtitles to the top of a plot, left, center, and right justified
def set_subtitles(wks, plot, left_string, center_string, right_string, font_height=0.01):
   ttres         = ngl.Resources()
   ttres.nglDraw = False

   # Use plot extent to call ngl.text(), otherwise you will see this error:
   # GKS ERROR NUMBER   51 ISSUED FROM SUBROUTINE GSVP  : --RECTANGLE DEFINITION IS INVALID
   strx = ngl.get_float(plot,"trXMinF")
   stry = ngl.get_float(plot,"trYMinF")

   ttres.txFontHeightF = font_height

   # Set annotation resources to describe how close text is to be attached to plot
   amres = ngl.Resources()
   if not hasattr(ttres,"amOrthogonalPosF"):
      amres.amOrthogonalPosF = -0.52   # Top of plot plus a little extra to stay off the border
   else:
      amres.amOrthogonalPosF = ttres.amOrthogonalPosF

   # if hasattr(amres,'tmEqualizeXYSizes') : del amres.tmEqualizeXYSizes
   # if hasattr(ttres,'tmEqualizeXYSizes') : del ttres.tmEqualizeXYSizes

   # Add strings to the top of the plot
   if left_string != '':
      # txidl = ngl.add_text(wks, plot, left_string, 0., 0., ttres)
      tx_id_l = ngl.text(wks, plot, left_string, strx, stry, ttres)
      amres.amJust         = "BottomLeft"
      amres.amParallelPosF = -0.5   # Left-justified
      anno_id_l            = ngl.add_annotation(plot, tx_id_l, amres)

   if center_string != '':
      tx_id_c = ngl.text(wks, plot, center_string, strx, stry, ttres)
      amres.amJust         = "BottomCenter"
      amres.amParallelPosF = 0.0   # Centered
      anno_id_c            = ngl.add_annotation(plot, tx_id_c, amres)

   if right_string != '':
      tx_id_r = ngl.text(wks, plot, right_string, strx, stry, ttres)
      amres.amJust         = "BottomRight"
      amres.amParallelPosF = 0.5   # Right-justifed
      anno_id_r            = ngl.add_annotation(plot, tx_id_r, amres)

   return
#-------------------------------------------------------------------------------
# add letter labels to each plot
def set_plot_labels(wks, plot_list, font_height=0.02, justify='left', case='lower'):
   if case=='lower': label = list( 'abcdefghijklmnopqrstuvwxyz' )
   if case=='upper': label = list( 'ABCDEFGHIJKLMNOPQRSTUVWXYZ' )
   
   ttres               = ngl.Resources()
   ttres.nglDraw       = False
   ttres.txFontHeightF = font_height

   # Set annotation resources to describe how close text is to be attached to plot
   amres = ngl.Resources()
   if not hasattr(ttres,"amOrthogonalPosF"):
      amres.amOrthogonalPosF = -0.4   # Top of plot minus a little 
   else:
      amres.amOrthogonalPosF = ttres.amOrthogonalPosF

   l = 0
   for plot in plot_list:
      strx = ngl.get_float(plot,"trXMinF")
      stry = ngl.get_float(plot,"trYMinF")
      tx_id = ngl.text(wks, plot, '('+label[l]+')', strx, stry, ttres)
      amres.amJust = "BottomCenter"
      if justify=='left' : amres.amParallelPosF = -0.42   # Left-justified
      if justify=='right': amres.amParallelPosF =  0.42   # Right-justifed
      anno_id_l = ngl.add_annotation(plot, tx_id, amres)
      l = l+1

      # # Add box - how can we know where the text is?
      # xx = strx
      # yy = stry
      # dx = ( ngl.get_float(plot,"trXMaxF") - ngl.get_float(plot,"trXMinF") ) * 0.03
      # dy = ( ngl.get_float(plot,"trYMaxF") - ngl.get_float(plot,"trYMinF") ) * 0.03
      # bx = np.array([ xx-dx, xx+dx, xx+dx, xx-dx, xx-dx ])
      # by = np.array([ yy+dy, yy+dy, yy-dy, yy-dy, yy+dy ])

      # lres = res_xy()
      # lres.xyLineThicknessF = 1.
      # lres.xyLineColor = 'red'
      # ngl.overlay(plot, ngl.xy(wks,bx,by,lres) )

   return
#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------
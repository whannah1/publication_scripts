begin

wks  = gsn_open_wks("png","figs/test_map")    ; send graphics to PNG file
res = True
res@mpProjection = "Robinson"
res@mpGridAndLimbOn = True
res@mpPerimOn = False
plot = gsn_csm_map(wks,res)        ; draw global map

end
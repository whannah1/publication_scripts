load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
; load "$HOME/NCL/custom_functions_E3SM.ncl"
load "$HOME/NCL/custom_functions.ncl"
begin

   ; scratch_dir = "/gpfs/alpine/scratch/hannah6/cli115/e3sm_scratch"
   scratch_dir = "/global/cscratch1/sd/whannah/e3sm_scratch/cori-knl/"

   ; case = (/"E3SM.PGVAL-XTRA.ne30_r05_oECv3.F2010SC5-CMIP6.master-32921c"/)
   ; name = (/"ne30np4"/)
   ; clr = (/"black"/)

   case = (/ \
            "E3SM.PGVAL-XTRA.ne30_r05_oECv3.F2010SC5-CMIP6.master-32921c" \
            ,"E3SM.PGVAL-XTRA.ne30pg4_r05_oECv3.F2010SC5-CMIP6.master-6a7c9b" \
            ,"E3SM.PGVAL-XTRA.ne30pg3_r05_oECv3.F2010SC5-CMIP6.master-6a7c9b" \
            ,"E3SM.PGVAL-XTRA.ne30pg2_r05_oECv3.F2010SC5-CMIP6.master-6a7c9b" \
            ; ,"E3SM.PGVAL-XTRA.ne30_r05_oECv3.F2010SC5-CMIP6.master-32921c" \
         /)
   name = (/"ne30np4","ne30pg4","ne30pg3","ne30pg2","ne30np4"/)
   clr = (/"black","blue","green","red","black"/)

   var = (/"DYN_PTTEND","PTTEND"/)
   ; var = (/"DYN_PTTEND"/)

   ; k_lev = 26 ; ~100mb
   k_lev = 45 ; ~500mb
   ; pnew = 500

   fig_type = "png"
   fig_file = "figs/F12-PT-spectra"

   write_file = False

   debug_mode = False

   temp_dir = "data/spectra_PT/"

;===============================================================================
;===============================================================================
   
   if debug_mode then
      print("")
      print("WARNING: DEBUG mode enabled!")
      print("")
   end if

   num_c = dimsizes(case)
   num_v = dimsizes(var)

   wks = gsn_open_wks(fig_type,fig_file)

   res = True
   res@gsnDraw                     = False
   res@gsnFrame                    = False
   res@tmXTOn                      = False
   res@tmXBMajorOutwardLengthF     = 0.
   res@tmXBMinorOutwardLengthF     = 0.
   res@tmYLMajorOutwardLengthF     = 0.
   res@tmYLMinorOutwardLengthF     = 0.
   res@tmYLLabelFontHeightF        = 0.015
   res@tmXBLabelFontHeightF        = 0.015
   res@tiXAxisFontHeightF          = 0.02
   res@tiYAxisFontHeightF          = 0.02
   ; res@tmYLMinorOn                 = False
   ; res@tmXBMinorOn                 = False
   

   ; res@vpHeightF = 0.4
   res@xyLineThicknessF            = 14.

   ; res@tiYAxisString = "Kinetic Energy [m~S~2~N~/s~S~2~N~]"
   ; res@tiYAxisString = "Physics Temperature Tendency [K/s]"
   res@tiYAxisString = "Physics Temperature Tendency [K/day]"
   res@tiXAxisString = "Spherical Wavenumber"
   

   res@xyXStyle ="Log"
   res@xyYStyle ="Log"
   res@tmXBMinorPerMajor = 8   
   res@tmXTMinorPerMajor = 8   
   res@tmYRMinorPerMajor = 8   
   res@tmYLMinorPerMajor = 8   
   ; res@tmYLLabelFont = 21    ;  21 = helvetica
   ; res@tmXBLabelFont = 21    ;  22 = helvetica-bold 

   ; res@xyLineThicknessF = 1.0
    
   ; res@xyLineColors = (/"turquoise","greenyellow","DarkGoldenRod1","red3","blue"/)
   ; res@xyDashPattern = 0

   res@tmGridDrawOrder = "PreDraw"

   res@tmXMajorGrid = True
   res@tmXMinorGrid = True
   res@tmXMajorGridLineColor = "gray"
   res@tmXMinorGridLineColor = "gray"
   res@tmXMajorGridThicknessF = 4
   res@tmXMinorGridThicknessF = 1

   res@tmYMajorGrid = True
   res@tmYMinorGrid = True
   res@tmYMajorGridLineColor = "gray"
   res@tmYMinorGridLineColor = "gray"
   res@tmYMajorGridThicknessF = 4
   res@tmYMinorGridThicknessF = 1

   ;;; use for PT in K/s
   ; res@trYMinF = 10e-13
   ; res@trYMaxF = 4e-11

   ;;; use for PT in K/day
   res@trYMinF = 3e-3
   ; res@trYMinF = 1e-2
   res@trYMaxF = 5e-1

   res@trXMinF = 1e1
   ; res@trXMinF = 4e1
   ; res@trXMaxF = 180  ; approx 2dx for ne30
   res@trXMaxF = 120  ; approx 3dx for ne30

;===============================================================================
;===============================================================================

   do c = 0,num_c-1
      
      tcase = case(c)

      print("")
      print("    "+tcase)

      do v = 0,num_v-1

         

         tfile = temp_dir+"spectra.phys_tend.v1."    +tcase+"."+var(v)+".klev_"+k_lev+".nc"
         ; tfile = temp_dir+"spectra.phys_tend.v1."    +tcase+".plev_"+pnew+".nc"

         print("")
         print("      "+var(v))
         ; print("      "+tfile)

         if write_file then
            ;----------------------------------------------------------------------
            ;----------------------------------------------------------------------
            sub_folder = "data_remap_721x1440"
            ; sub_folder = "data_remap_721x1440_dyn"
            ; sub_folder = "run"


            comp = "eam"

            ; ps_files   := systemfunc( "ls "+scratch_dir+"/"+tcase+"/"+sub_folder+"/"+tcase+"."+comp+".h1.*.nc" )
            data_files := systemfunc( "ls "+scratch_dir+"/"+tcase+"/"+sub_folder+"/"+tcase+"."+comp+".h2.*."+var(v)+".nc" )

            num_f = dimsizes(data_files)

            if debug_mode then num_f = 2 end if
            
            
            f_cnt = 0
            ; do f = num_f-1-num_ff,num_f-1
            do f = 0,num_f-1

               ; if debug_mode then f = 1 end if

               print("      "+data_files(f))
               ; exit()

               infile = addfile(data_files(f),"r")
               ; psfile = addfile(ps_files(f),"r")

               ;;; check levels
               ; lev = infile->lev
               ; nlev = dimsizes(lev)
               ; do k = 0,nlev-1
               ;    print(k+"  "+lev(k))
               ; end do
               ; exit

               num_t = dimsizes(infile->time)
               nlat  = dimsizes(infile->lat)
               ; nlev  = dimsizes(infile->lev)
               ; nlon  = dimsizes(infile->lon)
               

               ; if f.eq.0 then print("      nlat: "+nlat)  end if

               ; define spectra wavenumber
               spc := new ( (/nlat/), "double", 1d-99 )
               spc = spc@_FillValue
               spc!0 = "wavenumber"

               if f_cnt.eq.0 then
                  wavenum := spc     ; wavenumber for plotting
                  avg_spc := spc     ; average of multiple spectra
                  avg_spc(:) = 0
                  spc_cnt = 0
               end if

               ;-------------------------------------------------------------------
               ;-------------------------------------------------------------------

               if debug_mode then num_t = 2 end if

               do t=0,num_t-1

                  X := infile->$var(v)$(t,k_lev,:,:)

                  ; convert K/s to K/day
                  X = X*86400.

                  ;----------------------------------------------------------------
                  ; COMPUTE SPECTRA FROM U,V with vector spherical harmonics
                  ;-----------------------------------------------------------------
                  ; no longer needed, but save this code to verify it gives the
                  ; same answer as computing 'spc' from vor/div spectra
                  
                  dims = dimsizes(X)
                  nmax = dims(0)
                  mmax = dims(1)

                  print ("      computing spherical harmonic transform...")
                  ab = shagC(X(:,:)) ; scalar on Gauss grid

                     
                  ; compute energy in u,v components
                  cr = ab(0,:,:)                ; real coef  (nlat,nlat)
                  ci = ab(1,:,:)                ; imaginary  (nlat,nlat)
                  pwr = (cr^2 + ci^2)/2.        ; (nlat,nlat)  array
                  ab=0
                     
                  ; for clarity use do loops
                  do n1 = 0,nlat-1
                     wavenum(n1) = n1
                     spc(n1) = (/pwr(n1,0)/)  ; so we preserve attributes in spc
                     do m = 1,min( (/n1,mmax-1/) )
                        spc(n1) = spc(n1) + 2.*pwr(n1,m)
                     end do
                     spc(n1) = 0.25*spc(n1)
                  end do

                  avg_spc = ( spc_cnt*avg_spc + spc )/(spc_cnt+1)
                  spc_cnt = spc_cnt+1

               end do ; t

               f_cnt = f_cnt + 1

            end do ; f

            avg_spc@long_name = "avg_spc"
            printMAM(avg_spc)

            print("")
            print("f = "+f)
            do n = 0,10
               print("  "+sprintf("%5.1f", wavenum(n))+"  "+sprintf("%15.10f", avg_spc(n)) )
            end do

            ;----------------------------------------------------------------------
            ;----------------------------------------------------------------------
            print("")
            print("      Writing tmp file...  ("+tfile+")")
            if fileexists(tfile) then system("rm "+tfile) end if
            outfile = addfile(tfile,"c")
            outfile->avg_spc = avg_spc
            outfile->wavenum = wavenum
            outfile->cnt     = spc_cnt
            print("      done.")
            ;----------------------------------------------------------------------
            ;----------------------------------------------------------------------
         else

            tmpfile = addfile(tfile,"r")
            ; printVarSummary(tmpfile->avg_spc)
            avg_spc = tmpfile->avg_spc
            wavenum = tmpfile->wavenum
            spc_cnt = tmpfile->cnt

            ; printMAM(avg_spc)

         end if ; write_file
         ;-------------------------------------------------------------------------
         ;-------------------------------------------------------------------------      
         ; res@tiMainString     = "Climatology"
         res@xyLineColor = clr(c)
         ; res@xyDashPattern = dsh(c)

         if isStrSubset(var(v),"DYN") then
            res@xyDashPattern = 1
            if isStrSubset(case(c),"ne30_") then
               res@xyDashPattern = 0
            end if
         else
            res@xyDashPattern = 0
         end if

         ; if dsh(c).eq.16 then
         ;    newdash = NhlNewDashPattern(wks,(/"$$$$______$$$$______$$$$______$$$$______$$$$______"/))
         ;    res@xyDashPatterns = newdash
         ;    dsh(c) = newdash
         ; end if

         tplot = gsn_csm_xy(wks,wavenum,avg_spc,res)

         if c.eq.0 .and. v.eq.0 then
            plot = tplot
         else
            overlay( plot, tplot)
         end if

      end do ; v
   end do ; c

   polyres                  = True
   polyres@gsLineThicknessF = 1.0
   polyres@gsLineColor = "black"

   xref = (/10.,300./)
   spc_ref = 28*(xref/10)^(-3)
   gsn_polyline(wks,plot,xref,spc_ref,polyres)

   ;;; plot multiple reference lines
   ; xref = (/10.,300./)
   ; do i = -6,1
   ;    spc_ref = 2^i*15*(xref/10)^(-3)
   ;    gsn_polyline(wks,plot,xref,spc_ref,polyres)
   ; end do

;===============================================================================
;===============================================================================
   pres = True;setres_panel()
   ; pres = setres_panel_labeled()
   pres@gsnFrame = False
   pres@gsnPanelYWhiteSpacePercent = 5
   pres@gsnPanelBottom = 0.1
   ; pres@txString = season
   ; pres@gsLineThicknessF = 1.0
   ; pres@gsLineColor = "black"


   draw(plot)

   legend = create "Legend" legendClass wks
      "lgAutoManage"             : False
      "vpXF"                     : 0.25
      "vpYF"                     : 0.48
      "vpWidthF"                 : 0.15
      "vpHeightF"                : 0.2
      "lgPerimOn"                : False
      "lgTitleString"            : ""
      "lgTitleFontHeightF"       : 0.015
      "lgTitleOffsetF"           : 0.01 ; default = 0.03
      "lgItemCount"              : num_c
      "lgLabelStrings"           : " "+name(:num_c-1:-1)
      "lgLabelsOn"               : True
      "lgLineLabelsOn"           : False
      "lgLabelFontHeightF"       : 0.020
      ; "lgDashIndexes"            : dsh(:num_c-1:-1)
      "lgDashIndex"              : 0
      "lgMonoDashIndex"          : True
      "lgLineThicknessF"         : 30.
      "lgLineColors"             : clr(:num_c-1:-1)
      "lgLabelJust"              : "CenterLeft"
   end create
   draw(legend)

   frame(wks)

   trimPNG(fig_file)
;===============================================================================
;===============================================================================
end


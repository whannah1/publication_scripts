load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
; load "$HOME/NCL/custom_functions_E3SM.ncl"
load "$HOME/NCL/custom_functions.ncl"
begin

   ; scratch_dir = "/gpfs/alpine/scratch/hannah6/cli115/e3sm_scratch"
   scratch_dir = "/global/cscratch1/sd/whannah/e3sm_scratch/cori-knl/"

   ; case = (/ \
   ;          "E3SM.PGVAL.ne30_r05_oECv3.F2010SC5-CMIP6.master-cbe53b" \
   ;          ,"E3SM.PGVAL.ne30pg2_r05_oECv3.F2010SC5-CMIP6.master-cbe53b" \
   ;          ,"E3SM.PGVAL.ne30pg3_r05_oECv3.F2010SC5-CMIP6.master-cbe53b" \
   ;          ,"E3SM.PGVAL.ne30pg4_r05_oECv3.F2010SC5-CMIP6.master-cbe53b" \
   ;          ,"E3SM.PGVAL.ne45pg2_r05_oECv3.F2010SC5-CMIP6.master-327537" \
   ;          ,"E3SM.PGVAL.ne30pg2_r05_oECv3.F2010SC5-CMIP6.master-cbe53b.dyn" \
   ;          ; ,"E3SM.PGVAL.ne30pg3_r05_oECv3.F2010SC5-CMIP6.master-cbe53b.dyn" \
   ;          ; ,"E3SM.PGVAL.ne30pg4_r05_oECv3.F2010SC5-CMIP6.master-cbe53b.dyn" \
   ;          ; ,"E3SM.PGVAL.conusx4v1_r05_oECv3.F2010SC5-CMIP6.master-cbe53b" \
   ;          ; ,"E3SM.PGVAL.conusx4v1pg2_r05_oECv3.F2010SC5-CMIP6.master-cbe53b" \
   ;       /)
   ; name = (/"ne30np4","ne30pg2","ne30pg3","ne30pg4","ne45pg2","ne30pg2 (dyn)"/)
   ; clr = (/"black","red","green","blue","purple","red"/)
   ; dsh = (/0,0,0,0,2,16/)

   case = (/ \
            "E3SM.PGVAL.ne45pg2_r05_oECv3.F2010SC5-CMIP6.master-327537" \
            ,"E3SM.PGVAL.ne30_r05_oECv3.F2010SC5-CMIP6.master-cbe53b" \
            ,"E3SM.PGVAL.ne30pg2_r05_oECv3.F2010SC5-CMIP6.master-cbe53b" \
            ,"E3SM.PGVAL.ne30pg2_r05_oECv3.F2010SC5-CMIP6.master-cbe53b.dyn" \
         /)
   name = (/"ne45pg2","ne30np4","ne30pg2","ne30pg2 (dyn)"/)
   clr = (/"blue","black","red","red"/)
   dsh = (/0,0,0,16/)    


   ; define pressure level for interpolation
   ; pnew = 250
   pnew = 200  ; use this for model level
   ; pnew = 300
   ; pnew = 500

   fig_type = "png"
   fig_file = "figs/F13-KE-spectra"

   write_file = False

   debug_mode = False


   temp_dir = "data/spectra_KE/"

;===============================================================================
;===============================================================================
   
   if debug_mode then
      print("")
      print("WARNING: DEBUG mode enabled!")
      print("")
   end if

   num_c = dimsizes(case)

   wks = gsn_open_wks(fig_type,fig_file)

   ; res = setres_xy()

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
   res@xyLineThicknessF            = 16.

   res@tiYAxisString = "Kinetic Energy [m~S~2~N~/s~S~2~N~]"
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

   res@trYMinF = 1e-4
   res@trYMaxF = 5e1

   res@trXMinF = 1e1
   ; res@trXMaxF = 180  ; approx 2dx for ne30
   res@trXMaxF = 120  ; approx 3dx for ne30

;===============================================================================
;===============================================================================

   do c = 0,num_c-1
      
      tcase = case(c)

      if isStrSubset(tcase,".dyn") then
         tcase = str_sub_str(tcase,".dyn","")
         tfile = temp_dir+"spectra.KE.v1.dyn."+tcase+".plev_"+pnew+".nc"
      else
         tfile = temp_dir+"spectra.KE.v1."    +tcase+".plev_"+pnew+".nc"
      end if

      print("")
      print("    "+tcase)

      if write_file then
         ;----------------------------------------------------------------------
         ;----------------------------------------------------------------------
         sub_folder = "data_remap_721x1440"
         ; sub_folder = "run"

         if isStrSubset(tfile,"dyn")
            sub_folder = "data_remap_721x1440_dyn"
         end if

         comp = "cam"
         if isStrSubset(tcase,"ne45") then comp = "eam" end if

         ; ps_files   := systemfunc( "ls "+scratch_dir+"/"+tcase+"/"+sub_folder+"/"+tcase+"."+comp+".h1.*.nc" )
         data_files := systemfunc( "ls "+scratch_dir+"/"+tcase+"/"+sub_folder+"/"+tcase+"."+comp+".h2.*.nc" )

         ; if dimsizes(ps_files) .ne. dimsizes(data_files) then
         ;    print("ERROR: file number do not match!")
         ;    print(ps_files)
         ;    print("")
         ;    print(data_files)
         ;    print("")
         ;    exit
         ; end if

         num_f = dimsizes(data_files)

         if debug_mode then num_f = 2 end if
         
         
         f_cnt = 0
         ; do f = num_f-1-num_ff,num_f-1
         do f = 0,num_f-1

            print("      "+data_files(f))

            infile = addfile(data_files(f),"r")
            ; psfile = addfile(ps_files(f),"r")

            ; ntimes = dimsizes(infile->time)
            ; nlev   = dimsizes(infile->lev)
            nlat   = dimsizes(infile->lat)
            ; nlon   = dimsizes(infile->lon)

            ; if f.eq.0 then
            ;    print("      nlat: "+nlat)
            ; end if

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

            ;----------------------------------------------------------------------
            ; interpolate to pressure level
            ;----------------------------------------------------------------------
            ntime = 0      ; time index to use
            interp = 2     ; type of interpolation: 1 = linear, 2 = log, 3 = loglog
            extrap = True  ; is extrapolation desired if data is outside the range of PS
            P0mb = 1000    ; ps in Pa, but this argument must be in mb
            
            hyam = infile->hyam
            hybm = infile->hybm

            ; ps   = dim_avg_n_Wrap(psfile->PS(ntime:ntime+2,:,:),0)
            ; U := vinth2p( infile->U(ntime,:,:,:), hyam, hybm, pnew, ps, interp, P0mb, 1 ,extrap)
            ; V := vinth2p( infile->V(ntime,:,:,:), hyam, hybm, pnew, ps, interp, P0mb, 1 ,extrap)

            ; ; get rid of lev dimension
            ; U := U(0,:,:)
            ; V := V(0,:,:)

            ; KE = sqrt( U^2 + V^2 )

            num_t = dimsizes(infile->time)

            if debug_mode then num_t = 2 end if

            do t=0,num_t-1

               k_lev = 35 ; ~200mb

               if isStrSubset(tfile,"dyn") then
                  U := infile->DYN_U(t,k_lev,:,:)
                  V := infile->DYN_V(t,k_lev,:,:)
               else
                  U := infile->U(t,k_lev,:,:)
                  V := infile->V(t,k_lev,:,:)
                  ; ps := dim_avg_n_Wrap(psfile->PS(t:t+2,:,:),0)
                  ; U := vinth2p( infile->U(t,:,:,:), hyam, hybm, pnew, ps, interp, P0mb, 1 ,extrap)
                  ; V := vinth2p( infile->V(t,:,:,:), hyam, hybm, pnew, ps, interp, P0mb, 1 ,extrap)
                  ; ; get rid of lev dimension
                  ; U := U(0,:,:)
                  ; V := V(0,:,:)
               end if
               
               

               ;----------------------------------------------------------------------
               ; COMPUTE SPECTRA FROM U,V with vector spherical harmonics
               ;----------------------------------------------------------------------
               ; no longer needed, but save this code to verify it gives the
               ; same answer as computing 'spc' from vor/div spectra
               
               dims = dimsizes(U)
               nmax = dims(0)
               mmax = dims(1)

               print ("      computing vector spherical harmonic transform...")
               ab = vhagC(U(:,:),V(:,:));   ; Gauss grid
                  
               ; compute energy in u,v components
               cr = ab(0,:,:)                ; real coef  (nlat,nlat)
               ci = ab(1,:,:)                ; imaginary  (nlat,nlat)
               pwr = (cr^2 + ci^2)/2.        ; (nlat,nlat)  array
               cr = ab(2,:,:)                ; real coef  (nlat,nlat)
               ci = ab(3,:,:)                ; imaginary  (nlat,nlat)
               pwr = pwr + (cr^2 + ci^2)/2.        ; (nlat,nlat)  array
               ab=0

               ; printVarSummary(ps)
               ; printVarSummary(pwr)
                  
               ; for clarity use do loops
               do n1 = 0,nlat-1
                  wavenum(n1) = n1
                  spc(n1) = (/pwr(n1,0)/)  ; so we preserve attributes in spc
                  do m = 1,min( (/n1,mmax-1/) )
                     spc(n1) = spc(n1) + 2.*pwr(n1,m)
                  end do
                  spc(n1) = 0.25*spc(n1)
               end do

               ; spc@long_name = "spc"
               ; printMAM(spc)

               ; print("")
               ; print("f = "+f)
               ; do n = 0,10
               ;    print("  "+sprintf("%5.1f", wavenum(n))+"  "+sprintf("%15.10f", spc(n)) )
               ; end do

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

      end if ; write_file
      ;-------------------------------------------------------------------------
      ;-------------------------------------------------------------------------      
      ; res@tiMainString     = "Climatology"
      res@xyLineColor = clr(c)
      res@xyDashPattern = dsh(c)

      if dsh(c).eq.16 then
         newdash = NhlNewDashPattern(wks,(/"$$$$______$$$$______$$$$______$$$$______$$$$______"/))
         res@xyDashPatterns = newdash
         dsh(c) = newdash
      end if

      tplot = gsn_csm_xy(wks,wavenum,avg_spc,res)

      if c.eq.0 then
         plot = tplot
      else
         overlay( plot, tplot)
      end if

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
      "lgDashIndexes"            : dsh(:num_c-1:-1)
      ; "lgDashIndex"              : 0
      ; "lgMonoDashIndex"          : True
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


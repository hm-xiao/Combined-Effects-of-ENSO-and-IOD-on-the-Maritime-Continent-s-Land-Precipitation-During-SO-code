undef("ReadData")
function ReadData(fi, ccID)
begin
    f = addfile("Data/plot_data/"+fi,"r")

    u1 = f->plot_u(ccID,:,:)
    w1 = f->plot_w(ccID,:,:)
    u1p= f->prob_u(ccID,:,:)
    w1p= f->prob_w(ccID,:,:)

    u1@_FillValue = default_fillvalue("float")
    w1@_FillValue = default_fillvalue("float")
    return([/u1, w1, u1p, w1p/])
end 

begin
    sig =  0.05
    YAxisString = (/"Only El Nino", "Only PIOD", "El Nino+PIOD"/)
    choose_name = (/"Only_Nino", "Only_PIOD", "Nino_PIOD"/)

    scale = -100.0 ; w = w*-100

    ;figtype = "oldps";"eps"
    figtype = "eps"
    figname = "/work4/L.xhm/1091/paper/Picture/Fig2"
 
    f = addfile("Data/plot_data/ERA20C.nc","r")
    cc = f->choose_yy_name
    ccID = new(dimsizes(choose_name), integer)
    do i = 0,  dimsizes(choose_name)-1
        ccID(i) = ind(choose_name(i) .eq. cc)
    end do
   
    U1 =  ReadData("ERA20C.nc", ccID)

    area = (/1000., 100., 90., 150./)
    ylat = (/area(0), area(1), area(1), area(0), area(0)/)
    xlon = (/area(2), area(2), area(3), area(3), area(2)/)

   
   
  ; plot ************************************************************
    wks = gsn_open_wks(figtype, figname)

    res                      = True             ; plot mods desired
    res@gsnDraw              = False
    res@gsnFrame             = False
    res@gsnAddCyclic         = False    ; data already has cyclic point

    res@gsnPresHgtHeightLabelOn = False
    res@gsnLeftStringFontHeightF = 0.039
    res@gsnCenterStringFontHeightF = 0.042

    res@cnFillOn             = True               ; turn on color for contours
    res@cnFillMode           = "CellFill"
    res@cnLinesOn            = False              ; turn off contour lines
;    res@cnLineThicknessF     = 0.8
    res@cnLineLabelsOn       = False              ; turn off contour line labels
    cmap = read_colormap_file("ViBlGrWhYeOrRe")
    res@cnFillPalette        = cmap(:79,:)

    res@lbLabelBarOn         = False

    res@vpWidthF             = 0.88
    res@vpHeightF            = 0.68

    res@tmXTOn               = False
    res@tmYROn               = False
    res@trYMinF              = 100.
    res@trXMinF              = 30. ;60.
    res@trXMaxF              = 270.

    res@gsnRightString        = "hPa"
    res@tmXBLabelFontHeightF = 0.126 
    res@tmYLLabelFontHeightF = 0.052 
    res@cnInfoLabelOn        = False
    res@cnLevelSelectionMode = "ManualLevels"     ; set manual contour levels
    res@cnMinLevelValF       = -0.027
    res@cnMaxLevelValF       = 0.015
    res@cnLevelSpacingF      = 0.003

    resvc                    = True   
    resvc@gsnDraw            = False
    resvc@gsnFrame           = False
    resvc@gsnAddCyclic       = False    ; data already has cyclic point
    resvc@gsnMaximize        = True

    resvc@vcRefMagnitudeF    = 2.0                ; define vector ref mag
    resvc@vcRefLengthF       = 0.046              ; define length of vec ref
    resvc@vcGlyphStyle       = "CurlyVector"      ; turn on curly vectors
    resvc@vcMinDistanceF     = 0.082              ; thin out vectors
    resvc@vcMapDirection     = False
    resvc@vcLineArrowThicknessF = 1.62
    resvc@vcLineArrowHeadMaxSizeF = 0.022
    resvc@vcLineArrowHeadMinSizeF = 0.020
    resvc@vcLineArrowColor   = "gray58"

    resvc2 = resvc
    resvc2@vcLineArrowColor  = "red"


  ; Plot ________________________________________________________
    ERA = (/"(a)", "(b)", "(c)"/)
    Num1= (/11, 8, 15/) 
   
  ;  U1u = U1[0]
  ;  U1w = U1[1]
  ;  U2u = U2[0]
  ;  U2w = U2[1]

 ;   U1omg = U1[1] 
 ;   U2omg = U2[1] 

    U1u = where(U1[2].lt.sig .or. U1[3].lt.sig, U1[0], 0.)
    U1w = where(U1[2].lt.sig .or. U1[3].lt.sig, U1[1], 0.)

    U1omg = where(U1[3].lt.sig, U1[1], default_fillvalue("float"))

    copy_VarCoords(U1[0], U1u)
    copy_VarCoords(U1[1], U1w)

    copy_VarCoords(U1[1], U1omg)

    plot = new(3, graphic)
    text = (/"2a","2b","2c"/) 
    do i = 0,2
        res@tiYAxisString    = ""
        res@gsnLeftString    = "("+text(i)+") "+YAxisString(i)+" (n="+Num1(i)+")"
        res@tiMainString  = "";"ERA20C (1900-2010)"

        if(i.ne.0)then
            delete(res@tiMainString)
        end if

        plot(i)   = gsn_csm_pres_hgt(wks, U1omg(i,:,:),res)
    ;    plot2       = gsn_csm_vector(  wks, U1[0](i,:,:), U1[1](i,:,:)*scale, resvc)
        plot3       = gsn_csm_vector(  wks,   U1u(i,:,:),   U1w(i,:,:)*scale, resvc2)
        overlay(plot(i),  plot3) 
        delete(plot3)  
    end do


 ; Box _________________________________
    resline                  = True
    resline@gsLineThicknessF = 2.58
    resline@gsLineColor      = "black"

    dum = new((/3,4/),graphic)
    do i = 0,2
        do k = 0,3
            dum(i,k) = gsn_add_polyline(wks, plot(i), xlon, ylat, resline)
        end do
    end do

;  Map Only_________________________
    resline                  = True
    resline@gsLineThicknessF = 1.8
    resline@gsLineColor      = "blue"

    ResMap = True
    ResMap@gsnDraw            = False
    ResMap@gsnFrame           = False
    ResMap@mpOutlineBoundarySets       = "Geophysical"   ; turn on country boundaries
    ResMap@mpMaxLatF                   = 30           ; choose subregion
    ResMap@mpMinLatF                   = -30
    ResMap@mpMaxLonF                   = 280
    ResMap@mpMinLonF                   = 20
    ResMap@mpCenterLonF                = 150
    ResMap@vpWidthF             = 0.88
    ResMap@vpHeightF            = 0.68
    ResMap@tmXBLabelFontHeightF = 0.03
    ResMap@tmYLLabelFontHeightF = 0.03
    ResMap@gsnLeftString       = "(2d)"
    plot_map = gsn_csm_map_ce(wks,ResMap)                 ;
    DumLine =  gsn_add_polyline(wks, plot_map, (/30,30,180,270,270,180,30/), (/-10,10,10,10,-10,-10,-10/), resline)


 ; Panel ________________________________

   resP                     = True     ; plot mods desired
   resP@gsnAddCyclic        = False    ; data already has cyclic point
   resP@gsnMaximize         = True
   resP@gsnPanelLabelBar    = True
   resP@gsnPanelMainString  = "Vertical Profile" 
   resP@lbOrientation       = "Vertical"
   resP@lbLabelFontHeightF  = 0.014
   resP@lbTitleOn           = True
   resP@lbTitleString       = "Pa/s"
   resP@lbTitleFontHeightF  = 0.018


   resP@pmLabelBarOrthogonalPosF = 0.01
   resP@pmLabelBarParallelPosF   = 0.15
   resP@lbLabelFontHeightF       = 0.008
   resP@lbTitleFontHeightF       = 0.016
   gsn_panel(wks,(/plot(0),plot(1),plot(2),plot_map/), (/4,1/),resP)
;   gsn_panel(wks,(/plot(0),plot(1),plot(2),plot_map/), (/1,4/),resP)


 ;========================Transform Image File========================

 if (figtype.eq."oldps") then
     figtype = "jpg"
     system("convert -density 600 -trim "+figname+".ps  "+figname+"."+figtype)
     system("/bin/rm "+figname+".ps ")
 end if
 

end   

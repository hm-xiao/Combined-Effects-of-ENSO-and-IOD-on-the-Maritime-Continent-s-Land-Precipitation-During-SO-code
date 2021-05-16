begin
    yyB = 1900
    yyE = 2019

    neof   = 1        ; number of EOFs
    optEOF = True       
    optEOF@jopt = 0   ; This is the default; most commonly used; no need to specify.
    optEOF@pcrit = 50
    optETS = False

    Colors = (/"plum", "darkturquoise", "coral1", "gray88"/)
    ;Colors = (/"darkseagreen1", "darkturquoise", "orangered", "gray88"/)
    YAxisString = (/"Only PIOD", "Only El Nino", "El Nino+PIOD"/)

  ; ==============================================================
    f = addfile("/work4/L.xhm/1091/paper/Data/index_txt/ENSO_IOD_year.nc","r")
   var0 = f->Only_PIOD
   var1 = f->Only_Nino
   var2 = f->Nino_PIOD
     cp = f->cp_nino
     ep = f->ep_nino

    bar = new(yyE-yyB+1, string)
    bar =  "gray58"
    bar(var0(ind(var0.le.yyE))-yyB) = Colors(0)
    bar(var1(ind(var1.le.yyE))-yyB) = Colors(1)
    bar(var2(ind(var2.le.yyE))-yyB) = Colors(2)

    bar2 = new(yyE-yyB+1, string)
    bar2 = ""
    bar2(cp-yyB) = "C"
    bar2(ep-yyB) = "E"

    fn = addfile ("/work4/L.xhm/Data/Rain/precip.mon.total.v2020.nc", "r")

     tt = cd_calendar( (fn->time), 0)
    ttID = ind(yyB.le.tt(:,0) .and. tt(:,0).le.yyE)
     lat = tofloat(fn->lat({-10.:10.}))
     lon = tofloat(fn->lon({90.:150.}))
    Data = fn->precip(ttID,{-10.:10.},{90.:150.})

    Data@_FillValue    = default_fillvalue("float")
    Data@missing_value = default_fillvalue("float")

    xx  = ispan(0,dimsizes(Data(:,0,0))-1,1)
    Data = dtrend_msg_n(xx, Data, False, False, 0)
    Data = calcMonAnomTLL(  Data, clmMonTLL(Data) )

  wgtlat = cos(lat*4.*atan(1.)/180.)
    rain = wgt_areaave(Data, wgtlat, 1.0 ,0)
  rain_id = (rain(9-1::12)+rain(10-1::12))/2.

; WINDS
    f = addfile("/work4/L.xhm/1091/paper/Data/plot_data/ERA20C.nc","r")
    U = f->CompU
    W = f->CompU
    u = (U(9-1::12,:,:)+U(10-1::12,:,:))/2.
    w = (W(9-1::12,:,:)+W(10-1::12,:,:))/2.

  ; ==============================================================
    SLP = (Data(9-1::12,:,:)+Data(10-1::12,:,:))/2.
    SLP!0="time"
    SLP!1="lat"
    SLP!2="lon"
    SLP&lat = lat
    SLP&lon = lon
    SLP&lat@units = "degrees_north"
    SLP&lon@units = "degrees_east"

    rad    = 4.*atan(1.)/180.
    clat   = lat 
    clat   = sqrt( cos(rad*clat) )            
  
    wSLP   = SLP   ; Don't need ?
    wSLP   = SLP*conform(SLP, clat, 1)

    xw     = wSLP(lat|:,lon|:,time|:)
    x      = wSLP(time|:,lat|:,lon|:)
  
    eof      = eofunc_Wrap(xw, neof, optEOF)      
    eof_ts   = eofunc_ts_Wrap (xw, eof, optETS)

  ; =================================================================
  ; Normalize time series: Sum spatial weights over the area of used
  ; =================================================================
    NN = SLP(0,:,:)
    
    dimxw  = dimsizes( xw )
    mln    = dimxw(1)
    sumWgt = mln*sum( clat )
;    eof_ts = eof_ts/sumWgt

  ;============================================================
  ; Regression
  ;============================================================
    u_reg = regCoef_n(eof_ts(0,:110), u, 0, 0)
    w_reg = regCoef_n(eof_ts(0,:110), w, 0, 0)

  ;============================================================
  ; PLOTS
  ;============================================================
    figtype = "oldps"
    ;figtype = "eps" 
    figname = "Picture/Fig3"

    wks = gsn_open_wks(figtype, figname)
;    plot  = new(neof,graphic)                ; create graphic array
;    plot2 = new(neof,graphic)                ; create graphic array
                                          ; only needed if paneling
  ; EOF patterns
    res                      = True         
    res@gsnDraw              = False        ; don't draw yet
    res@gsnFrame             = False        ; don't advance frame yet
    res@gsnAddCyclic         = False        ; plotted dataa are not cyclic

    res@vpHeightF        = 0.58            ; Changes the aspect ratio
    res@vpWidthF         = 0.68
 
    res@mpFillOn             = False        ; turn off map fill
    res@mpMinLatF            = -10.
    res@mpMaxLatF            = 10.
    res@mpMinLonF            = 90.
    res@mpMaxLonF            = 150.

    res@cnFillOn             = True         ; turn on color fill
    res@cnLinesOn            = False        ; True is default
    res@cnLineLabelsOn       = False        ; True is default
    res@cnFillPalette        = "BlRe"
    res@cnFillMode           = "CellFill"
                                          ; set symmetric plot min/max
    symMinMaxPlt(eof, 16, False, res)       ; contributed.ncl


  ;*******************************************
  ; first plot
  ;*******************************************
    n = 0
    res@gsnLeftString  = "(3a)"
    res@gsnRightString = "EOF "+(n+1)+" "+sprintf("%5.1f", eof@pcvar(n)) +"%"
    plot = gsn_csm_contour_map_ce(wks,eof(n,:,:),res)

  ;*******************************************
  ; second plot
  ;*******************************************
  ; EOF time series  [bar form]
    xaxis = ispan(yyB,yyE,1)
    rts                  = True             ; plot mods desired
    rts@gsnScale         = True
    rts@gsnFrame         = False
    rts@gsnDraw          = False
    ;rts@gsnMaximize      = True

    rts@gsnXYBarChart         = True
    rts@gsnXYBarChartColors2  = bar
    rts@gsnXYBarChartBarWidth = .96
    rts@gsnXYBarChartOutlineThicknessF = .01
    rts@xyLineColor = "white"  
 
    rts@gsnYRefLine           = 0. 
    rts@gsnYRefLineColor      = "black"
    rts@gsnYRefLineThicknessF = 2.
  
    rts@tmEqualizeXYSizes    = True
    rts@tmXBMode         = "Explicit"
    rts@tmXBValues       = xaxis(::10)
    rts@tmXBLabels       = ispan(yyB,yyE,10)
    rts@tmXBLabelFontHeightF = 0.01

    rts@tmXTOn              = False
    rts@tmYROn              = False
    rts@tmXBMajorLengthF    = 0.009
    rts@trYMaxF             = 4400
    rts@trYMinF             = -4400

    rts@vpHeightF        = 0.26            ; Changes the aspect ratio
    rts@vpWidthF         = 0.628

    rts@tiYAxisString    = "Amplitude"
    rts@tiXAxisString    = "Year"
    rts@tiXAxisFontHeightF = 0.01

    rts2 = True
    rts2@tiYAxisSide      = "Right"
    rts2@tiYAxisString    = "Precipitation Anomaly (mm/mon)"
    rts2@tiYAxisFontColor = "purple"
    rts2@xyLineThicknesses  = 2.18
    rts2@xyLineColors       = "purple" 
    rts2@tmYRMajorLineColor = "purple"
    rts2@tmYRLabelFontColor = "purple"
    rts2@trYMaxF            = 100
    rts2@trYMinF            = -100

;    do n=0,neof-1
      n = 0
       rts@gsnLeftString  = "(3b)"
;       rts@gsnRightString = "EOF "+(n+1)+" "+sprintf("%5.1f", eof@pcvar(n)) +"%"
       plot2 = gsn_csm_xy2 (wks, xaxis, eof_ts(n,:), rain_id, rts, rts2)
;    end do

    txres               = True
    txres@txFontHeightF = .0072

    yy = where(eof_ts(0,:).lt.0, eof_ts(0,:)-200., eof_ts(0,:)+200.)
    text = gsn_add_text(wks, plot2, bar2, xaxis, yy, txres)

    resD                     = True
    resD@gsMarkerOpacityF    = 1.
    resD@gsMarkerThicknessF  = 1.58

    resD@gsMarkerIndex       = 1
    resD@gsMarkerColor       = "yellow"
    dum = gsn_add_polymarker(wks, plot2@xy2, xaxis, rain_id, resD)


    getvalues plot                ; get plot size for use in creating labelbar
    "vpXF"      : vpx
    "vpYF"      : vpy
    "vpHeightF" : vph
    "vpWidthF"  : vpw
    end getvalues

    lbw    = 0.4 * vpw           ; Make labelbar size a fraction of the plot.
    nboxes = 1 ;dimsizes(VarName)

    lbres                    = True          ; labelbar only resources
    lbres@vpWidthF           = 0.1          ; labelbar width
    lbres@vpHeightF          = 0.026          ; labelbar height
    lbres@lbBoxMajorExtentF  = 0.58          ; puts space between color boxes
    lbres@lbMonoFillPattern  = True          ; Solid fill pattern
    lbres@lbLabelFontHeightF = 0.0128          ; font height. default is small

    lbres@lbFillColors       = Colors(0)
    gsn_labelbar_ndc(wks,nboxes,YAxisString(0), 0.2, 0.5, lbres)

    lbres@lbFillColors       = Colors(1)
    gsn_labelbar_ndc(wks,nboxes,YAxisString(1), 0.38, 0.5, lbres)

    lbres@lbFillColors       = Colors(2)
    gsn_labelbar_ndc(wks,nboxes,YAxisString(2), 0.56, 0.5, lbres)

    lbres@lbFillColors       = Colors(3)
    gsn_labelbar_ndc(wks,nboxes,"Other Year", 0.76, 0.5, lbres)
  
    resP                     = True         ; modify the panel plot
    resP@gsnMaximize         = True         ; large format
    resP@gsnPanelLeft        = 0.12
    gsn_panel(wks,(/plot,plot2/),(/2,1/),resP)     ; now draw as one plot
 


  ;*******************************************
  ; third plot
  ;*******************************************
  ; Regression Wind
    rsv                      = True             ; plot mods desired
    rsv@gsnDraw              = False
    rsv@gsnFrame             = False
    rsv@gsnAddCyclic         = False    ; data already has cyclic point

    rsv@cnFillOn             = True               ; turn on color for contours
    rsv@cnFillMode           = "CellFill"
    rsv@cnLinesOn            = False              ; turn off contour lines
;    rsv@cnLineThicknessF     = 0.8
    rsv@cnLineLabelsOn       = False              ; turn off contour line labels
    cmap = read_colormap_file("ViBlGrWhYeOrRe")
    rsv@cnFillPalette        = cmap(:79,:)

    rsv@vpWidthF             = 0.7
    rsv@vpHeightF            = 0.58

    rsv@tmXTOn               = False
    rsv@tmYROn               = False
    rsv@trYMinF              = 100.
    rsv@trXMinF              = 30. ;60.
    rsv@trXMaxF              = 270.

    rsv@tmXBLabelFontHeightF = 0.1
    rsv@tmYLLabelFontHeightF = 0.04
    rsv@cnInfoLabelOn        = False
;    rsv@cnLevelSelectionMode = "ManualLevels"     ; set manual contour levels
;    rsv@cnMinLevelValF       = -0.027
;    rsv@cnMaxLevelValF       = 0.015
;    rsv@cnLevelSpacingF      = 0.003

    rsvvc                    = True
    rsvvc@gsnDraw            = False
    rsvvc@gsnFrame           = False
    rsvvc@gsnAddCyclic       = False    ; data already has cyclic point
    rsvvc@gsnMaximize        = True

;    rsvvc@vcRefMagnitudeF    = 2.0                ; define vector ref mag
    rsvvc@vcRefLengthF       = 0.046              ; define length of vec ref
    rsvvc@vcGlyphStyle       = "CurlyVector"      ; turn on curly vectors
    rsvvc@vcMinDistanceF     = 0.082              ; thin out vectors
    rsvvc@vcMapDirection     = False
    rsvvc@vcLineArrowThicknessF = 1.62
    rsvvc@vcLineArrowHeadMaxSizeF = 0.022
    rsvvc@vcLineArrowHeadMinSizeF = 0.020
    rsvvc@vcLineArrowColor   = "gray58"

    w_reg!0 = W!1
    w_reg&lev = W&lev
    w_reg&lev@units = W&lev@units

    w_reg!1 = W!2
    w_reg&lon = W&lon
    w_reg&lon@units = W&lon@units
    copy_VarCoords(w_reg, u_reg)

    plot3a = gsn_csm_pres_hgt(wks, w_reg, rsv)
    plot3b = gsn_csm_vector(  wks, u_reg, w_reg*-100., rsvvc)
;    overlay(plot3a,  plot3b)
;    draw(plot3a)
;    frame(wks)

;========================Transform Image File========================

 if (figtype.eq."oldps") then
     figtype = "jpg"
     system("convert -density 600 -trim "+figname+".ps  "+figname+"."+figtype)
     system("/bin/rm "+figname+".ps ")
 end if


end


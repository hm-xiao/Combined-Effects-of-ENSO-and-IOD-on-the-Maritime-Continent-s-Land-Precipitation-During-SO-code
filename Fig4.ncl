begin
    yyB = 1900
    yyE = 2019
    N   = yyE-yyB+1

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


;************************************
; compute wavelet (Missing values not allowed due to use of FFT)
;************************************  
    mother  = 0   ; Morlet
    param   = 6.0 ; Morlet
    dt      = 1.              ;for NAO (timesteps in units of years)
    s0      = dt
    dj      = 0.25
    jtot    = 1+floattointeger(((log10(N*dt/s0))/dj)/log10(2.)) 
    npad    = N
    nadof   = new(2,float)
    noise   = 1 
    siglvl  = .05
    isigtest= 1
  
   w = wavelet(rain_id,mother,dt,param,s0,dj,jtot,npad,noise,isigtest,siglvl,nadof)

;************************************
; create coodinate arrays for plot
;************************************
  power            = onedtond(w@power,(/jtot,N/))
  power!0          = "period"                        ; Y axis
  power&period     = w@period

  power!1          = "time"                          ; X axis
  power&time       = ispan(yyB,yyE,1)

  power@long_name  = "Power Spectrum"
  power@units      = "mm/mon^2"

  coi = w@coi

; compute significance ( >= 1 is significant)
  SIG              = power                            ; transfer meta data
  SIG              = power/conform (power,w@signif,0)
  SIG@long_name    = "Significance"
  SIG@units        = " "
;************************************
; initial resource settings
;************************************
  ;figtype = "eps" ;"oldps"
  figtype = "oldps"
  figname = "Picture/Fig4"
  wks = gsn_open_wks(figtype, figname)

  res                     = True                  ; plot mods desired 
  res@gsnFrame            = False                 ; don't advance frame yet
  res@gsnDraw             = False                 ; don't draw yet

  res@cnFillOn            = True                  ; turn on color
  res@cnFillPalette       = "WhBlGrYeRe" 
  res@cnFillMode          = "RasterFill"          ; turn on raster mode
  res@cnRasterSmoothingOn = True                  ; turn on raster smoothing
  res@cnLinesOn           = False                 ; turn off contour lines
  res@cnInfoLabelOn       = False
  res@lbOrientation       = "Vertical"            ; vertical label bar

  res@tmYLMode            = "Explicit"
  res@tmYLValues          = (/1,2,4,8,16,32,64,128/)
  res@tmYLLabels          = "" + res@tmYLValues
  res@tmLabelAutoStride   = True

  res@trYReverse          = True                  ; reverse y-axis
  res@tmLabelAutoStride   = True
  res@trYMaxF             = max(coi)

  res@vpHeightF           = .4                    ; height and width of plot
  res@vpWidthF            = .55

  res@gsnLeftString       =  "(4a)"
  res@tiYAxisString       = "Period (years)"
  plot = gsn_csm_contour(wks,power,res)  

  res2 = True                            ; res2 probability plots
  res2@gsnDraw             = False       ; Do not draw plot
  res2@gsnFrame            = False       ; Do not advance frome 
  res2@trYReverse          = True

  res2@cnLevelSelectionMode = "ManualLevels" ; set manual contour levels
  res2@cnMinLevelValF      = 0.00        ; set min contour level
  res2@cnMaxLevelValF      = 2.00        ; set max contour level
  res2@cnLevelSpacingF     = 1.00        ; set contour spacing
  res2@cnInfoLabelOn       = False
  res2@cnLinesOn           = False       ; do not draw contour lines
  res2@cnLineLabelsOn      = False       ; do not draw contour labels

  res2@tmYLMode            = "Explicit"
  res2@tmYLValues          = (/1,2,4,8,16,32,64,128/)
  res2@tmYLLabels          = "" + res2@tmYLValues

  res2@gsnLeftString       = ""
  res2@gsnRightString      = ""

  iplot_stipple   = gsn_csm_contour(wks,SIG,res2)
  opt                      = True
  opt@gsnShadeHigh         = 17          ; stipple
  opt@gsnShadeFillType     = "pattern"
  opt@gsnShadeFillScaleF   = 0.5         ; add extra density
  opt@gsnShadeFillDotSizeF = 0.002      ; make dots larger
  iplot_stipple = gsn_contour_shade(iplot_stipple,0, 0.8, opt)

  overlay(plot,iplot_stipple)

  ress = True
  ress@gsLineColor      = "purple"
  ress@gsLineThicknessF = 4.0
  aa = gsn_add_polyline(wks, plot, power&time, coi, ress)

  ;pres = True
  ;pres@gsnMaximize = True
  ;pres@gsnPaperOrientation = "portrait"
  ;gsn_panel(wks,(/plot,plot2/),(/2,1/),pres)

; ==================================
; PLOT2
; ==================================
  f2 = addfile("/work4/L.xhm/1091/CMIP5/per_comb_box_data.nc","r")
  plot_data = f2->plot_data
  xaxis = (/1., 2.5, 4., 5.5/)

  Res            = True                         ; plot mods desired
  Res@tiYAxisString        = "Combining Percentage"
  Res@trYMaxF              = .8
  Res@trYMinF              = 0.
  Res@tmBorderThicknessF   = 4.
  Res@tmYLLabelFontHeightF = 0.022
  Res@tmXBLabelFontHeightF = 0.022
  Res@vpHeightF  = .3    
  Res@vpWidthF   = .58

  opti          = True
  opti@boxWidth = .8               ; Width of box (x units)

  llRes                   = True
  llRes@gsLineThicknessF  = 3.2

  opti@boxColors = (/"blue","brown"/)
  Res@tmXBLabels = (/"1905-1944","1965-2004"/)
  plot2  = boxplot(wks, (/1,2/), plot_data(0:1,:), opti, Res, llRes)

  opti@boxColors = (/"green","red"/)
  Res@tmXBLabels = (/"piControl","rcp85"/)
  plot3  = boxplot(wks, (/1,2/), plot_data(2:3,:), opti, Res, llRes)

  mRes               = True                     ; marker mods desired
  mRes@gsMarkerIndex = 16                        ; polymarker style
  mRes@gsMarkerSizeF = 10.                   ; polymarker size
  mRes@gsMarkerThicknessF = 7.6

  mRes@gsMarkerColor = "blue"                    ; polymarker color
  dum1 = gsn_add_polymarker(wks,plot2, 1, 0.30769,mRes)
  mRes@gsMarkerColor = "brown"                    ; polymarker color
  dum2 = gsn_add_polymarker(wks,plot2, 2, 0.57143,mRes)

  txRes               = True
  txRes@txFontHeightF = .017
  text  = gsn_add_text(wks, plot2, "(4b)",  0.22, 0.8, txRes)
  text2 = gsn_add_text(wks, plot3, "(4c)",  0.22, 0.8, txRes)
    
  resP                     = True         ; modify the panel plot
  resP@gsnMaximize         = True         ; large format
  gsn_panel(wks,(/plot,plot2,plot3/),(/2,2/),resP)     ; now draw as one plot

;========================Transform Image File========================

 if (figtype.eq."oldps") then
     figtype = "jpg"
     system("convert -density 600 -trim "+figname+".ps  "+figname+"."+figtype)
     system("/bin/rm "+figname+".ps ")
 end if




end

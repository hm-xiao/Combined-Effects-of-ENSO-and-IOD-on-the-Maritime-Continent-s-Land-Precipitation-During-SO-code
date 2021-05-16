function COMP(fi, var, yyB, yyE)
begin

Num = yyE-yyB+1
   tt_array = ispan(6,  (yyE-yyB+1)*12-7, 1) -1 ; 1st Jun to 2nd  May

; Year
       year = ispan(yyB,yyE,1)
        rad = 4.*atan(1.)/180.

; Data__________________________________________________________
    f = addfile("/work4/L.xhm/1091/paper/Data/index_txt/ENSO_IOD_year.nc","r")
   OnlyPIOD = f->Only_PIOD
   OnlyNino = f->Only_Nino
   NinoPIOD = f->Nino_PIOD

   OP = OnlyPIOD(ind(yyB.le.OnlyPIOD .and. OnlyPIOD.le.yyE-1)) - yyB
   ON = OnlyNino(ind(yyB.le.OnlyNino .and. OnlyNino.le.yyE-1)) - yyB
   NP = NinoPIOD(ind(yyB.le.NinoPIOD .and. NinoPIOD.le.yyE-1)) - yyB

   if(var.eq."precip")then
       fn = addfile("/work4/L.xhm/Data/Rain/"+fi,"r")
       tt = cd_calendar( (fn->time), 0)
       ttID = ind(yyB.le.tt(:,0) .and. tt(:,0).le.yyE)
        lat = fn->lat({-10.:10.})
       Data = fn->$var$(ttID,{-10.:10.},{90.:150.})
   else
       fn = addfile("/work4/L.xhm/1091/TWS/Data/"+fi,"r")
       tt = cd_calendar( (fn->time), 0)
       ttID = ind(yyB.le.tt(:,0) .and. tt(:,0).le.yyE)
       Data = fn->$var$(0,ttID,:,:)
   end if
   
   Data@_FillValue    = default_fillvalue("float")
   Data@missing_value = default_fillvalue("float")

   xx  = ispan(0,dimsizes(Data(:,0,0))-1,1)
   Data = dtrend_msg_n(xx, Data, False, False, 0)
   
   if(var.eq."TWS")then
       Data_avg = clmMonTLL(Data)
       do i = 0,11
           Data_avg(i,:,:) = dim_avg_n(Data,0)
       end do
       Data = calcMonAnomTLL(Data, Data_avg)
   end if

  data    = dim_avg_n(dim_avg_n(Data,1),1)
  data_cc = data(tt_array)

  res_data = reshape(data_cc, (/dimsizes(data_cc)/12, 12/))
  cli_data = dim_avg_n(res_data,0)
  var_data = dim_stddev_n(res_data,0)

    lev_95 = new(12,float)
    lev_99 = new(12,float)

  do i = 0,11
      lev_95(i) = (cdft_t(0.975, Num-2))*var_data(i)/(Num^0.5)     
      lev_99(i) = (cdft_t(0.995, Num-2))*var_data(i)/(Num^0.5)     
  end do

   OP_data = dim_avg_n(res_data(OP,:), 0)
   ON_data = dim_avg_n(res_data(ON,:), 0)
   NP_data = dim_avg_n(res_data(NP,:), 0)

   OP_var  = dim_variance_n(res_data(OP,:), 0)
   ON_var  = dim_variance_n(res_data(ON,:), 0)
   NP_var  = dim_variance_n(res_data(NP,:), 0)
   
print(dimsizes(OP))
print(dimsizes(ON))
print(dimsizes(NP))
   return([/cli_data, OP_data, ON_data, NP_data, lev_95, lev_99,  var_data, OP_var, ON_var, NP_var, \
            dimsizes(data_cc)/12, dimsizes(OP), dimsizes(ON), dimsizes(NP)/])
end

;  =========================================================

function STA(xaxis, cli_data, lev_95, lev_99)
begin
  xp    = new( 2*12, float )
  yp_95 = new( 2*12, float )
  yp_99 = new( 2*12, float )

  do k = 0, 12-1
      xp(k)     = xaxis(k)
      yp_95(k)  = cli_data(k)+lev_95(k) 
      yp_99(k)  = cli_data(k)+lev_99(k) 

      xp(2*12-1-k) = xaxis(k)
      yp_95(2*12-1-k)  = cli_data(k)-lev_95(k) 
      yp_99(2*12-1-k)  = cli_data(k)-lev_99(k) 
  end do
  return([/xp, yp_95, yp_99/])
end

;  =========================================================

begin
  colors = (/"black", "green", "blue", "red"/)
   xaxis = ispan(1,12,1) 

   yyB = 1900
   yyE = 2019

  gpcc = COMP("precip.mon.total.v2020.nc",  "precip",  1900,  2019)
   tws = COMP("box_tws_clm.nc",  "TWS",  1911,  2015)

  
; plot ************************************************************
;  figtype = "oldps"
  figtype = "eps"
;figtype@wkWidth = 6000
;figtype@wkHeight = 6000

  figname = "Picture/Fig1"

  wks = gsn_open_wks(figtype, figname)

  res                      = True             ; plot mods desired
  res@gsnMaximize          = True
  res@gsnFrame             = False
  res@gsnDraw              = False

  res@vpWidthF             = 1.86
  res@vpHeightF            = 1.06

  res@tmXTOn               = False
  res@tmYROn               = False
  res@tmBorderThicknessF   = 5.8
  res@tmXBLabelFontHeightF = 0.049
  res@tmYLLabelFontHeightF = 0.049
  res@tmXBMode             = "Explicit"
  res@tmXBValues           = xaxis
  res@tmXBLabels           = (/6,7,8,9,10,11,12,1,2,3,4,5/)

  res@trXMaxF              =  12.4
  res@trXMinF              =  0.6

  res@xyLineColors         = colors 
  res@xyLineThicknesses    = (/4.2, 4.2, 4.2, 4.2, 4.2/)
  res@xyLineOpacities      = (/.9, .6, .6, .6, .6/)
  res@xyDashPatterns       = (/0, 0, 0, 0, 0/)

  res@tiXAxisString        = "Month" 
  res@tiXAxisFontHeightF   = 0.058
  res@tiYAxisFontHeightF   = 0.058

  ;add a legend
  res@pmLegendDisplayMode    = "Always"              ; turn on legend
  res@pmLegendSide           = "Top"                 ; Change location of
  res@pmLegendParallelPosF   = .82                   ; move units right
  res@pmLegendOrthogonalPosF =  -1.06   ; more neg = down

  res@pmLegendWidthF         = 0.2                  ; Change width and
  res@pmLegendHeightF        = 0.3                  ; height of legend.
  res@lgLabelFontHeightF     = .036                   ; change font height
  res@lgPerimOn              = False                 ; no box aroun
  res@tiMainFontHeightF    = 0.066


  gsres                   = True                        ; poly res
  gsres@tfPolyDrawOrder   = "Predraw"                   ; draw this first
  gsres@gsFillColor       = "grey68"                 ; color chosen

  resD                     = True
  resD@gsMarkerOpacityF    = 3.2
  resD@gsMarkerThicknessF  = 6.8

  resD@gsMarkerIndex       =  16 
  resD@gsMarkerColor       =  colors(0)

; ===========================================
; Precipitation
; ===========================================
  res@trYMaxF              =  300.
  res@trYMinF              =  120. 
  res@gsnLeftString        = "(1a)"
  res@tiYAxisString        = "precipitation (mm/mon)"
  res@xyExplicitLegendLabels = (/" Climatology", " Only PIOD (n=9)", " Only El Nino (n=10)", " El Nino+PIOD (n=17)"/)

  ; cli_data, OP_data, ON_data, NP_data, lev_95, lev_99, var_data, OP_var, ON_var, NP_var
  ; dimsizes(data_cc)/12, dimsizes(OP), dimsizes(ON), dimsizes(NP)/])
  res@tiMainString        = "MC's Land Precipitation Climatology";("+yyB+"-"+yyE+")"
  plot_gpcc  = gsn_csm_xy(wks, xaxis, (/gpcc[0], gpcc[1], gpcc[2], gpcc[3]/), res)


  ; STA(xaxis, cli_data, lev_95, lev_99)
  sts_gpcc = STA(xaxis, gpcc[0], gpcc[4], gpcc[5])

  ; xp, yp_95, yp_99
  dummy_gpcc = gsn_add_polygon (wks,plot_gpcc, sts_gpcc[0], sts_gpcc[1], gsres)


; Test
  sigl = 0.05 
  gpcc_test = new((/3,12/),float) ; OP, ON, NP
  do i = 0,2
      test = ttest(gpcc[0],gpcc[6],gpcc[10],gpcc[1+i],gpcc[7+i],gpcc[11+i],True,False)
      gpcc_test(i,:) = where(test.lt.sigl, 0, 1)
  end do

  id2a = ind(gpcc_test(0,:).eq.0)
  if(any(.not.ismissing(id2a)))then
    dum2a = gsn_add_polymarker(wks, plot_gpcc, xaxis(id2a), gpcc[1](id2a), resD)
  end if

  id3a = ind(gpcc_test(1,:).eq.0)
  if(any(.not.ismissing(id3a)))then
    dum3a = gsn_add_polymarker(wks, plot_gpcc, xaxis(id3a), gpcc[2](id3a), resD)
  end if

  id4a = ind(gpcc_test(2,:).eq.0)
  if(any(.not.ismissing(id4a)))then
    dum4a = gsn_add_polymarker(wks, plot_gpcc, xaxis(id4a), gpcc[3](id4a), resD)
  end if

; ===========================================
; TWS
; ===========================================
  res@trYMaxF              = 22.; 5540.
  res@trYMinF              = -62. 
  res@gsnLeftString        = "(1b)"
  res@tiYAxisString        = "TWS Anomaly(mm)"
  res@xyExplicitLegendLabels = (/" Climatology", " Only PIOD (n=9)", " Only El Nino (n=9)", " El Nino+PIOD (n=14)"/)

  ; cli_data, OP_data, ON_data, NP_data, lev_95, lev_99, var_data, OP_var, ON_var, NP_var
  ; dimsizes(data_cc)/12, dimsizes(OP), dimsizes(ON), dimsizes(NP)/])
  res@tiMainString        = "MC's TWSA Climatology"; ("+yyB+"-"+yyE+")"
  plot_tws  = gsn_csm_xy(wks, xaxis, (/tws[0], tws[1], tws[2], tws[3]/), res)

  ; STA(xaxis, cli_data, lev_95, lev_99)
  sts_tws = STA(xaxis, tws[0], tws[4], tws[5])

  ; xp, yp_95, yp_99
  dummy_tws = gsn_add_polygon (wks,plot_tws, sts_tws[0], sts_tws[1], gsres)

; Test
  sigl = 0.05
  tws_test = new((/3,12/),float) ; OP, ON, NP
  do i = 0,2
      test = ttest(tws[0],tws[6],tws[10],tws[1+i],tws[7+i],tws[11+i],True,False)
      tws_test(i,:) = where(test.lt.sigl, 0, 1)
  end do

  id2b = ind(tws_test(0,:).eq.0)
  if(any(.not.ismissing(id2b)))then
    dum2b = gsn_add_polymarker(wks, plot_tws, xaxis(id2b), tws[1](id2b), resD)
  end if

  resD@gsMarkerColor       =  colors(0)
  id3b = ind(tws_test(1,:).eq.0)
  if(any(.not.ismissing(id3b)))then
    dum3b = gsn_add_polymarker(wks, plot_tws, xaxis(id3b), tws[2](id3b), resD)
  end if

  resD@gsMarkerColor       =  colors(0)
  id4b = ind(tws_test(2,:).eq.0)
  if(any(.not.ismissing(id4b)))then
    dum4b = gsn_add_polymarker(wks, plot_tws, xaxis(id4b), tws[3](id4b), resD)
  end if


; Panel ________________________________

  resP                     = True     ; plot mods desired
  resP@gsnMaximize         = True
  resP@gsnPanelSave        = True

;  gsn_panel(wks, (/plot_gpcc/), (/1,1/), resP)
  gsn_panel(wks, (/plot_gpcc, plot_tws/), (/2,1/), resP)


 ;========================Transform Image File========================

 if (figtype.eq."oldps") then
     figtype = "jpg"
     system("convert -density 600 -trim "+figname+".ps  "+figname+"."+figtype)
     system("/bin/rm "+figname+".ps ")
 end if

end

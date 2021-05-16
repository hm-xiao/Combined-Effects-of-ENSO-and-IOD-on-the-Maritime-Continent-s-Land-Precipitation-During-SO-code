load  "year_Func.ncl"
;load  "/work4/L.xhm/1091/paper/Function/dtrend_Func.ncl"

begin

  yyB = 1900
  yyE = 2020

 year = ispan(yyB,yyE,1)
 rad  = 4.*atan(1.)/180.


; Read Data___________________________________________________________________
; SST
    f = addfile("/work4/L.xhm/Data/SST/HadISST_sst.nc","r")               ; data from:1870
 time = f->time
 tt   = cd_calendar(time,0)
 ttID = ind(yyB.le.tt(:,0) .and. tt(:,0).le.yyE)
  
  lat = f->latitude(::-1)
  lon = f->longitude
  sst = f->sst( ttID, ::-1, :)

  sst@_FillValue    = default_fillvalue("float")
  sst@missing_value = default_fillvalue("float")

  xx  = ispan(0,dimsizes(sst(:,0,0))-1,1)
  SST = dtrend_msg_n(xx, sst, False, False, 0)
  sst_5 = calcMonAnomTLL(  SST, clmMonTLL(SST) )
  copy_VarCoords(sst, sst_5)
 

; ENSO______________ 
; ONI
enso_sst  = sst_5(:,{-5.:5.},{-170.:-120.})  ; ONI (5N-5S, 170W-120W)
  wgtlat_enso = cos(lat({-5.:5.})*rad)
enso_id   = wgt_areaave(enso_sst, wgtlat_enso, 1.0 ,0)
ENSO_id   = runave(enso_id, 3, 0)
enso_id   = (enso_id-avg(enso_id))/stddev(enso_id)  ; ENSO index
enso_id   = runave(enso_id, 3, 0)

; EMI_______________
    lonID = ind(165..le.lon .or. lon.le.-140.)
    EMI   =  wgt_areaave(sst_5(:,{-10.:10.},lonID), cos(lat({-10.:10.})*rad), 1.0 ,0)-\ 
         0.5*wgt_areaave(sst_5(:,{-15.:5.}, {-110.:-70.}), cos(lat({-15.:5.})*rad), 1.0 ,0)-\ 
         0.5*wgt_areaave(sst_5(:,{-10.:20.},{125.:145.}),  cos(lat({-10.:20.})*rad), 1.0 ,0) 

    EMI   =  runave(EMI,  3, 0)
    EMI_NDJ = EMI(12-1::12)

    
; IOD_______________
; IODE 
  iode    = sst_5(:,{-10.:0.}, {90.:110.})
  wgtlate = cos(lat({-10.:0.})*rad)
  iode_id = wgt_areaave(iode, wgtlate, 1.0 ,0)
 nor_iode = (iode_id-avg(iode_id)) /stddev(iode_id)

; IODW
  iodw    = sst_5(:,{-10.:10.}, {50.:70.})
  wgtlatw = cos(lat({-10.:10.})*rad)
  iodw_id = wgt_areaave(iodw, wgtlatw, 1.0 ,0)
 nor_iodw = (iodw_id-avg(iodw_id)) /stddev(iodw_id)

  iod_ano = nor_iodw-nor_iode ; IOD index
  iod_ano = (iod_ano-avg(iod_ano)) /stddev(iod_ano)

  iod_cc = iod_ano
 enso_cc = enso_id


; SON & JJA
 dmiSON   = iod_cc(9-1::12)+iod_cc(10-1::12)+iod_cc(11-1::12)
 dmiSON   = dmiSON/3.
 dmiJJA   = iod_cc(6-1::12)+iod_cc(7-1::12)+iod_cc(8-1::12) 
 dmiJJA   = dmiJJA/3.
 
;Function  ***********************
PIOD = year(ind(dmiJJA.gt.1  .or. dmiSON.gt.1  )) ;positive IOD
NIOD = year(ind(dmiJJA.lt.-1 .or. dmiSON.lt.-1 )) ;negative IOD
Nino = year(ind( enso_cc(12-1::12).gt.1  ))  ;nino  NDJ > 1 std
Nina = year(ind( enso_cc(12-1::12).lt.-1 ))  ;nina  NDJ <-1 std

  cp = year( ind(EMI_NDJ.ge.(0.7*stddev(EMI_NDJ))) )
  ff = year_Func(year, Nino, Nina, PIOD, NIOD, cp)
;exit
; ***********************


; Output ________________________________
  OutFile = "index.nc"
  if (fileexists(OutFile)) then
      system("rm "+OutFile)  ; do this in ncl
  end if

; Output ************
    fout = addfile(OutFile,"c")
    fout->time=time(ttID)
    fout->enso_id=enso_cc
    fout->ENSO_id=ENSO_id

    fout->nor_iode=nor_iode;(ttID)
    fout->nor_iodw=nor_iodw;(ttID)
    fout->iod_id=iod_cc

    fout->dmiJJA=dmiJJA
    fout->dmiSON=dmiSON

    fout->EMI=EMI

    att               = True            ; assign file attributes
    att@creation_date = systemfunc ("date")
    att@contact       = "r07229022@ntu.edu.tw"
    fileattdef( fout, att )

print("end")

end

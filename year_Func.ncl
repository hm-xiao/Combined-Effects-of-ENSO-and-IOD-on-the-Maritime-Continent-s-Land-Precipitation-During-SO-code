undef("year_Func")
function year_Func(year, Nino, Nina, PIOD, NIOD, cp)
begin

;AMO
   f = addfile("/work4/L.xhm/1091/paper/Data/amo/amo_index.nc", "r")
 AMO = f->amo_dtrend_reg
 amo = ( AMO(9-1::12)+AMO(10-1::12) )/2.
 ;amo = (amo-avg(amo))/stddev(amo)
PAMO = year(ind(amo.gt.0))
NAMO = year(ind(amo.lt.0))

Norm      = new(dimsizes(year), typeof(year))
Only_PIOD = new(dimsizes(year), typeof(year))
Only_NIOD = new(dimsizes(year), typeof(year))
Only_Nino = new(dimsizes(year), typeof(year))
Only_Nina = new(dimsizes(year), typeof(year))

Nino_NIOD = new(dimsizes(year), typeof(year))
Nina_PIOD = new(dimsizes(year), typeof(year))
Nino_PIOD = new(dimsizes(year), typeof(year))
Nina_NIOD = new(dimsizes(year), typeof(year))

Nino_PAMO = new(dimsizes(year), typeof(year))
Nino_NAMO = new(dimsizes(year), typeof(year))
Nino_PIOD_PAMO = new(dimsizes(year), typeof(year))
Nino_PIOD_NAMO = new(dimsizes(year), typeof(year))

cp_nino      = new(dimsizes(year), typeof(year))
ep_nino      = new(dimsizes(year), typeof(year))

do i = 0,dimsizes(year)-1
    if(any(year(i).eq.Nino))then
        if(any(year(i).eq.cp))then
            cp_nino(i) = year(i)
        else
            ep_nino(i) = year(i)
        end if

        if(any(year(i).eq.PAMO))then
            Nino_PAMO(i) = year(i)
        end if
        if(any(year(i).eq.NAMO))then
            Nino_NAMO(i) = year(i)
        end if

        if(any(year(i).eq.NIOD))then
            Nino_NIOD(i) = year(i)
        else if(any(year(i).eq.PIOD))then
            Nino_PIOD(i) = year(i)

            if(any(year(i).eq.PAMO))then
                Nino_PIOD_PAMO(i) = year(i)
            end if
            if(any(year(i).eq.NAMO))then
                Nino_PIOD_NAMO(i) = year(i)
            end if

        else 
            Only_Nino(i) = year(i)
        end if
        end if


    else if(any(year(i).eq.Nina))then
        if(any(year(i).eq.NIOD))then
            Nina_NIOD(i) = year(i)
        else if(any(year(i).eq.PIOD))then
            Nina_PIOD(i) = year(i)
        else
            Only_Nina(i) = year(i)
        end if
        end if
    else if(any(year(i).eq.PIOD))then
        Only_PIOD(i) = year(i)
    else if(any(year(i).eq.NIOD))then
        Only_NIOD(i) = year(i)
    else 
        Norm(i) = year(i)
    end if
    end if
    end if
    end if
end do


; *******************
  OutFile = "ENSO_IOD_year.nc"
  if (fileexists(OutFile)) then
    system("rm "+OutFile)  ; do this in ncl
  end if

; Output ************
    fout = addfile(OutFile,"c")
    fout->PIOD=PIOD
    fout->NIOD=NIOD

    fout->Nino=Nino

if(any(.not.ismissing(Nino_PIOD)))
    fout->Nino_PIOD = Nino_PIOD(ind(.not.ismissing(Nino_PIOD)))
end if

if(any(.not.ismissing(Nino_NIOD)))
    fout->Nino_NIOD = Nino_NIOD(ind(.not.ismissing(Nino_NIOD)))
end if

if(any(.not.ismissing(Only_Nino)))
    fout->Only_Nino = Only_Nino(ind(.not.ismissing(Only_Nino)))
end if

if(any(.not.ismissing(Only_PIOD)))
    fout->Only_PIOD = Only_PIOD(ind(.not.ismissing(Only_PIOD)))    
end if


    fout->Nina=Nina
if(any(.not.ismissing(Nina_NIOD)))
    fout->Nina_NIOD = Nina_NIOD(ind(.not.ismissing(Nina_NIOD)))
end if

if(any(.not.ismissing(Nina_PIOD)))
    fout->Nina_PIOD = Nina_PIOD(ind(.not.ismissing(Nina_PIOD)))
end if

if(any(.not.ismissing(Only_Nina)))
    fout->Only_Nina = Only_Nina(ind(.not.ismissing(Only_Nina)))
end if

if(any(.not.ismissing(Only_NIOD)))
    fout->Only_NIOD= Only_NIOD(ind(.not.ismissing(Only_NIOD)))
end if

if(any(.not.ismissing(Norm)))
    fout->Norm = Norm(ind(.not.ismissing(Norm)))
end if

if(any(.not.ismissing(Nino_PAMO)))
    fout->Nino_PAMO = Nino_PAMO(ind(.not.ismissing(Nino_PAMO)))
end if

if(any(.not.ismissing(Nino_NAMO)))
    fout->Nino_NAMO = Nino_NAMO(ind(.not.ismissing(Nino_NAMO)))
end if

if(any(.not.ismissing(Nino_PIOD_PAMO)))
    fout->Nino_PIOD_PAMO = Nino_PIOD_PAMO(ind(.not.ismissing(Nino_PIOD_PAMO)))
end if

if(any(.not.ismissing(Nino_PIOD_NAMO)))
    fout->Nino_PIOD_NAMO = Nino_PIOD_NAMO(ind(.not.ismissing(Nino_PIOD_NAMO)))
end if

if(any(.not.ismissing(cp_nino)))
    fout->cp_nino = cp_nino(ind(.not.ismissing(cp_nino)))
end if
if(any(.not.ismissing(ep_nino)))
    fout->ep_nino = ep_nino(ind(.not.ismissing(ep_nino)))
end if


    att               = True            ; assign file attributes
    att@creation_date = systemfunc ("date")
    att@contact       = "r07229022@ntu.edu.tw"
    fileattdef( fout, att )

    return(fout)


end

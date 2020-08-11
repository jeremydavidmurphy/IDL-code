; This code is used to generate the plots of the final model mass. It
; requires that the dl.dat, dM.dat, dMhalo.dat, bindemo_r.dat and
; bindemo_v.dat files exist in the calling directory.

PRO Pradmass

; Modified on Jan 3,2012: This code now reads in a list called
; mass.list. This list is a list of the directories you want to pull
; dL.dat, dM.dat and dMhalo.dat files from in order to determine the
; uncertainties. THE BEST VALUE MUST BE THE FIRST ONE IN THIS LIST!

;*********************************************************
bhmass = 1.8e9
bhmin = 1.8e9
bhmax = 1.8e9
angleRmax = 2250.0 ;from galaxy.params file
xr = [0.1,angleRmax] ;xrange for plotting
yr = [1e7,1e14] ;yrange for plotting
scale = 0 ;this is the scale of the galaxy (pc/arcsec). This will plot a top axis in kpc. If set to 0 then no top axis is plotted.
mederr = 'yes' ;set to 'yes' if you want the best value to come from the median of all values. set to 'no' and the code will take the first model as the best model.
;*********************************************************

readcol,'mass.list',f='a,f',silent=1,directory,ML
ndir = n_elements(directory)
if (ndir eq 1) then begin
    print,'No uncertainies to be calculated!'
    print,'Changing to directory '+directory
    cd,directory
    readcol,'bindemo_r.out',silent=1,format='i,x,f,f,f,x,x',$
      skipline=4,rad,radl,radm,radh
;the initial zeros are dumped as they don't make it into dL.dat,
;dM.dat, etc
    rad = rad[1:*]
    radl = radl[1:*]
    radm = radm[1:*]
    radh = radh[1:*]
    readcol,'bindemo_v.out',silent=1,format='i,x,f,f,f',$
      skipline=4,angle,angl,angm,angh
    
    readcol,'dL.dat',silent=1,format='i,i,f',Lr,La,Lf ;mass from stars
    readcol,'dM.dat',silent=1,format='i,i,f',Mr,Ma,Mf ;total mass
    readcol,'dMhalo.dat',silent=1,format='i,i,f',Hr,Ha,Hf ;mass from DM halo
    
    Lf = Lf * ML[0]
    
    radl = radl * angleRmax
    radm = radm * angleRmax ;only this one is used
    radh = radh * angleRmax
;    angl = asin(angl) * 360.0/(2*!pi)
;    angm = asin(angm) * 360.0/(2*!pi)
;    angh = asin(angh) * 360.0/(2*!pi)
    
    na = n_elements(angle)
    nr = n_elements(rad)
    
    Ltot = dblarr(nr)
    Mtot = dblarr(nr)
    Htot = dblarr(nr)
    
    for j=0,nr-1 do begin       ;a step through each radial bin
        ilr = where(Lr eq rad[j])
        Lp = total(Lf[ilr])  ;the total stellar mass in one radial bin
        if (j ne 0) then Ltot[j] = Ltot[j-1] + Lp else Ltot[j] = Lp
        imr = where(Mr eq rad[j])
        Mp = total(Mf[imr])     ;the total mass in one radial bin
        if (j ne 0) then Mtot[j] = Mtot[j-1] + Mp else Mtot[j] = Mp + bhmass
        ihr = where(Hr eq rad[j])
        Hp = total(Hf[ihr])  ;the total DM halo mass in one radial bin
        if (j ne 0) then Htot[j] = Htot[j-1] + Hp else Htot[j] = Hp
    endfor
    
    if (scale gt 0.0) then begin
        rkpc = (radh * scale) / 1000.0
        xr2 = [min(rkpc),max(rkpc)]
        scalin = 'on'
    endif else scalin = 'off'
    
    set_plot,'ps'
    device,file='total_mass.ps',/color
    
    loadct,0
    if (scalin eq 'off') then begin
        plot,radm,htot,/xlog,/ylog,/nodata,yrange=[yr[0],yr[1]],$
          xrange=[xr[0],xr[1]],/xstyle,/ystyle,xthick=3,ythick=3,$
          charthick=3,thick=3,ytitle='Mass(R) (M!dsolar!n)',$
          xtitle='R (arcsec)',position=[0.1,0.11,0.98,0.91]
        oplot,radm,mtot,thick=3,linestyle=0
        loadct,4,/silent
        oplot,radm,htot,color=110,thick=3,linestyle=5
        oplot,radm,ltot,color=150,thick=3,linestyle=2
        device,/close
    endif
    if (scalin eq 'on') then begin
        plot,radh,htot,/xlog,/ylog,/nodata,yrange=[yr[0],yr[1]],$
          xrange=[xr[0],xr[1]],xstyle=9,/ystyle,xthick=3,ythick=3,$
          charthick=3,thick=3,ytitle='Mass(R) (M!dsolar!n)',$
          xtitle='R (arcsec)',position=[0.1,0.11,0.98,0.91]
        oplot,radm,mtot,thick=3,linestyle=0
        loadct,4,/silent
        oplot,radm,htot,color=110,thick=3,linestyle=5
        oplot,radm,ltot,color=150,thick=3,linestyle=2
        loadct,0,/silent
        axis,xaxis=1,xrange=[xr2[0],xr2[1]],/save,xtitle='R (kpc)',$
          charthick=3,/xstyle,xthick=3
        device,/close
    endif
    
    set_plot,'x'
    cd,'..'
endif

if (ndir gt 1) then begin
    print,'Uncertainies will be calculated!'
    for d=0,ndir-1 do begin ;a loop through each directory
        cd,directory[d]
        print,'Changing to directory '+directory[d]
        if (d eq 0) then begin
            readcol,'bindemo_r.out',silent=1,format='i,x,f,f,f,x,x',$
              skipline=4,rad,radl,radm,radh
            rad = rad[1:*]
            radl = radl[1:*]
            radm = radm[1:*]
            radh = radh[1:*]
            readcol,'bindemo_v.out',silent=1,format='i,x,f,f,f',$
              skipline=4,angle,angl,angm,angh
            na = n_elements(angle)
            nr = n_elements(rad)
            Ltot = dblarr(nr) ;the best values
            Mtot = dblarr(nr) ;the best values
            Htot = dblarr(nr) ;the best values
            radl = radl * angleRmax
            radm = radm * angleRmax
            radh = radh * angleRmax
            errarray = dblarr(nr,3,ndir) ;radial bins, Ltot/Mtot/Htot, # of directories
        endif

        readcol,'dL.dat',silent=1,format='i,i,f',Lr,La,Lf ;mass from stars
        readcol,'dM.dat',silent=1,format='i,i,f',Mr,Ma,Mf ;total mass
        readcol,'dMhalo.dat',silent=1,format='i,i,f',Hr,Ha,Hf ;mass from DM halo
        Lf = Lf * ML[d]
    
        for j=0,nr-1 do begin   ;a step through each radial bin
            ilr = where(Lr eq rad[j])
            Lp = total(Lf[ilr]) ;the total stellar mass in one radial bin
            if (j ne 0 and d eq 0) then Ltot[j] = Ltot[j-1] + Lp else if (d eq 0) then Ltot[j] = Lp
            if (j ne 0) then errarray[j,0,d] = errarray[j-1,0,d] + Lp else errarray[j,0,d] = Lp
            imr = where(Mr eq rad[j])
            Mp = total(Mf[imr]) ;the total mass in one radial bin
            if (j ne 0 and d eq 0) then Mtot[j] = Mtot[j-1] + Mp else if (d eq 0) then Mtot[j] = Mp + bhmass
            if (j ne 0) then errarray[j,1,d] = errarray[j-1,1,d] + Mp else errarray[j,1,d] = Mp + bhmass
            ihr = where(Hr eq rad[j])
            Hp = total(Hf[ihr]) ;the total DM halo mass in one radial bin
            if (j ne 0 and d eq 0) then Htot[j] = Htot[j-1] + Hp else if (d eq 0) then Htot[j] = Hp
            if (j ne 0) then errarray[j,2,d] = errarray[j-1,2,d] + Hp else errarray[j,2,d] = Hp
        endfor

        cd,'..'
    endfor

    error = dblarr(nr,6)
    for j=0,nr-1 do begin
        error[j,0] = min(errarray[j,0,*])
        error[j,1] = max(errarray[j,0,*])
        error[j,2] = min(errarray[j,1,*])-bhmin
        error[j,3] = max(errarray[j,1,*])+bhmax
        error[j,4] = min(errarray[j,2,*])
        error[j,5] = max(errarray[j,2,*])
    endfor

    if (scale gt 0.0) then begin ;to plot the scale in kpc along the top axis
        rkpc = (radh * scale) / 1000.0
        xr2 = [min(rkpc),max(rkpc)]
        scalin = 'on'
    endif else scalin = 'off'
    
    if (mederr eq 'yes') then begin ;to take the median of all the models as the best model (i.e. when you don't have a X^2 constraint)
        for j=0,nr-1 do begin
            Mtot[j] = median(errarray[j,1,*],/even)
            Ltot[j] = median(errarray[j,0,*],/even)
            Htot[j] = median(errarray[j,2,*],/even)
        endfor
    endif

    set_plot,'ps'
    device,file='total_mass.ps',/color

    if (ndir gt 2) then begin    
        loadct,0
        if (scalin eq 'off') then begin
            plot,radh,mtot,/xlog,/ylog,/nodata,yrange=[yr[0],yr[1]],$
              xrange=[xr[0],xr[1]],/xstyle,/ystyle,xthick=3,ythick=3,$
              charthick=3,thick=3,ytitle='Mass(R) (M!dsolar!n)',$
              xtitle='R (arcsec)',position=[0.1,0.11,0.98,0.91]
            oplot,radm,mtot,thick=3,linestyle=0
            oplot,radm,error[*,2],thick=2,linestyle=2
            oplot,radm,error[*,3],thick=2,linestyle=2
            loadct,4,/silent
            oplot,radm,htot,color=110,thick=3,linestyle=0
            oplot,radm,error[*,4],thick=2,linestyle=2,color=110
            oplot,radm,error[*,5],thick=2,linestyle=2,color=110
            oplot,radm,ltot,color=150,thick=3,linestyle=0
            oplot,radm,error[*,0],thick=2,linestyle=2,color=150
            oplot,radm,error[*,1],thick=2,linestyle=2,color=150
            device,/close
        endif
        if (scalin eq 'on') then begin
            plot,radh,htot,/xlog,/ylog,/nodata,yrange=[yr[0],yr[1]],$
              xrange=[xr[0],xr[1]],xstyle=9,/ystyle,xthick=3,ythick=3,$
              charthick=3,thick=3,ytitle='Mass(R) (M!dsolar!n)',$
              xtitle='R (arcsec)',position=[0.1,0.11,0.98,0.91]
            oplot,radm,mtot,thick=3,linestyle=0
            oplot,radm,error[*,2],thick=2,linestyle=2
            oplot,radm,error[*,3],thick=2,linestyle=2
            loadct,4,/silent
            oplot,radm,htot,color=110,thick=3,linestyle=0
            oplot,radm,error[*,4],thick=2,linestyle=2,color=110
            oplot,radm,error[*,5],thick=2,linestyle=2,color=110
            oplot,radm,ltot,color=150,thick=3,linestyle=0
            oplot,radm,error[*,0],thick=2,linestyle=2,color=150
            oplot,radm,error[*,1],thick=2,linestyle=2,color=150
            loadct,0,/silent
            axis,xaxis=1,xrange=[xr2[0],xr2[1]],/save,xtitle='R (kpc)',$
              charthick=3,/xstyle,xthick=3
            device,/close
        endif
    endif

    if (ndir eq 2) then begin    
        loadct,0
        if (scalin eq 'off') then begin
            plot,radh,mtot,/xlog,/ylog,/nodata,yrange=[yr[0],yr[1]],$
              xrange=[xr[0],xr[1]],/xstyle,/ystyle,xthick=3,ythick=3,$
              charthick=3,thick=3,ytitle='Mass(R) (M!dsolar!n)',$
              xtitle='R (arcsec)',position=[0.1,0.11,0.98,0.91]
            oplot,radm,errarray[*,1,0],thick=2
            oplot,radm,errarray[*,1,1],thick=2
            loadct,4,/silent
            oplot,radm,errarray[*,2,0],thick=2,color=110
            oplot,radm,errarray[*,2,1],thick=2,color=110
            oplot,radm,errarray[*,0,0],thick=2,color=150
            oplot,radm,errarray[*,0,1],thick=2,color=150
            device,/close
        endif
        if (scalin eq 'on') then begin
            plot,radh,htot,/xlog,/ylog,/nodata,yrange=[yr[0],yr[1]],$
              xrange=[xr[0],xr[1]],xstyle=9,/ystyle,xthick=3,ythick=3,$
              charthick=3,thick=3,ytitle='Mass(R) (M!dsolar!n)',$
              xtitle='R (arcsec)',position=[0.1,0.11,0.98,0.91]
            oplot,radm,errarray[*,1,0],thick=2
            oplot,radm,errarray[*,1,1],thick=2
            loadct,4,/silent
            oplot,radm,errarray[*,2,0],thick=2,color=110
            oplot,radm,errarray[*,2,1],thick=2,color=110
            oplot,radm,errarray[*,0,0],thick=2,color=150
            oplot,radm,errarray[*,0,1],thick=2,color=150
            loadct,0,/silent

            axis,xaxis=1,xrange=[xr2[0],xr2[1]],/save,xtitle='R (kpc)',$
              charthick=3,/xstyle,xthick=3
            device,/close
        endif
    endif

    set_plot,'x'
endif

stop
END

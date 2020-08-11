; Modified on Feb 7, 2012: This code now reads in a list of
; directories. Each directory is assumed to contain the intmom files
; you want to use for the uncertainties. IF ONLY 1 DIRECTORY IS IN THE
; LIST, THE CODE WILL HANDLE THIS BY NOT PLOTTING THE UNCERTAINIES.

PRO pintmom, dirlist

; DIRLIST: a list of the directories to cd into to pull the
; intmomSTARS and intmomGC files.
; SCALE: This is the size scale of the galaxy. EX: M87 has a scale of
; 86.5 pc/arcsec. This allows a second x-axis plotted showing kpc.
; N2832:475, N3842:440, N4889:450, n6086:655
; GALAXY: This is just for naming conventions. If used, the output
; file will be named "galaxy_intmom.ps'. Also, if the scale keyword IS
; NOT used, then the galaxy name is put into the title. (If 'scale' is
; used, there's no good space for this to be included.)

galaxy = 'M49'
scale = 86.5 ;set to 0 if you don't want to plot the top x-axis in kpc
x1 = 0.1 ;the ps plotting ranges
x2 = 1500
xs = [0.1,500] ;the range over which to plot the stars
xg = [0.1,1100] ;the range over which to plot the GCs
y1 = 0.4
y2 = 3.0
first = 'no' ;if set to 'yes' then it's assumed the first directory in the list is the best-fit model and is used. If set to 'no' then the median of all the values are taken and plotted as the best.
;******************************************************************
angles = [5.77,17.56,30.22,45.00,71.57]

readcol,dirlist,f='a',directory
n00 = n_elements(directory)

if (n00 eq 1) then begin
    cd,directory    
    print,'No uncertainies to be calculated!'
    print,'Changing to directory '+directory
    readcol,'intmomSTARS.out',format='f,f,f,f,x,f,f,f',R,theta,sigR,St,Sp,vp,beta,skipline=1,silent=1
    n0 = n_elements(R)
    n1 = n_elements(angles)
    n2 = n0 / n1
    sigarrS = dblarr(n2,3,n1)
    radiusS = fltarr(n2,3)

    cntr1 = 0
    cntr2 = 0
    for j=0,n0-1 do begin
        t1 = theta[j]
        t2 = angles[cntr1]
        if (t1 eq t2) then begin
            sigarrS[cntr2,0,cntr1] = sigR[j]
            sigarrS[cntr2,1,cntr1] = sqrt((St[j]^2 + Sp[j]^2 + vp[j]^2)*0.5)
            sigarrS[cntr2,2,cntr1] = sqrt((vp[j]^2)*0.5)
            cntr2 = cntr2 + 1
            if (cntr1 eq 0) then radiusS[j,0] = R[j]
        endif else begin
            cntr2 = 0
            cntr1 = cntr1 + 1
            sigarrS[cntr2,0,cntr1] = sigR[j]
            sigarrS[cntr2,1,cntr1] = sqrt((St[j]^2 + Sp[j]^2 + vp[j]^2)*0.5)
            sigarrS[cntr2,2,cntr1] = sqrt((vp[j]^2)*0.5)
            cntr2 = 1
        endelse
    endfor

    for j=0,n2-1 do begin
        RR = sigarrS[j,0,*]
        TH = sigarrS[j,1,*]
        div = RR / TH
        radiusS[j,1] = mean(div)
        TH = sigarrS[j,2,*]
        div = RR / TH
        radiusS[j,2] = mean(div)
    endfor

    radkpc1 = (x1 * scale) / 1000.0
    radkpc2 = (x2 * scale) / 1000.0

;and the same thing for the GCs
    readcol,'intmomGC.out',format='f,f,f,f,x,f,f,f',R,theta,sigR,St,Sp,vp,beta,skipline=1,silent=1
    n0 = n_elements(R)
    n1 = n_elements(angles)
    n2 = n0 / n1
    sigarrG = dblarr(n2,3,n1)
    radiusG = fltarr(n2,3)
    
    cntr1 = 0
    cntr2 = 0
    for j=0,n0-1 do begin
        t1 = theta[j]
        t2 = angles[cntr1]
        if (t1 eq t2) then begin
            sigarrG[cntr2,0,cntr1] = sigR[j]
            sigarrG[cntr2,1,cntr1] = sqrt((St[j]^2 + Sp[j]^2 + vp[j]^2)*0.5)
            sigarrG[cntr2,2,cntr1] = sqrt((vp[j]^2)*0.5)
            cntr2 = cntr2 + 1
            if (cntr1 eq 0) then radiusG[j,0] = R[j]
        endif else begin
            cntr2 = 0
            cntr1 = cntr1 + 1
            sigarrG[cntr2,0,cntr1] = sigR[j]
            sigarrG[cntr2,1,cntr1] = sqrt((St[j]^2 + Sp[j]^2 + vp[j]^2)*0.5)
            sigarrG[cntr2,2,cntr1] = sqrt((vp[j]^2)*0.5)
            cntr2 = 1
        endelse
    endfor

    for j=0,n2-1 do begin
        RR = sigarrG[j,0,*]
        TH = sigarrG[j,1,*]
        div = RR / TH
        radiusG[j,1] = mean(div)
        TH = sigarrG[j,2,*]
        div = RR / TH
        radiusG[j,2] = mean(div)
    endfor

;now the values that are not to be plotted (as determined by 'xs' and
;'xg') are driven out of the plotting range by being set to -1000
    i1 = where(radiusS[*,0] lt xs[0],c)
    if (c gt 0) then radiusS[i1,0] = -1000 
    i1 = where(radiusS[*,0] gt xs[1],c)
    if (c gt 0) then radiusS[i1,0] = -1000 
    i1 = where(radiusG[*,0] lt xg[0],c)
    if (c gt 0) then radiusG[i1,0] = -1000 
    i1 = where(radiusG[*,0] gt xg[1],c)
    if (c gt 0) then radiusG[i1,0] = -1000 

    set_plot,'ps'
    device,file=galaxy+'_intmom.ps',/color
    loadct,0,/silent
    plot,radiusS[*,0],radiusS[*,1],/xlog,xrange=[x1,x2],xstyle=9,$
      xtitle='!6Radius (arcsec)',ytitle='!7r!3!dr!n/!7r!3!dt',$
      charsize=1.2,xthick=3,ythick=3,thick=5,charthick=3,/nodata,$
      yrange=[y1,y2],/ys,position=[0.13,0.13,0.93,0.89]
    loadct,4,/silent
    oplot,radiusS[*,0],radiusS[*,1],thick=5,color=150
    oplot,radiusS[*,0],radiusS[*,2],thick=5,color=150,linestyle=2
    oplot,radiusG[*,0],radiusG[*,1],thick=5,color=60
    oplot,radiusG[*,0],radiusG[*,2],thick=5,color=60,linestyle=2
    loadct,0,/silent
    axis,xaxis=1,xrange=[radkpc1,radkpc2],/save,xtitle='!6R (kpc)',$
      charthick=3,charsize=1.2,xstyle=1,xthick=3
    device,/close_file
    set_plot,'x'
    cd,'../'
endif

; THE SECOND LOOP: where there are more than 1 directory and
; uncertainies are calculated
if (n00 gt 1) then begin
    print,'Uncertainies will be calculated!'
    for jj=0,n00-1 do begin
        cd,directory[jj]
        print,'Changing to directory '+directory[jj]
        readcol,'intmomSTARS.out',format='f,f,f,f,x,f,f,f',R,theta,sigR,St,Sp,vp,beta,skipline=1,silent=1
        readcol,'intmomGC.out',format='f,f,f,f,x,f,f,f',R,theta,sigR2,St2,Sp2,vp2,beta2,skipline=1,silent=1
        if (jj eq 0) then begin
            n0 = n_elements(R)
            n1 = n_elements(angles)
            n2 = n0 / n1
            sigarrS = dblarr(n2,2,n1) ;radial position,sig_r & sig_t,
            sigarrG = dblarr(n2,2,n1)
            radiusS = fltarr(n2,2,n00) 
            radiusG = fltarr(n2,2,n00)
        endif
        
        cntr1 = 0
        cntr2 = 0
        for k=0,n0-1 do begin ;a loop through the radial positions
            t1 = theta[k]
            t2 = angles[cntr1]
            if (t1 eq t2) then begin ;we are at a new angular position
                sigarrS[cntr2,0,cntr1] = sigR[k]^2
                sigarrS[cntr2,1,cntr1] = (St[k]^2 + Sp[k]^2 + vp[k]^2)*0.5
                sigarrG[cntr2,0,cntr1] = sigR2[k]^2
                sigarrG[cntr2,1,cntr1] = (St2[k]^2 + Sp2[k]^2 + vp2[k]^2)*0.5
                cntr2 = cntr2 + 1
                if (cntr1 eq 0) then radiusS[k,0,jj] = R[k]
                if (cntr1 eq 0) then radiusG[k,0,jj] = R[k]
            endif else begin
                cntr2 = 0
                cntr1 = cntr1 + 1
                sigarrS[cntr2,0,cntr1] = sigR[k]^2
                sigarrS[cntr2,1,cntr1] = (St[k]^2 + Sp[k]^2 + vp[k]^2)*0.5
                sigarrG[cntr2,0,cntr1] = sigR2[k]^2
                sigarrG[cntr2,1,cntr1] = (St2[k]^2 + Sp2[k]^2 + vp2[k]^2)*0.5
                cntr2 = 1
            endelse
        endfor
        for k=0,n2-1 do begin
            RR = sigarrS[k,0,*]
            TH = sigarrS[k,1,*]
            div = RR / TH
            radiusS[k,1,jj] = sqrt(mean(div))
            RR = sigarrG[k,0,*]
            TH = sigarrG[k,1,*]
            div = RR / TH
            radiusG[k,1,jj] = sqrt(mean(div))
        endfor
        cd,'../'
    endfor ;a loop over each directory
    ;now all the values are read in
    if (first eq 'no') then stars = median(radiusS[*,1,*],dim=3) else stars = radiusS[*,1,0]
    if (first eq 'no') then gc = median(radiusG[*,1,*],dim=3) else gc = radiusR[*,1,0]
;now the values that are not to be plotted (as determined by 'xs' and
;'xg') are driven out of the plotting range by being set to -1000

    i1 = where(radiusS[*,0] lt xs[0],c)
    if (c gt 0) then radiusS[i1,0] = -1000 
    i1 = where(radiusS[*,0] gt xs[1],c)
    if (c gt 0) then radiusS[i1,0] = -1000 
    i1 = where(radiusG[*,0] lt xg[0],c)
    if (c gt 0) then radiusG[i1,0] = -1000 
    i1 = where(radiusG[*,0] gt xg[1],c)
    if (c gt 0) then radiusG[i1,0] = -1000 

    radS = radiusS[*,0,0]
    radG = radiusG[*,0,0]

    starsU = fltarr(n2)
    starsL = fltarr(n2)
    gcU = fltarr(n2)
    gcL = fltarr(n2)

    for j=0,n2-1 do begin
        starsU[j] = max(radiusS[j,1,*])
        gcU[j] = max(radiusG[j,1,*])
        starsL[j] = min(radiusS[j,1,*])
        gcL[j] = min(radiusG[j,1,*])
    endfor

    if (scale ne 0) then begin
        radkpc1 = (float(x1) * float(scale)) / 1000.0
        radkpc2 = (float(x2) * float(scale)) / 1000.0
    endif

    set_plot,'ps'
    device,file=galaxy+'_intmom.ps',/color
    loadct,0,/silent
    if (scale gt 0) then begin
        plot,radS,stars,/nodata,/xlog,xrange=[x1,x2],xstyle=9,xtitle='!3Radius (arcsec)',$
          ytitle='!7r!3!dr!n/!7r!3!dt!3',yrange=[y1,y2],/ystyle,xthick=4,ythick=4,$
          charthick=4,charsize=1.1,position=[0.13,0.13,0.93,0.89],font=-1
    endif
    if (scale eq 0) then begin
        plot,radS,stars,/nodata,/xlog,xrange=[x1,x2],xstyle=1,xtitle='!3Radius (arcsec)',$
          ytitle='!7r!3!dr!n/!7r!3!dt!3',yrange=[y1,y2],/ystyle,xthick=4,ythick=4,$
          charthick=4,charsize=1.0,position=[0.13,0.13,0.93,0.89],title=galaxy,font=-1
    endif
    
    oplot,[x1,x2],[1.0,1.0],thick=3,color=150
    loadct,4,/silent
    oplot,radS[0:*],stars[0:*],color=150,thick=10
    oplot,radS[0:*],starsU[0:*],color=150,thick=5,linestyle=2
    oplot,radS[0:*],starsL[0:*],color=150,thick=5,linestyle=2
    oplot,radG[0:*],gc[0:*],color=60,thick=10
    oplot,radG[0:*],gcU[0:*],color=60,thick=5,linestyle=2
    oplot,radG[0:*],gcL[0:*],color=60,thick=5,linestyle=2

    if (scale ne 0) then begin
        loadct,0
        axis,xaxis=1,xtitle='!3R (kpc)',xrange=[radkpc1,radkpc2],/save,$
          charthick=4,charsize=1.1,xstyle=1,xthick=4,font=-1
    endif

    device,/close_file
    set_plot,'x'
endif

stop
END


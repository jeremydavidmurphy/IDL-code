;This is a modification of radplot.pro. There is little change from
;the original expect now it runs like a function rather than a procedure.

;The purpose of this file is to pick off the bad fibers in the data
;due to foreground stars, etc. It allows for an iterative procedure
;to eliminate fibers one at a time. The output is the frame_nam.R
;file, which has 6 columns:

;fiber #   Radius   deltaRA    deltaDec    Flux (normalized)   Flux (median)

;THE "NAME.R" FILES ARE WRITTEN WITHIN THE ROUTINE, SO THIS IS IN SOME
;SENSE A FAUX FUNCTION IN THAT IT RETURNS JUST THE NAME GIVEN TO THE
;.R FILE. This is necessary as the fibers that are being kept are
;continually being written and reread back at each iteration of the rejection.

;******************************************************************
function radplotf, r1, f1, name1, torm, SKIP=skipreject
;******************************************************************
; R1 is the name of the name.rad file
; F1 is the name of the data file
; NAME1 is the name of the .R output file.
; TORM is either 't' or 'm', depending on whether you want to plot the
; total for each fiber or the median.
; SKIP: Allows for the rejection of fiber step to be skipped.

;******************************************************************
;The os is the radial offset, in arcseconds, used to align the fiber
;values with their numbers in the display. THIS IS JUST FOR PLOTTING PURPOSES!
os = 10
;******************************************************************

;plotscreen = 'no' ;will suppress the radial plots. good when the galaxy is small on sky as these plots give little information.
plotscreen = 'yes'

readcol,r1,format='I,F,F,F',fibers,radius,dRA,dDec

nn1 = n_elements(f1[*,0]) ;the number of wavelength elements
nn2 = n_elements(f1[0,*]) ;the number of fibers
flux = dblarr(nn2)

check = n_elements(where(f1 eq -666))
if(check ne 1) then f1(where(f1 eq -666)) = 0.0

;The totals are determined and normalized to one
if (torm eq 't') then begin
    FOR j=0,nn2-1 DO BEGIN
        flux[j] = total(f1[*,j])
    ENDFOR
    fluxn = flux/max(flux)
    whytitle = 'Normalized Flux (based on total Fiber counts)'
endif
if (torm eq 'm') then begin
    for j=0,nn2-1 do begin
        flux[j] = median(f1[*,j],/even)
    endfor
    fluxn = flux/max(flux)
    whytitle = 'Normalized Flux (based on median fiber counts)'
endif

;Here, the dead fibers are tossed out of the array
badindex = where(fluxn eq 0.0) 
if (badindex[0] ne -1) then begin
    radius[badindex] = -1.0
    flux[badindex] = -1.0
    fluxn[badindex] = -1.0
    dRA[badindex] = -1.0
    dDec[badindex] = -1.0
    fibers[badindex] = -1.0
endif

outfib = transpose(fibers)
outrad = transpose(radius)
outflux = transpose(flux)
outfluxn = transpose(fluxn)
outdRA = transpose(dRA)
outdDec = transpose(dDec)

form1 = '(i4,3x,f9.3,2x,f9.3,3x,f9.3,4x,f10.7,3x,f11.5)'
form2 = '(a4,3x,a9,a11,a12,a14,a14)'
name = ['Fib#','Radius','dRA','dDec','Nflux','flux']
openw, 6, name1
printf,6,name,format=form2
for j=0,n_elements(outfib)-1 do begin
printf,6,outfib[j],outrad[j],outdRA[j],outdDec[j],outfluxn[j],$
  outflux[j],format=form1
endfor
free_lun, 6

ans='c'

if (plotscreen eq 'yes' and skipreject eq 'no') then begin
    window,5,retain=2,xsize=360,ysize=400,xpos=860,ypos=50
    window,4,retain=2,xsize=360,ysize=400,xpos=480,ypos=50
    window,3,retain=2,xsize=360,ysize=400,xpos=100,ypos=50
endif

REPEAT BEGIN
    readcol,name1,format='I,F,F,F,F,F',skipline=1,fiber,radius,dRA,dDec,fluxn,flux,silent=1
    xmin = min(radius[where(radius ne -1.0)])-15
    xmax = max(radius[where(radius ne -1.0)])+15
    yup = max(fluxn[where(fluxn ne -1.0)])+(max(fluxn[where(fluxn ne -1.0)])*0.1)
    ydown = min(fluxn[where(fluxn ne -1.0)])-(min(fluxn[where(fluxn ne -1.0)])*0.1)
    if (ydown lt 0.00001) then ydown = 0.00001

    if (plotscreen eq 'yes' and skipreject eq 'no') then begin
        wset,5
        loadct,0,/silent
        plot,radius,fluxn,psym=2,title=file,xtitle='Radius (arcsec)',xrange=[xmin,xmax],xstyle=1,$
          /nodata,ytitle=name1,charsize=1.3,/ylog,yrange=[ydown,yup],ystyle=1
;        plot,radius,fluxn,psym=2,title=file,xtitle='Radius (arcsec)',xrange=[xmin,xmax],xstyle=1,$
;          yrange=[ydown,yup],ystyle=1,/nodata,ytitle=name1,charsize=1.3,/ylog
        loadct,4,/silent
        oplot,radius,fluxn,psym=2,symsize=0.5,color=255
        loadct,0
        for j=2,n_elements(fluxn)-3,3 do xyouts,radius[j]-os,fluxn[j],fiber[j],charsize=1.5
        
        wset,4
        loadct,0,/silent
        plot,radius,fluxn,psym=2,title=file,xtitle='Radius (arcsec)',xrange=[xmin,xmax],xstyle=1,$
          /nodata,ytitle=name1,charsize=1.3,/ylog,yrange=[ydown,yup],ystyle=1
;        plot,radius,fluxn,psym=2,title=file,xtitle='Radius (arcsec)',xrange=[xmin,xmax],xstyle=1,$
;          yrange=[y,down,yup],ystyle=1,/nodata,ytitle=name1,charsize=1.3,/ylog
        loadct,4,/silent
        oplot,radius,fluxn,psym=2,symsize=0.5,color=255
        loadct,0,/silent
        for j=1,n_elements(fluxn)-2,3 do xyouts,radius[j]-os,fluxn[j],fiber[j],charsize=1.5
        
        wset,3
        loadct,0,/silent
        plot,radius,fluxn,psym=2,title=file,xtitle='Radius (arcsec)',xrange=[xmin,xmax],xstyle=1,$
          /nodata,ytitle=name1,charsize=1.3,/ylog,yrange=[ydown,yup],ystyle=1
;        plot,radius,fluxn,psym=2,title=file,xtitle='Radius (arcsec)',xrange=[xmin,xmax],xstyle=1,$
;          yrange=[ydown,yup],ystyle=1,/nodata,ytitle=name1,charsize=1.3,/ylog
        loadct,4,/silent
        oplot,radius,fluxn,psym=2,symsize=0.5,color=255
        loadct,0,/silent
        for j=0,n_elements(fluxn)-1,3 do xyouts,radius[j]-os,fluxn[j],fiber[j],charsize=1.5
    endif

    if (skipreject eq 'yes') then begin ;this now allows you to quickly run through this, keeping all fibers.
        ans = 'd'
        wait,1.0
        goto,jumpend
    endif

    print,'Enter a fiber to reject'
    print,'(Enter "d" if you are done.)'
    read,ans
    IF (ans NE 'd') AND (ans NE '') THEN BEGIN
        killfib = uint(ans)
        indy = where(outfib EQ killfib)
        WHILE (indy EQ -1) DO BEGIN
            print,'That fiber is no longer in the mix!  Try another.'
            read,ans
            if (ans eq 'd') then goto,jumpend else killfib = uint(ans)
            indy = where(outfib EQ killfib)
        ENDWHILE
        outfib[indy] = -1
        outrad[indy] = -1
        outfluxn[indy] = -1
        outflux[indy] = -1
        outdRA[indy] = -1
        outdDec[indy] = -1

        openw, 6, name1
        printf,6,name,format=form2
        for j=0,n_elements(outfib)-1 do begin
            printf,6,outfib[j],outrad[j],outdRA[j],outdDec[j],outfluxn[j],outflux[j],format=form1
        endfor
        free_lun, 6
    ENDIF

jumpend:
ENDREP UNTIL (ans EQ 'd')

if (plotscreen eq 'yes' and skipreject eq 'no') then begin
    window,5,retain=2,xsize=420,ysize=450,xpos=860,ypos=50
    
    loadct,0,/silent
    plot,outrad,outfluxn,psym=2,title=file,xtitle='Radius (arcsec)',xrange=[xmin,xmax],xstyle=1,$
      yrange=[ydown,yup],ystyle=1,ytitle='Normalized Flux',charsize=1.2,/ylog
    loadct,4,/silent
    xyouts,0.45,0.85,'This is the final profile.',/normal,charsize=1.5,color=180
    
    print,'The next ENTER deletes the plot...'
    pause
    wdelete,5
endif

return,name1
END

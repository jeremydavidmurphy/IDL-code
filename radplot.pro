;This purpose of this file is to pick off the bad fibers in the data
;due to foreground stars, etc. It allows for an iterative procedure
;to eliminate fibers one at a time. The output is the frame_nam.R
;file, which is three columns:

;  fiber #      sorted radial position       estimate of flux

;The code needs the name.rad output from radius.pro in the calling directory.
;This works for either fiber bundle.

PRO radplot, radius_file, frame_num, Percent=p
; Ex: 'NGC1600a', 'jm0056'

ans=''
IF (n_elements(frame_num) EQ 0) THEN BEGIN
    print,'Enter a collapsed .fits file (w/o the pefsmc.fits):'
    read,ans
    f1 = ans+'pefsmc.fits'
ENDIF ELSE f1 = frame_num+'pefsmc.fits'

if (n_elements(p) eq 0) then p = 0.5

f2 = radius_file+'.rad' ;This is the radial array (unsorted) that comes from radius.pro

readcol,f2,format='I,F',fibers,radius

data = readfits(f1)
n1 = n_elements(data[*,0]) ;the number of wavelength elements
n2 = n_elements(data[0,*]) ;the number of fibers
flux = dblarr(n2)

check = min(data)
IF(check EQ -666) THEN data(Where(data EQ -666))=0.0

;IF (n1 EQ 1024) THEN BEGIN
;    r1 = 150
;    r2 = 800
;ENDIF
;IF (n1 EQ 2048) THEN BEGIN
;    r1 = 200
;    r2 = 1650
;ENDIF

;A SWITCH BETWEEN SELECTING THE 'TOTAL' AND THE MEDIAN
;(added on Oct 2,2008)
jumpback:
anstm=''
print,'Use the total or the median for plotting and rejection? (t or m):'
read,anstm
if (anstm eq 't') then begin
    FOR j=0,n2-1 DO BEGIN
        flux[j] = total(data[*,j])
        print,'The rough total for fiber '+strn(j+1)+' is: '+strn(flux[j])
    ENDFOR
    whytitle = 'Total Fiber Counts'
endif
if (anstm eq 'm') then begin
    for j=0,n2-1 do begin
        flux[j] = median(data[*,j],/even)
        print,'The median for fiber '+strn(j+1)+' is: '+strn(flux[j])
    endfor
    whytitle = 'Median of Fiber'
endif
if (anstm ne 't') and (anstm ne 'm') then begin
    print,'Not an option. Try again'
    goto, jumpback
endif
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;Here, the dead fibers are tossed out of the array
bigindex = where(flux GT 0) 
flux = flux[bigindex] 
fibers = fibers[bigindex]
radius = radius[bigindex]
n = n_elements(bigindex)
xlow = min(radius)-10
xhigh = max(radius)+10
cuthigh = max(flux)

;Now the fibers are sorted by radial position.
index = bsort(radius)
sfibers = fibers[index]
sradius = radius[index]
sflux = flux[index]

cutman = fltarr(1)
cutlow = min(flux)

REPEAT BEGIN
    rejected = fltarr(4,n)
    pp=1.+p
    ctr=-1
;    IF (cutman EQ 0) THEN BEGIN
        FOR j=2,n-3 DO BEGIN
            av = (sflux[j-2]+sflux[j-1]+sflux[j+1]+sflux[j+2])/4.
            ct = pp*av
            IF (sflux[j] GT ct) THEN BEGIN ;The potential stars are pinned down here...
                ctr = ctr+1
                rejected[0,ctr] = sfibers[j]
                rejected[1,ctr] = sradius[j]
                rejected[2,ctr] = sflux[j]
                rejected[3,ctr] = j ;The index for the positions in the SORTED arrays is kept.
            ENDIF
        ENDFOR
;    ENDIF ELSE BEGIN
;        FOR j = 0,n_elements(sflux)-1 DO BEGIN
;            indo = where(sflux GE cutman)
;            rejected[0,*] = sfibers[indo]
;            rejected[1,*] = sradius[indo]
;            rejected[2,*] = sflux[indo]
;            rejected[3,*] = indo
;        ENDELSE

    nn = n_elements(where(rejected[0,*]) NE 0.)
    IF (nn NE 0) THEN rejected = rejected[*,0:nn-1]
    IF (nn EQ 0) THEN print,'No fibers were rejected!?'
    indy = uint(rejected[2,*]) ;The index to the SORTED ARRAYS for the potential rejects.

    window,2, retain=2
    device,decomposed=0
    loadct,0
    IF (n_elements(cuthigh) NE 0) THEN BEGIN 
        plot,sradius,sflux,charsize=1.3,psym=7,xtitle='Radius (arcsec)',$
          xrange=[xlow,xhigh],xstyle=1,yrange=[cutlow,cuthigh],ytitle=whytitle,$
          title='Fiber Flux for frame '+frame_num+'  ('+radius_file+')'
    ENDIF ELSE BEGIN
        plot,sradius,sflux,charsize=1.3,psym=7,xtitle='Radius (arcsec)',$
          xrange=[xlow,xhigh],xstyle=1,ytitle=whytitle,$
          title='Fiber Flux for frame '+frame_num+'  ('+radius_file+')',/ynozero
    ENDELSE

    loadct,4
    xyouts,0.15,0.13,'Marked fibers will be rejected from the final list.',charsize=1.5,/normal
    
    loadct,4
    IF (nn NE 0) THEN BEGIN
        FOR j=0,nn-1 DO plots,rejected[1,j],rejected[2,j],psym=7,color=60
        FOR j=0,nn-1 DO xyouts,rejected[1,j],rejected[2,j],strn(uint(rejected[0,j])),charsize=1.5,color=110
    ENDIF
    
    print,'Enter P, Y, M or D:'
    print,'P: percent change, Y: limit y-axis, M: manual cut, D:done'
    read,ans
    IF (ans EQ 'p') THEN BEGIN
        print,'The current deviation is set at: '+strn(p)
        print,'Enter the new deviation value:'
        read,p
    ENDIF
    IF (ans EQ 'y') THEN BEGIN
        print,'Enter a new Y upper limit:'
        read,cuthigh
    ENDIF
    IF (ans EQ 'm') THEN BEGIN
        print,'Enter a manual cut:'
        read,cutman
    ENDIF

ENDREP UNTIL (ans EQ 'd')

;The stars are now actually tossed from the fiber list

outfib = fltarr(n)
outflux = fltarr(n)
outrad = fltarr(n)
ctr=0
FOR j=0,n-1 DO BEGIN
    IF(j NE rejected[3,ctr]) THEN BEGIN
        outfib[j] = sfibers[j]
        outflux[j] = sflux[j]
        outrad[j] = sradius[j]
    ENDIF ELSE BEGIN
        ctr = ctr+1
        IF (ctr EQ nn) THEN ctr = nn-1
        outfib[j] = -1
        outflux[j] = -1
        outrad[j] = -1
    ENDELSE
ENDFOR

outfib = outfib[where(outfib NE -1)] ;Outfib is now sorted fibers WITHOUT stars.
outfib = transpose(outfib)
outfib = uint(outfib)
outflux= outflux[where(outflux NE -1)]
outflux = transpose(outflux)
outrad = outrad[where(outrad NE -1)]
outrad = transpose(outrad)

openw, 6, frame_num+'.R'
FOR j=0,n_elements(outfib)-1 DO printf, 6, outfib[j], outrad[j], outflux[j]
free_lun, 6

;print,'Save the REJECTED file?'
;read,ans
;IF (ans Eq 'y') THEN BEGIN
;    outr = frame_num+'.rejected'
;    openw,5,outr
;    printf,5,rejected
;    free_lun,5
;ENDIF

ans='c'
file = frame_num+'.R'
REPEAT BEGIN
    readcol,file,fiber,radius,flux
    fiber = uint(fiber)
    xmin = min(radius)-15
    xmax = max(radius)+15
    ymax = max(flux)+(0.10*max(flux))
    ymin = min(flux)-(0.10*min(flux))

    loadct,0
    plot,radius,flux,psym=2,title=file,xtitle='Radius (arcsec)',xrange=[xmin,xmax],xstyle=1,$
      yrange=[ymin,ymax],ystyle=1
    xyouts,0.6,0.85,'These are the final fibers.',/normal,charsize=1.5
    loadct,4
    xyouts,radius-5,flux,fiber,charsize=1.2,color=110
    print,'Enter a fiber to reject'
    print,'(Enter "d" if you are done.)'
    read,ans
    IF (ans NE 'd') AND (ans NE '') THEN BEGIN
        killfib = uint(ans)
        indy = where(outfib EQ killfib)
        WHILE (indy EQ -1) DO BEGIN
            print,'That fiber is no longer in the mix!  Try another.'
            read,killfib
            indy = where(outfib EQ killfib)
        ENDWHILE
        outfib[indy] = -1
        outrad[indy] = -1
        outflux[indy] = -1
        outfib = outfib[where(outfib NE -1)]
        outflux= outflux[where(outflux NE -1)]
        outrad = outrad[where(outrad NE -1)]
        openw, 6, frame_num+'.R'
        FOR j=0,n_elements(outfib)-1 DO printf, 6, outfib[j], outrad[j], outflux[j]
        free_lun, 6
    ENDIF

    ENDREP UNTIL (ans EQ 'd')

ans='y'
;print,'Save a hardcopy of the final profile? (y or n):'
;read,ans
if (ans eq 'y') then begin
    set_plot,'ps'
    device,file=frame_num+'.ps',/color
    loadct,0
    plot,radius,flux,psym=2,title='Flux for all fibers KEPT for frame '+frame_num,xtitle='Radius (arcsec)',$
      ytitle=whytitle,charsize=1.2,xrange=[xmin,xmax],xstyle=1,yrange=[ymin,ymax],ystyle=1,$
      xthick=2,ythick=2,charthick=2
    xyouts,0.2,0.2,'These are the final fibers.',/normal
    loadct,4
    xyouts,radius-5,flux,fiber,charsize=0.8,color=110
    device,/close
    set_plot,'x'
    print,''
    print,'The plot was saved as '+frame_num+'.ps'
    print,''
endif else print,'No plot will be saved.'

print,'Next ENTER delete the plot:'
pause
wdelete,2

END

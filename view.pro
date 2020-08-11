; This routine is used to visually inspect which fibers are being
; combined into a given bin.
pro view, pa

; the PA (position angle) is input in degrees. It's then converted to
; radians

pa = pa/57.29578

readcol,'bin.list',format='a',bins

n1 = n_elements(bins)

window,2,retain=2
device,decomposed=0

for k=0,n1-1 do begin
    readcol,bins[k],format='a',fibers
    
    n1 = n_elements(fibers)
    fib = intarr(n1)
    let = strarr(n1)
    
    for j=0,n1-1 do begin
        f = strsplit(fibers[j],'_',/extract)
        fib[j] = uint(f[0])
        let[j] = f[1]
    endfor
    
    readcol,'letter.coord',format='a,a',letters,coords
    
    ia = where(let eq 'a')
    ib = where(let eq 'b')
    ic = where(let eq 'c')
    id = where(let eq 'd')
    ie = where(let eq 'e')
    
    if (ia[0] ne -1) then fa = fib[ia]-1
    if (ib[0] ne -1) then fb = fib[ib]-1
    if (ic[0] ne -1) then fc = fib[ic]-1
    if (id[0] ne -1) then fd = fib[id]-1
    if (ie[0] ne -1) then fe = fib[ie]-1
    
    allra = fltarr(1)
    alldec = fltarr(1)
    
    if (ia[0] ne -1) then begin
        readcol,coords[0],format='i,x,f,f',fiber,dra,ddec,silent=1
        Adra = dra[fa]
        Addec = ddec[fa]
        allra = [allra,Adra]
        alldec = [alldec,Addec]
        sa = 1
    endif else sa = 0
    
    if (ib[0] ne -1) then begin
        readcol,coords[1],format='i,x,f,f',fiber,dra,ddec,silent=1
        Bdra = dra[fb]
        Bddec = ddec[fb]
        allra = [allra,Bdra]
        alldec = [alldec,Bddec]
        sb = 1
    endif else sb = 0
    
    if (ic[0] ne -1) then begin
        readcol,coords[2],format='i,x,f,f',fiber,dra,ddec,silent=1
        Cdra = dra[fc]
        Cddec = ddec[fc]
        allra = [allra,Cdra]
        alldec = [alldec,Cddec]
        sc = 1
    endif else sc = 0
    
    if (id[0] ne -1) then begin
        readcol,coords[3],format='i,x,f,f',fiber,dra,ddec,silent=1
        Ddra = dra[fd]
        Dddec = ddec[fd]
        allra = [allra,Ddra]
        alldec = [alldec,Dddec]
        sd = 1
    endif else sd = 0
    
    if (ie[0] ne -1) then begin
        readcol,coords[4],format='i,x,f,f',fiber,dra,ddec,silent=1
        Edra = dra[fe]
        Eddec = ddec[fe]
        allra = [allra,Edra]
        alldec = [alldec,Eddec]
        se = 1
    endif else se = 0
    
    allra = allra[1:*]
    alldec = alldec[1:*]
    
    yup = max(alldec)+5
    ydown = min(alldec)-5
    xup = max(allra)+5
    xdown = min(allra)-5
    
    
    xup = 150.0
    xdown = -150.0
    yup = 250.0
    ydown= -215
    
    pax1 = [yup*tan(pa),0]
    pay1 = [yup,0]
    pax2 = [0,ydown*tan(pa)]
    pay2 = [0,ydown]
    
    
    

    loadct,0
    plot,allra,alldec,psym=sym(1),xtitle='RA Offsets (arcsec)',$
      ytitle='Dec Offsets (arcsec)',charsize=1.2,$
      title='Fiber positions for '+bins[k],/nodata,$
      xrange=[xdown,xup],yrange=[ydown,yup],xstyle=1,ystyle=1
    plots,0,0,psym=2,symsize=3
    oplot,pax1,pay1,thick=2
    oplot,pax2,pay2,thick=2
    
    loadct,4
    if(ia[0] ne -1) then oplot,Adra,Addec,psym=sym(1),color=60
    if(ib[0] ne -1) then oplot,Bdra,Bddec,psym=sym(1),color=110
    if(ic[0] ne -1) then oplot,Cdra,Cddec,psym=sym(1),color=150
    if(id[0] ne -1) then oplot,Ddra,Dddec,psym=sym(1),color=180
    if(ie[0] ne -1) then oplot,Edra,Eddec,psym=sym(1),color=255
    
    xyouts,0.15,0.15,'Pointing A',color=60,/normal,charsize=1.5
    xyouts,0.15,0.185,'Pointing B',color=110,/normal,charsize=1.5
    xyouts,0.15,0.22,'Pointing C',color=150,/normal,charsize=1.5
    xyouts,0.15,0.255,'Pointing D',color=180,/normal,charsize=1.5
    xyouts,0.15,0.29,'Pointing E',color=255,/normal,charsize=1.5
    pause
endfor

print,'Next enter deletes the plot...'
pause

;set_plot,'ps'
;device,file='fiberview.ps',/color
;loadct,0
;plot,allra,alldec,psym=sym(1),xtitle='RA Offsets (arcsec)',$
;  ytitle='Dec Offsets (arcsec)',charsize=1.2,$
;  title='Fiber positions for '+bins[k],/nodata,$
;  xrange=[xdown,xup],yrange=[ydown,yup],xstyle=1,ystyle=1,$
;  xthick=3,ythick=3,charthick=3
;plots,0,0,psym=2,symsize=3,thick=3
;oplot,pax1,pay1,thick=2
;oplot,pax2,pay2,thick=2

;loadct,4
;if(ia[0] ne -1) then oplot,Adra,Addec,psym=sym(1),color=60
;if(ib[0] ne -1) then oplot,Bdra,Bddec,psym=sym(1),color=110
;if(ic[0] ne -1) then oplot,Cdra,Cddec,psym=sym(1),color=150
;if(id[0] ne -1) then oplot,Ddra,Dddec,psym=sym(1),color=180
;loadct,0
;if(ie[0] ne -1) then oplot,Edra,Eddec,psym=sym(1)
;loadct,4
;xyouts,0.18,0.15,'Pointing A',color=60,/normal,charsize=1.2,charthick=3
;xyouts,0.18,0.185,'Pointing B',color=110,/normal,charsize=1.2,charthick=3
;xyouts,0.18,0.22,'Pointing C',color=150,/normal,charsize=1.2,charthick=3
;xyouts,0.18,0.255,'Pointing D',color=180,/normal,charsize=1.2,charthick=3
;loadct,0
;xyouts,0.18,0.29,'Pointing E',/normal,charsize=1.2,charthick=3
;device,/close_file
;set_plot,'x'

wdelete,2
stop
end

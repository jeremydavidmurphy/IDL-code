PRO Ploghalo

type = ''
print,'What type of halo (nfw/log)?'
read,type

if (type eq 'log') then begin
    vc = fltarr(1)
    print,'Enter a circular velocity (km/sec):'
    read,vc
    
    rs = fltarr(1)
    print,'Enter a scale radius (kpc):'
    read,rc
    
    rad = findgen(2000)-0.9
    halo = fltarr(2000)
    
    for j=0,1999 do begin
        halo[j] = vc^2 * (3 * rc^2 + rad[j]^2) / (rc^2 + rad[j]^2)^2
    endfor
endif

if (type eq 'nfw') then begin
    rs = fltarr(1)
    print,'Enter a scale radius (kpc):'
    read,rc
    
    rad = findgen(2000)-0.9
    halo = fltarr(2000)
    
    for j=0,1999 do begin
        halo[j] = 1.0 / ((rad[j]/rs) * (1 + (rad[j]/rs))^2)
    endfor
endif
    
window,1,retain=2,ypos=50
device,decomposed=0
loadct,0,/silent
plot,rad,halo,/xlog,xrange=[0.1,2000],/xstyle

ans = ''
print,'Overplot another halo (y/n)?:'
read,ans

if (ans eq 'n') then begin
wdelete
stop
endif 
if (ans eq 'y') then begin
    colors = [60,80,100,120,140,160,180,200,220,240,255]
    cntr = 0
    repeat begin
        if (type eq 'log') then begin
            print,'Enter a circular velocity (km/sec):'
            read,vc
            print,'Enter a scale radius (kpc):'
            read,rc
            for j=0,1999 do begin
                halo[j] = vc^2 * (3 * rc^2 + rad[j]^2) / (rc^2 + rad[j]^2)^2
            endfor
        endif
        if (type eq 'nfw') then begin
            print,'Enter a scale radius (kpc):'
            read,rc
            for j=0,1999 do begin
                halo[j] = 1.0 / ((rad[j]/rs) * (1 + (rad[j]/rs))^2)
            endfor
        endif
        loadct,4,/silent
        oplot,rad,halo,color=colors[cntr]
        print,'Plot another? (y/n):'
        read,ans
        cntr = cntr + 1
        if (cntr eq 11) then begin
            print,'That is too many! I am full!'
            pause
            wdelete
            stop
        endif
    endrep until (ans eq 'n')
    print,'The next ENTER deletes the plot...'
    pause
    wdelete
endif

stop
END

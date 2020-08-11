PRO Psfit

; This routine plots an LOSVD profile, with errors. It works from a
; list, called 'list.sfit' and allows for several to be plotted to
; screen in a row, with pause between to visually inspect.

readcol,'list.sfit',format='a',files

n0 = n_elements(files)

ans = ''
print,'Plot with errors? (y/n)'
read,ans

if (ans eq 'y') then begin
    ans2 = ''
    print,'Plot the output to screen? (y/n)'
    read,ans2
    data = fltarr(4,1000,n0)
    window,0,retain=2,
    device,decomposed=0
    for j=0,n0-1 do begin
        readcol,files[j],format='f,f,f,f',vel,i,id,iu
        data[0,*,j] = vel
        data[1,*,j] = i
        data[2,*,j] = id
        data[3,*,j] = iu
        loadct,0
        if (ans2 eq 'y') then begin
            plot,vel,i,xtitle='Velocity (km/sec)',thick=2,title=files[j]
            loadct,4
            oplot,vel,id,linestyle=2,color=150,thick=2
            pause
        endif
    endfor
endif

END

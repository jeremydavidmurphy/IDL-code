PRO Pprofiles, list

; This routine is used to generate plots of the X and Y profiles of a
; list of fiber spots.

readcol,list,f='a',files
n0 = n_elements(files)
colors = intarr(n0)
for j=0,n0-1 do colors[j] = (j * floor(255.0/n0))


for j=0,n0-1 do begin
    frame = readfits(files[j],/silent)
    index = where(frame gt 5000.0)
    temp = Fit_Ellipse(index,CENTER=cent)
    x = round(cent[0])
    y = round(cent[1])
    sliceX = frame[*,y-10:y+10]
    sliceY = transpose(frame[x-10:x+10,*])
    profileX = median(sliceX,dim=2)
    profileY = median(sliceY,dim=2)
    if (j eq 0) then begin
        xsize = n_elements(profileX)
        ysize = n_elements(profileY)
        window,0,retain=2,xsize=500,ysize=500
        device,decomposed=0
        loadct,0
        plot,profileX,xtitle='Pixels',ytitle='CCD Counts',xrange=[0,xsize],$
          /xstyle,title='X-PROFILE'
        loadct,33
        window,2,retain=2,xsize=500,ysize=500
        device,decomposed=0
        loadct,0
        plot,profiley,xtitle='Pixels',ytitle='CCD Counts',xrange=[0,xsize],$
          /xstyle,title='Y-PROFILE'
        loadct,33

        set_plot,'ps'
        device,file='profiles.ps',/color
        !p.multi = [0,2,1,0,1]
        loadct,0
        plot,profileX,xtitle='Pixels',ytitle='CCD Counts',xrange=[0,xsize],$
          /xstyle,title='X-PROFILE',position=[0.1,0.1,0.45,0.9]
        plot,profiley,xtitle='Pixels',ytitle='CCD Counts',xrange=[0,xsize],$
          /xstyle,title='Y-PROFILE',position=[0.45,0.1,0.9,0.9]
        loadct,33
        set_plot,'x'
    endif else begin
        wset,0
        oplot,profileX,color=colors[j]
        wset,2
        oplot,profileY,color=colors[j]
        set_plot,'ps'
        oplot,profileX,color=colors[j]
        oplot,profileY,color=colors[j]
        set_plot,'x'
    endelse

endfor

stop
END

pro plotfield, frame

; this routine takes a collapsed frame and plots the fiber flux. you
; know. the poor man's photometry
;******************************************************************

readcol,'/home/grad78/murphy/finder_code_delta/scripts/IFUcen_vp2.txt',silent=1,f='x,f,f',dra,ddec
;readcol,'IFUcent_vp2.txt',silent=1,f='x,f,f',dra,ddec
n0 = n_elements(dra)

data = readfits(frame,h,/silent)

if (n_elements(data[0,*]) gt 300) then begin
    print,'Collapsing your file...'
    mask = ''
    wgtf = ''
    print,'Enter the name of the weight file:'
    read,wgtf
    print,'Enter the name of the mask file:'
    read,mask
    wgt = readfits(wgtf,/silent)
    datac = wgtcollapsef(data,wgt,mask,5)
    data = datac
endif

mdata = median(data,dim=1)
ilow = where(mdata lt 0.0,ctlow)
if (ctlow gt 0) then print,'YOU HAVE '+strn(ctlow)+' fibers below zero!!!'
ibad = where(mdata lt -660,badcount)
if (badcount gt 0) then begin
    mdata[ibad] = 0.0
    dra[ibad] = 1000.0 ;the bad fibers are just driven out of the plotting range.
    ddec[ibad] = 1000.0
endif
datamax = max(mdata)
pscale = 255.0 / datamax
ans = 'y'
window,0,retain=2,xsize=650,ysize=650
device,decomposed=0
plotposition = aspect(1.0/1.0)
window,2,retain=2,xsize=650,ysize=650
device,decomposed=0
plotposition = aspect(1.0/1.0)
plotposition[0] = plotposition[0]*0.9
plotposition[1] = plotposition[1]*0.9
plotposition[2] = plotposition[2]*1.1
plotposition[3] = plotposition[3]*1.1

repeat begin
    jumpback:
    colors = floor(pscale * mdata)
    iover = where(colors gt 255.0,ctover)
    iunder = where(colors le 255.0,ctunder)

;the fiber # figure...
    wset,2
    loadct,0,/silent
    plot,dra,ddec,psym=sym(1),/nodata,xtitle='RA offset (arcsec)',$
      ytitle='Dec offset (arcsec)',xrange=[-60,60],yrange=[-60,60],/xs,/ys,$
      position=plotposition,charsize=1.5,title='Fiber Number'
    loadct,33,/silent
    if (ctover gt 0) then begin
        pmdata = mdata
        pmdata[iover] = 0.0
        colors = floor(pscale * pmdata)
        for j=0,n0-1 do plots,dra[j],ddec[j],psym=sym(1),color=colors[j],symsize=4
        loadct,0,/silent
        for j=0,ctover-1 do plots,dra[iover[j]],ddec[iover[j]],psym=sym(1),color=0,symsize=4
    endif else $
      for j=0,n0-1 do plots,dra[j],ddec[j],psym=sym(1),color=colors[j],symsize=4
    loadct,0,/silent
    for j=0,n0-1 do xyouts,dra[j]-2.5,ddec[j]-1,strn(j+1),charsize=1.2,color=0

;the intensity figure...
    wset,0
    loadct,0,/silent
    plot,dra,ddec,psym=sym(1),/nodata,xtitle='RA offset (arcsec)',$
      ytitle='Dec offset (arcsec)',xrange=[-60,60],yrange=[-60,60],/xs,/ys,$
      position=plotposition,charsize=1.5,title='Median Fiber Flux'
    loadct,33,/silent
    if (ctover gt 0) then begin
        for j=0,n0-1 do plots,dra[j],ddec[j],psym=sym(1),color=colors[j],symsize=4
        loadct,0,/silent
        for j=0,ctover-1 do plots,dra[iover[j]],ddec[iover[j]],psym=sym(1),color=0,symsize=4
    endif else $
      for j=0,n0-1 do plots,dra[j],ddec[j],psym=sym(1),color=colors[j],symsize=4
    loadct,0,/silent
    for j=0,n0-1 do xyouts,dra[j]-2.5,ddec[j]-1,strn(round(mdata[j])),charsize=1.2,color=0

    print,'Stretch the image? (y/n):'
    read,ans
    if (ans eq 'n') then goto,jump1
    if (ans eq 'y') then begin
        print,'The current stretch value is '+strn(pscale)
        print,'The max of the current image (color) is '+strn(max(mdata[iunder]))+' ('+strn(max(colors[iunder]))+')'
        print,'The min of the current image (color) is '+strn(min(mdata[iunder]))+' ('+strn(min(colors[iunder]))+')'
        print,'Enter a new stretch:'
        read,pscale
        goto,jumpback
    endif
    jump1:
endrep until (ans eq 'n')    

print,'Want a postscript plot of this? (y/n)'
read,ans
if (ans eq 'y') then begin
    jumpback2:
    print,'Fiber #, fiber flux or both? (n/f/b)'
    read,ans
    if (ans eq 'n' or ans eq 'b') then begin
        set_plot,'ps'
        device,file='plotfield_number.ps',/color
        loadct,0,/silent
        plot,dra,ddec,psym=sym(1),/nodata,xtitle='RA offset (arcsec)',$
          ytitle='Dec offset (arcsec)',xrange=[-60,60],yrange=[-60,60],/xs,/ys,$
          position=plotposition,charsize=1.0,title='Fiber Number',$
          xthick=3,ythick=3,charthick=3
        loadct,33,/silent
        if (ctover gt 0) then begin
            for j=0,n0-1 do plots,dra[j],ddec[j],psym=sym(1),color=colors[j],symsize=3
            loadct,0,/silent
            for j=0,ctover-1 do plots,dra[iover[j]],ddec[iover[j]],psym=sym(1),color=0,symsize=3
        endif else $
          for j=0,n0-1 do plots,dra[j],ddec[j],psym=sym(1),color=colors[j],symsize=3
        loadct,0,/silent
        for j=0,n0-1 do xyouts,dra[j]-2.5,ddec[j]-1,strn(j+1),charsize=0.8,color=0
        device,/close_file
        set_plot,'x'
    endif

    if (ans eq 'f' or ans eq 'b') then begin
        set_plot,'ps'
        device,file='plotfield_flux.ps',/color
        loadct,0,/silent
        plot,dra,ddec,psym=sym(1),/nodata,xtitle='RA offset (arcsec)',$
          ytitle='Dec offset (arcsec)',xrange=[-60,60],yrange=[-60,60],/xs,/ys,$
          position=plotposition,charsize=1.0,title='Median Fiber Flux',$
          xthick=3,ythick=3,charthick=3
        loadct,33,/silent
        if (ctover gt 0) then begin
            for j=0,n0-1 do plots,dra[j],ddec[j],psym=sym(1),color=colors[j],symsize=3
            loadct,0,/silent
            for j=0,ctover-1 do plots,dra[iover[j]],ddec[iover[j]],psym=sym(1),color=0,symsize=3
        endif else $
          for j=0,n0-1 do plots,dra[j],ddec[j],psym=sym(1),color=colors[j],symsize=3
        loadct,0,/silent
        for j=0,n0-1 do xyouts,dra[j]-2.5,ddec[j]-1,strn(round(mdata[j])),charsize=0.8,color=0
        device,/close_file
        set_plot,'x'
    endif
    if (ans ne 'n' and ans ne 'f' and ans ne 'b') then begin
        print,'Oops! Try again!'
        goto,jumpback2
    endif
endif

wdelete,0
wdelete,2
stop
END

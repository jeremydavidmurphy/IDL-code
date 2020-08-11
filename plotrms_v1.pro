;This routine is used to plot the RMS output of the rfitlov
;routine. The screen output is piped into a file, then modified to
;include only the rms values.

;THIS ROUTINE IS USED TO LOOK AT ONE RMS.OUT FILE THAT MAY HAVE
;SEVERAL DIFFERENT REGIONS FIT WITHIN IT. THEY NEED TO BE IN ORDERED
;BLOCKS, EACH OF THE SAME SIZE.

;If the name keyword call isn't used, this routine will look for a
;file entitled "rms.out" in the calling directory

;----------------------------------------------------------------
pro plotrms_v1, NAME=rmsfile, BNAME=binnamelist 
;----------------------------------------------------------------

;if the binlist is included then these are included in the plot.

;A FILE NAMED RMS.LIST IS REQUIRED TO EXIST IN THE CALLING
;DIRECTORY. This is a list of the regions within the rms file you're
;interested in plotting. The number of names in this file forms the
;denominator from which the rms input file is split upon.

if (n_elements(binnamelist) ne 0) then begin
    readcol,binnamelist,format='a',blist
    swchname = 'on'
endif else swchname = 'off'

if (n_elements(rmsfile) eq 0) then begin
    readcol,'rms.out',format='d,x,x,x',rms
endif else readcol,rmsfile,format='d,x,x,x',rms

readcol,'rms.list',format='a',names
div = n_elements(names)

if (div le 4) then begin
    color = [60,110,180,150]
    swch = 'on'
endif else swch = 'off'

if (div gt 4) and (div le 7) then cs=30
if (div ge 8) then cs = 20
symbol = [1,2,4,5,6,7,1,2,4,5,6,7,1,2,4,5,6,7]

n0 = n_elements(rms)
step = n0/div
n1 = 0
n2 = step
x0 = findgen(n2)+1
xs = x0
xar = intarr(n2,div)
rmsarr = dblarr(1,step,div)
rmsmin = strarr(3,step)

for j=0,div-1 do begin
    rmsarr[0,*,j] = rms[n1:n2-1]
    n1 = n1+step
    n2 = n2+step
    xar[*,j] = xs
    xs = xs + step
endfor

for j=0,step-1 do begin
    mn = min(rmsarr[0,j,*])
    mni = where(rmsarr[0,j,*] eq mn)
    mn = string(mn)
    nm = names[mni]
    rmsmin[0,j] = blist[j]
    rmsmin[1,j] = nm
    rmsmin[2,j] = mn
endfor

f1 = '(A10,2x,A10,A16)'
openw,5,'rmsmin.txt'
for j=0,step-1 do begin
    printf,5,rmsmin[0,j],rmsmin[1,j],rmsmin[2,j],format=f1
endfor
free_lun,5

xupo = max(x0)+2
xdown = -2
xups = max(xar)+2
yup = max(rmsarr)
ydown = min(rmsarr)-(min(rmsarr)*0.2)

window,0,retain=2
device,decomposed=0
loadct,0

plot,x0,rmsarr[0,*,0],psym=3,xtitle='Files',ytitle='RMS',$
  yrange=[ydown,yup],charsize=1.5,xrange=[xdown,xupo],xstyle=1,$
  title='RMS plot taken from the output of FITLOV',/nodata,$
  ystyle=1,/ynozero

loadct,4
y = 0.86
if (swch eq 'on') then begin
    for j=0,div-1 do begin
        oplot,x0,rmsarr[0,*,j],psym=symbol[j],color=color[j]
        xyouts,0.02,y,names[j],charthick=1.5,charsize=1.7,$
          color=color[j],/normal
        y = y-0.04
        pause
    endfor
    loadct,0
endif else begin
    for j=0,div-1 do begin
        oplot,x0,rmsarr[0,*,j],psym=symbol[j],color=60+(j*cs)
        xyouts,0.02,y,names[j],charthick=1.5,charsize=1.7,color=60+(j*cs),/normal
        y = y-0.04
        pause
    endfor
    loadct,0
endelse

if (swchname eq 'on') then begin
    loadct,0
    for j=0,step-1 do begin
        xyouts,x0[j],ydown+0.0002,blist[j],orientation=90,charsize=1.0
    endfor
endif

ans=''
print,'Save the plot? (y or n):'
read,ans

if (ans eq 'y') then begin
    outname = ''
    print,'Name the output file (w/o the .ps):'
    read,outname
    outname = outname+'.ps'
    set_plot,'ps'
    device,file=outname,/color
    plot,x0,rmsarr[0,*,0],psym=3,xtitle='Files',ytitle='RMS',$
      yrange=[ydown,yup],charsize=1.2,xthick=2,ythick=2,charthick=2,$
      title='RMS plot taken from the output of FITLOV',/nodata,$
      xrange=[xdown,xupo],xstyle=1,ystyle=1
    
    loadct,4
    y = 0.86
    if (swch eq 'on') then begin
        for j=0,div-1 do begin
            oplot,x0,rmsarr[0,*,j],psym=symbol[j],color=color[j],thick=3
            xyouts,0.02,y,names[j],charthick=2,charsize=1.5,color=color[j],$
              /normal
            y = y-0.04
        endfor
    endif else begin
        for j=0,div-1 do begin
            oplot,x0,rmsarr[0,*,j],psym=symbol[j],color=60+(j*cs),thick=3
            xyouts,0.02,y,names[j],charthick=2,charsize=1.5,color=60+(j*cs),$
              /normal
            y = y-0.04
        endfor
    endelse
    if (swchname eq 'on') then begin
        loadct,0
        for j=0,step-1 do begin
            xyouts,x0[j],ydown+0.0002,blist[j],orientation=90,charsize=1.0
        endfor
    endif
    device,/close
    set_plot,'x'
    loadct,0
endif else print,'No plot for you!'

print,'The next ENTER deletes the plot...'
pause
wdelete,0

stop
end

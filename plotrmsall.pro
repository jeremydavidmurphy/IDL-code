;This routine is used to plot the rms values coming from the output of
;rfitlov. This is now an obsolete code and has been supplanted by a
;more flexible version. This code accepts the rmsall.out file and
;plots it, yet it requires that each fitting region is the same
;length.

pro plotrmsall

readcol,'rms.out',format='X,X,D,X',rms

div = intarr(1)
print,'The number of regions in the rms.out file:'
read,div
div = div[0]
for j=0,div-1 do begin
    names = strarr(div)
    print,'Enter the name of region #'+strn(j+1)
    read,names[j]
endfor

n0 = n_elements(rms)
step = n0/div
n1 = 0
n2 = step
x0 = findgen(n2)+1
xs = x0
xar = intarr(n2,div)
rmsarr = dblarr(1,step,div)

for j=0,div-1 do begin
    rmsarr[0,*,j] = rms[n1:n2-1]
    n1 = n1+step
    n2 = n2+step
    xar[*,j] = xs
    xs = xs + n1
endfor

window,0,retain=2
device,decomposed=0

loadct,0
plot,x1,rms,psym=3,xtitle='Files',ytitle='RMS',xrange=[-10,n1+10],$
  xstyle=1,charsize=1.5,title='RMS plot taken from the output of FITLOV'

loadct,4
oplot,x1,rmshk,psym=1,color=60
oplot,x2,rmsgband,psym=1,color=110
oplot,x3,rmshbeta,psym=1,color=150
oplot,x4,rmsmg,psym=1,color=180
oplot,x5,rmsall,psym=1,color=255

xyouts,0.45,0.86,'Blue: H+K',charthick=1.5,charsize=2,color=60,/normal
xyouts,0.45,0.82,'Green: Gband',charthick=1.5,charsize=2,color=110,/normal
xyouts,0.45,0.78,'Red: H-Beta',charthick=1.5,charsize=2,color=150,/normal
xyouts,0.45,0.74,'Orange: Mg',charthick=1.5,charsize=2,color=180,/normal
xyouts,0.45,0.70,'Yellow: All',charthick=1.5,charsize=2,color=255,/normal

ans=''
print,'Save the plot? (y or n):'
read,ans

if (ans eq 'y') then begin
    set_plot,'ps'
    device,file='rms.ps',/color
    
    loadct,0
    plot,x1,rms,psym=3,xtitle='Files',ytitle='RMS',xrange=[-10,n1+10],$
      xstyle=1,xthick=2,ythick=2,charthick=2,title='RMS plot taken from the output of FITLOV'
    oplot,x5,rmsall,psym=1,thick=2
    xyouts,0.45,0.70,'Black: All',charsize=1.5,charthick=2,/normal
    
    loadct,4
    oplot,x1,rmshk,psym=1,color=60,thick=2
    oplot,x2,rmsgband,psym=1,color=110,thick=2
    oplot,x3,rmshbeta,psym=1,color=150,thick=2
    oplot,x4,rmsmg,psym=1,color=180,thick=2

    xyouts,0.45,0.86,'Blue: H+K',charsize=1.5,charthick=2,color=60,/normal
    xyouts,0.45,0.82,'Green: Gband',charsize=1.5,charthick=2,color=110,/normal
    xyouts,0.45,0.78,'Red: H-Beta',charsize=1.5,charthick=2,color=150,/normal
    xyouts,0.45,0.74,'Orange: Mg',charsize=1.5,charthick=2,color=180,/normal

    loadct,0
    device,/close
    set_plot,'x'
endif else print,'Word. No plot will be generated...'

stop
end

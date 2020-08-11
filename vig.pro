pro vig, list

readcol,silent=1,list,format='a',files

n1=n_elements(files)
test = readfits(files[0])
n2 = n_elements(test)
sum1 = dblarr(n1)
sum2 = dblarr(n1)
sum3 = dblarr(n1)
sum4 = dblarr(n1)
sum5 = dblarr(n1)

bias = readfits('bias.fits')
print,min(bias)
print,median(bias)

for j=0,n1-1 do begin
    data = readfits(files[j],/silent)
    print,'working on '+files[j]
;    print,median(data)
;    print,min(data)
;    print,imax
    mn = 1800
    m = max(data)
    imax = where(data eq m)
    data[imax] = 0.0
    m = max(data)
    imax = where(data eq m)
    data[imax] = 0.0
    m = max(data)
    imax = where(data eq m)
    data[imax] = 0.0
    sum1[j] = total(data)
    datam = data - mn
    sum3[j] = total(datam)
    datab = data - bias
    sum2[j] = total(datab)
    
    sorted = datam[bsort(datam)]
;    loadct,0
;    window,0,retain=2
;    plot,sorted,charsize=2.0,psym=3
;    print,'enter an x-axis cutoff:'
;    read,cutoff
    cutoff = 8.1e6
;    i = where(sorted lt cutoff)
;    trimmed = sorted[i]
    trimmed = sorted[cutoff:*]
    xaxis = findgen(n2-cutoff)+cutoff
;    loadct,4
;    oplot,xaxis,trimmed,color=150,thick=3
    remains = sorted[0:cutoff]
;    window,2,retain=2
;    plot,remains,psym=3
    sum4[j] = total(trimmed)
    sum5[j] = total(remains)
    if (j eq 0) then begin
        set_plot,'ps'
        device,file='halo.ps',/color
        loadct,0
        plot,remains,psym=3,/nodata,xthick=3,ythick=3,charthick=3,$
          title='Spot Halo Estimates',charsize=1.0,ytitle='Counts',$
          yrange=[0,400],ystyle=1,xrange=[0,8.1e6],xstyle=1
        loadct,4
        oplot,remains,psym=3,color=60
;        set_plot,'x'
    endif else begin
;        set_plot,'ps'
;        device,file='halo.ps',/color
        if (j le 5) then begin
;            loadct,4
            oplot,remains,psym=3,color=60
;            set_plot,'x'
        endif
        if (j gt 5) then begin
            oplot,remains,psym=3,color=150
;            set_plot,'x'
        endif
    endelse
endfor
device,/close_file
;loadct,0
;window,0,retain=2
;loadct,0
set_plot,'ps'
device,file='totals.ps',/color

plot,sum1,xrange=[-2,n1+2],xstyle=1,psym=2,yrange=[6e9,1.15e10],ystyle=1,$
  title='Fiber Total Estimates',charsize=1.0,xthick=3,ythick=3,charthick=3
oplot,sum2,psym=1,thick=3
oplot,sum3,psym=4,thick=3
oplot,sum4,psym=5,thick=3

device,/close_file

;print,'The percent of light tossed from sum #4 (which was the visual clipping method) is plotted now...'
;pause
percent = sum5/sum4

;window,1,retain=2

set_plot,'ps'
device,file='percent_loss.ps'
plot,percent,/ynozero,psym=2,symsize=2,title='percent of light tossed',$
  xrange=[-2,n1+2],xstyle=1,ytitle='Percentage of light OUTSIDE of the spot',$
  xthick=3,ythick=3,charthick=3,thick=3
device,/close_file

stop
end

pro wavecheck_v1,data,ptow,mask,LINE=line,PIXEL=pixel,LIST=linelist

;this routine reads in an uncollapsed data file, the mask file and
;ptow file. it then returns where it thinks the 5577.34A Oxygen line
;is located.
step = 2

if (n_elements(line) eq 0) then begin
;line = 5577.34
;line = 4359.0
;line = 5460.75  ;mercury....?
;line = 5769.6
    
;line = 3611.2867
;line = 4046.5539
;line = 4077.8298
;line = 4678.1474
;line = 4799.9080
;line = 4916.0661
;line = 5085.8173
;line = 5460.7366
;line = 5769.6
endif


data = readfits(data)
readcol,ptow,format='d,d,d,d,d',c0,c1,c2,c3,c4
readcol,mask,format='i,x',fib

n1 = n_elements(data[*,0])
n2 = n_elements(data[0,*])
n3 = n_elements(fib)

wave = dblarr(n1,n2)

if (n_elements(pixel) ne 0) then begin
    for k=0,n3-1 do begin
        for j=0,n1-1 do begin
            wave[j,k] = c0[k] + c1[k]*j + c2[k]*j*j + c3[k]*j*j*j + c4[k]*j*j*j*j
            if (j eq pixel and k eq 0) then wavecenter = wave[j,k]
        endfor
    endfor
    print,'The estimated wavelength of the given pixel position is '+strn(wavecenter)    
stop
endif

for k=0,n3-1 do begin
    for j=0,n1-1 do begin
        wave[j,k] = c0[k] + c1[k]*j + c2[k]*j*j + c3[k]*j*j*j + c4[k]*j*j*j*j
    endfor
endfor

maski = fib-1

cold = dblarr(n1,n3)

for j=0,n3-1 do begin
    fiber = data[*,maski[j]-step:maski[j]+step]
    if (median(fiber) eq -666) then goto, jump1 ;the -666 fibers are left at zeros
    wgts = median(fiber,dimension=1)
    wgts = wgts/max(wgts)
    for k=0,n1-1 do begin
        fiber[k,*] = fiber[k,*]/wgts
    endfor
    cold[*,j] = total(fiber,2)
    jump1:
endfor

ans = ''
ans2 = ''
;print,'Plot the fits to screen?'
;read,ans

aa = dblarr(4,n3)
wdelta = dblarr(n3)
vdelta = dblarr(n3)
;the gaussian loop
ans='y'
if (ans eq 'y') then window,0,retain=2

for j=0,n3-1 do begin
    spec = cold[*,j]
    if (median(spec) ne 0.0) then begin
        w = wave[*,j]
        i1 = where(w gt floor(line-15))
        i1 = i1[0]
        i2 = where(w gt floor(line+15))
        i2 = i2[0]
        ps = spec[i1:i2]
        pw = w[i1:i2]
        g = gaussfit(pw,ps,a,nterms=4)
        if (ans eq 'y') then begin
            plot,pw,ps,psym=-2,title='Fiber '+strn(j+1),/ynozero
            oplot,pw,g,thick=3
            read,ans2
            if(ans2 eq 'q') then ans = 'n'
        endif
        aa[*,j] = a
        wdelta[j] = a[1] - line
        temp = (wdelta[j]/line)*299792.458
        if (temp lt 150 and temp gt -150) then vdelta[j] = temp 
    endif
endfor

i = where(vdelta ne 0.0)

aa = aa[*,i]
wdelta = wdelta[i]
vdelta = vdelta[i]
wcent = mean(wdelta)
vcent = mean(vdelta)
scatter = meanabsdev(wdelta[1,*])

fibernumber = indgen(n_elements(vdelta))+1

;print,'Plot VS wavelength or velocity? (enter "w" or "v")
;read,ans
ans = 'w'
set_plot,'ps'
;device,file='wavecheck.ps'
device,file=strn(line)+'.ps'

if (ans eq 'v') then begin
    plot,fibernumber,vdelta,psym=sym(1),xtitle='Fiber Number',ytitle='Velocity Error (!ukm!n/!dsec!n)',$
      xthick=3,ythick=3,charthick=3,xrange=[-5,250],xstyle=1,position=[0.12,0.12,0.9,0.9]
    xyouts,0.2,0.85,'Assumed line center: '+strn(line)+' A',charsize=1.1,charthick=3,/normal
    xyouts,0.2,0.81,'Mean calculated wavelenth: '+strn(wcent)+' A',charsize=1.1,charthick=3,/normal
    xyouts,0.2,0.77,'Mean velocity error: '+strn(vcent)+' !ukm!n/!dsec!n',charsize=1.1,charthick=3,/normal
endif
if (ans eq 'w') then begin
    plot,fibernumber,wdelta,psym=sym(1),xtitle='Fiber Number',ytitle='Wavelength Error (A)',$
      xthick=3,ythick=3,charthick=3,xrange=[-5,250],xstyle=1,position=[0.12,0.12,0.9,0.9];,yrange=[-.2,0.6],ystyle=1
    xyouts,0.2,0.85,'Assumed line center: '+strn(line)+' A',charsize=1.1,charthick=3,/normal
    xyouts,0.2,0.81,'Mean calculated wavelenth: '+strn(wcent)+' A',charsize=1.1,charthick=3,/normal
    xyouts,0.2,0.77,'Mean wavelength error: '+strn(wcent)+' A',charsize=1.1,charthick=3,/normal
endif

;axis,yaxis=1,ytitle='Wavelength Error (A)',charthick=3
;oplot,fibernumber,wdelta,psym=3
device,/close_file

device,file='wavecheck.ps'
plot,fibernumber,wdelta,psym=sym(1),xtitle='Fiber Number',ytitle='Wavelength Error (A)',$
  xthick=3,ythick=3,charthick=3,xrange=[-5,250],xstyle=1,position=[0.12,0.12,0.9,0.9] ;,yrange=[-.2,0.6],ystyle=1
xyouts,0.2,0.85,'Assumed line center: '+strn(line)+' A',charsize=1.1,charthick=3,/normal
xyouts,0.2,0.81,'Mean calculated wavelenth: '+strn(wcent)+' A',charsize=1.1,charthick=3,/normal
xyouts,0.2,0.77,'Mean wavelength error: '+strn(wcent)+' A',charsize=1.1,charthick=3,/normal
device,/close_file


set_plot,'x'
stop
end

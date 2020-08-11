pro wavecheck,data,ptow,mask,list,vorw
;this routine reads in an uncollapsed data file, the mask file and
;ptow file. it then reads in a list of lines that it will attempt to
;fit to a gaussian, one fiber at a time, and generate .ps files of the
;offsets between the assumed line center (the values in the list) and
;the calculated values. The abcissa is fiber number and the ordinate
;either delta wavelength or delta velocity. The choice of which is set
;by the "vorw" parameter.

step = 0 ;this dictates the number of rows of the fiber to bin.

;DATA: The file you are running the gaussian fitting routine on
;PTOW: The ptow.dat file
;MASK: The mask.dat file
;LIST: The list of lines you want to attempt to fit. THIS HAS TWO
;COLUMNS. THE FIRST IS THE LINE CENTRIOD WHILE THE SECOND SETS THE
;RANGE OVER WHICH THE GAUSSIAN FUNCTION IS DEFINED (bigger lines want
;bigger ranges...)
;VORW: Type either 'v' or 'w'. This just sets the ordinate values to
;plot either delta wavelength (in Angstroms) or velocity (km/sec)

; GAUSSFIT:
; result = gaussfit(x,y,A,nterms=4)
; where 'A' are coefficients of the form
; A_0*e^(-(x-A_1/2*A_2)) + A_3

readcol,list,silent=1,format='a,i',lines,range
n0 = n_elements(lines)

data = readfits(data)
readcol,ptow,format='d,d,d,d,d',c0,c1,c2,c3,c4
readcol,mask,format='i,x',fib

n1 = n_elements(data[*,0])
n2 = n_elements(data[0,*])
n3 = n_elements(fib)

wave = dblarr(n1,n2)

for k=0,n3-1 do begin
    for j=0,n1-1 do begin
        wave[j,k] = c0[k] + c1[k]*j + c2[k]*j*j + c3[k]*j*j*j + c4[k]*j*j*j*j
    endfor
endfor

maski = fib-1

cold = dblarr(n1,n3)

;the data is collapsed
for j=0,n3-1 do begin
    if (step eq 0) then fiber = data[*,maski[j]] else begin
        fiber = data[*,maski[j]-step:maski[j]+step]
    endelse

    if (median(fiber) eq -666) then goto, jump1 ;the -666 fibers are left at zeros
    wgts = median(fiber,dimension=1)
    wgts = wgts/max(wgts)
    for k=0,n1-1 do fiber[k,*] = fiber[k,*]/wgts
    if (step eq 0) then cold[*,j] = fiber else cold[*,j] = total(fiber,2)
    jump1:
endfor

out = dblarr(6,n0)

for l=0,n0-1 do begin
    window,0,retain=2
    loadct,0
    ans = ''
    
    aa = dblarr(4,n3)
    wdelta = dblarr(n3)
    vdelta = dblarr(n3)
    varr = dblarr(n3)
    for j=0,n3-1 do begin
        spec = cold[*,j]
        if (median(spec) ne 0.0) then begin
            w = wave[*,j]
            i1 = where(w gt floor(lines[l]-range[l]))
            i1 = i1[0]
            i2 = where(w gt floor(lines[l]+range[l]))
            i2 = i2[0]
            ps = spec[i1:i2]
            pw = w[i1:i2]
            g = gaussfit(pw,ps,a,nterms=4)
            if (ans ne 'q') then begin
                plot,pw,ps,psym=-2,title='Fiber '+strn(j+1),/ynozero,xtitle=strn(lines[l]),charsize=1.2
                oplot,pw,g,thick=3
                read,ans
            endif
        endif
        aa[*,j] = a
        wdelta[j] = a[1] - lines[l]
        varr[j] = (a[1]/lines[l])*299792.458
        temp = (wdelta[j]/lines[l])*299792.458
        if (temp lt 150 and temp gt -150) then vdelta[j] = temp 
    endfor

    i = where(vdelta ne 0.0)
    
    aa = aa[*,i]
    wdelta = wdelta[i]
    vdelta = vdelta[i]
    varr = varr[i]

    wcent = median(aa[1,*])
    vcent = (wcent/lines[l])*299792.458
    wsd = stddev(aa[1,*])
    vsd = stddev(varr)

    fibernumber = indgen(n_elements(vdelta))+1
    
    out[0,l] = lines[l]
    out[1,l] = wcent
    out[2,l] = mean(wdelta)
    out[3,l] = wsd
    out[4,l] = mean(vdelta)
    out[5,l] = vsd
    
    if (vorw eq 'v') then begin
        set_plot,'ps'
        device,file=strn(lines[l])+'v.ps',/color
        loadct,0
        plot,fibernumber,vdelta,psym=sym(1),xtitle='Fiber Number',ytitle='Velocity Error (!ukm!n/!dsec!n)',$
          xthick=3,ythick=3,charthick=3,xrange=[-5,250],xstyle=1,position=[0.12,0.12,0.9,0.9]
        loadct,4
        xyouts,0.2,0.85,'Assumed line center: '+strn(lines[l])+' A',charsize=1.1,charthick=3,/normal,color=150
        xyouts,0.2,0.81,'Mean velocity error: '+strn(vcent)+' !ukm!n/!dsec!n',charsize=1.1,charthick=3,/normal,color=150
        xyouts,0.2,0.77,'Standard Deviation: '+strn(vsd)+' !ukm!n/!dsec!n',charsize=1.1,charthick=3,/normal,color=150
        device,/close_file

        set_plot,'x'
        window,0,retain=2,xsize=800,ysize=600
        device,decomposed=0
        
        loadct,0
        plot,fibernumber,vdelta,psym=sym(1),xtitle='Fiber Number',ytitle='Velocity Error (!ukm!n/!dsec!n)',$
          xrange=[-5,250],xstyle=1,position=[0.12,0.12,0.9,0.9],charsize=1.5
        loadct,4
        xyouts,0.2,0.85,'Assumed line center: '+strn(lines[l])+' A',charsize=1.3,/normal,color=255
        xyouts,0.2,0.81,'Mean velocity error: '+strn(vcent)+' !ukm!n/!dsec!n',charsize=1.3,/normal,color=255
        xyouts,0.2,0.77,'Standard Deviation: '+strn(vsd)+' !ukm!n/!dsec!n',charsize=1.3,/normal,color=255

        pause
        wdelete,0
    endif

    if (vorw eq 'w') then begin
        set_plot,'ps'
        device,file=strn(lines[l])+'w.ps',/color
        loadct,0
        plot,fibernumber,wdelta,psym=sym(1),xtitle='Fiber Number',ytitle='Wavelength Error (A)',$
          xthick=3,ythick=3,charthick=3,xrange=[-5,250],xstyle=1,position=[0.12,0.12,0.9,0.9]
        loadct,4
        xyouts,0.2,0.85,'Assumed line center: '+strn(lines[l])+' A',charsize=1.1,charthick=3,/normal,color=150
        xyouts,0.2,0.81,'Mean estimated wavelenth: '+strn(wcent)+' A',charsize=1.1,charthick=3,/normal,color=150
        xyouts,0.2,0.77,'Standard Deviation: '+strn(wsd)+' A',charsize=1.1,charthick=3,/normal,color=150
        device,/close_file

        set_plot,'x'
        window,0,retain=2,xsize=800,ysize=600
        device,decomposed=0

        loadct,0
        plot,fibernumber,wdelta,psym=sym(1),xtitle='Fiber Number',ytitle='Wavelength Error (A)',$
          xrange=[-5,250],xstyle=1,position=[0.12,0.12,0.9,0.9],charsize=1.5
        loadct,4
        xyouts,0.2,0.85,'Assumed line center: '+strn(lines[l])+' A',/normal,charsize=1.3,color=255
        xyouts,0.2,0.81,'Mean estimated wavelenth: '+strn(wcent)+' A',/normal,charsize=1.3,color=255
        xyouts,0.2,0.77,'Standard Deviation: '+strn(wsd)+' A',/normal,charsize=1.3,color=255
    
        pause
        wdelete,0
    endif
    
    
endfor

f1 = '(f9.4,2x,f9.4,2x,f9.5,2x,f9.5,2x,f9.5,2x,f9.5)'
openw,5,'wavecheck.dat'
for j=0,n0-1 do printf,5,out[*,j],format=f1
free_lun,5

end

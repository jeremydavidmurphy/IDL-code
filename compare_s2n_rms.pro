; this routine generates a plot comparing the rms values (as
; calculated from the signal-2-noise of your spectra via s2ncalc.pro)
; and that output during the fitlov routine (i.e. the rms of the
; template star-to-data fit)

; Modified on March 6, 2012: This code now reads in the bin_S2N.fits
; file and the rmsall.out files from Hoku and makes a comparison. This
; can then generate the numbers you feed into the Monte Carlo fitting
; routines.

PRO compare_s2n_rms,s2nlist,rmslist

galname = 'NGC 307'
plotPS = 'yes'
range = [0,0.6] ;controls the plotting range for both axes
hk = [3600,4050]
gb = [4215,4575]
hb = [4445,4945]
mg = [4900,5420]
fe = [5060,5720]
al = [3600,5500]

;****************************************************************

readcol,s2nlist,f='a',s2nfiles
readcol,rmslist,f='a',rmsfiles
n0 = n_elements(s2nfiles)
if (n_elements(rmsfiles) ne n0) then stop

for j=0,n0-1 do begin ;a loop over each bin
    s2nfile = s2nfiles[j]
    rmsfile = rmsfiles[j]
    s2n = readfits(s2nfile,header,/silent)
    if (j eq 0) then n1 = n_elements(s2n)
    wi = sxpar(header,'CRVAL1',count=cnt)
    wd = sxpar(header,'CDELT1',count=cnt)
    wave = findgen(n1) * wd + wi

    readcol,rmsfile,silent=1,format='a,f,x,x,x',name,rms
    if (j eq 0) then n2 = n_elements(name)
    if (j eq 0) then rmsarray = fltarr(2,n2,n0) ;rms and 1/sn, regions, spatial bins
    regions = strarr(n2)
    for k=0,n2-1 do begin
        s = strsplit(name[k],'_.',/extract) ;NOTE: This assumes a certain naming convention
        regions[k] = s[1]
    endfor
    ihk = where(regions eq 'HKi',chk)
    igb = where(regions eq 'GBi',cgb)
    ihb = where(regions eq 'HBi',chb)
    img = where(regions eq 'MGi',cmg)
    ife = where(regions eq 'FEi',cfe)
    ial = where(regions eq 'ALi',cal)

    if (chk gt 0) then begin
        rmsarray[0,0,j] = rms[ihk]
        i1 = where(wave ge hk[0] and wave le hk[1])
        rmsarray[1,0,j] = 1.0 / median(s2n[i1],/even)
    endif
    if (cgb gt 0) then begin
        rmsarray[0,1,j] = rms[igb]
        i1 = where(wave ge gb[0] and wave le gb[1])
        rmsarray[1,1,j] = 1.0 / median(s2n[i1],/even)
    endif
    if (chb gt 0) then begin
        rmsarray[0,2,j] = rms[ihb]
        i1 = where(wave ge hb[0] and wave le hb[1])
        rmsarray[1,2,j] = 1.0 / median(s2n[i1],/even)
    endif
    if (cmg gt 0) then begin
        rmsarray[0,3,j] = rms[img]
        i1 = where(wave ge mg[0] and wave le mg[1])
        rmsarray[1,3,j] = 1.0 / median(s2n[i1],/even)
    endif
    if (cfe gt 0) then begin
        rmsarray[0,4,j] = rms[ife]
        i1 = where(wave ge fe[0] and wave le fe[1])
        rmsarray[1,4,j] = 1.0 / median(s2n[i1],/even)
    endif
    if (cal gt 0) then begin
        rmsarray[0,5,j] = rms[ial]
        i1 = where(wave ge al[0] and wave le al[1])
        rmsarray[1,5,j] = 1.0 / median(s2n[i1],/even)
    endif
endfor

set_plot,'x'
window,0,retain=2,xsize=700,ysize=700
device,decomposed=0
loadct,0,/silent
plot,rmsarray[0,0,*],rmsarray[1,0,*],psym=1,/iso,xtitle='RMS from template fit',$
  ytitle='RMS from S2N calculation',xr=[range[0],range[1]],$
  yr=[range[0],range[1]],/xs,/ys,title=galname
oplot,[0,1],[0,1]
loadct,33,/silent
oplot,rmsarray[0,ihk,*],rmsarray[1,ihk,*],psym=1,color=30
oplot,rmsarray[0,igb,*],rmsarray[1,igb,*],psym=1,color=70
oplot,rmsarray[0,ihb,*],rmsarray[1,ihb,*],psym=1,color=110
oplot,rmsarray[0,img,*],rmsarray[1,img,*],psym=1,color=150
oplot,rmsarray[0,ife,*],rmsarray[1,ife,*],psym=1,color=190
oplot,rmsarray[0,ial,*],rmsarray[1,ial,*],psym=1,color=230
    
set_plot,'ps'
device,file='s2n_v_rms.ps',/color
loadct,0,/silent
plot,rmsarray[0,0,*],rmsarray[1,0,*],psym=1,/iso,xtitle='RMS from template fit',$
  ytitle='RMS from S2N calculation',xr=[range[0],range[1]],$
  yr=[range[0],range[1]],/xs,/ys,title=galname,position=[0.13,0.13,0.93,0.89],$
  xthick=4,ythick=4,charthick=4
loadct,33,/silent
oplot,rmsarray[0,ihk,*],rmsarray[1,ihk,*],psym=1,color=30,thick=3
oplot,rmsarray[0,igb,*],rmsarray[1,igb,*],psym=1,color=70,thick=3
oplot,rmsarray[0,ihb,*],rmsarray[1,ihb,*],psym=1,color=110,thick=3
oplot,rmsarray[0,img,*],rmsarray[1,img,*],psym=1,color=150,thick=3
oplot,rmsarray[0,ife,*],rmsarray[1,ife,*],psym=1,color=190,thick=3
oplot,rmsarray[0,ial,*],rmsarray[1,ial,*],psym=1,color=230,thick=3
loadct,0
oplot,[0,1],[0,1],thick=2
device,/close
set_plot,'x'

stop
END


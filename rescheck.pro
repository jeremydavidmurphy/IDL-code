PRO rescheck, month, night

; This routine is used to calculate the FWHM of the reduced arcs from
; an observing run. It generates a PS file showing FWHM, in angstroms,
; for both X and Y on the chip.

;MONTH and NIGHT: The code will look for the extracted arc, ptow and
;mask files for a given month and night.

; This code will look for the HgCd.dat in the calling directory. 

; There is also an option to generate statistics on a given
; bin. Returned as bin####.[BGR]stats, a text file is generated giving
; the median, mean and variance for FWHM in angstroms and km/sec for
; the given bin. The user is prompted to enter a list of bins.

;*********************************************************************
nprofile = 5 ;the number of rows in a fiber profile
step = (nprofile-1)/2

swch = 'ff' ;change to anything other than 'on' to stop the plotting routine.
;*********************************************************************
arc = 'arc_'+month+'_'+night+'e.fits'
ptow = 'ptow_'+month+'_'+night+'.dat'
mask = 'mask_'+month+'_'+night+'.dat'

data = readfits(arc,/silent)
n0 = n_elements(data[*,0])
n1 = n_elements(data[0,*])

readcol,mask,silent=1,format='i,i',row,fib
n2 = n_elements(fib)
index = row - 1

readcol,ptow,silent=1,format='f,f,f,f,f',c0,c1,c2,c3,c4
if (n_elements(c4) ne n2) then stop
wave = fltarr(n0,n2)
for j=0,n2-1 do begin
    for k=0,n0-1 do begin
        wave[k,j] = c0[j]+c1[j]*k+c2[j]*k*k+c3[j]*k*k*k+c4[j]*k*k*k*k
    endfor
endfor
disp = median(c1)

readcol,'HgCd.dat',silent=1,format='f,a',hgcd1,descript
n3 = n_elements(hgcd1)
hgcd0 = hgcd1-12
hgcd2 = hgcd1+12

collapsed = fltarr(n0,n2)
values = fltarr(n3,n2,4) ; wavelength,fiber #,(

print,''
print,'Working...'
print,''
if (swch eq 'on') then window,2,retain=2
for j=0,n2-1 do begin
    fiber = total(data[*,index[j]-step:index[j]+step],2)
    onew = wave[*,j]
    if (total(onew) eq 0) then goto, jump1
    for k=0,n3-1 do begin
        i = where(onew gt hgcd0[k] and onew lt hgcd2[k])
        pd = fiber[i]
        pw = onew[i]
        temp = gaussfit(pw,pd,a,nterms=4)
        fwhmA = 2*sqrt(2*alog(2))*a[2]
        values[k,j,0] = a[1]
        values[k,j,1] = a[1]-hgcd1[k]
        values[k,j,2] = fwhmA
        values[k,j,3] = 299298.0*(a[2]*disp/a[1])
        if (swch eq 'on') then begin
            plot,pw,pd,title=arc,charsize=1.5,xtitle='Angstroms'
            xyouts,0.18,0.88,'Calculated wavelength: '+strn(a[1]),charsize=1.2,/normal
            xyouts,0.18,0.84,'Theoretical wavelength: '+strn(hgcd1[k]),charsize=1.2,/normal
            xyouts,0.18,0.80,'Offset: '+strn(values[k,j,1]),charsize=1.2,/normal
            xyouts,0.18,0.76,'FWHM in angstroms: '+strn(fwhmA),charsize=1.2,/normal
        endif
    endfor
    jump1:
endfor

temp = strsplit(arc,'.',/extract)
outname1 = temp[0]+'_FWHM.ps'
outname2 = temp[0]+'_IDISP.ps'
;colors = [60,75,90,105,120,255,240,225,210,195,180,165,150]
colors = [60,75,90,105,120,225,210,195,180,165,150]
fibern = indgen(n2)+1

i = where(values[0,*,0] eq 0)
if (i[0] ne -1) then values[*,i,*] = !Values.F_NAN
ydn = min(values[*,*,2])-0.2
yup = max(values[*,*,2])+0.2

;FWHM PLOT
set_plot,'ps'
device,filename=outname1,/color

!p.multi = [0,1,2]
!y.omargin = [2,4]

loadct,0
plot,values[0,*,0],values[0,*,2],title='FWHM for '+temp[0],psym=1,/nodata,$
  xrange=[min(hgcd1)-50,max(hgcd1)+50],yrange=[ydn,yup],$
  xstyle=1,ystyle=1,xtitle='Wavelength (A)',ytitle='FWMH (A)',charsize=1.0,$
  xthick=2,ythick=2,charthick=2
loadct,4
for j=0,n3-1 do begin
    oplot,values[j,*,0],values[j,*,2],psym=1,color=colors[j],symsize=0.5,thick=1.5
endfor

loadct,0
plot,fibern,values[0,*,2],title='FWHM for '+temp[0],psym=1,/nodata,xrange=[-5,n2+5],$
  xstyle=1,yrange=[ydn,yup],$
  ystyle=1,xtitle='Fiber Number',ytitle='FWMH (A)',charsize=1.0,$
  xthick=2,ythick=2,charthick=2,symsize=0.5,thick=1.5
loadct,4
for j=0,n3-1 do begin
    oplot,fibern,values[j,*,2],psym=1,color=colors[j],symsize=0.5,thick=1.5
endfor
device,/close_file

; VELOCITY PLOT
ydn = min(values[*,*,3])-20
yup = max(values[*,*,3])+20
set_plot,'ps'
device,filename=outname2,/color

!p.multi = [0,1,2]
!y.omargin = [2,4]

loadct,0
plot,values[0,*,0],values[0,*,3],title='Instrumental Dispersion for '+temp[0],psym=1,/nodata,$
  xrange=[min(hgcd1)-50,max(hgcd1)+50],yrange=[ydn,yup],$
  xstyle=1,ystyle=1,xtitle='Wavelength (A)',ytitle='Velocity (km/sec)',charsize=1.0,$
  xthick=2,ythick=2,charthick=2,symsize=0.5,thick=2
loadct,4
for j=0,n3-1 do begin
    oplot,values[j,*,0],values[j,*,3],psym=1,color=colors[j],symsize=0.5,thick=1.5
endfor

loadct,0
plot,fibern,values[0,*,3],title='Instrumental Dispersion for '+temp[0],psym=1,/nodata,xrange=[-5,n2+5],$
  xstyle=1,yrange=[ydn,yup],$
  ystyle=1,xtitle='Fiber Number',ytitle='Velocity (km/sec)',charsize=1.0,$
  xthick=2,ythick=2,charthick=2
loadct,4
for j=0,n3-1 do begin
    oplot,fibern,values[j,*,3],psym=1,color=colors[j],symsize=0.5,thick=1.5
endfor
device,/close_file

ans = ''
print,'Enter a list of bins to do statistics on? (y/n)'
read,ans
ans = 'n'
if (ans eq 'y') then begin
    binlist=''
    print,'Enter the name of the bin list:'
    read,binlist
    readcol,binlist,silent=1,format='a',list
    n4 = n_elements(list)
    for k=0,n4-1 do begin ;a loop through each bin
        onelist = list[k]
        readcol,onelist,format='a',binlist
        n5 = n_elements(binlist)
        temp = strarr(2,n5)
        for l=0,n5-1 do temp[*,l] = strsplit(binlist[l],'_',/extract)
        fibers = temp[0,*]
        ifibers = fibers - 1
        fvalues = values[*,ifibers,*]
        igood = where(finite(fvalues[0,*,3]) eq 1)
        n6 = n_elements(igood)
        outval = fltarr(n3,n6,4)
        for l=0,n3-1 do outval[l,*,*] = fvalues[l,igood,*]
        blue = fltarr(6)
        blue[0] = median(outval[0:2,*,2])
        blue[1] = mean(outval[0:2,*,2])
        blue[2] = variance(outval[0:2,*,2])
        blue[3] = median(outval[0:2,*,3])
        blue[4] = mean(outval[0:2,*,3])
        blue[5] = variance(outval[0:2,*,3])
        openw,5,onelist+'.Bstats'
        printf,5,blue
        free_lun,5
        green = fltarr(6)
        green[0] = median(outval[3:5,*,2])
        green[1] = mean(outval[3:5,*,2])
        green[2] = variance(outval[3:5,*,2])
        green[3] = median(outval[3:5,*,3])
        green[4] = mean(outval[3:5,*,3])
        green[5] = variance(outval[3:5,*,3])
        openw,5,onelist+'.Gstats'
        printf,5,green
        free_lun,5
        red = fltarr(6)
        red[0] = median(outval[6:8,*,2])
        red[1] = mean(outval[6:8,*,2])
        red[2] = variance(outval[6:8,*,2])
        red[3] = median(outval[6:8,*,3])
        red[4] = mean(outval[6:8,*,3])
        red[5] = variance(outval[6:8,*,3])
        if (red[5] gt 50) then stop
        openw,5,onelist+'.Rstats'
        printf,5,red
        free_lun,5
    endfor
endif

stop
END

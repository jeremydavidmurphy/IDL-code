; This routine is used AFTER the pfitlov.out file has been generated
; and the bin????.coord files exist. It prepares the files for input
; into Michele's PLOT_VELFIELD.PRO routine.

pro preplotvelfield, REGION=region, PFIT=pfitfile, COORD=coord

; VALUE is a string of either 'vel' 'disp' 'h3' or 'h4'. Two lists are
; output. A master, containing all values, and a version ready for
; feeding into plot_velfield.pro.

; The REGION is a string, matching the pfitlov.out naming convention,
; for the region you want to plot. IF OMITTED, THE MEDIAN OF ALL THE
; RELEVANT REGIONS IS RETURNED.

; The COORD list is a list of all the b[in]n#### files.

if (n_elements(region) eq 0) then med = 'on' else med = 'off'
if (n_elements(pfitfile) eq 0) then pfitfile = 'pfitlov.out'
if (n_elements(coord) eq 0) then coord = 'coord.list'

readcol,pfitfile,silent=1,format='a,f,f,x,f,f,x,x,x',name,vel,disp,h3,h4
readcol,coord,silent=1,format='a',coordfiles

n0 = n_elements(name)
bins = strarr(n0)
regions = strarr(n0)

;ans1=''
;ans2=''
;test = strsplit(name[0],'_.',/extract)
;print,''
;print,transpose(test)
;print,''
;print,'Enter the element that gives the bin:'
;read,ans1
;print,'Enter the element that gives the spectral region:'
;read,ans2
;ans1=uint(ans1)-1
;ans2=uint(ans2)-1
ans1=0
ans2=2
for j=0,n0-1 do begin
    temp = strsplit(name[j],'_.',/extract)
    bins[j] = temp[ans1]
    regions[j] = temp[ans2]
endfor

n1 = n_elements(coordfiles)


for j=0,n1-1 do begin
;    nRA = intarr(1)
;    nDec = intarr(1)
    file = coordfiles[j]
    readcol,file,silent=1,skipline=1,format='i,f,f',fiber,dRA,dDec
    ind = where(dRA ne -1.0) ;bad fibers are tossed
    fiber = fiber[ind]
    dRA = dRA[ind]
    dDec = dDec[ind]
    n3 = n_elements(dRA)
    print,file
;    for k=0,n3-2 do begin
;        t1 = dRA[k]
;        t2 = dRA[k+1]
;        if (n3-2 eq 1) then begin
;            dRA = dRA[k]
;            nDec = dDec[k]
;            goto,jump1
;        endif
;        if (t1 ne t2) then begin
;            nRA = [nRA,dRA[k]]
;            nDec = [nDec,dDec[k]]
;        endif
;    endfor
;    nRA = [nRA[1:*],dRA[k]]
;    nDec = [nDec[1:*],dDec[k]]
;    jump1:
;    n4 = n_elements(nRA)
    file = strsplit(file,'.',/extract)
    onebin = file[0]
    ibin = where(bins eq onebin)
    if (med eq 'on') then begin
        piece = transpose([[vel[ibin]],[disp[ibin]],[h3[ibin]],[h4[ibin]]])
        medpiece = median(piece,dimension=2)
;        print,''
;        print,piece
;        print,'And the median is...'
;        print,medpiece
;        pause
        if (j eq 0) then openw,5,'plotfield_median.dat'
        for k=0,n3-1 do printf,5,onebin,'median',dRA[k],dDec[k],medpiece,$
          format='(a9,2x,a6,2x,f8.3,2x,f8.3,2x,f10.3,2x,f10.3,2x,f7.3,2x,f7.3)'
;        printf,5,onebin,'median',dRA[k],dDec[k],medpiece,$
;          format='(a9,2x,a6,2x,f8.3,2x,f8.3,2x,f10.3,2x,f10.3,2x,f7.3,2x,f7.3)'
;         printf,5,onebin,'median',dRA,dDec,medpiece,$
;          format='(a9,2x,a6,2x,f8.3,2x,f8.3,2x,f10.3,2x,f10.3,2x,f7.3,2x,f7.3)'
    endif
    if (med eq 'off') then begin
        if (j eq 0) then openw,5,'plotfield_'+region+'.dat'
        ispec = where(regions eq region)
        indall = intarr(n_elements(ibin))
        for k=0,n_elements(ibin)-1 do indall[k] = where(ispec eq ibin[k])
        itemp = indall[where(indall ne -1)]
        if (itemp[0] eq -1) then stop
        ione = ispec[itemp]
;       ione = ibin[where(indall ne -1)]
        piece = [vel[ione],disp[ione],h3[ione],h4[ione]]
;        print,piece & pause
        for k=0,n3-1 do printf,5,onebin,region,dRA[k],dDec[k],piece,$
;        printf,5,onebin,region,dRA,dDec,piece,$
          format='(a9,2x,a4,2x,f8.3,2x,f8.3,2x,f10.3,2x,f10.3,2x,f7.3,2x,f7.3)'
    endif
endfor
free_lun,5    

;readcol,'temp',format='x,x,f,f,f,f,f,f',ra,dec,vel,disp,h3,h4,silent=1
if (med eq 'on') then readcol,'plotfield_median.dat',format='x,x,f,f,f,f,f,f',ra,dec,vel,disp,h3,h4,silent=1
if (med eq 'off') then readcol,'plotfield_'+region+'.dat',format='x,x,f,f,f,f,f,f',ra,dec,vel,disp,h3,h4,silent=1
;window,0,retain=2,xsize=900,ysize=900
;device,decomposed=0
;ans=''
;print,'Plot which thing?'
;print,'(vel/disp/h3/h4)'
;read,ans
nnn = n_elements(ra)
nn = floor(nnn/2.0)
ra1 = ra[0:nn]
ra2 = ra[nn+1:*]
dec1 = dec[0:nn]
dec2 = dec[nn+1:*]
vel1 = vel[0:nn]
vel2 = vel[nn+1:*]

rup = max(vel)
rdown = min(vel)

set_plot,'ps'
device,file='M87velocity.ps'
loadct,0
levels = (rup-rdown)/255.*findgen(255)+rdown
loadct,0
;contour,vel1,ra1,dec1,/fill,iso=1,levels=levels,/irregular,position=[0.1,0.1,0.5,0.5],$
;  yrange=[min(dec),max(dec)],/ystyle
contour,vel,ra,dec,iso=1,levels=levels,/irregular,position=[0.33,0.09,0.85,0.85],$
  yrange=[min(dec),max(dec)],/ystyle,/fill,/nodata,xthick=3,ythick=3,thick=3,$
  xtitle='RA (arcsec)',ytitle='Dec (arcsec)',charthick=3

loadct,33
;contour,vel1,ra1,dec1,/fill,/overplot,levels=levels,/irregular,iso=1,color=0
;contour,vel2,ra2,dec2,/fill,/overplot,levels=levels,/irregular,iso=1,color=0
contour,vel,ra,dec,/fill,/overplot,levels=levels,/irregular,iso=1,color=0

loadct,0
oplot,ra,dec,psym=3

;the colorbar
loadct,33
ncolors = n_elements(levels)
loc = [0.22,0.90,0.69,0.95]
bar = bindgen(256) # replicate(1b,10)
xsize = (loc[2] - loc[0]) * !d.x_vsize
ysize = (loc[3] - loc[1]) * !d.y_vsize
xstart = loc[0] * !d.x_vsize
ystart = loc[1] * !d.y_vsize

;bar = bytscl(bar, top=ncolors-1)
;bar = BYTSCL(bar, TOP=(ncolors-1) < (255-bottom)) + bottom
;bar = CONGRID(bar, CEIL(xsize*!D.X_VSize), CEIL(ysize*!D.Y_VSize), /INTERP)

tv, bar, xstart, ystart, xsize=xsize, ysize=ysize
loadct,0

plots, [loc[0], loc[0], loc[2], loc[2], loc[0]], [loc[1], loc[3], loc[3], loc[1], loc[1]],/normal

middle = (max(vel)-min(vel))/2
vup = middle
vdown = -middle
;vup = max(vel)
;vdown = min(vel)
xyouts,0.22,0.86,strn(vdown),/normal,charthick=3,charsize=1.5
xyouts,0.60,0.86,strn(vup),/normal,charthick=3,charsize=1.5
xyouts,0.42,0.86,'km/sec',/normal,charthick=3,charsize=1.5
device,/close_file



rup = max(disp)
rdown = min(disp)
set_plot,'ps'
device,file='M87dispersion.ps'
loadct,0
levels = (rup-rdown)/255.*findgen(255)+rdown
loadct,0
;contour,vel1,ra1,dec1,/fill,iso=1,levels=levels,/irregular,position=[0.1,0.1,0.5,0.5],$
;  yrange=[min(dec),max(dec)],/ystyle
contour,disp,ra,dec,iso=1,levels=levels,/irregular,position=[0.33,0.09,0.85,0.85],$
  yrange=[min(dec),max(dec)],/ystyle,/fill,/nodata,xthick=3,ythick=3,thick=3,$
  xtitle='RA (arcsec)',ytitle='Dec (arcsec)',charthick=3

loadct,33
;contour,vel1,ra1,dec1,/fill,/overplot,levels=levels,/irregular,iso=1,color=0
;contour,vel2,ra2,dec2,/fill,/overplot,levels=levels,/irregular,iso=1,color=0
contour,disp,ra,dec,/fill,/overplot,levels=levels,/irregular,iso=1,color=0

loadct,0
oplot,ra,dec,psym=3

;the colorbar
loadct,33
ncolors = n_elements(levels)
loc = [0.22,0.90,0.69,0.95]
bar = bindgen(256) # replicate(1b,10)
xsize = (loc[2] - loc[0]) * !d.x_vsize
ysize = (loc[3] - loc[1]) * !d.y_vsize
xstart = loc[0] * !d.x_vsize
ystart = loc[1] * !d.y_vsize

;bar = bytscl(bar, top=ncolors-1)
;bar = BYTSCL(bar, TOP=(ncolors-1) < (255-bottom)) + bottom
;bar = CONGRID(bar, CEIL(xsize*!D.X_VSize), CEIL(ysize*!D.Y_VSize), /INTERP)

tv, bar, xstart, ystart, xsize=xsize, ysize=ysize
loadct,0

plots, [loc[0], loc[0], loc[2], loc[2], loc[0]], [loc[1], loc[3], loc[3], loc[1], loc[1]],/normal

;middle = (max(disp)-min(disp))/2
;vup = middle
;vdown = -middle
xyouts,0.22,0.86,strn(rdown),/normal,charthick=3,charsize=1.5
xyouts,0.60,0.86,strn(rup),/normal,charthick=3,charsize=1.5
xyouts,0.42,0.86,'km/sec',/normal,charthick=3,charsize=1.5
device,/close_file
set_plot,'x'

stop
end

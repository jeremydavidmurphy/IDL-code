PRO Ppallmc_v1, file, COORD=coordfile

; This routine is used to plot the output of the pallmc routine. This
; is an estimate of vel, disp, h3 and h4 from the monte carlo
; simulations, thus includes errors.

readcol,file,silent=1,format='x,f,f,f,f,f,f,f,f,a,i',$
  v,verr,disp,disperr,h3,h3err,h4,h4err,name,ang
n0 = n_elements(v)
rflip = intarr(n0)
angs = [01,02,03,04,05]
colors = [0,60,110,180,150]

for j=0,n0-1 do begin
    temp = strsplit(name[j],'_.',/extract)
    name[j] = temp[0]
    region = temp[2]
    temp = strsplit(name[j],'0123456789',/extract)
    temp = temp[0]
    if (temp eq 'bin') then rflip[j] = 1
    if (temp eq 'bnn') then rflip[j] = -1
    if (temp ne 'bin') and (temp ne 'bnn') then stop
endfor

if (n_elements(coordfile) eq 0) then $
readcol,'coord.list',silent=1,format='a',coordfiles else $
readcol,coordfile,silent=1,format='a',coordfiles

n1 = n_elements(coordfiles)

if (n0 ne n1) then stop
radius = fltarr(n1)

for j=0,n1-1 do begin
    readcol,coordfiles[j],silent=1,skipline=1,$
      format='x,f,f',dra,ddec
    n2 = n_elements(dra)
    radarr = fltarr(n2)
    for k=0,n2-1 do radarr[k] = sqrt(dra[k]^2 + ddec[k]^2)
    radius[j] = mean(radarr)
endfor

free_lun,5
openw,5,'bin.radius'
for j = 0,n1-1 do begin
    printf,5,coordfiles[j],radius[j],format='(a16,2x,f7.3)'
endfor
free_lun,5

radius = radius * rflip
xup = max(radius)+10
xdown = min(radius)-10

;-------------------------------------------------------------
;                    VELOCITY
;-------------------------------------------------------------

ans=''
print,'Fix the y-axis values for velocity? (y/n)'
read,ans
if (ans eq 'y') then begin
    yup=intarr(1)
    ydown=intarr(1)
    print,'Enter a lower limit:'
    read,ydown
    print,'Enter an upper limit:'
    read,yup
endif
    
set_plot,'ps'
device,file='velocity_'+region+'.ps',/color
loadct,0
ind = where(ang eq angs[0])
if (ans eq 'y') then ploterror,radius[ind],v[ind],verr[ind],psym=1,$
  title='Velocity values for the '+region+' region',yrange=[ydown,yup],$
  xtitle='Radius (arcsec)',ytitle='Velocity (km/sec)',xrange=[xdown,xup],$
  xthick=3,ythick=3,thick=3,xstyle=1,charthick=3,errthick=3,ystyle=1 else $
ploterror,radius[ind],v[ind],verr[ind],psym=1,$
  title='Velocity values for the '+region+' region',$
  xtitle='Radius (arcsec)',ytitle='Velocity (km/sec)',xrange=[xdown,xup],$
  xthick=3,ythick=3,thick=3,xstyle=1,charthick=3,errthick=3

plots,0.15,0.91,psym=2,/normal,thick=3
xyouts,0.18,0.90,'01',/normal,charthick=3

loadct,4
for j=1,4 do begin
    ind = where(ang eq angs[j])
    oploterror,radius[ind],v[ind],verr[ind],psym=1,color=colors[j],$
      errcolor=colors[j],errthick=3,thick=3
    plots,0.15,0.91-(j*0.03),psym=2,/normal,thick=3,color=colors[j]
    xyouts,0.18,0.90-(j*0.03),'0'+strn(angs[j]),/normal,$
      charthick=3,color=colors[j]
endfor
device,/close_file

;-------------------------------------------------------------
;                    VELOCITY DISPERION
;-------------------------------------------------------------

print,'Fix the y-axis values for dispersion? (y/n)'
read,ans
if (ans eq 'y') then begin
    yup=intarr(1)
    ydown=intarr(1)
    print,'Enter a lower limit:'
    read,ydown
    print,'Enter an upper limit:'
    read,yup
endif

set_plot,'ps'
device,file='dispersion_'+region+'.ps',/color
loadct,0
ind = where(ang eq angs[0])
if (ans eq 'y') then ploterror,radius[ind],disp[ind],disperr[ind],psym=1,$
  title='Velocity dispersion values for the '+region+' region',$
  xtitle='Radius (arcsec)',ytitle='Velocity Dispersion (km/sec)',$
  xrange=[xdown,xup],yrange=[ydown,yup],ystyle=1,$
  xthick=3,ythick=3,thick=3,xstyle=1,charthick=3,errthick=3 else $
ploterror,radius[ind],disp[ind],disperr[ind],psym=1,$
  title='Velocity dispersion values for the '+region+' region',$
  xtitle='Radius (arcsec)',ytitle='Velocity Dispersion (km/sec)',$
  xrange=[xdown,xup],$
  xthick=3,ythick=3,thick=3,xstyle=1,charthick=3,errthick=3

plots,0.15,0.91,psym=2,/normal,thick=3
xyouts,0.18,0.90,'01',/normal,charthick=3

loadct,4
for j=1,4 do begin
    ind = where(ang eq angs[j])
    oploterror,radius[ind],disp[ind],disperr[ind],psym=1,color=colors[j],$
      errcolor=colors[j],errthick=3,thick=3
    plots,0.15,0.91-(j*0.03),psym=2,/normal,thick=3,color=colors[j]
    xyouts,0.18,0.90-(j*0.03),'0'+strn(angs[j]),/normal,$
      charthick=3,color=colors[j]
endfor
device,/close_file

;-------------------------------------------------------------
;                           H3
;-------------------------------------------------------------
set_plot,'ps'
device,file='H3_'+region+'.ps',/color
loadct,0
ind = where(ang eq angs[0])
ploterror,radius[ind],h3[ind],h3err[ind],psym=1,$
  title='H3 values for the '+region+' region',$
  xtitle='Radius (arcsec)',ytitle='H3',xrange=[xdown,xup],$
  xthick=3,ythick=3,thick=3,xstyle=1,charthick=3,errthick=3
plots,0.15,0.91,psym=2,/normal,thick=3
xyouts,0.18,0.90,'01',/normal,charthick=3

loadct,4
for j=1,4 do begin
    ind = where(ang eq angs[j])
    oploterror,radius[ind],h3[ind],h3err[ind],psym=1,color=colors[j],$
      errcolor=colors[j],errthick=3,thick=3
    plots,0.15,0.91-(j*0.03),psym=2,/normal,thick=3,color=colors[j]
    xyouts,0.18,0.90-(j*0.03),'0'+strn(angs[j]),/normal,$
      charthick=3,color=colors[j]
endfor
device,/close_file

;-------------------------------------------------------------
;                            H4
;-------------------------------------------------------------
set_plot,'ps'
device,file='H4_'+region+'.ps',/color
loadct,0
ind = where(ang eq angs[0])
ploterror,radius[ind],h4[ind],h4err[ind],psym=1,$
  title='H4 values for the '+region+' region',$
  xtitle='Radius (arcsec)',ytitle='H4',xrange=[xdown,xup],$
  xthick=3,ythick=3,thick=3,xstyle=1,charthick=3,errthick=3
plots,0.15,0.91,psym=2,/normal,thick=3
xyouts,0.18,0.90,'01',/normal,charthick=3

loadct,4
for j=1,4 do begin
    ind = where(ang eq angs[j])
    oploterror,radius[ind],h4[ind],h4err[ind],psym=1,color=colors[j],$
      errcolor=colors[j],errthick=3,thick=3
    plots,0.15,0.91-(j*0.03),psym=2,/normal,thick=3,color=colors[j]
    xyouts,0.18,0.90-(j*0.03),'0'+strn(angs[j]),/normal,$
      charthick=3,color=colors[j]
endfor
device,/close_file

stop
END

PRO Ppallmc

;*************************************************************************
;plotting parameters
y1 = [-50,70];velocity
y2 = [270,370];disp
y3 = [-0.09,0.09];h3
y4 = [-0.09,0.07];h4

xup = 300
xdn = -300

galname = 'M87'
pmjr = 'no';change to anything other than 'no' if you want to plot the major axis only. (this supresses all off-major-axis values)
pw = 0.8 ;plot width. set this to a smaller number and your figure gets narrower. this is in normalized coordinates.
;error = 1 ;standard deviation of the 4 spectral region values
error = 2 ;the mean of the error of the pallmc output.

;readcol,'bin.radiusmajor',format='a,f',bin2,rad
readcol,'bin.radius',format='a,f',bin2,rad
outta = 'outta'
outta = 'inna'
;*************************************************************************

readcol,'pfitlov.Mlosvd',format='a,f,f,x,f,f,x,x,x',bin1,vel,disp,h3,h4
n0 = n_elements(bin1)

n1 = n_elements(bin2)

bin1t = strarr(n0)
bin2t = strarr(n1)

for j=0,n0-1 do begin
    t = strsplit(bin1[j],'_.',/extract)
    bin1t[j] = t[0]
endfor
for j=0,n1-1 do begin
    t = strsplit(bin2[j],'.',/extract)
    bin2t[j] = t[0]
endfor

i = where(bin1t eq bin2t)
binp = bin1t[i]
n2 = n_elements(binp)

radt = rad[i]

radlin = fltarr(n2)
i1 = 0
i2 = 0
i3 = 0
i4 = 0
i5 = 0

;**********************************************************
readcol,'pallmc.gb.out',f='f,f,f,f,f,f,f,f,f,x',r,v1,ve1,s1,se1,h31,h3e1,h41,h4e1
readcol,'pallmc.hb.out',f='f,f,f,f,f,f,f,f,f,x',r,v2,ve2,s2,se2,h32,h3e2,h42,h4e2
readcol,'pallmc.mgw.out',f='f,f,f,f,f,f,f,f,f,x',r,v3,ve3,s3,se3,h33,h3e3,h43,h4e3
readcol,'pallmc.mgwo.out',f='f,f,f,f,f,f,f,f,f,x',r,v4,ve4,s4,se4,h34,h3e4,h44,h4e4
readcol,'pallmc.fe.out',f='f,f,f,f,f,f,f,f,f,x',r,v5,ve5,s5,se5,h35,h3e5,h45,h4e5

if (error eq 1) then begin ;SD of the values of each spectral region
    Evel = [[v1],[v2],[v4],[v5]]
    Edisp = [[s1],[s2],[s4],[s5]]
    Eh3 = [[h31],[h32],[h34],[h35]]
    Eh4 = [[h41],[h42],[h44],[h45]]
endif
if (error eq 2) then begin ;mean of the error from pallmc.out
    Evel = [[ve1],[ve2],[ve4],[ve5]]
    Edisp = [[se1],[se2],[se4],[se5]]
    Eh3 = [[h3e1],[h3e2],[h3e4],[h3e5]]
    Eh4 = [[h4e1],[h4e2],[h4e4],[h4e5]]
endif

ihigh = where(disp gt 420,count) ;rejection based on dispersion outliers
if (count gt 0) then Evel[ihigh] = !values.f_nan
if (count gt 0) then Edisp[ihigh] = !values.f_nan
if (count gt 0) then Eh3[ihigh] = !values.f_nan
if (count gt 0) then Eh4[ihigh] = !values.f_nan

ihigh = where(abs(vel) gt 100,count) ;rejection based on velocity outliers
if (count gt 0) then Evel[ihigh] = !values.f_nan
if (count gt 0) then Edisp[ihigh] = !values.f_nan
if (count gt 0) then Eh3[ihigh] = !values.f_nan
if (count gt 0) then Eh4[ihigh] = !values.f_nan

sdV = fltarr(n1)
sdS = fltarr(n1)
sdH3 = fltarr(n1)
sdH4 = fltarr(n1)
if (error eq 1) then begin
    for j=0,n1-1 do sdV[j] =  stddev(Evel[j,*],/nan)
    for j=0,n1-1 do sdS[j] =  stddev(Edisp[j,*],/nan)
    for j=0,n1-1 do sdH3[j] = stddev(EH3[j,*],/nan)
    for j=0,n1-1 do sdH4[j] = stddev(EH4[j,*],/nan)
endif
if (error eq 2) then begin
    for j=0,n1-1 do sdV[j] =  max(Evel[j,*],/nan)
    for j=0,n1-1 do sdS[j] =  max(Edisp[j,*],/nan)
    for j=0,n1-1 do sdH3[j] = max(EH3[j,*],/nan)
    for j=0,n1-1 do sdH4[j] = max(EH4[j,*],/nan)
endif

    EV = sdV;*0.5
    ES = sdS;*0.5
    EH3 = sdH3;*0.5
    EH4 = sdH4;*0.5

if (outta eq 'outta') then begin
f1 = '(a7,f10.5,f10.5,f10.5,f11.5,f9.5,f10.5,f10.5,f10.5,f10.5)'
openw,5,'Moments_+_errors.txt'
for j=0,n0-1 do begin
    t = strsplit(bin1[j],'.',/extract)
    t = t[0]
    printf,5,t,rad[j],vel[j],EV[j],disp[j],ES[j],h3[j],Eh3[j],h4[j],Eh4[j],format=f1
endfor
free_lun,5
stop
endif

;**********************************************************

for j=0,n2-1 do begin
    t = strsplit(binp[j],'0123456789',/extract)
    tt = strsplit(binp[j],'n0',/extract)
    if (t[0] eq 'bin') then radlin[j] = radt[j]
    if (t[0] eq 'bnn') then radlin[j] = radt[j] * (-1.0)
    if (tt[2] eq 1) then i1 = [i1,j]
    if (tt[2] eq 2) then i2 = [i2,j]
    if (tt[2] eq 3) then i3 = [i3,j]
    if (tt[2] eq 4) then i4 = [i4,j]
    if (tt[2] eq 5) then i5 = [i5,j]
endfor
i1 = i1[1:*]
i2 = i2[1:*]
i3 = i3[1:*]
i4 = i4[1:*]
i5 = i5[1:*]
icolor = [60,110,180,150]

ibad = where(vel gt y1[1] or vel lt y1[0],count)
if (count ne 0) then vel[ibad] = !Values.F_NAN
ibad = where(disp gt y2[1] or disp lt y2[0],count)
if (count ne 0) then disp[ibad] = !Values.F_NAN
ibad = where(h3 gt y3[1] or h3 lt y3[0],count)
if (count ne 0) then h3[ibad] = !Values.F_NAN
ibad = where(h4 gt y4[1] or h4 lt y4[0],count)
if (count ne 0) then h4[ibad] = !Values.F_NAN

set_plot,'ps'
device,file=galname+'_moments.ps',/color

!p.region = [0.1,0.1,0.9,0.9]
!p.multi = [0,1,4,0,1]
!y.omargin = [2,4]

loadct,0
plot,radlin,vel[i],psym=sym(4),charsize=1.0,xthick=2,ythick=2,thick=3,$
  charthick=1.5,xrange=[xdn,xup],xstyle=1,yrange=[y1],ystyle=1,$
  position=[0.15,0.675,pw,0.85],symsize=0.5,/nodata,xtickformat='(A1)',$
  xticklen=0.07,yminor=2
xyouts,0.115,0.7,'Velocity (km/sec)',/normal,orientation=90,charsize=0.5,$
  charthick=1.5
if (pmjr eq 'no') then begin
    loadct,4
    plots,radlin[i2],vel[i2],psym=sym(9),color=60,thick=3,symsize=0.5
    plots,radlin[i3],vel[i3],psym=sym(9),color=110,thick=3,symsize=0.5
    plots,radlin[i4],vel[i4],psym=sym(9),color=180,thick=3,symsize=0.5
    plots,radlin[i5],vel[i5],psym=sym(9),color=150,thick=3,symsize=0.5
    loadct,0
endif
plots,radlin[i1],vel[i1],psym=sym(4),thick=3,symsize=0.5
errplot,radlin[i1],vel[i1]-EV,vel[i1]+EV

loadct,0
plot,radlin,disp[i],psym=sym(4),charsize=1.0,xthick=2,ythick=2,thick=3,$
  charthick=1.5,xrange=[xdn,xup],xstyle=1+8,yrange=[y2],ystyle=1,$
  position=[0.15,0.50,pw,0.675],symsize=0.5,/nodata,xtickformat='(A1)',$
  xticklen=0.07,yminor=2
xyouts,0.115,0.51,'Dispersion (km/sec)',/normal,orientation=90,charsize=0.5,$
  charthick=1.5
if (pmjr eq 'no') then begin
    loadct,4
    plots,radlin[i2],disp[i2],psym=sym(9),color=60,thick=3,symsize=0.5
    plots,radlin[i3],disp[i3],psym=sym(9),color=110,thick=3,symsize=0.5
    plots,radlin[i4],disp[i4],psym=sym(9),color=180,thick=3,symsize=0.5
    plots,radlin[i5],disp[i5],psym=sym(9),color=150,thick=3,symsize=0.5
    loadct,0
endif
plots,radlin[i1],disp[i1],psym=sym(4),thick=3,symsize=0.5
errplot,radlin[i1],disp[i1]-ES,disp[i1]+ES

loadct,0
plot,radlin,h3[i],psym=sym(4),charsize=1.0,xthick=2,ythick=2,thick=3,$
  charthick=1.5,xrange=[xdn,xup],xstyle=1+8,yrange=[y3],ystyle=1,$
  position=[0.15,0.325,pw,0.5],symsize=0.5,/nodata,xtickformat='(A1)',$
  xticklen=0.07,yminor=2
xyouts,0.115,0.412,'H3',/normal,orientation=90,charsize=0.5,charthick=1.5
if (pmjr eq 'no') then begin
    loadct,4
    plots,radlin[i2],h3[i2],psym=sym(9),color=60,thick=3,symsize=0.5
    plots,radlin[i3],h3[i3],psym=sym(9),color=110,thick=3,symsize=0.5
    plots,radlin[i4],h3[i4],psym=sym(9),color=180,thick=3,symsize=0.5
    plots,radlin[i5],h3[i5],psym=sym(9),color=150,thick=3,symsize=0.5
    loadct,0
endif
plots,radlin[i1],h3[i1],psym=sym(4),thick=3,symsize=0.5
errplot,radlin[i1],h3[i1]-EH3,h3[i1]+EH3

loadct,0
plot,radlin,h4[i],psym=sym(4),charsize=1.0,xthick=2,ythick=2,thick=3,$
  charthick=1.5,xrange=[xdn,xup],xstyle=1+8,yrange=[y4],ystyle=1,$
  position=[0.15,0.15,pw,0.325],symsize=0.5,/nodata,xticklen=0.07,$
  yminor=2,xtitle='Radius (arcsec)'
xyouts,0.115,0.237,'H4',/normal,orientation=90,charsize=0.5,charthick=1.5
if (pmjr eq 'no') then begin
    loadct,4
    plots,radlin[i2],h4[i2],psym=sym(9),color=60,thick=3,symsize=0.5
    plots,radlin[i3],h4[i3],psym=sym(9),color=110,thick=3,symsize=0.5
    plots,radlin[i4],h4[i4],psym=sym(9),color=180,thick=3,symsize=0.5
    plots,radlin[i5],h4[i5],psym=sym(9),color=150,thick=3,symsize=0.5
    loadct,0
endif
plots,radlin[i1],h4[i1],psym=sym(4),thick=3,symsize=0.5
errplot,radlin[i1],h4[i1]-EH4,h4[i1]+EH4

device,/close_file
set_plot,'x'

stop
END

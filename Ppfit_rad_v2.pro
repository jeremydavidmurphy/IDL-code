PRO Ppfit_rad_v2, file

; This code is just Ppfit_rad.pro, only instead of coloring by angular
; bin, it colors by spectral region. The bin.radius.5 text file is
; just the bin.radius file, only doubled up until it matches the
; length of the concatenated pfitlov.out file
;*************************************************************************
;plotting parameters
ydn1 = -100;velocity
yup1 = 100;velocity
ydn3 = 100;disp
yup3 = 300;disp
ydn2 = -0.1;h3
yup2 = 0.1;h3
ydn4 = -0.1;h4
yup4 = 0.1;h4

xup = 100
xdn = -100

pltmjr = 'no' ;set to 'yes' if you want to emphasize the major axis (i.e. all bin##01's) 

;galname = 'M87'

pw = 0.9 ;plot width. set this to a smaller number and your figure gets narrower. this is in normalized coordinates.
;*************************************************************************

; FILE: The name of the pfitlov.out file

; GALNAME: If you don't enter this, you will be prompted to do so by
; the code. This is just for titles in the figures and naming of the
; output files.

if (n_elements(galname) eq 0) then begin
    galname = ''
    print,'Enter the name of the galaxy:'
    read,galname
endif

readcol,file,format='a,f,f,x,f,f,x,x,x',bin1,vel,disp,h3,h4
n0 = n_elements(bin1)

readcol,'bin.radius.4',format='a,i,f',bin2,index,rad
n1 = n_elements(bin2)

if n0 ne n1 then begin
    print,'Your bin list is not the same length as your pfitlov list'
    stop
endif

ineg = where(bin2 eq 'bin')
rad[ineg] = rad[ineg] * (-1.0)

i1 = where(index eq 01)
i2 = where(index eq 02)
i3 = where(index eq 03)
i4 = where(index eq 04)
i5 = where(index eq 05)

icolor = [60,110,180,150]

ibad = where(vel gt yup1 or vel lt ydn1,count)
if (count ne 0) then vel[ibad] = !Values.F_NAN
ibad = where(disp gt yup3 or disp lt ydn3,count)
if (count ne 0) then disp[ibad] = !Values.F_NAN
ibad = where(h3 gt yup2 or h3 lt ydn2,count)
if (count ne 0) then h3[ibad] = !Values.F_NAN
ibad = where(h4 gt yup4 or h4 lt ydn4,count)
if (count ne 0) then h4[ibad] = !Values.F_NAN

set_plot,'ps'
device,file=galname+'_moments.ps',/color

!p.region = [0.1,0.1,0.9,0.9]
!p.multi = [0,1,4,0,1]
!y.omargin = [2,4]

loadct,0
plot,rad,vel,psym=sym(4),charsize=1.0,xthick=2,ythick=2,thick=3,$
  charthick=1.5,xrange=[xdn,xup],xstyle=1,yrange=[ydn1,yup1],ystyle=1,$
  position=[0.15,0.675,pw,0.85],symsize=0.5,/nodata,xtickformat='(A1)',$
  xticklen=0.07,yminor=2,title=galname
xyouts,0.115,0.7,'Velocity (km/sec)',/normal,orientation=90,charsize=0.5,$
  charthick=1.5
loadct,4
plots,rad[i2],vel[i2],psym=sym(9),color=60,thick=3,symsize=0.5
plots,rad[i3],vel[i3],psym=sym(9),color=110,thick=3,symsize=0.5
plots,rad[i4],vel[i4],psym=sym(9),color=180,thick=3,symsize=0.5
plots,rad[i5],vel[i5],psym=sym(9),color=150,thick=3,symsize=0.5
loadct,0
plots,rad[i1],vel[i1],psym=sym(4),thick=3,symsize=0.5

loadct,0
plot,rad,disp,psym=sym(4),charsize=1.0,xthick=2,ythick=2,thick=3,$
  charthick=1.5,xrange=[xdn,xup],xstyle=1+8,yrange=[ydn3,yup3],ystyle=1,$
  position=[0.15,0.50,pw,0.675],symsize=0.5,/nodata,xtickformat='(A1)',$
  xticklen=0.07,yminor=2
xyouts,0.115,0.51,'Dispersion (km/sec)',/normal,orientation=90,charsize=0.5,$
  charthick=1.5
loadct,4
plots,rad[i2],disp[i2],psym=sym(9),color=60,thick=3,symsize=0.5
plots,rad[i3],disp[i3],psym=sym(9),color=110,thick=3,symsize=0.5
plots,rad[i4],disp[i4],psym=sym(9),color=180,thick=3,symsize=0.5
plots,rad[i5],disp[i5],psym=sym(9),color=150,thick=3,symsize=0.5
loadct,0
plots,rad[i1],disp[i1],psym=sym(4),thick=3,symsize=0.5

loadct,0
plot,rad,h3,psym=sym(4),charsize=1.0,xthick=2,ythick=2,thick=3,$
  charthick=1.5,xrange=[xdn,xup],xstyle=1+8,yrange=[ydn2,yup2],ystyle=1,$
  position=[0.15,0.325,pw,0.5],symsize=0.5,/nodata,xtickformat='(A1)',$
  xticklen=0.07,yminor=2
xyouts,0.115,0.412,'H3',/normal,orientation=90,charsize=0.5,charthick=1.5
loadct,4
plots,rad[i2],h3[i2],psym=sym(9),color=60,thick=3,symsize=0.5
plots,rad[i3],h3[i3],psym=sym(9),color=110,thick=3,symsize=0.5
plots,rad[i4],h3[i4],psym=sym(9),color=180,thick=3,symsize=0.5
plots,rad[i5],h3[i5],psym=sym(9),color=150,thick=3,symsize=0.5
loadct,0
plots,rad[i1],h3[i1],psym=sym(4),thick=3,symsize=0.5

loadct,0
plot,rad,h4,psym=sym(4),charsize=1.0,xthick=2,ythick=2,thick=3,$
  charthick=1.5,xrange=[xdn,xup],xstyle=1+8,yrange=[ydn4,yup4],ystyle=1,$
  position=[0.15,0.15,pw,0.325],symsize=0.5,/nodata,xticklen=0.07,$
  yminor=2,xtitle='Radius (arcsec)'
xyouts,0.115,0.237,'H4',/normal,orientation=90,charsize=0.5,charthick=1.5
loadct,4
plots,rad[i2],h4[i2],psym=sym(9),color=60,thick=3,symsize=0.5
plots,rad[i3],h4[i3],psym=sym(9),color=110,thick=3,symsize=0.5
plots,rad[i4],h4[i4],psym=sym(9),color=180,thick=3,symsize=0.5
plots,rad[i5],h4[i5],psym=sym(9),color=150,thick=3,symsize=0.5
loadct,0
plots,rad[i1],h4[i1],psym=sym(4),thick=3,symsize=0.5

device,/close_file
set_plot,'x'

stop
END

PRO frplot, fbrnbr


ans=''
ftemp = dblarr(30,8)
files = findfile('fib'+fbrnbr+'*.frd',count=num)
fvals = dblarr(30,num)
frange = dblarr(30,num,6)
ee = dblarr(30,num)
FOR k=0,num-1 DO BEGIN
    openr,5,files[k]
    readf,5,ftemp
    ee[*,k] = ftemp[*,0]
    fvals[*,k] = ftemp[*,7]
    FOR j=0,5 DO frange[*,k,j] = ftemp[*,j+1]
    free_lun,5
ENDFOR


xulim = max(fvals)+0.1
xllim = min(fvals)-0.1

lx1=[3.65,3.65]
ly1=[0.6,1.1]
lx337=[fvals[25,0],fvals[25,0]]
lx400=[fvals[25,1],fvals[25,1]]
lx600=[fvals[25,2],fvals[25,2]]

window,0,retain=2
device,decomposed=0
loadct,0
plot,fvals[*,0],ee[*,0],xrange=[xulim,xllim],yrange=[0.6,1.1],title='FRD Wavelength Dependence for Fiber '+fbrnbr,xstyle=1, $
  xtitle='Output F/#',ytitle='Encircled Energy',charsize=1.5,/nodata
oplot,lx1,ly1

loadct,4
red=150
blue=60
green=110
aqua=90
orange=185
clr = [blue,green,red]

FOR k=0,num-1 DO BEGIN
    oplot,fvals[*,k],ee[*,k],psym=4,color=clr[k]
ENDFOR
oplot,lx337,ly1,color=blue
plots,xulim-0.1,1.05,psym=4,symsize=2,color=blue
xyouts,xulim-0.12,1.045,': 3370A',color=blue,charsize=1.5
oplot,lx400,ly1,color=green
plots,xulim-0.1,1.03,psym=4,symsize=2,color=green
xyouts,xulim-0.12,1.025,': 4000A',color=green,charsize=1.5
oplot,lx600,ly1,color=red
plots,xulim-0.1,1.01,psym=4,symsize=2,color=red
xyouts,xulim-0.12,1.005,': 6000A',color=red,charsize=1.5

xyouts,fvals[25,0]-0.05,0.62,'95% ENCIRCLED ENERGY FOR 3 BANDS',orientation=90,charsize=1.3,color=blue

loadct,0

xyouts,3.6,0.62,'INPUT BEAM F/#',orientation=90,charsize=1.3

;print,'Save the file?'
;read,ans
ans='y'
IF (ans EQ 'y') THEN BEGIN
    set_plot,'ps'
    device,filename='F'+fbrnbr+'eeplot.ps',/color
    loadct,0
    plot,fvals[*,0],ee[*,0],psym=3,xrange=[4.25,2.7],yrange=[0.6,1.1], $
      title='FRD Wavelength Dependence for Fiber '+fbrnbr+' (WITH cover plates)',xstyle=1,xtitle='Output F/#', $
      ytitle='Encircled Energy',charthick=2.3,/nodata
    oplot,lx1,ly1,thick=2

    loadct,4
    FOR k=0,num-1 DO BEGIN
        oplot,fvals[*,k],ee[*,k],psym=4,color=clr[k],thick=2.5
    ENDFOR
    oplot,lx337,ly1,color=blue,thick=3
    plots,xulim-0.1,1.05,psym=4,symsize=2,color=blue,thick=2
    xyouts,xulim-0.12,1.042,' 3370A',color=blue,charsize=1.5,charthick=2
    oplot,lx400,ly1,color=green,thick=3
    plots,xulim-0.1,1.02,psym=4,symsize=2,color=green,thick=2
    xyouts,xulim-0.12,1.012,' 4000A',color=green,charsize=1.5,charthick=2
    oplot,lx600,ly1,color=red,thick=3
    plots,xulim-0.1,0.99,psym=4,symsize=2,color=red,thick=2
    xyouts,xulim-0.12,0.982,' 6000A',color=red,charsize=1.5,charthick=2
    
;    xyouts,fvals[25,0]-0.04,0.62,'95% ENCIRCLED ENERGY FOR 3370A',orientation=90,charsize=1.5,color=blue,charthick=2
    xyouts,fvals[25,0]-0.055,0.62,'F/#: '+strn(fvals[25,0]),orientation=90,charsize=1.5,color=blue,charthick=2.6
    
    loadct,0
    xyouts,3.61,0.62,'INPUT BEAM F/3.65',orientation=90,charsize=0.9,charthick=2

    device,/close_file
    set_plot,'x'
ENDIF ELSE print,'Word.  No plot for you!'
pause
wdelete,0
stop
END

    
    

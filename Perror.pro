pro Perror

restore,'kin.idl'

device,file='M87Vel.ps',/color
loadct,4
ploterror,radF,velocity,velocitye,xtitle='Radius (arcsec)',ytitle='Velocity (km/sec)',$
  psym=sym(1),thick=3,xthick=3,ythick=3,charthick=3,type=2,errthick=3,title='M87 Velocity',$
/nodata,xrange=[0,300],style=1
oploterror,rad1,vel1,vel1e,psym=sym(1),errthick=3
oploterror,rad2,vel2,vel2e,psym=sym(1),errthick=3,color=150
device,/close_file

device,file='M87Disp.ps',/color
loadct,4
ploterror,radF,dispersion,dispersione,xtitle='Radius (arcsec)',ytitle='Dispersion (km/sec)',$
  psym=sym(1),thick=3,xthick=3,ythick=3,charthick=3,type=2,errthick=3,title='M87 Dispersion',$
/nodata,xrange=[0,300],xtyle=1
oploterror,rad1,disp1,disp1e,psym=sym(1),errthick=3
oploterror,rad2,disp2,disp2e,psym=sym(1),errthick=3,color=150
device,/close_file

stop
end

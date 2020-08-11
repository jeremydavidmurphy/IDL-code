set_plot,'ps'
device,filename='Allcolor_ew.ps',/color
loadct,0
plot,radius,m87data[0,*],psym=2,title=lines[0],xtitle='!3Radius (arcsec)',ytitle='EW (A)',charsize=1.2,yrange=[6,11],symsize=2,ystyle=1,thick=2,charthick=1.5,color=60,/nodata
loadct,4
plot,radius,m87data[0,*],psym=2,title=lines[0],xtitle='!3Radius (arcsec)',ytitle='EW (A)',charsize=1.2,yrange=[6,11],symsize=2,ystyle=1,thick=2,charthick=1.5,color=60
loadct,0
plot,radius,m87data[1,*],psym=2,title=lines[1],xtitle='!3Radius (arcsec)',ytitle='EW (A)',charsize=1.2,yrange=[6,11],symsize=2,ystyle=1,thick=2,charthick=1.5,color=60,/nodata
loadct,4
plot,radius,m87data[1,*],psym=2,title=lines[1],xtitle='!3Radius (arcsec)',ytitle='EW (A)',charsize=1.2,yrange=[6,11],symsize=2,ystyle=1,thick=2,charthick=1.5,color=60
loadct,0
plot,radius,m87data[5,*],psym=2,title=lines[5],xtitle='!3Radius (arcsec)',ytitle='EW (A)',charsize=1.2,yrange=[3,6],symsize=2,ystyle=1,thick=2,charthick=1.5,color=110,/nodata
loadct,4
plot,radius,m87data[5,*],psym=2,title=lines[5],xtitle='!3Radius (arcsec)',ytitle='EW (A)',charsize=1.2,yrange=[3,6],symsize=2,ystyle=1,thick=2,charthick=1.5,color=110
loadct,0
plot,radius,m87data[7,*],psym=2,title=lines[7],xtitle='!3Radius (arcsec)',ytitle='EW (A)',charsize=1.2,yrange=[3,5],symsize=2,ystyle=1,thick=2,charthick=1.5,color=140,/nodata
loadct,4
plot,radius,m87data[7,*],psym=2,title=lines[7],xtitle='!3Radius (arcsec)',ytitle='EW (A)',charsize=1.2,yrange=[3,5],symsize=2,ystyle=1,thick=2,charthick=1.5,color=140
loadct,0
plot,radius,m87data[8,*],psym=2,title=lines[8],xtitle='!3Radius (arcsec)',ytitle='EW (A)',charsize=1.2,yrange=[0.5,2.5],symsize=2,ystyle=1,thick=2,charthick=1.5,color=200,/nodata
loadct,4
plot,radius,m87data[8,*],psym=2,title=lines[8],xtitle='!3Radius (arcsec)',ytitle='EW (A)',charsize=1.2,yrange=[0.5,2.5],symsize=2,ystyle=1,thick=2,charthick=1.5,color=200
loadct,0
plot,radius,m87data[3,*],psym=2,title=lines[3],xtitle='!3Radius (arcsec)',ytitle='EW (A)',charsize=1.2,yrange=[0.5,2.5],symsize=2,ystyle=1,thick=2,charthick=1.5,color=140,/nodata
loadct,4
plot,radius,m87data[3,*],psym=2,title=lines[3],xtitle='!3Radius (arcsec)',ytitle='EW (A)',charsize=1.2,yrange=[0.5,2.5],symsize=2,ystyle=1,thick=2,charthick=1.5,color=140
device,/close

set_plot,'ps'
device,filename='Allcolor_ew.ps',/color
loadct,4
plot,radius,m87data[0,*],psym=2,title=lines[0],xtitle='!3Radius (arcsec)',ytitle='EW (A)',charsize=1.2,yrange=[6,11],symsize=2,ystyle=1,thick=2,charthick=1.5,color=60
plot,radius,m87data[1,*],psym=2,title=lines[1],xtitle='!3Radius (arcsec)',ytitle='EW (A)',charsize=1.2,yrange=[6,11],symsize=2,ystyle=1,thick=2,charthick=1.5,color=60
plot,radius,m87data[5,*],psym=2,title=lines[5],xtitle='!3Radius (arcsec)',ytitle='EW (A)',charsize=1.2,yrange=[3,6],symsize=2,ystyle=1,thick=2,charthick=1.5,color=110
plot,radius,m87data[7,*],psym=2,title=lines[7],xtitle='!3Radius (arcsec)',ytitle='EW (A)',charsize=1.2,yrange=[3,5],symsize=2,ystyle=1,thick=2,charthick=1.5,color=170
plot,radius,m87data[8,*],psym=2,title=lines[8],xtitle='!3Radius (arcsec)',ytitle='EW (A)',charsize=1.2,yrange=[0.5,2.5],symsize=2,ystyle=1,thick=2,charthick=1.5,color=200
plot,radius,m87data[3,*],psym=2,title=lines[3],xtitle='!3Radius (arcsec)',ytitle='EW (A)',charsize=1.2,yrange=[0.5,2.5],symsize=2,ystyle=1,thick=2,charthick=1.5,color=140
device,/close

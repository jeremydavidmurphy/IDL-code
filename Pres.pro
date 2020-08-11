; This routine was used to generate an analysis of the VIRUS-P IR from
; a number of different sources. It's primarily a wrapper for
; rescheck2F.pro


PRO Pres

arc1 = 'iron_arcpe.fits'
ptow1 = 'ptow_sept10_n1.dat'
mask1 = 'mask_sept10_n1.dat'

;out1 = rescheck3F(arc1,ptow1,mask1,60)
out1 = rescheck3F(arc1,ptow1,mask1,60,splot='on')

arc2 = 'arc_sept10_n1e.fits'
;arc2 = 'hgcd_arcpe.fits'
ptow2 = 'ptow_sept10_n1.dat'
mask2 = 'mask_sept10_n1.dat'

out2 = rescheck3F(arc2,ptow2,mask2,15)
;out2 = rescheck3F(arc2,ptow2,mask2,12,splot='on')

arc3 = 'arc_M29_n1e.fits'
ptow3 = 'ptow_M29_n1.dat'
mask3 = 'mask_M29_n1.dat'

out3 = rescheck3F(arc3,ptow3,mask3,15)
;out3 = rescheck3F(arc3,ptow3,mask3,15,splot='on')

arc4 = 'arc_may10_n1e.fits'
ptow4 = 'ptow_may10_n1.dat'
mask4 = 'mask_may10_n1.dat'

out4 = rescheck3F(arc4,ptow4,mask4,15)
;out4 = rescheck3F(arc4,ptow4,mask4,15,splot='on')

arc5 = 'pr0136pe.fits'
ptow5 = 'ptow_M29_n1.dat'
mask5 = 'mask_M29_n1.dat'

out5 = rescheck3F(arc5,ptow5,mask5,3)
;out5 = rescheck3F(arc5,ptow5,mask5,5,splot='on')



readcol,'ptop_may10_n5s.dat',silent=1,a,b,c,d
ptop = fltarr(2048,246)
for j=0,2047 do for k=0,245 do ptop[j,k]= a[k] + b[k]*j + c[k]*j*j + d[k]*j*j*j
readcol,'ptow_may10_n5s.dat',silent=1,a,b,c,d,e
pwave = fltarr(2048,246)
for j=0,2047 do for k=0,245 do pwave[j,k]= a[k] + b[k]*j + c[k]*j*j + d[k]*j*j*j + e[k]*j*j*j*j

mptop = median(ptop,dim=2)
mpwave = median(pwave,dim=2)
i = where(ptop[0,*] eq max(ptop[0,*]))
upy = ptop[*,i]
upx = pwave[*,i]
i = where(ptop[0,*] eq min(ptop[0,*]))
dny = ptop[*,i]
dnx = pwave[*,i]

window,2,retain=2
device,decomposed=0

loadct,0,/silent
plot,out1[*,0,0],out1[*,0,3],psym=2,yrange=[50,250],/nodata,xrange=[3500,6000]
loadct,4,/silent

for j=0,245 do oplot,out1[*,j,0],out1[*,j,4],psym=sym(1),symsize=0.5,color=60
for j=0,245 do oplot,out2[*,j,0]+12,out2[*,j,4],psym=sym(1),symsize=0.5,color=110
for j=0,245 do oplot,out3[*,j,0],out3[*,j,4],psym=sym(1),symsize=0.5,color=150
for j=0,245 do oplot,out4[*,j,0]-12,out4[*,j,4],psym=sym(1),symsize=0.5,color=180
for j=0,245 do plots,out5[*,j,0],out5[*,j,4],psym=sym(1),symsize=0.5,color=255

oplot,mpwave,mptop,thick=3,color=255
oplot,dnx,dny,linestyle=1,color=255
oplot,upx,upy,linestyle=1,color=255
loadct,0,/silent
oplot,mpwave,299792.548*2.3/mpwave,thick=1.5,linestyle=2
oplot,mpwave,299792.548*2.3/mpwave+12,thick=1.5,linestyle=2
oplot,mpwave,299792.548*2.3/mpwave-12,thick=1.5,linestyle=2
pause

loadct,0,/silent
plot,out1[*,0,0],out1[*,0,3],psym=2,yrange=[3,7],/nodata,xrange=[3500,6000]
loadct,4,/silent

for j=0,245 do oplot,out1[*,j,0],out1[*,j,3],psym=sym(1),symsize=0.5,color=60
pause
for j=0,245 do oplot,out2[*,j,0]+25,out2[*,j,3],psym=sym(1),symsize=0.5,color=110
pause
for j=0,245 do oplot,out3[*,j,0],out3[*,j,3],psym=sym(1),symsize=0.5,color=150
pause
for j=0,245 do oplot,out4[*,j,0]-25,out4[*,j,3],psym=sym(1),symsize=0.5,color=180
pause
for j=0,245 do plots,out5[*,j,0],out5[*,j,3],psym=sym(1),symsize=0.5,color=255

;oplot,mpwave,mptop,thick=3,color=255
;oplot,dnx,dny,linestyle=1,color=255
;oplot,upx,upy,linestyle=1,color=255
;loadct,0,/silent
;oplot,mpwave,299792.548*2.3/mpwave,thick=1.5,linestyle=2
;oplot,mpwave,299792.548*2.3/mpwave+12,thick=1.5,linestyle=2
;oplot,mpwave,299792.548*2.3/mpwave-12,thick=1.5,linestyle=2
pause

set_plot,'ps'
device,file='VIRUS-P_IR.ps',/color

loadct,0,/silent
plot,out1[*,0,0],out1[*,0,3],psym=2,yrange=[3,7],xthick=2,/ystyle,$
  ythick=2,charthick=2,xtitle='Wavelength (A)',$
  title='VIRUS-P Instrumental Resolution',$
  ytitle='IR (FWHM in Angstroms)',/nodata,xrange=[3500,6000]
loadct,4,/silent

for j=0,245 do oplot,out1[*,j,0],out1[*,j,3],psym=sym(1),color=60,symsize=0.3
for j=0,245 do oplot,out2[*,j,0]+12,out2[*,j,3],psym=sym(1),color=110,symsize=0.3
for j=0,245 do oplot,out3[*,j,0],out3[*,j,3],psym=sym(1),color=150,symsize=0.3
for j=0,245 do oplot,out4[*,j,0]-12,out4[*,j,3],psym=sym(1),color=180,symsize=0.3
loadct,0,/silent
for j=0,245 do plots,out5[*,j,0],out5[*,j,3],psym=sym(1),symsize=0.3
loadct,4,/silent

;oplot,mpwave,mptop,thick=5,color=180
;oplot,dnx,dny,linestyle=2,color=180,thick=2
;oplot,upx,upy,linestyle=2,color=180,thick=2
;loadct,0,/silent
;oplot,mpwave,299792.548*2.3/mpwave,thick=5,linestyle=2
device,/close_file

set_plot,'x'

;x = [out4[0,*,0],out4[2:3,*,0],out4[5,*,0],out4[7:14,*,0],out5[*,*,0]]
;y = [out4[0,*,3],out4[2:3,*,3],out4[5,*,3],out4[7:14,*,3],out5[*,*,3]]
x = [out1[*,*,0],out2[*,*,0]+25,out3[*,*,0],out4[*,*,0]-25,out5[*,*,0]]
y = [out1[*,*,3],out2[*,*,3],out3[*,*,3],out4[*,*,3],out5[*,*,3]]


window,0,retain=2
device,decomposed=0

loadct,0,/silent
plot,out4[*,0,0],out4[*,0,3],psym=2,yrange=[3,7],/nodata,xrange=[3500,6000]
loadct,4,/silent

for j=0,245 do begin
    oplot,x[*,j],y[*,j],color=j,psym=sym(1),symsize=0.5
;    wait,0.1
endfor

yy = median(y,dim=2,/even)

pause
while !d.window ne -1 do wdelete, !d.window

stop
END

; This routine is used to generate the FRD plots for many fibers. It's
; perhaps going obsolete as I am re-writing my code to produce these
; plots in another way.

PRO frplot_all

;This routine works with all the vp1 and vp2 data
restore,'vp1data.idl'
restore,'vp2data.idl'

on1='vp2_fibers.ps'
on2='Comparison of VP1 and VP2 Polymicro Fibers'

colour = [60,110,180,150]

xulim1 = max(fnumbervp1)+0.05
xllim1 = min(fnumbervp1)-0.05
xulim2 = max(fnumbervp2)+0.05
xllim2 = min(fnumbervp2)-0.05

lx1=[3.65,3.65]
lx2=[3.36,3.36]
ly1=[0.65,1.105]

ulim951 = min(datavp1[25,7,*])
llim951 = max(datavp1[25,7,*])
ulim952 = min(datavp2[25,7,*])
llim952 = max(datavp2[25,7,*])

ly2 = [0.88,1.02]
lx2a1 = [llim951, llim951]
lx2b1 = [ulim951, ulim951]
lx2a2 = [llim952, llim952]
lx2b2 = [ulim952, ulim952]

l11 = [xulim1-0.045,xulim1-0.17]
l12 = [xulim2-0.045,xulim2-0.17]
l2 = [1.03,1.03]

;==========================================================================
;BEGINNING OF THE PLOT
;==========================================================================

;this first part should be the same for all, except for possibily the
;limits set by using vp2 data over vp1 data

set_plot,'ps'
device,filename=on1,/color
;loadct,0
;plot,datavp2[*,7,0],datavp2[*,0,0],xthick=2,ythick=2,thick=2,charthick=2,$
;  xtitle='Output F/#',ytitle='Encircled Energy (EE)',title=on2,$
;  yrange=[0.65,1.05],ystyle=1,xrange=[xulim2,xllim2],xstyle=1,/NODATA
;loadct,4

;================================================================
;this section gets turned on for the fiber comparison plots only
;the plot part above gets turned off
;================================================================
!P.MULTI=[0,2,2,0,0]
loadct,0
plot,m337pm,mee,xthick=2,ythick=2,charthick=2,xtitle='Output F/#',$
  ytitle='Encircled Energy (EE)',title='3370A F/#!dout!n for 4 Fiber Types',yrange=[0.65,1.05],ystyle=1,$
  xrange=[xulim2,xllim2],xstyle=1,charsize=0.7,/NODATA
loadct,4
oplot,m337pm,mee,psym=-3,thick=3,symsize=0.7,color=colour[0]
oplot,m337co,mee,psym=-3,thick=3,symsize=0.7,color=colour[1]
oplot,m337ft16,mee,psym=-3,thick=3,symsize=0.7,color=colour[2]
oplot,m337ft22,mee,psym=-3,thick=3,symsize=0.7,color=colour[3]
xyouts,xllim2+0.65,0.795,'Polymicro  NA: 0.22',color=colour[0],charsize=0.8,charthick=3
xyouts,xllim2+0.65,0.765,'CeramOptec  NA: 0.22',color=colour[1],charsize=0.8,charthick=3
xyouts,xllim2+0.65,0.735,'FiberTech  NA: 0.16',color=colour[2],charsize=0.8,charthick=3
xyouts,xllim2+0.65,0.705,'Fibertech  NA: 0.22',color=colour[3],charsize=0.8,charthick=3

loadct,0
plot,m365pm,mee,xthick=2,ythick=2,charthick=2,xtitle='Output F/#',$
  ytitle='Encircled Energy (EE)',title='3650A F/#!dout!n for 4 Fiber Types',yrange=[0.65,1.05],ystyle=1,$
  xrange=[xulim2,xllim2],xstyle=1,charsize=0.7,/NODATA
loadct,4
oplot,m365pm,mee,psym=-3,thick=3,symsize=0.7,color=colour[0]
oplot,m365co,mee,psym=-3,thick=3,symsize=0.7,color=colour[1]
oplot,m365ft16,mee,psym=-3,thick=3,symsize=0.7,color=colour[2]
oplot,m365ft22,mee,psym=-3,thick=3,symsize=0.7,color=colour[3]
xyouts,xllim2+0.65,0.795,'Polymicro  NA: 0.22',color=colour[0],charsize=0.8,charthick=3
xyouts,xllim2+0.65,0.765,'CeramOptec  NA: 0.22',color=colour[1],charsize=0.8,charthick=3
xyouts,xllim2+0.65,0.735,'FiberTech  NA: 0.16',color=colour[2],charsize=0.8,charthick=3
xyouts,xllim2+0.65,0.705,'Fibertech  NA: 0.22',color=colour[3],charsize=0.8,charthick=3

loadct,0
plot,m400pm,mee,xthick=2,ythick=2,charthick=2,xtitle='Output F/#',$
  ytitle='Encircled Energy (EE)',title='4000A F/#!dout!n for 4 Fiber Types',yrange=[0.65,1.05],ystyle=1,$
  xrange=[xulim2,xllim2],xstyle=1,charsize=0.7,/NODATA
loadct,4
oplot,m400pm,mee,psym=-3,thick=3,symsize=0.7,color=colour[0]
oplot,m400co,mee,psym=-3,thick=3,symsize=0.7,color=colour[1]
oplot,m400ft16,mee,psym=-3,thick=3,symsize=0.7,color=colour[2]
oplot,m400ft22,mee,psym=-3,thick=3,symsize=0.7,color=colour[3]
xyouts,xllim2+0.65,0.795,'Polymicro  NA: 0.22',color=colour[0],charsize=0.8,charthick=3
xyouts,xllim2+0.65,0.765,'CeramOptec  NA: 0.22',color=colour[1],charsize=0.8,charthick=3
xyouts,xllim2+0.65,0.735,'FiberTech  NA: 0.16',color=colour[2],charsize=0.8,charthick=3
xyouts,xllim2+0.65,0.705,'Fibertech  NA: 0.22',color=colour[3],charsize=0.8,charthick=3

loadct,0
plot,m600pm,mee,xthick=2,ythick=2,charthick=2,xtitle='Output F/#',$
  ytitle='Encircled Energy (EE)',title='6000A F/#!dout!n for 4 Fiber Types',yrange=[0.65,1.05],ystyle=1,$
  xrange=[xulim2,xllim2],xstyle=1,charsize=0.7,/NODATA
loadct,4
oplot,m600pm,mee,psym=-3,thick=3,symsize=0.7,color=colour[0]
oplot,m600co,mee,psym=-3,thick=3,symsize=0.7,color=colour[1]
oplot,m600ft16,mee,psym=-3,thick=3,symsize=0.7,color=colour[2]
oplot,m600ft22,mee,psym=-3,thick=3,symsize=0.7,color=colour[3]
xyouts,xllim2+0.65,0.795,'Polymicro  NA: 0.22',color=colour[0],charsize=0.8,charthick=3
xyouts,xllim2+0.65,0.765,'CeramOptec  NA: 0.22',color=colour[1],charsize=0.8,charthick=3
xyouts,xllim2+0.65,0.735,'FiberTech  NA: 0.16',color=colour[2],charsize=0.8,charthick=3
xyouts,xllim2+0.65,0.705,'Fibertech  NA: 0.22',color=colour[3],charsize=0.8,charthick=3
!P.MULTI=[0]

;this next part just puts in the input line and 95% requirement line
;and won't need to change.

;oplot,lx1,ly1,color=60,thick=3
;xyouts,3.67,0.68,'INPUT BEAM F/#',orientation=90,color=60,charsize=1,charthick=3
;oplot,lx2,ly1,color=150,thick=3
;xyouts,3.38,0.68,'VIRUS 95% EE REQUIREMENT',orientation=90,color=150,charsize=1,charthick=3

;===========================================================================
;this chunk gets turned ON for all.ps and OFF for mean.ps
;oplot,lx2a,ly2,thick=6,color=110,linestyle=2
;oplot,lx2b,ly2,thick=6,color=110,linestyle=2
;oplot,l1,l2,thick=6,color=110,linestyle=2
;xyouts,xulim-0.20,1.025,': 95% EE Range',color=110,charsize=1.5,charthick=3
;===========================================================================

;Plots points...
;=======================================================================================
;FOR k=0,numvp1-1,4 DO BEGIN
;    oplot,datavp1[*,7,k+1],datavp1[*,0,k+1],psym=1,symsize=0.7,thick=2,color=colour[1]
;    oplot,datavp1[*,7,k],datavp1[*,0,k],psym=1,symsize=0.7,thick=2,color=colour[0]
;    oplot,datavp1[*,7,k+2],datavp1[*,0,k+2],psym=1,symsize=0.7,thick=2,color=colour[2]
;    oplot,datavp1[*,7,k+3],datavp1[*,0,k+3],psym=1,symsize=0.7,thick=2,color=colour[3]
;ENDFOR
;FOR k=0,numvp2-1,4 DO BEGIN
;    oplot,datavp2[*,7,k+1],datavp2[*,0,k+1],psym=4,symsize=0.7,thick=2,color=colour[1]
;    oplot,datavp2[*,7,k],datavp2[*,0,k],psym=4,symsize=0.7,thick=2,color=colour[0]
;    oplot,datavp2[*,7,k+2],datavp2[*,0,k+2],psym=4,symsize=0.7,thick=2,color=colour[2]
;    oplot,datavp2[*,7,k+3],datavp2[*,0,k+3],psym=4,symsize=0.7,thick=2,color=colour[3]
;ENDFOR
;=======================================================================================

;Plots the mean values off ALL THE FIBERS...
;=======================================================================================
;oplot,m337vp1,mee,psym=1,symsize=1.0,thick=2,color=colour[0]
;oplot,m365vp1,mee,psym=1,symsize=1.0,thick=2,color=colour[1]
;oplot,m400vp1,mee,psym=1,symsize=1.0,thick=2,color=colour[2]
;oplot,m600vp1,mee,psym=1,symsize=1.0,thick=2,color=colour[3]

;oplot,m337vp2,mee,psym=4,symsize=1.0,thick=2,color=colour[0]
;oplot,m365vp2,mee,psym=4,symsize=1.0,thick=2,color=colour[1]
;oplot,m400vp2,mee,psym=4,symsize=1.0,thick=2,color=colour[2]
;oplot,m600vp2,mee,psym=4,symsize=1.0,thick=2,color=colour[3]
;=======================================================================================



;versions of the LEGEND
;this gives the symbols with wavelengths...
;YOU'LL NEED TO CHANGE THE 'XULIM' TO EITHER XULIM1 OR XULIM2
;plots,xulim2-0.17,1.0,psym=4,symsize=1.8,thick=3,color=colour[0]
;xyouts,xulim2-0.20,0.995,': 3370A',color=colour[0],charsize=1.5,charthick=3
;plots,xulim2-0.17,0.97,psym=4,symsize=1.8,thick=3,color=colour[1]
;xyouts,xulim2-0.20,0.965,': 3650A',color=colour[1],charsize=1.5,charthick=3
;plots,xulim2-0.17,0.94,psym=4,symsize=1.8,thick=3,color=colour[2]
;xyouts,xulim2-0.20,0.935,': 4000A',color=colour[2],charsize=1.5,charthick=3
;plots,xulim2-0.17,0.91,psym=4,symsize=1.8,thick=3,color=colour[3]
;xyouts,xulim2-0.20,0.905,': 6000A',color=colour[3],charsize=1.5,charthick=3

;===================
;this one is used to show different symbols for different types of
;fibers
;loadct,0
;plots,xulim-0.12,1.0,psym=1,symsize=1.8,thick=3
;xyouts,xulim-0.15,0.995,': VP1 Polymicro Fibers',charsize=1.5,charthick=3
;plots,xulim-0.12,0.97,psym=4,symsize=1.8,thick=3
;xyouts,xulim-0.15,0.965,': VP2 Polymicro Fibers',charsize=1.5,charthick=3
;plots,xulim-0.12,0.94,psym=4,symsize=1.8,thick=3
;xyouts,xulim-0.15,0.94,'  (Colors as in other plots)',charsize=1.0,charthick=3
;plots,xulim-0.12,0.91,psym=4,symsize=1.8,thick=3
;xyouts,xulim-0.15,0.905,': 6000A',charsize=1.5,charthick=3

device,/close_file
set_plot,'x'

stop
END

PRO frplot_v3

;This routine is for plotting many different FRD measurements.  It
;works from a list of any length

;Unlike the first version, this code works from a list, CALLED 'LIST', which is all
;the 'name.frd' files FOR A GIVEN FIBER.

;The order of the list is important. It needs to be ordered first by
;fiber number, then by filter number:
;ex:
;fiber18_filter3370A.frd
;fiber18_filter3650A.frd
;fiber18_filter4000A.frd
;Fiber18_filter6000A.frd
;fiber19_filter3370A.frd
;fiber19_filter3650A.frd
;fiber19_filter4000A.frd
;Fiber19_filter6000A.frd

readcol,'listvp2',files2,FORMAT='A'
num2 = n_elements(files2)
readcol,'listvp1',files1,FORMAT='A'
num1 = n_elements(files1)

on1='vp1_vp2_pm.ps'
on2='Comparison of VP1 and VP2 Polymicro Fibers'

colour = [60,110,180,150]

data2 = dblarr(30,8,num2)
data1 = dblarr(30,8,num1)

temp = dblarr(8,30)
names2 = strarr(num2)
names1 = strarr(num1)

FOR k=0,num1-1 DO BEGIN
    name1 = strsplit(files1[k],'.',/EXTRACT)
    openr,5,files1[k]
    readf,5,temp
    free_lun,5
    trans = transpose(temp)
    data1[*,*,k] = trans
    names1[k] = name[0]
ENDFOR
FOR k=0,num2-1 DO BEGIN
    name2 = strsplit(files2[k],'.',/EXTRACT)
    openr,5,files2[k]
    readf,5,temp
    free_lun,5
    trans = transpose(temp)
    data2[*,*,k] = trans
    names2[k] = name[0]
ENDFOR

ee1= dblarr(30,num1)
ee2= dblarr(30,num2)
fnumber1 = ee1
fnumber2 = ee2

FOR j=0,num1-1 DO BEGIN
    ee1[*,j] = data1[*,0,j]
    temp = data1[*,7,j]
    fnumber1[*,j] = reverse(temp)
ENDFOR
FOR j=0,num2-1 DO BEGIN
    ee2[*,j] = data2[*,0,j]
    temp = data2[*,7,j]
    fnumber2[*,j] = reverse(temp)
ENDFOR

xulim1 = max(fnumber1)+0.1
xllim1 = min(fnumber1)-0.1
xulim2 = max(fnumber2)+0.1
xllim2 = min(fnumber2)-0.1

lx1=[3.65,3.65]
lx2=[3.36,3.36]
ly1=[0.65,1.105]

ulim951 = min(data1[25,7,*])
llim951 = max(data1[25,7,*])
ulim952 = min(data2[25,7,*])
llim952 = max(data2[25,7,*])

ly2 = [0.88,1.02]
lx2a1 = [llim951, llim951]
lx2b1 = [ulim951, ulim951]
lx2a2 = [llim952, llim952]
lx2b2 = [ulim952, ulim952]

l11 = [xulim1-0.045,xulim1-0.17]
l12 = [xulim2-0.045,xulim2-0.17]
l2 = [1.03,1.03]

;finds the mean of the list at each color
;this also sets up the means for each fiber type in VP2
;====================================================================
mee = dblarr(30)
m600pm = dblarr(30)
m600co = dblarr(30)
m600ft16 = dblarr(30)
m600ft22 = dblarr(30)
m600 = dblarr(30)
m400pm = dblarr(30)
m400co = dblarr(30)
m400ft16 = dblarr(30)
m400ft22 = dblarr(30)
m400 = dblarr(30)
m365pm = dblarr(30)
m365co = dblarr(30)
m365ft16 = dblarr(30)
m365ft22 = dblarr(30)
m365 = dblarr(30)
m337pm = dblarr(30)
m337co = dblarr(30)
m337ft16 = dblarr(30)
m337ft22 = dblarr(30)
m337 = dblarr(30)

slice = dblarr(30,num2/4)
qtr = num2/4

;====================================================================
;These individual fiber values get turned on only for VP2
;====================================================================

FOR k=0,qtr-1 DO slice[*,k] = data2[*,0,k*4]
FOR j=0,29 DO mee[j] = mean(slice[j,*])

FOR k=0,qtr-1 DO slice[*,k] = data2[*,7,k*4]
FOR j=0,29 DO m337[j] = mean(slice[j,*])
FOR j=0,29 DO m337pm[j] = mean(slice[j,0:9])
FOR j=0,29 DO m337co[j] = mean(slice[j,10:13])
FOR j=0,29 DO m337ft16[j] = mean(slice[j,14:18])
FOR j=0,29 DO m337ft22[j] = mean(slice[j,19:23])

FOR k=0,qtr-1 DO slice[*,k] = data2[*,7,k*4+1]
FOR j=0,29 DO m365[j] = mean(slice[j,*])
FOR j=0,29 DO m365pm[j] = mean(slice[j,0:9])
FOR j=0,29 DO m365co[j] = mean(slice[j,10:13])
FOR j=0,29 DO m365ft16[j] = mean(slice[j,14:18])
FOR j=0,29 DO m365ft22[j] = mean(slice[j,19:23])

FOR k=0,qtr-1 DO slice[*,k] = data2[*,7,k*4+2]
FOR j=0,29 DO m400[j] = mean(slice[j,*])
FOR j=0,29 DO m400pm[j] = mean(slice[j,0:9])
FOR j=0,29 DO m400co[j] = mean(slice[j,10:13])
FOR j=0,29 DO m400ft16[j] = mean(slice[j,14:18])
FOR j=0,29 DO m400ft22[j] = mean(slice[j,19:23])

FOR k=0,qtr-1 DO slice[*,k] = data2[*,7,k*4+3]
FOR j=0,29 DO m600[j] = mean(slice[j,*])
FOR j=0,29 DO m600pm[j] = mean(slice[j,0:9])
FOR j=0,29 DO m600co[j] = mean(slice[j,10:13])
FOR j=0,29 DO m600ft16[j] = mean(slice[j,14:18])
FOR j=0,29 DO m600ft22[j] = mean(slice[j,19:23])

;==========================================================================
;BEGINNING OF THE PLOT
;==========================================================================

set_plot,'ps'
device,filename=on1,/color
loadct,0
plot,data[*,7,0],data[*,0,0],xthick=2,ythick=2,thick=2,charthick=2,$
  xtitle='Output F/#',ytitle='Encircled Energy (EE)',title=on2,$
  yrange=[0.65,1.05],ystyle=1,xrange=[xulim,xllim],xstyle=1,/NODATA
loadct,4

oplot,lx1,ly1,color=60,thick=3
xyouts,3.67,0.68,'INPUT BEAM F/#',orientation=90,color=60,charsize=1,charthick=3
oplot,lx2,ly1,color=150,thick=3
xyouts,3.38,0.68,'VIRUS 95% EE REQUIREMENT',orientation=90,color=150,charsize=1,charthick=3

;===========================================================================
;this chunk gets turned ON for all.ps and OFF for mean.ps
;oplot,lx2a,ly2,thick=6,color=110,linestyle=2
;oplot,lx2b,ly2,thick=6,color=110,linestyle=2
;oplot,l1,l2,thick=6,color=110,linestyle=2
;xyouts,xulim-0.20,1.025,': 95% EE Range',color=110,charsize=1.5,charthick=3
;===========================================================================

;Plots points...
;FOR k=0,num-1,4 DO BEGIN
;    oplot,data[*,7,k+1],data[*,0,k+1],psym=4,symsize=0.7,thick=2,color=colour[1]
;    oplot,data[*,7,k],data[*,0,k],psym=4,symsize=0.7,thick=2,color=colour[0]
;    oplot,data[*,7,k+2],data[*,0,k+2],psym=4,symsize=0.7,thick=2,color=colour[2]
;    oplot,data[*,7,k+3],data[*,0,k+3],psym=4,symsize=0.7,thick=2,color=colour[3]
;ENDFOR

;Plots the mean values off ALL THE FIBERS...
;oplot,m337,mee,psym=1,symsize=1.0,thick=5,color=colour[0]
;oplot,m365,mee,psym=1,symsize=1.0,thick=5,color=colour[1]
;oplot,m400,mee,psym=1,symsize=1.0,thick=5,color=colour[2]
;oplot,m600,mee,psym=1,symsize=1.0,thick=5,color=colour[3]

;Plots the mean values FOR THE 4 DIFFERENT TYPES IN VP2
oplot,m337pm,mee,psym=1,thick=5,symsize=1.0,color=colour[0]
;oplot,m337co,mee,psym=1,thick=3,symsize=0.7,color=colour[1]
;oplot,m337ft16,mee,psym=2,thick=3,symsize=0.7,color=colour[2]
;oplot,m337ft22,mee,psym=7,thick=3,symsize=0.7,color=colour[3]

oplot,m365pm,mee,psym=1,thick=5,symsize=1.0,color=colour[1]
;oplot,m365co,mee,psym=1,thick=3,symsize=0.7,color=colour[1]
;oplot,m365ft16,mee,psym=2,thick=3,symsize=0.7,color=colour[2]
;oplot,m365ft22,mee,psym=7,thick=3,symsize=0.7,color=colour[3]

oplot,m400pm,mee,psym=1,thick=5,symsize=1.0,color=colour[2]
;oplot,m400co,mee,psym=1,thick=3,symsize=0.7,color=colour[1]
;oplot,m400ft16,mee,psym=2,thick=3,symsize=0.7,color=colour[2]
;oplot,m400ft22,mee,psym=7,thick=3,symsize=0.7,color=colour[3]

oplot,m600pm,mee,psym=1,thick=5,symsize=1.0,color=colour[3]
;oplot,m600co,mee,psym=1,thick=3,symsize=0.7,color=colour[1]
;oplot,m600ft16,mee,psym=2,thick=3,symsize=0.7,color=colour[2]
;oplot,m600ft22,mee,psym=7,thick=3,symsize=0.7,color=colour[3]


;THE LEGEND
plots,xulim-0.17,1.0,psym=4,symsize=1.8,thick=3,color=colour[0]
xyouts,xulim-0.20,0.995,': 3370A',color=colour[0],charsize=1.5,charthick=3
plots,xulim-0.17,0.97,psym=4,symsize=1.8,thick=3,color=colour[1]
xyouts,xulim-0.20,0.965,': 3650A',color=colour[1],charsize=1.5,charthick=3
plots,xulim-0.17,0.94,psym=4,symsize=1.8,thick=3,color=colour[2]
xyouts,xulim-0.20,0.935,': 4000A',color=colour[2],charsize=1.5,charthick=3
plots,xulim-0.17,0.91,psym=4,symsize=1.8,thick=3,color=colour[3]
xyouts,xulim-0.20,0.905,': 6000A',color=colour[3],charsize=1.5,charthick=3


;=========================================================================
;OBSOLETE
;=========================================================================

;Plots lines...
;FOR k=0,num-1,4 DO BEGIN
;    oplot,data[*,7,k+1],data[*,0,k+1],thick=2,color=colour[1]
;    oplot,data[*,7,k],data[*,0,k],thick=2,color=colour[0]
;    oplot,data[*,7,k+2],data[*,0,k+2],thick=2,color=colour[2]
;    oplot,data[*,7,k+3],data[*,0,k+3],thick=2,color=colour[3]
;    print,files[k]
;ENDFOR

;oplot,data[*,7,0],data[*,0,0],psym=4,symsize=0.5,thick=2,color=60
;oplot,data[*,7,1],data[*,0,1],psym=4,symsize=0.5,thick=2,color=110
;oplot,data[*,7,2],data[*,0,2],psym=4,symsize=0.5,thick=2,color=180
;oplot,data[*,7,3],data[*,0,3],psym=4,symsize=0.5,thick=2,color=150

;oplot,data[*,7,4],data[*,0,4],psym=4,symsize=0.5,thick=2,color=60
;oplot,data[*,7,5],data[*,0,5],psym=4,symsize=0.5,thick=2,color=110
;oplot,data[*,7,6],data[*,0,6],psym=4,symsize=0.5,thick=2,color=180
;oplot,data[*,7,7],data[*,0,7],psym=4,symsize=0.5,thick=2,color=150

;oplot,data[*,7,8],data[*,0,8],psym=7,symsize=0.5,thick=2,color=60
;oplot,data[*,7,9],data[*,0,9],psym=7,symsize=0.5,thick=2,color=110
;oplot,data[*,7,10],data[*,0,10],psym=7,symsize=0.5,thick=2,color=180
;oplot,data[*,7,11],data[*,0,11],psym=7,symsize=0.5,thick=2,color=150

;oplot,data[*,7,12],data[*,0,12],psym=7,symsize=0.5,thick=2,color=60
;oplot,data[*,7,13],data[*,0,13],psym=7,symsize=0.5,thick=2,color=110
;oplot,data[*,7,14],data[*,0,14],psym=7,symsize=0.5,thick=2,color=180
;oplot,data[*,7,15],data[*,0,15],psym=7,symsize=0.5,thick=2,color=150

;oplot,data[*,7,16],data[*,0,16],psym=7,symsize=0.5,thick=2,color=60
;oplot,data[*,7,17],data[*,0,17],psym=7,symsize=0.5,thick=2,color=110
;oplot,data[*,7,18],data[*,0,18],psym=7,symsize=0.5,thick=2,color=180
;oplot,data[*,7,19],data[*,0,19],psym=7,symsize=0.5,thick=2,color=150

;oplot,data[*,7,20],data[*,0,20],psym=7,symsize=0.5,thick=2,color=60
;oplot,data[*,7,21],data[*,0,21],psym=7,symsize=0.5,thick=2,color=110
;oplot,data[*,7,22],data[*,0,22],psym=7,symsize=0.5,thick=2,color=180
;oplot,data[*,7,23],data[*,0,23],psym=7,symsize=0.5,thick=2,color=150

;oplot,data[*,7,24],data[*,0,24],psym=7,symsize=0.5,thick=2,color=60
;oplot,data[*,7,25],data[*,0,25],psym=7,symsize=0.5,thick=2,color=110
;oplot,data[*,7,26],data[*,0,26],psym=7,symsize=0.5,thick=2,color=180
;oplot,data[*,7,27],data[*,0,27],psym=7,symsize=0.5,thick=2,color=150

;oplot,data[*,7,25],data[*,0,25],psym=7,symsize=0.5,thick=2,color=60
;oplot,data[*,7,26],data[*,0,26],psym=7,symsize=0.5,thick=2,color=110
;oplot,data[*,7,27],data[*,0,27],psym=7,symsize=0.5,thick=2,color=180
;oplot,data[*,7,28],data[*,0,28],psym=7,symsize=0.5,thick=2,color=150

;oplot,data[*,7,29],data[*,0,29],psym=7,symsize=0.5,thick=2,color=60
;oplot,data[*,7,30],data[*,0,30],psym=7,symsize=0.5,thick=2,color=110
;oplot,data[*,7,31],data[*,0,31],psym=7,symsize=0.5,thick=2,color=180
;oplot,data[*,7,32],data[*,0,32],psym=7,symsize=0.5,thick=2,color=150

;oplot,data[*,7,33],data[*,0,33],psym=7,symsize=0.5,thick=2,color=60
;oplot,data[*,7,34],data[*,0,34],psym=7,symsize=0.5,thick=2,color=110
;oplot,data[*,7,35],data[*,0,35],psym=7,symsize=0.5,thick=2,color=180
;oplot,data[*,7,36],data[*,0,36],psym=7,symsize=0.5,thick=2,color=150

;oplot,data[*,7,37],data[*,0,37],psym=7,symsize=0.5,thick=2,color=60
;oplot,data[*,7,38],data[*,0,38],psym=7,symsize=0.5,thick=2,color=110
;oplot,data[*,7,39],data[*,0,39],psym=7,symsize=0.5,thick=2,color=180
;oplot,data[*,7,47],data[*,0,40],psym=7,symsize=0.5,thick=2,color=150

device,/close_file
set_plot,'x'
stop
END

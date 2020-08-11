;UNDER CONSTRUCTION

;This routine is used to compare the best rms values for different
;fitting regions so as to determine the best level of sky
;subtraction. It runs AFTER plotrmsky.pro, which generates the output
;name_regionrms.min

;The list is of the form:
;framename_region_rms.min
;EX:   jm0043_Gband_rms.min
;      jm0043_H+K_rms.min
;      jm0043_Hbeta_rms.min
;      jm0043_Mg_rms.min     etc.

;THIS CODE REQUIRES THAT ALL OF THE RMS.MIN FILES ARE ORDERED SUCH
;THAT ALL THE FILE NAMES ARE IN ORDER, AS ABOVE.

;This requires 4 regions. Ff longer or shorter, change the region
;denoted by the ****

;---------------------------------------------------------------------
PRO plotrms_v2, list, step
;---------------------------------------------------------------------

yup = 0.15  ;This is the max RMS that will make it onto the plot

readcol,list,format='a',rmslist
n0 = n_elements(rmslist)
names = strarr(2,n0)
xaxis = findgen(n0)+1
yaxis = intarr(n0)
xup = max(xaxis)+1
xdown = min(xaxis)-1

for j=0,n0-1 do begin
   temp = rmslist[j]
   temp = strsplit(temp,'_',/extract)
   names[*,j] = [temp[0],temp[1]]
endfor

window,0,retain=2
device,decomposed=0

;*******************************************************
for j=0,n0-1,step do begin
   readcol,rmslist[j],format='a,f',f1,rms1
   readcol,rmslist[j+1],format='a,f',f2,rms2
   readcol,rmslist[j+2],format='a,f',f3,rms3
   readcol,rmslist[j+3],format='a,f',f4,rms4
   n1 = n_elements(f1)
   f11 = strarr(n1)
   n2 = n_elements(f2)
   f22 = strarr(n2)
   n3 = n_elements(f3)
   f33 = strarr(n3)
   n4 = n_elements(f4)
   f44 = strarr(n4)
   
   yup = max([rms1,rms2,rms3,rms4])
   ydown = min([rms1,rms2,rms3,rms4])

   for k=0,n1-1 do begin
       temp = strsplit(f1[k],'T',/extract)
       f11[k] = temp[0]
       f11[k] = strsplit(f11[k],'1234567890',/extract)
       f11 = f11[1]

   loadct,0
   plot,xaxis,yaxis,psym=3,charsize=1.2,xtitle='File Name',ytitle='RMS',$
        title='COMPARISON OF RMS VALUES',xrange = [xdown,xup],$
        yrange=[0,yup],xstyle=1,ystyle=1,/nodata
   loadct,4
;   plots,xaxis[j],rms1[0],psym=1,charsize=1.5,color=110 
   xyouts,xaxis[j],rms1[0],f11,charthick=1.5,color=110
;   plots,xaxis[j],rms2[0],psym=2,charsize=1.5,color=60
   xyouts,xaxis[j],rms2[0],f22,charthick=1.5,color=60
;   plots,xaxis[j],rms3[0],psym=4,charsize=1.5,color=180
   xyouts,xaxis[j],rms3[0],f33,charthick=1.5,color=180
;   plots,xaxis[j],rms4[0],psym=5,charsize=1.5,color=150
   xyouts,xaxis[j],rms4[0],f33,charthick=1.5,color=150

   for k=1,n1-1 do begin
      xyouts,xaxis[j],rms1[k],charsize=0.8,color=110
      xyouts,xaxis[j],rms2[k],charsize=0.8,color=60
      xyouts,xaxis[j],rms3[k],charsize=0.8,color=180
      xyouts,xaxis[j],rms4[k],charsize=0.8,color=150
   endfor
endfor


stop
end

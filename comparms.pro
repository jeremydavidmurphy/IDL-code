;This routine is used to compare the best rms values for different
;fitting regions so as to determine the best level of sky
;subtraction. It runs AFTER plotrmsky.pro, which generates the output
;name_region_rms.min AND AFTER you've PASTED all the four regions
;together into one file. The output of this gets called filename.rms
;EX: jm0077.rms. This is generally 8 columns of the 4 files below
;file_H+K  RMS   file_Gband   RMS  etc.

;list example
;jm0077.rms
;jm0079.rms
;jm0081.rms
;etc.
;with each of these lists in the form of
;jm0077xxTTf.fits  0.023481922  jm0077yxTTF.fits  0.023417182
;etc, where each column is sorted from lowest to highest in terms of
;RMS and the eight columns ARE ORDERED H+K, GBAND, HBETA AND MG

;THE NUMBER OF FILES INCLUDED IN THE FILE.RMS IS THE NUMBER FOR "STEP"
;---------------------------------------------------------------------
PRO comparms, list
;---------------------------------------------------------------------

step = 4
yup = 0.25  ;This is the max RMS that will make it onto the plot
h = 0.20  ;the 'h'eight of the file names
s = 0.01   ;the 's'tep between each

readcol,list,format='a',rmslist
n0 = n_elements(rmslist)
files = strarr(n0)
outfiles = strarr(4,n0)

xaxis = findgen(n0)+1
yaxis = intarr(n0)
xup = max(xaxis)+2
xdown = min(xaxis)-4

for j=0,n0-1 do begin
   temp = rmslist[j]
   temp = strsplit(temp,'_.',/extract)
   files[j] = temp[0] 
endfor

window,0,retain=2,xsize=1200,ysize=900
device,decomposed=0

;*******************************************************
for j=0,n0-1 do begin
    readcol,rmslist[j],format='a,d,a,d,a,d,a,d',f1,rms1,f2,rms2,f3,rms3,f4,rms4
    n1 = n_elements(f1)
    nm1 = strarr(n1)
    n2 = n_elements(f2)
    nm2 = strarr(n2)
    n3 = n_elements(f3)
    nm3 = strarr(n3)
    n4 = n_elements(f4)
    nm4 = strarr(n4)

    outfiles[0,j] = f1[0]
    outfiles[1,j] = f2[0]
    outfiles[2,j] = f3[0]
    outfiles[3,j] = f4[0]

    for k=0,n1-1 do begin
        temp = strsplit(f1[k],'T',/extract)
        temp = temp[0]
        temp = strsplit(temp,'1234567890',/extract)
        nm1[k] = temp[1]
    endfor
    for k=0,n2-1 do begin
        temp = strsplit(f2[k],'T',/extract)
        temp = temp[0]
        temp = strsplit(temp,'1234567890',/extract)
        nm2[k] = temp[1]
    endfor
    for k=0,n3-1 do begin
        temp = strsplit(f3[k],'T',/extract)
        temp = temp[0]
        temp = strsplit(temp,'1234567890',/extract)
        nm3[k] = temp[1]
    endfor
    for k=0,n4-1 do begin
        temp = strsplit(f4[k],'T',/extract)
        temp = temp[0]
        temp = strsplit(temp,'1234567890',/extract)
        nm4[k] = temp[1]
    endfor
    
    if (j eq 0) then begin
        loadct,0
        plot,xaxis,yaxis,charsize=1.2,xtitle='Files',ytitle='RMS',$
          title='COMPARISON OF RMS VALUES',xrange=[xdown,xup],$
          yrange=[0,yup],xstyle=1,ystyle=1,/nodata
    endif

    loadct,4
    xyouts,xaxis[j],h,nm1[0],charthick=1.5,charsize=1.5,color=60
    xyouts,-2,h,'H+K',charthick=1.5,charsize=1.5,color=60
    xyouts,xaxis[j],h+s,nm2[0],charthick=1.5,charsize=1.5,color=110
    xyouts,-2,h+s,'Gband',charthick=1.5,charsize=1.5,color=110
    xyouts,xaxis[j],h+s+s,nm3[0],charthick=1.5,charsize=1.5,color=180
    xyouts,-2,h+s+s,'Hbeta',charthick=1.5,charsize=1.5,color=180
    xyouts,xaxis[j],h+s+s+s,nm3[0],charthick=1.5,charsize=1.5,color=150
    xyouts,-2,h+s+s+s,'Mg',charthick=1.5,charsize=1.5,color=150
    xyouts,xaxis[j]+0.5,h-s-s-s,files[j],charthick=1.5,charsize=1.5,color=255,orientation=90
    
    for k=0,n1-1 do begin
        plots,xaxis[j],rms1[k],psym=4,symsize=1.0,color=60
        plots,xaxis[j],rms2[k],psym=4,symsize=1.0,color=110
        plots,xaxis[j],rms3[k],psym=4,symsize=1.0,color=180
        plots,xaxis[j],rms4[k],psym=4,symsize=1.0,color=150
    endfor

endfor

print,'Next ENTER deletes the plot...'
pause
wdelete,0

ans = ''
print,'Save the plot and lowest RMS values? (y or n)'
read,ans
if (ans eq 'y') then begin
    name=''
    print,'Enter the name of the files (w/o .ps):'
    read,name
    nameps = name+'.ps'
    nametxt = name+'.txt'

    openw,5,nametxt
    printf,5,outfiles
    free_lun,5

    set_plot,'ps'
    device,file=nameps,/color
    loadct,0
    
    for j=0,n0-1 do begin
        readcol,rmslist[j],format='a,d,a,d,a,d,a,d',f1,rms1,f2,rms2,f3,rms3,f4,rms4
        n1 = n_elements(f1)
        nm1 = strarr(n1)
        n2 = n_elements(f2)
        nm2 = strarr(n2)
        n3 = n_elements(f3)
        nm3 = strarr(n3)
        n4 = n_elements(f4)
        nm4 = strarr(n4)
        
        for k=0,n1-1 do begin
            temp = strsplit(f1[k],'T',/extract)
            temp = temp[0]
            temp = strsplit(temp,'1234567890',/extract)
            nm1[k] = temp[1]
        endfor
        for k=0,n2-1 do begin
            temp = strsplit(f2[k],'T',/extract)
            temp = temp[0]
            temp = strsplit(temp,'1234567890',/extract)
            nm2[k] = temp[1]
        endfor
        for k=0,n3-1 do begin
            temp = strsplit(f3[k],'T',/extract)
            temp = temp[0]
            temp = strsplit(temp,'1234567890',/extract)
            nm3[k] = temp[1]
        endfor
        for k=0,n4-1 do begin
            temp = strsplit(f4[k],'T',/extract)
            temp = temp[0]
            temp = strsplit(temp,'1234567890',/extract)
            nm4[k] = temp[1]
        endfor
        
        if (j eq 0) then begin
            loadct,0
            plot,xaxis,yaxis,charsize=1.2,xtitle='Files',ytitle='RMS',$
              title='COMPARISON OF RMS VALUES',xrange=[xdown,xup],$
              yrange=[0,yup],xstyle=1,ystyle=1,/nodata,xthick=2,ythick=2,charthick=2
        endif
        
        loadct,4
        xyouts,xaxis[j],h,nm1[0],charthick=1.5,charsize=0.7,color=60,orientation=90
        xyouts,-2.5,h,'H+K',charthick=1.5,charsize=0.7,color=60
        xyouts,xaxis[j],h+s,nm2[0],charthick=1.5,charsize=0.7,color=110,orientation=90
        xyouts,-2.5,h+s,'Gband',charthick=1.5,charsize=0.7,color=110
        xyouts,xaxis[j],h+s+s,nm3[0],charthick=1.5,charsize=0.7,color=180,orientation=90
        xyouts,-2.5,h+s+s,'Hbeta',charthick=1.5,charsize=0.7,color=180
        xyouts,xaxis[j],h+s+s+s,nm3[0],charthick=1.5,charsize=0.7,color=150,orientation=90
        xyouts,-2.5,h+s+s+s,'Mg',charthick=1.5,charsize=0.7,color=150
        loadct,0
        xyouts,xaxis[j],h-s-s-s,files[j],charthick=1.5,charsize=0.7,orientation=90
        loadct,4
        for k=0,n1-1 do begin
            plots,xaxis[j],rms1[k],psym=4,symsize=0.8,color=60
            plots,xaxis[j],rms2[k],psym=4,symsize=0.8,color=110
            plots,xaxis[j],rms3[k],psym=4,symsize=0.8,color=180
            plots,xaxis[j],rms4[k],psym=4,symsize=0.8,color=150
        endfor
    endfor
    device,/close
    set_plot,'x'
endif

stop
end

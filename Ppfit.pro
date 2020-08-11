;This code expects (and doesn't really work without) the pfitlov.out
;file in the form of containing discreet chunks, broken my wavelength
;regions. So, the lists fed into pfitlov needs to have the different
;wavelength regions separated by chunks and not mixed together.
;IT SHOULD BE ABLE TO HANDLE ANY NUMBER OF CHUNKS... 

;EX:
;bin1234_H+K.fit
;bin1235_H+K.fit
;bin1236_H+K.fit
;bin1234_Gband.fit
;bin1235_Gband.fit
;bin1236_Gband.fit
;bin1234_Hbeta.fit
;bin1235_Hbeta.fit
;bin1236_Hbeta.fit
;bin1234_Mg.fit
;bin1235_Mg.fit
;bin1236_Mg.fit

;It's also nice to have them in order from blue to red and there is
;color-coded that occurs, matching the first spectra with a blue color
;and red to red.

;The names printed on the plots will split the pfitlov.out names on
;THE FIRST '_'. As the files should be of the form
;filename_region.fits
;The filename gets included in the plots AND THE DATA IS SPLIT ON THE
;REGION. THERE MUST BE AN EVEN NUMBER OF FILES IN EACH REGION! If this
;criteria is not met, you will be misplotting by mixing up the
;wavelength regions.

;-------------------------------------------------------------------
pro Ppfit,  NAME=pfitname
;-------------------------------------------------------------------

;to analysize pfitlov.out files. it looks for the pfitlov.out file in
;the calling directory. It generates plots of velocity, vel. disp. h3
;and h4, both to the screen and as .ps files.

;This code requires that the fitting regions ARE ALL TOGETHER IN THE
;PFITLOV LIST.

ans = ''
ans1=''
note = ''

color = [60,110,180,150,220,255,80,130,170,200,240,60,110,180,150,220,255,80,130,170,200,240,60,110,180,150,220,255,80,130,170,200,240,60,110,180,150,220,255,80,130,170,200,240]
;symbol = [-1,-2,-4,-5,-6,-7,-1,-2,-4,-5,-6,-7,-1,-2,-4,-5,-6,-7,-1,-2,-4,-5,-6,-7,-1,-2,-4,-5,-6,-7,-1,-2,-4,-5,-6,-7]
symbol = [1,2,4,5,6,7,1,2,4,5,6,7,1,2,4,5,6,7,1,2,4,5,6,7,1,2,4,5,6,7,1,2,4,5,6,7,1,2,4,5,6,7,1,2,4,5,6,7,1,2,4,5,6,7,1,2,4,5,6,7,1,2,4,5,6,7,1,2,4,5,6,7,1,2,4,5,6,7,1,2,4,5,6,7,1,2,4,5,6,7]
;color = [60,60,110,110,180,180,150,150,220,220]
;symbol = [1,1,2,2,4,4,5,5,6,6,7,7]
;symbol = [-1,-1,-1,-1,-2,-2,-2,-2,-4,-4,-4,-4,-5,-5,-5,-5]
;symbol = [1,1,1,2,2,2,4,4,4,5,5,5,6,6,6,7,7,7]
;color = [60,110,180,150,60,110,180,150,60,110,180,150]
;symbol = [1,2,4,5,1,2,4,5,1,2,4,5]
;symbol = [-1,-2,-4,-5,-1,-2,-4,-5,-1,-2,-4,-5]
;color = [60,150,60,150,60,150,60,150]
;symbol = [1,1,1,1,2,2,2,2,4,4,4,4,4,4,4,5,5,5,5,5,5,5]
;color = [110,110,110,110,110,110,110,110,110,110,110,110,110,110,150,150,150,150,150,150,150,150,150,150,150,150,150,150]
;color = [60,60,60,60,60,60,60,110,110,110,110,110,110,110,150,150,150,150,150,150,150,180,180,180,180,180,180,180]
;color = [110,130,110,130,110,130,110,130,110,130,110,130,110,130,150,150,180,150,180,180,150,180,150,180,150,180,150,180]
;color = [110,110,110,110,110,110,110,150,150,150,150,150,150,150]

t = n_elements(pfitname)
if (t eq 0) then begin
    readcol,'pfitlov.out',wholename,vel,sigma,h3,h4,format='A,F,F,X,F,F,X,X,X'
    pfitname='pfitlov.out'
endif else begin
    readcol,pfitname,wholename,vel,sigma,h3,h4,format='A,F,F,X,F,F,X,X,X'
endelse

n1 = n_elements(wholename)
t = wholename[0]
t = strsplit(t,'_.',/extract)
n2 = n_elements(t)
names = strarr(n2,n1)
plotnames = strarr(n1)
for j=0,n1-1 do begin
    t = wholename[j]
    t1 = strsplit(t,'_.',/extract)
    plotnames[j] = t1[0]
    names[*,j] = t1
endfor

print,names[*,0]
element = intarr(1)
print,'Which array element to differentiate on?'
print,'Type 0 if you want to manually enter the length of a region.'
read,element
if (element eq 0) then begin
    print,'Enter the LENGTH of the region you want to split the file on'
    read,element
    n3 = element
    n3 = n3[0]
    n3 = fix(n3)
    n4 = n1/n3
    n4 = n4[0]
    n4 = fix(n4)
    print,'Which array element to split on for the naming?'
    read,element
    element = element[0]-1
endif else begin
    element = element[0]-1
    test = names[element,0]
    index = where(names[element,*] eq test)
    n3 = n_elements(index) ;the length of a fitting region.
    n4 = n1/n3
endelse

namearr = strarr(2,n3,n4)
velarr = dblarr(n3,n4)
sigmaarr = dblarr(n3,n4)
h3arr = dblarr(n3,n4)
h4arr = dblarr(n3,n4)

i0 = 0
i1 = n3-1
for j=0,n4-1 do begin
    namearr[0,*,j] = plotnames[i0:i1]
    namearr[1,0,j] = names[element,i0]
    velarr[*,j] = vel[i0:i1]
    sigmaarr[*,j] = sigma[i0:i1]
    h3arr[*,j] = h3[i0:i1]
    h4arr[*,j] = h4[i0:i1]
    i0 = i0+n3
    i1 = i1+n3
endfor

;the median for each set is calculated
medvelarr = dblarr(n3)
medsigarr = dblarr(n3)
medh3arr = dblarr(n3)
medh4arr = dblarr(n3)

for j=0,n3-1 do begin
    medvelarr[j] = median(velarr[j,*],/even)
    medsigarr[j] = median(sigmaarr[j,*],/even)
    medh3arr[j] = median(h3arr[j,*],/even)
    medh4arr[j] = median(h4arr[j,*],/even)
endfor
medvelarr = round(medvelarr)
medsigarr = round(medsigarr)
ttt = ['velocity', 'dispersion']

openw,5,pfitname+'_medvalues.txt'
;printf,ttt
for j=0,n3-1 do printf,5,medvelarr[j],medsigarr[j]
free_lun,5

xaxis = findgen(n3)+1.1
xup = max(xaxis)+3
xdown = min(xaxis)-3
;*********************************************************************************************
;                                    VELOCITY
;*********************************************************************************************
shift = 0.035
window,0,retain=2,xpos=50,ypos=500,xsize=1200,ysize=500
device,decomposed=0

ymed = median(velarr)
;yup = ymed+200
;ydown = ymed-200
yup=100
ydown=-200

loadct,0
plot,xaxis,velarr[*,0],ytitle='VELOCITY (km/sec)',$
  yrange=[ydown,yup],xrange=[xdown,xup],psym=3,/ynozero,/nodata,$
  ystyle=1,xstyle=1,title='Velocity data for '+pfitname

loadct,4
for j=0,n4-1 do begin
    oplot,xaxis,velarr[*,j],psym=symbol[j],color=color[j]
    xyouts,0.06,0.9-(j*shift),namearr[1,0,j],color=color[j],/normal,$
      charsize=1.2,charthick=1.5
endfor
for k=0,n3-1 do begin
    xyouts,xaxis[k],ydown+5,namearr[0,k,0],orientation=90,charsize=1.2
    xyouts,xaxis[k],yup-50,strn(medvelarr[k]),orientation=90,charsize=1.2
endfor

;print,'Save a hardcopy?'
;read,ans
ans='y'
if(ans eq 'y') then begin
    print,'The plot will be saved as "velocity.ps".'
    print,'Add a note to all the plots?'
    read,ans1
    if (ans1 eq 'y') then begin
        note = ''
        print,'Type a note...'
        read,note
    endif
    
    set_plot,'ps'
    device,filename='velocity.ps',/color,/portrait
    loadct,0
    plot,xaxis,velarr[*,0],ytitle='VELOCITY (km/sec)',$
      yrange=[ydown,yup],psym=3,/ynozero,/nodata,xthick=2,ythick=2,charthick=2,$
      ystyle=1,title='Velocity data for '+pfitname,xrange=[xdown,xup],xstyle=1
    loadct,4
    for j=0,n4-1 do begin
        oplot,xaxis,velarr[*,j],psym=symbol[j],color=color[j],symsize=0.8,thick=1.5
        xyouts,0.02,0.9-(j*0.03),namearr[1,0,j],color=color[j],/normal,charsize=1.0,charthick=2.0
    endfor
    for k=0,n3-1 do begin
        xyouts,xaxis[k],ydown+10,namearr[0,k,0],orientation=90,charsize=0.8,charthick=2.0
        xyouts,xaxis[k],yup-50,strn(medvelarr[k]),orientation=90,charsize=0.8,charthick=2.0
    endfor
    if (ans1 eq 'y') then begin
        xyouts,0.02,0.02,note,/normal,charthick=2.0,charsize=0.8
    endif
    device,/close_file
    set_plot,'x'
    
endif else print,'No velocity plot will be generates.'
print,'The next ENTER moves to the dispersion plot.'
pause

;*********************************************************************************************
;                                    DISPERSION
;*********************************************************************************************
shift = 0.035
window,1,retain=2,xpos=50,ypos=450,xsize=1200,ysize=500
device,decomposed=0

;ymed = median(sigmaarr)
;yup = max(sigmaarr)+30
;ydown = min(sigmaarr)-10
;if (yup gt 600) then yup = ymed+200
ydown=150
yup=500

loadct,0
plot,xaxis,sigmaarr[*,0],ytitle='VELOCITY DISPERSION (km/sec)',$
  yrange=[ydown,yup],psym=3,/ynozero,ystyle=1,/nodata,xrange=[xdown,xup],$
  xstyle=1,title='Velocity Dispersion data for '+pfitname

loadct,4
for j=0,n4-1 do begin
    oplot,xaxis,sigmaarr[*,j],psym=symbol[j],color=color[j]
    xyouts,0.06,0.9-(j*shift),namearr[1,0,j],color=color[j],/normal,$
      charsize=1.2,charthick=1.5
endfor
for k=0,n3-1 do begin
    xyouts,xaxis[k],ydown+5,namearr[0,k,0],orientation=90,charsize=1.2
    xyouts,xaxis[k],yup-50,strn(medsigarr[k]),orientation=90,charsize=1.2
endfor

;print,'Save a hardcopy?'
;read,ans
ans='y'
if(ans eq 'y') then begin
    print,'The plot will be saved as "dispersion.ps".'
;    print,'Add a note to the plot?'
;    read,ans1
;    if (ans1 eq 'y') then begin
;        note = ''
;        print,'Type a note...'
;        read,note
;    endif
    
    set_plot,'ps'
    device,filename='dispersion.ps',/color,/portrait
    plot,xaxis,sigmaarr[*,0],ytitle='VELOCITY DISPERSION (km/sec)',$
      yrange=[ydown,yup],psym=3,xthick=2,ythick=2,charthick=2,/ynozero,$
      ystyle=1,/nodata,title='Velocity Dispersion data for '+pfitname,$
      xrange=[xdown,xup],xstyle=1
    
    loadct,4
    for j=0,n4-1 do begin
        oplot,xaxis,sigmaarr[*,j],psym=symbol[j],color=color[j],symsize=0.8,thick=1.5
        xyouts,0.02,0.9-(j*0.03),namearr[1,0,j],color=color[j],/normal,charsize=1.0,charthick=2.0
    endfor
    for k=0,n3-1 do begin
        xyouts,xaxis[k],ydown+10,namearr[0,k,0],orientation=90,charsize=0.8,charthick=2.0
        xyouts,xaxis[k],yup-50,strn(medsigarr[k]),orientation=90,charsize=0.8,charthick=2.0
    endfor
    if (ans1 eq 'y') then begin
        xyouts,0.02,0.02,note,/normal,charthick=2.0,charsize=0.8
    endif
    device,/close_file
    set_plot,'x'
    
endif else print,'No velocity dispersion plot will be generated.'
print,'The next ENTER moves to the h3 plot.'
pause

;*********************************************************************************************
;                                    H3
;*********************************************************************************************
shift = 0.035
window,2,retain=2,xpos=50,ypos=400,xsize=1200,ysize=500
device,decomposed=0

;yup = max(h3arr)
;ydown = min(h3arr)
;if (yup gt 0.5) then yup = 0.5
;if (ydown lt -0.5) then ydown = -0.5
yup = 0.1
ydown = -0.2

loadct,0
plot,xaxis,h3arr[*,0],ytitle='H_3',yrange=[ydown,yup],psym=3,$
  /ynozero,/nodata,title='H_3 for '+pfitname,ystyle=1,xrange=[xdown,xup],xstyle=1

loadct,4
for j=0,n4-1 do begin
    oplot,xaxis,h3arr[*,j],psym=symbol[j],color=color[j]
    xyouts,0.06,0.9-(j*shift),namearr[1,0,j],color=color[j],/normal,$
      charsize=1.2,charthick=1.5
endfor
for k=0,n3-1 do begin
    xyouts,xaxis[k],ydown+0.02,namearr[0,k,0],orientation=90,charsize=1.2
;    xyouts,xaxis[k],yup-0.1,strn(medh3arr[k]),orientation=90,charsize=1.2
endfor

;print,'Save a hardcopy?'
;read,ans
ans='y'
if(ans eq 'y') then begin
;    print,'Add a note to the plot?'
;    read,ans1
;    if (ans1 eq 'y') then begin
;        note = ''
;        print,'Type a note...'
;        read,note
;    endif
    
    set_plot,'ps'
    device,filename='h3.ps',/color,/portrait
    loadct,0
    plot,xaxis,h3arr[*,0],ytitle='H_3',yrange=[ydown,yup],psym=3,$
      /ynozero,/nodata,title='H_3 for '+pfitname,ystyle=1,xrange=[xdown,xup],xstyle=1,$
      xthick=2,ythick=2,charthick=2
    loadct,4
    for j=0,n4-1 do begin
        oplot,xaxis,h3arr[*,j],psym=symbol[j],color=color[j],symsize=0.8,thick=1.5
        xyouts,0.02,0.9-(j*shift),namearr[1,0,j],color=color[j],/normal,$
          charsize=1.0,charthick=2.0
    endfor
    for k=0,n3-1 do begin
        xyouts,xaxis[k],ydown+0.02,namearr[0,k,0],orientation=90,charsize=0.8,charthick=2.0
;        xyouts,xaxis[k],yup-0.1,strn(medh3arr[k]),orientation=90,charsize=0.8,charthick=2.0
    endfor
    if (ans1 eq 'y') then begin
        xyouts,0.02,0.02,note,/normal,charthick=2.0,charsize=0.8
    endif
    device,/close_file
    set_plot,'x'

    print,'The plot was saved as "h3.ps".'
endif else print,'No h3 plot will be generated.'
print,'The next ENTER moves to the h4 plot.'
pause

;*********************************************************************************************
;                                    H4
;*********************************************************************************************
shift = 0.035
window,3,retain=2,xpos=50,ypos=350,xsize=1200,ysize=500
device,decomposed=0

;yup = max(h4arr)
;ydown = min(h4arr)
;if (yup gt 1) then yup = 1
;if (ydown lt -0.2) then ydown = -0.2
yup = 0.1
down = -0.2

loadct,0
plot,xaxis,h4arr[*,0],ytitle='H_4',yrange=[ydown,yup],psym=3,$
  /ynozero,/nodata,title='H_4 for '+pfitname,ystyle=1,xrange=[xdown,xup],xstyle=1

loadct,4
for j=0,n4-1 do begin
    oplot,xaxis,h4arr[*,j],psym=symbol[j],color=color[j]
    xyouts,0.06,0.9-(j*shift),namearr[1,0,j],color=color[j],/normal,$
      charsize=1.2,charthick=1.5
endfor
for k=0,n3-1 do begin
    xyouts,xaxis[k],ydown+0.02,namearr[0,k,0],orientation=90,charsize=1.2
;    xyouts,xaxis[k],yup-0.1,strn(medh4arr[k]),orientation=90,charsize=1.2
endfor

;print,'Save a hardcopy?'
;read,ans
ans='y'
if(ans eq 'y') then begin
;    print,'Add a note to the plot?'
;    read,ans1
;    if (ans1 eq 'y') then begin
;        note = ''
;        print,'Type a note...'
;        read,note
;    endif
    
    set_plot,'ps'
    device,filename='h4.ps',/color,/portrait
    loadct,0
    plot,xaxis,h4arr[*,0],ytitle='H_4',yrange=[ydown,yup],psym=3,$
      /ynozero,/nodata,title='H_4 for '+pfitname,ystyle=1,xrange=[xdown,xup],xstyle=1,$
      xthick=2,ythick=2,charthick=2
    loadct,4
    for j=0,n4-1 do begin
        oplot,xaxis,h4arr[*,j],psym=symbol[j],color=color[j],symsize=0.8,thick=1.5
        xyouts,0.02,0.9-(j*shift),namearr[1,0,j],color=color[j],/normal,$
          charsize=1.0,charthick=2.0
    endfor
    for k=0,n3-1 do begin
        xyouts,xaxis[k],ydown+0.02,namearr[0,k,0],orientation=90,charsize=0.8,charthick=2.0
;        xyouts,xaxis[k],yup-0.1,strn(medh4arr[k]),orientation=90,charsize=0.8,charthick=2.0
    endfor
    if (ans1 eq 'y') then begin
        xyouts,0.02,0.02,note,/normal,charthick=2.0,charsize=0.8
    endif
    device,/close_file
    set_plot,'x'

    print,'The plot was saved as "h4.ps".'
endif else print,'No h4 plot will be generated.'

print,'Next ENTER deletes the plots.'
pause
wdelete,0,1,2,3
stop
end


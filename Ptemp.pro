PRO Ptemp

;*************************************************************************
;plotting parameters
y1 = [-50,70];velocity
y2 = [270,370];disp
y3 = [-0.09,0.09];h3
y4 = [-0.09,0.07];h4

xup = 300
xdn = -300

;galname = 'M87'
;pmjr = 'no';change to anything other than 'no' if you want to plot the major axis only. (this supresses all off-major-axis values)
pw = 0.8 ;plot width. set this to a smaller number and your figure gets narrower. this is in normalized coordinates.
;error = 1 ;standard deviation of the 4 spectral region values
;error = 2 ;the mean of the error of the pallmc output.

;readcol,'bin.radiusmajor',format='a,f',bin2,rad
readcol,'bin.radius',format='a,f',bin2,rad
;outta = 'outta'
;outta = 'inna'
;*************************************************************************

readcol,'pfit.list',f='a',list
n0 = n_elements(list)

window,2,retain=2,xsize=800,ysize=350
device,decomposed=0
;window,6,retain=2,xsize=800,ysize=350
;device,decomposed=0
;window,10,retain=2,xsize=800,ysize=350
;device,decomposed=0
;window,14,retain=2,xsize=800,ysize=350
loadct,0

colors = [60,110,180,150]

print,list

for j=0,n0-1 do begin
    file = list[j]
    readcol,file,format='a,f,f,x,f,f,x,x,x',bin1,vel,disp,h3,h4
    n1 = n_elements(bin1)
    
    if (j eq 0) then begin
        plot,rad,vel,psym=2,title='Velocity'
        loadct,4
    endif else oplot,rad,vel,psym=2,color=colors[j-1]
    pause
endfor

for j=0,n0-1 do begin
    file = list[j]
    readcol,file,format='a,f,f,x,f,f,x,x,x',bin1,vel,disp,h3,h4
    n1 = n_elements(bin1)
    
    if (j eq 0) then begin
        plot,rad,disp,psym=2,title='Velocity',/ynozero
        loadct,4
    endif else oplot,rad,disp,psym=2,color=colors[j-1]
    pause
endfor

for j=0,n0-1 do begin
    file = list[j]
    readcol,file,format='a,f,f,x,f,f,x,x,x',bin1,vel,disp,h3,h4
    n1 = n_elements(bin1)
    
    if (j eq 0) then begin
        plot,rad,h3,psym=2,title='Velocity'
        loadct,4
    endif else oplot,rad,h3,psym=2,color=colors[j-1]
    pause
endfor

for j=0,n0-1 do begin
    file = list[j]
    readcol,file,format='a,f,f,x,f,f,x,x,x',bin1,vel,disp,h3,h4
    n1 = n_elements(bin1)
    
    if (j eq 0) then begin
        plot,rad,h4,psym=2,title='Velocity'
        loadct,4
    endif else oplot,rad,h4,psym=2,color=colors[j-1]
    pause
endfor

stop
END

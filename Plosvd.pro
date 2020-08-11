pro Plosvd, list, ERROR=error

; The code overplots LOSVDs for visual inspection. The list is a
; list of lists. So, if you just have one, your 'list' will contain
; the name of a single list, within which are all the losvds you want
; to overplot. If the ERROR keyword is set to anything, then the code
; assumes you are plotting the final losvds (with uncertainties). If
; left off, the code assumes you are plotting the output of fitlov

; Plotting to a ps file has been killed for now...

;***********************************************************
v1 = -1500 ;v1 and v2 are the upper and lower plotting range
v2 =  1500
;***********************************************************

if (n_elements(error) ne 0) then error = 'on' else error = 'off'

readcol,list,silent=1,format='a',listolists
n0 = n_elements(listolists)
    
for l=0,n0-1 do begin ;a loop through each list of losvds
    list = listolists[l]
    name = strsplit(list,'_.',/extract)
    name = name[0]
    readcol,silent=1,list,format='a',files

    n1 = n_elements(files)
    regions = strarr(n1)
    colors = (floor(255/n1)*indgen(n1)) + floor(255/n1)
    for j=0,n1-1 do begin
        temp = strsplit(files[j],'_.',/extract)
        regions[j] = temp[1]
    endfor
    if (error eq 'on') then $
      readcol,silent=1,files[0],f='f,f,f,f',vel,int,down,up
    if (error eq 'off') then $
      readcol,silent=1,files[0],f='x,f,f',vel,int,skipline=1
    
    if (l eq 0) then begin
        set_plot,'x'
        window,0,retain=2
        device,decomposed=0
    endif else wset,0

    loadct,0,silent=1
    plot,vel,int,title='LOSVDs for '+name,xrange=[v1,v2],$
      xstyle=1,yrange=[0,0.20],ystyle=1,xtitle='Velocity (km/sec)',$
      charsize=1.2,/nodata
    loadct,33,silent=1

    for j=0,n1-1 do begin
        if (error eq 'off') then begin
            readcol,silent=1,skipline=1,files[j],format='x,f,f',vel,int
            oplot,vel,int,color=colors[j]
            xyouts,0.1,0.9-(0.03*j),regions[j],color=colors[j],/normal
        endif
        if (error eq 'on') then begin
            readcol,silent=1,files[j],f='f,f,f,f',vel,int,down,up
            oplot,vel,int,color=colors[j]
            oplot,vel,down,color=colors[j],linestyle=2
            oplot,vel,up,color=colors[j],linestyle=2
            xyouts,0.1,0.9-(0.03*j),regions[j],color=colors[j],/normal
        endif
    endfor
    pause
endfor

stop
END

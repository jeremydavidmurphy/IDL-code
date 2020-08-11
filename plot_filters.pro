
pro plot_filters, list

;This routine is used to read in and generate .ps files for the VP
;test bench interference filters tested out at the observatory on the
;Cary.

;LIST: A list of the .csv files

readcol,list,silent=1,format='a',files

n1 = n_elements(files)

xup = intarr(1)
xdown = intarr(1)

for j=0,n1-1 do begin
    file = files[j]
    outname = strsplit(file,'.',/extract)
    outname = outname[0]
    readcol,file,silent=1,format='f,f',w,t
    yup = max(t)+3
    set_plot,'x'
    window,0,retain=2
    plot,w,t,yrange=[0,yup],ystyle=1
    print,'Fix the x-axes (lower/upper)'
    read,xdown
    read,xup
    
    set_plot,'ps'
    device,file=outname+'.ps',/color
    loadct,0
    plot,w,t,xrange=[xdown,xup],yrange=[0,yup],xstyle=1,$
      xtitle='Wavelength (A)',ytitle='% Transmission',$
      title='Transmission for filter'+outname,$
      xthick=3,ythick=3,charthick=3,/nodata
    loadct,4
    oplot,w,t,thick=3,color=150
    device,/close_file
endfor

stop
end

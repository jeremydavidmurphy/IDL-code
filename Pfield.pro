; MODIFIED on Nov 20, 2009: To include options to plot dots rather
; than contours. To output .eps files for inclusion into Latex. (This
; effectively supercedes the justplot.pro routine.)

; This routine is used to make colored contour plots of the output
; from pfitlov.out. It uses the coord files (often of the form
; name.coord). These are in the form given below and are the 

; The REGION is a string, matching the pfitlov.out naming convention,
; for the spectral region you want to plot. This is will be HK, GB,
; HB, etc. IF OMITTED, THE ROUTINE ASSUMES YOU HAVE JUST ONE REGION IN
; THE PFITLOV.OUT FILE (i.e. the final values, or a single region). If
; 'all' is used, the median of all the regions will be used.

; The COORD list is a list of all the b[in]n#### files. If not used
; the file "coord.list" is searched for. The format of these files is:
; 13 % murphy@grad78 work> more bnn1504.coord
; F #             dRA          dDec
;  81        10.49200      21.79800
;  24         0.26400
; 106       -19.55600       9.59900
;  21         0.23430
; where, because the format for readcol is format = 'i,f,f', the first
; line, and the other lines with 2 elements are skipped.

; PFIT is the name of the pfitlov.out file. If not used, the file
; "pfitlov.out" is searched for,

pro Pfield, REGION=region, PFIT=pfitfile, COORD=coord

;*************************************************************************
dotsize = 0.6
encap = 'n' ;set to "n" if you want to generate regular .ps files
plotH = 'yes' ;set to something else to suppress plotting H3 and H4
;*************************************************************************

if (n_elements(region) eq 0) then begin
    med = 'one'
;    outname = ''
    outname = 'Pfield.txt'
;    print,'Enter the name of the output file:'
;    read,outname
endif else begin
    if (region eq 'all') then med = 'on' else med = 'off'
endelse

if (n_elements(pfitfile) eq 0) then pfitfile = 'pfitlov.out'
if (n_elements(coord) eq 0) then coord = 'coord.list'

readcol,pfitfile,silent=1,format='a,f,f,x,f,f,x,x,x',name,vel,disp,h3,h4
readcol,coord,silent=1,format='a',coordfiles

n0 = n_elements(name)
bins = strarr(n0)
regions = strarr(n0)

ans1=''
test = strsplit(name[0],'_.',/extract)
print,''
print,transpose(test)
print,''
print,'Enter the element that gives the bin:'
read,ans1
ans1=uint(ans1)-1
if (med eq 'off') then begin
    ans2=''
    print,'Enter the element that gives the spectral region'
    print,'Or enter "0" if all the files are to be plotted.'
    read,ans2
    ans2=uint(ans2)-1
endif 

for j=0,n0-1 do begin
    temp = strsplit(name[j],'_.',/extract)
    bins[j] = temp[ans1]
    if (med eq 'off') then regions[j] = temp[ans2]
    print,bins[j],regions[j]
endfor

n1 = n_elements(coordfiles)

for j=0,n1-1 do begin           ;a loop through each bin
    file = coordfiles[j]
    readcol,file,silent=1,skipline=1,format='i,f,f',fiber,dRA,dDec
    ind = where(dRA ne -1.0,count) ;bad fibers are tossed
    if (count ne 0) then begin
        fiber = fiber[ind]	
        dRA = dRA[ind]
        dDec = dDec[ind]
    endif
    n3 = n_elements(dRA)
    file = strsplit(file,'.',/extract)
    onebin = file[0]
    ibin = where(bins eq onebin,count)
    if (count eq 0) then begin ;if there is no information in the pfitlov.out file
        print,'No information on '+onebin+' found in your pfitlov.out file!'
        goto, jump1
    endif
    
    if (med eq 'on') then begin
        piece = transpose([[vel[ibin]],[disp[ibin]],[h3[ibin]],[h4[ibin]]])
        medpiece = median(piece,dimension=2)
        if (j eq 0) then openw,5,'plotfield_median.dat'
        for k=0,n3-1 do printf,5,onebin,'median',dRA[k],dDec[k],medpiece,$
          format='(a9,2x,a6,2x,f8.3,2x,f8.3,2x,f10.3,2x,f10.3,2x,f7.3,2x,f7.3)'
        if (j eq 0) then openw,6,'temp.txt'
        for k=0,n3-1 do printf,6,onebin,'median',dRA[k],dDec[k],medpiece,$
          format='(a9,2x,a6,2x,f8.3,2x,f8.3,2x,f10.3,2x,f10.3,2x,f7.3,2x,f7.3)'
    endif
    
    if (med eq 'off') then begin
        ispec = where(regions eq region)
        indall = intarr(n_elements(ibin))
        for k=0,n_elements(ibin)-1 do indall[k] = where(ispec eq ibin[k])
        itemp = indall[where(indall ne -1)]
        if (itemp[0] eq -1) then stop
        ione = ispec[itemp]
        piece = [vel[ione],disp[ione],h3[ione],h4[ione]]
        if (j eq 0) then openw,5,'plotfield_'+region+'.dat'
        for k=0,n3-1 do printf,5,onebin,region,dRA[k],dDec[k],piece,$
        format='(a9,2x,a4,2x,f8.3,2x,f8.3,2x,f10.3,2x,f10.3,2x,f7.3,2x,f7.3)'
        if (j eq 0) then openw,6,'temp.txt'
        for k=0,n3-1 do printf,6,onebin,region,dRA[k],dDec[k],piece,$
        format='(a9,2x,a4,2x,f8.3,2x,f8.3,2x,f10.3,2x,f10.3,2x,f7.3,2x,f7.3)'
    endif

    if (med eq 'one') then begin
        piece = transpose([[vel[ibin]],[disp[ibin]],[h3[ibin]],[h4[ibin]]])
        if (j eq 0) then openw,5,outname
        for k=0,n3-1 do printf,5,onebin,dRA[k],dDec[k],piece,$
          format='(a9,2x,f8.3,2x,f8.3,2x,f10.3,2x,f10.3,2x,f7.3,2x,f7.3)'
        if (j eq 0) then openw,6,'temp.txt'
        for k=0,n3-1 do printf,6,onebin,dRA[k],dDec[k],piece,$
          format='(a9,2x,f8.3,2x,f8.3,2x,f10.3,2x,f10.3,2x,f7.3,2x,f7.3)'
    endif
jump1:
endfor

free_lun,5
free_lun,6

readcol,'temp.txt',format='a,f,f,f,f,f,f',binz,ra,dec,vel,disp,h3,h4,silent=1

ans3 = ''
print,'Plot dots or contours? (d/c):'
read,ans3

rup = max(vel)
rdown = min(vel)

if (ans3 eq 'c') then begin

    set_plot,'ps'

    if (encap eq 'y') then device,/encapsulated,file='VELcont.eps',$
    /color else device,file='VELcont.ps',/color
    loadct,0
    levels = (rup-rdown)/255.*findgen(255)+rdown
    loadct,0
    contour,vel,ra,dec,iso=1,levels=levels,/irregular,position=[0.33,0.09,0.85,0.85],$
      yrange=[min(dec),max(dec)],/ystyle,/fill,/nodata,xthick=3,ythick=3,thick=3,$
      xtitle='RA (arcsec)',ytitle='Dec (arcsec)',charthick=3,closed=1

    loadct,33
    contour,vel,ra,dec,/fill,/overplot,levels=levels,/irregular,iso=1,color=0
    
    loadct,0
    oplot,ra,dec,psym=sym(6),symsize=0.4

    loadct,33
    ncolors = n_elements(levels)
    loc = [0.22,0.90,0.69,0.95]
    bar = bindgen(256) # replicate(1b,10)
    xsize = (loc[2] - loc[0]) * !d.x_vsize
    ysize = (loc[3] - loc[1]) * !d.y_vsize
    xstart = loc[0] * !d.x_vsize
    ystart = loc[1] * !d.y_vsize

    tv, bar, xstart, ystart, xsize=xsize, ysize=ysize
    loadct,0

    plots,[loc[0],loc[0],loc[2],loc[2],loc[0]],[loc[1],loc[3],loc[3],loc[1],loc[1]],/normal

    middle = (max(vel)-min(vel))/2
    vup = middle
    vdown = -middle

    xyouts,0.22,0.86,strn(round(vdown)),/normal,charthick=3,charsize=1.5
    xyouts,0.63,0.86,strn(round(vup)),/normal,charthick=3,charsize=1.5
    xyouts,0.42,0.868,'km/sec',/normal,charthick=3,charsize=1.0
    device,/close_file

    rup = max(disp)
    rdown = min(disp)
    if (encap eq 'y') then device,/encapsulated,file='DISPcont.eps',$
    /color else device,file='VELcont.ps',/color
    loadct,0
    levels = (rup-rdown)/255.*findgen(255)+rdown
    loadct,0

    contour,disp,ra,dec,iso=1,levels=levels,/irregular,position=[0.33,0.09,0.85,0.85],$
      yrange=[min(dec),max(dec)],/ystyle,/fill,/nodata,xthick=3,ythick=3,thick=3,$
      xtitle='RA (arcsec)',ytitle='Dec (arcsec)',charthick=3

    loadct,33
    contour,disp,ra,dec,/fill,/overplot,levels=levels,/irregular,iso=1,color=0

    loadct,0
    oplot,ra,dec,psym=3
    
    loadct,33
    ncolors = n_elements(levels)
    loc = [0.22,0.90,0.69,0.95]
    bar = bindgen(256) # replicate(1b,10)
    xsize = (loc[2] - loc[0]) * !d.x_vsize
    ysize = (loc[3] - loc[1]) * !d.y_vsize
    xstart = loc[0] * !d.x_vsize
    ystart = loc[1] * !d.y_vsize

    tv, bar, xstart, ystart, xsize=xsize, ysize=ysize
    loadct,0
    
    plots,[loc[0],loc[0],loc[2],loc[2],loc[0]],[loc[1],loc[3],loc[3],loc[1],loc[1]],/normal

    xyouts,0.22,0.86,strn(round(rdown)),/normal,charthick=3,charsize=1.5
    xyouts,0.63,0.86,strn(round(rup)),/normal,charthick=3,charsize=1.5
    xyouts,0.42,0.868,'km/sec',/normal,charthick=3,charsize=1.0
    device,/close_file
    set_plot,'x'
endif

goto, jump2
; VELOCITY
;*********************************************************************************
if (ans3 eq 'd') then begin
    nnn = n_elements(ra)
    rup = max(vel)
    rdown = min(vel)

    set_plot,'ps'

    if (encap eq 'y') then device,/encapsulated,file='VELdots.eps',$
    /color else device,file='VELdots.ps',/color
    loadct,0
    levels = (rup-rdown)/255.*findgen(255)+rdown
    loadct,0
    
    spotcolor = intarr(nnn)
    si = bsort(vel)
    sv = vel[si]
    sra = ra[si]
    sdec = dec[si]

    slope = 255./(max(sv)-min(sv))
    for j=0,nnn-1 do spotcolor[j] = slope * (sv[j]-min(sv))
    
    loadct,0
    plot,sra,sdec,psym=3,/nodata,/isotropic,position=[0.33,0.10,0.85,0.88],$
      xthick=3,ythick=3,xrange=[min(ra)-10,max(ra)+10],/xstyle,xtitle='arcsec',$
      ytitle='arcsec',charthick=2,yrange=[min(dec)-10,max(dec)+10],/ystyle
    loadct,33
    for j=0,nnn-1 do plots,sra[j],sdec[j],psym=sym(1),color=spotcolor[j],$
      symsize=dotsize
    
    ncolors = n_elements(spotcolor)
    loc = [0.27,0.915,0.63,0.945]
    bar = bindgen(256) # replicate(1b,10)
    xsize = (loc[2] - loc[0]) * !d.x_vsize
    ysize = (loc[3] - loc[1]) * !d.y_vsize
    xstart = loc[0] * !d.x_vsize
    ystart = loc[1] * !d.y_vsize
    tv, bar, xstart, ystart, xsize=xsize, ysize=ysize
    loadct,0
    plots,[loc[0],loc[0],loc[2],loc[2],loc[0]],[loc[1],loc[3],loc[3],loc[1],loc[1]],/normal

    middle = (max(vel)-min(vel))/2
    vup = middle
    vdown = -middle
    xyouts,0.27,0.885,strn(round(vdown)),/normal,charthick=2,charsize=1.1
    xyouts,0.60,0.885,strn(round(vup)),/normal,charthick=2,charsize=1.1
    xyouts,0.42,0.89,'km/sec',/normal,charthick=2,charsize=0.8
;    xyouts,0.40,0.955,'M87 Velocity',/normal,charthick=2,charsize=1.0
    device,/close_file

; DISPERSION
;*********************************************************************************
    rup = max(disp)
    rdown = min(disp)
    if (encap eq 'y') then device,/encapsulated,file='DISPdots.eps',$
    /color else device,file='DISPdots.ps',/color
    loadct,0
    
    spotcolor = intarr(nnn)
    si = bsort(disp)
    sv = disp[si]
    sra = ra[si]
    sdec = dec[si]
    
    slope = 230./(max(sv)-min(sv))
    for j=0,nnn-1 do spotcolor[j] = slope * (sv[j]-min(sv))
    
    loadct,0
    plot,sra,sdec,psym=3,/nodata,/isotropic,position=[0.33,0.10,0.85,0.88],$
      xthick=3,ythick=3,xrange=[min(ra)-10,max(ra)+10],/xstyle,xtitle='arcsec',$
      ytitle='arcsec',charthick=2,yrange=[min(dec)-10,max(dec)+10],/ystyle
    loadct,33
    for j=0,nnn-1 do plots,sra[j],sdec[j],psym=sym(1),color=spotcolor[j],$
      symsize=dotsize
    
    loc = [0.27,0.915,0.63,0.945]
    bar = bindgen(240) # replicate(1b,10)
    xsize = (loc[2] - loc[0]) * !d.x_vsize
    ysize = (loc[3] - loc[1]) * !d.y_vsize
    xstart = loc[0] * !d.x_vsize
    ystart = loc[1] * !d.y_vsize
    tv, bar, xstart, ystart, xsize=xsize, ysize=ysize
    loadct,0
    plots,[loc[0],loc[0],loc[2],loc[2],loc[0]],[loc[1],loc[3],loc[3],loc[1],loc[1]],/normal

    dup = max(disp)
    ddown = min(disp)
    xyouts,0.27,0.885,strn(round(ddown)),/normal,charthick=2,charsize=1.1
    xyouts,0.59,0.885,strn(round(dup)),/normal,charthick=2,charsize=1.1
    xyouts,0.42,0.89,'km/sec',/normal,charthick=2,charsize=0.8
;    xyouts,0.39,0.955,'M87 Dispersion',/normal,charthick=2,charsize=1.0
    device,/close_file

    if (plotH eq 'yes') then begin
; H3
;*********************************************************************************
        rup = max(h3)
        rdown = min(h3)
        if (encap eq 'y') then device,/encapsulated,file='H3dots.eps',$
          /color else device,file='H3dots.ps',/color
        loadct,0
        
        spotcolor = intarr(nnn)
        si = bsort(h3)
        sv = h3[si]
        sra = ra[si]
        sdec = dec[si]
        
        slope = 230./(max(sv)-min(sv))
        for j=0,nnn-1 do spotcolor[j] = slope * (sv[j]-min(sv))
        
        loadct,0
        plot,sra,sdec,psym=3,/nodata,/isotropic,position=[0.33,0.10,0.85,0.88],$
          xthick=3,ythick=3,xrange=[min(ra)-10,max(ra)+10],/xstyle,xtitle='arcsec',$
          ytitle='arcsec',charthick=2,yrange=[min(dec)-10,max(dec)+10],/ystyle
        loadct,33
        for j=0,nnn-1 do plots,sra[j],sdec[j],psym=sym(1),color=spotcolor[j],$
          symsize=dotsize
        
        loc = [0.27,0.915,0.63,0.945]
        bar = bindgen(240) # replicate(1b,10)
        xsize = (loc[2] - loc[0]) * !d.x_vsize
        ysize = (loc[3] - loc[1]) * !d.y_vsize
        xstart = loc[0] * !d.x_vsize
        ystart = loc[1] * !d.y_vsize
        tv, bar, xstart, ystart, xsize=xsize, ysize=ysize
        loadct,0
        plots,[loc[0],loc[0],loc[2],loc[2],loc[0]],[loc[1],loc[3],loc[3],loc[1],loc[1]],/normal
        
        dup = max(h3)
        ddown = min(h3)
        print,dup
        dup = ''
        print,'Enter the plot range for this number:'
        read,dup
        print,ddown
        ddown = ''
        print,'Enter the plot range for this number:'
        read,ddown
        xyouts,0.27,0.885,ddown,/normal,charthick=2,charsize=1.1
        xyouts,0.565,0.885,dup,/normal,charthick=2,charsize=1.1
        xyouts,0.45,0.89,'H3',/normal,charthick=2,charsize=0.8
;    xyouts,0.39,0.955,'M87 H3',/normal,charthick=2,charsize=1.0
        device,/close_file
        
;H4
;*********************************************************************************
        rup = max(h4)
        rdown = min(h4)
        if (encap eq 'y') then device,/encapsulated,file='H4dots.eps',$
          /color else device,file='H4dots.ps',/color
        loadct,0
        
        spotcolor = intarr(nnn)
        si = bsort(h4)
        sv = h4[si]
        sra = ra[si]
        sdec = dec[si]
        
        slope = 230./(max(sv)-min(sv))
        for j=0,nnn-1 do spotcolor[j] = slope * (sv[j]-min(sv))
        
        loadct,0
        plot,sra,sdec,psym=3,/nodata,/isotropic,position=[0.33,0.10,0.85,0.88],$
          xthick=3,ythick=3,xrange=[min(ra)-10,max(ra)+10],/xstyle,xtitle='arcsec',$
          ytitle='arcsec',charthick=2,yrange=[min(dec)-10,max(dec)+10],/ystyle
        loadct,33
        for j=0,nnn-1 do plots,sra[j],sdec[j],psym=sym(1),color=spotcolor[j],$
          symsize=dotsize
        
        loc = [0.27,0.915,0.63,0.945]
        bar = bindgen(240) # replicate(1b,10)
        xsize = (loc[2] - loc[0]) * !d.x_vsize
        ysize = (loc[3] - loc[1]) * !d.y_vsize
        xstart = loc[0] * !d.x_vsize
        ystart = loc[1] * !d.y_vsize
        tv, bar, xstart, ystart, xsize=xsize, ysize=ysize
        loadct,0
        plots,[loc[0],loc[0],loc[2],loc[2],loc[0]],[loc[1],loc[3],loc[3],loc[1],loc[1]],/normal
        
        dup = max(h4)
        ddown = min(h4)
        print,dup
        dup = ''
        print,'Enter the plot range for this number:'
        read,dup
        print,ddown
        ddown = ''
        print,'Enter the plot range for this number:'
        read,ddown
        
        xyouts,0.27,0.885,ddown,/normal,charthick=2,charsize=1.1
        xyouts,0.565,0.885,dup,/normal,charthick=2,charsize=1.1
        xyouts,0.45,0.89,'H4',/normal,charthick=2,charsize=0.8
;    xyouts,0.39,0.955,'M87 H4',/normal,charthick=2,charsize=1.0
        device,/close_file

jump2:
set_plot,'ps'
;BIN LOCATION
;*********************************************************************************
        if (encap eq 'y') then device,/encapsulated,file='BinLoc.eps',$
          /color else device,file='BinLoc.ps',/color
        loadct,0
        
        nnn = n_elements(ra)
        spotcolor = intarr(nnn)
        slope = 195./n0
        for mm=0,n0-1 do begin
            ibin = where(binz eq bins[mm]) 
            spotcolor[ibin] = 60 + slope * mm
        endfor

        loadct,0
        plot,ra,dec,psym=3,/nodata,/isotropic,position=[0.33,0.10,0.85,0.88],$
          xthick=3,ythick=3,xrange=[min(ra)-10,max(ra)+10],/xstyle,xtitle='arcsec',$
          ytitle='arcsec',charthick=2,yrange=[min(dec)-10,max(dec)+10],/ystyle
        loadct,4
        for j=0,nnn-1 do plots,ra[j],dec[j],psym=sym(1),color=spotcolor[j],$
          symsize=dotsize
        
;        loc = [0.27,0.915,0.63,0.945]
;        bar = bindgen(240) # replicate(1b,10)
;        xsize = (loc[2] - loc[0]) * !d.x_vsize
;        ysize = (loc[3] - loc[1]) * !d.y_vsize
;        xstart = loc[0] * !d.x_vsize
;        ystart = loc[1] * !d.y_vsize
;        tv, bar, xstart, ystart, xsize=xsize, ysize=ysize
;        loadct,0
;        plots,[loc[0],loc[0],loc[2],loc[2],loc[0]],[loc[1],loc[3],loc[3],loc[1],loc[1]],/normal
        
;        dup = max(h4)
;        ddown = min(h4)
;        print,dup
;        dup = ''
;        print,'Enter the plot range for this number:'
;        read,dup
;        print,ddown
;        ddown = ''
;        print,'Enter the plot range for this number:'
;        read,ddown
        
;        xyouts,0.27,0.885,ddown,/normal,charthick=2,charsize=1.1
;        xyouts,0.565,0.885,dup,/normal,charthick=2,charsize=1.1
;        xyouts,0.45,0.89,'H4',/normal,charthick=2,charsize=0.8
;    xyouts,0.39,0.955,'M87 H4',/normal,charthick=2,charsize=1.0
        device,/close_file
    endif


endif

set_plot,'x'

stop
end

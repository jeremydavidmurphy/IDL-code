PRO LOSVD_combine, list

; This routine is used to combine different LOSVDs (i.e. the
; bin####.fit files). The LIST is assumed to be a list of lists, with
; each individual list being all the name.fit files you want to
; combine.

;COMPILE: biweightEF.pro (this now calcuates errors via SD and
;uncertainty of the biweight)

; MODIFIED on May 26th. This now outputs the text files in the same
; format as the mcfit2 files that fall into the rest of Karl's
; routines. The format for these files is now:
; VELOCITY   LOSVD   LOSVD-UNCERTAINTY   LOSVD+UNCERTAINTY

;***********************************************************
v1 = -1500 ;v1 and v2 are the upper and lower plotting range
v2 =  1500
plotting = 'off'
writerror = 'yes' ;set this to 'no' to supress writing out the uncertainty, which is calculated in two ways- the standard deviation and the biweight uncertainty.
slow = 'on'
scale = 0.85
writenorm = 'no'
;***********************************************************

readcol,list,silent=1,format='a',listolists
n0 = n_elements(listolists)
    
for j=0,n0-1 do begin ;the loop through the different lists
    list = listolists[j]
    name = strsplit(list,'_.',/extract)
    binname = name[0]

    readcol,silent=1,list,format='a',files
    n1 = n_elements(files)
    readcol,files[0],f='d,d,d,d',a,b,c,d,silent=1
    if (j eq 0) then begin
        aa = a & bb = b & cc = c & dd = d
    endif else begin
        aa = [aa,a]
        bb = [bb,b]
        cc = [cc,c]
        dd = [dd,d]
    endelse
    readcol,files[0],f='x,d,d',t1,t2,silent=1
    n2 = n_elements(t1)-1
    array = dblarr(3,n2,n1)
    colors1 = indgen(n1) * floor(255.0/n1)
    for k=0,n1-1 do begin ;the loop for an individual list
        readcol,silent=1,files[k],f='x,d,d',vel,int,skipline=1
;        offset = total(int[3:25])
        offset = total(int[5:23])
        array[0,*,k] = vel
        array[1,*,k] = int
        array[2,*,k] = int/offset
    endfor

    piece = array[1,*,*]
    piece = transpose(piece)
    piece = transpose(piece)
    LOSVD = biweightEF(piece)
    LOSVDMED = median(piece,dim=2)
    ERROR = LOSVD[*,1] * scale
    izero = where(ERROR lt 0.01,c)
    if (c gt 0) then ERROR[izero] = 0.005
    LOSVD = LOSVD[*,0]
    SD = dblarr(n2)
    for l=0,n2-1 do SD[l] = stddev(piece[l,*]) * scale

    piece = array[2,*,*]
    piece = transpose(piece)
    piece = transpose(piece)
    LOSVDn = biweightEF(piece)
    LOSVDMEDn = median(piece,dim=2)
    ERRORn = LOSVDn[*,1] * scale
    LOSVDn = LOSVDn[*,0]
    SDn = dblarr(n2)
    for l=0,n2-1 do SDn[l] = stddev(piece[l,*]) * scale

;    ibig = where(abs(vel) gt 800) ;fixing the larger velocity uncertainty to large values
;    SD[ibig] = 0.03

    if (plotting eq 'on' and j eq 0) then begin
        set_plot,'x'
        window,0,retain=2
        device,decomposed=0
        window,2,retain=2
        device,decomposed=0
    endif

    if (plotting eq 'on') then begin
        wset,0
        loadct,0,/silent
        plot,array[0,*,0],array[1,*,0],xrange=[v1,v2],/xstyle,/nodata,$
          yrange=[0,0.16],/ystyle,title=binname,xtitle='Velocity'
        loadct,27,/silent
        for k=0,n1-1 do begin
            oplot,array[0,*,k],array[1,*,k],color=colors1[k]
        endfor
        loadct,0,/silent
        oplot,vel,LOSVD,thick=2.5

        wset,2
        loadct,0,/silent
        plot,vel,LOSVD,thick=1.5,title=binname+': Uncertainty: Biwt in white, Median and S.D. in yellow',$
          xrange=[v1,v2],/xstyle,xtitle='Velocity',yrange=[0,0.16],/ystyle
        oploterr,vel,LOSVD,ERROR
        loadct,4
        oplot,vel,LOSVDMED,color=255,thick=1.5
        oploterr,vel+10.0,LOSVDMED,SD
        if (slow eq 'on') then wait,1.0

        wset,0
        loadct,0,/silent
        plot,array[0,*,0],array[2,*,0],xrange=[v1,v2],/xstyle,/nodata,$
          yrange=[0,0.16],/ystyle,title='Normalized LOSVD: '+binname,xtitle='Velocity'
        loadct,27,/silent
        for k=0,n1-1 do begin
            oplot,array[0,*,k],array[2,*,k],color=colors1[k]
        endfor
        loadct,0,/silent
        oplot,vel,LOSVDn,thick=2.5

        wset,2
        loadct,0,/silent
        plot,vel,LOSVDn,thick=1.5,title=binname+': Uncertainty: Biwt in white, Median and S.D. in yellow',$
          xrange=[v1,v2],/xstyle,xtitle='Velocity',yrange=[0,0.16],/ystyle
        oploterr,vel,LOSVDn,ERRORn
        loadct,4
        oplot,vel,LOSVDMEDn,color=255,thick=1.5
        oploterr,vel+10.0,LOSVDMEDn,SDn
        if (slow eq 'on') then wait,1.0
        
    endif

if (writenorm eq 'yes') then begin
;*************************************************
    LOSVD = LOSVDn
    ERROR = ERRORn
    SD = SDn
    LOSVDMED = LOSVDMEDn
;*************************************************
endif

    free_lun,5
    openw,5,binname+'_LOSVDbi.txt'
    for k=0,28 do begin
        if (writerror eq 'yes') then $
        printf,5,strn(vel[k]),' ',strn(LOSVD[k]),' ',strn(LOSVD[k]-ERROR[k]),' ',strn(LOSVD[k]+ERROR[k]) $
          else printf,5,strn(k+1),' ',strn(vel[k]),' ',strn(LOSVD[k])
    endfor
    free_lun,5
    openw,5,binname+'_LOSVDmed.txt'
    for k=0,28 do begin
        if (writerror eq 'yes') then $
        printf,5,strn(vel[k]),' ',strn(LOSVDMED[k]),' ',strn(LOSVDMED[k]-SD[k]),' ',strn(LOSVDMED[k]+SD[k]) $
          else printf,5,strn(k+1),' ',strn(vel[k]),' ',strn(LOSVDMED[k])
    endfor

goto, jumpend    

temp = strsplit(listolists[j],'.',/extract)
psname = temp[0]+'_LOSVD.ps'
set_plot,'ps'
device,file=psname,/color
loadct,0,/silent
plot,array[0,*,0],array[1,*,0],xrange=[v1,v2],/xstyle,/nodata,$
  yrange=[0,0.16],/ystyle,title=binname,xtitle='Velocity'
loadct,27,/silent
for k=0,n1-1 do begin
    oplot,array[0,*,k],array[1,*,k],color=colors1[k]
endfor
loadct,0,/silent
oplot,vel,LOSVD,thick=2.5
device,/close_file
set_plot,'x'

jumpend:

endfor
if (plotting eq 'on') then begin
    print,'The next ENTER deletes the figures.'
    pause
    while !d.window ne -1 do wdelete, !d.window
endif

stop
END

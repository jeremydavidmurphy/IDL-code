pro Plosvd_over, list, LOOP=loop

; This routine generates overplotted LOSVD's to either the screen or a
; postscript file. If LOOP is set to loop then it assumes the list is
; a list of lists

; 

;***********************************************************
v1 = -1500 ;v1 and v2 are the upper and lower plotting range
v2 =  1500
;***********************************************************

if (n_elements(loop) ne 0) then loop = 'on' else loop = 'off'

skip = 'on'
if (loop eq 'on') then begin
    ans=''
    print,'Plot to screen or postscript file? (s/f)'
    read,ans

    readcol,list,silent=1,format='a',listolists
    n0 = n_elements(listolists)
    
    for l=0,n0-1 do begin
        list = listolists[l]
        name = strsplit(list,'_.',/extract)
        name = name[0]
        readcol,silent=1,list,format='a',files

        n1 = n_elements(files)
        bins = strarr(n1)
        regions = bins
        
        colors = (floor(255/n1)*indgen(n1)) + floor(255/n1)
;        colors = [60,110,180,150]
        for j=0,n1-1 do begin
            temp = strsplit(files[j],'_.',/extract)
            bins[j] = temp[0]
            regions[j] = temp[1]
        endfor
        ttt = bins[0]+'.losvd'
        readcol,ttt,f='f,f,f,f',silent=1,Mvel,Mint,Mdown,Mup
        print,total(Mint)
        Mscale = total(Mint)
        Mint = Mint / Mscale
        Mdown = Mdown / Mscale
        Mup = Mup / Mscale
;        print,total(Mint) / (1000.0/29.0)
;        Mscale = (1000.0/29.0)/total(Mint)
;        print,mscale
;        print,total(mint * mscale)
;        Mint = Mint * Mscale
;        Mdown = Mdown * Mscale
;        Mup = Mup * Mscale
;        print,total(Mint) / (1000.0/29.0)
;       print,mscale
;        readcol,silent=1,skipline=1,files[0],format='x,f,f',vel,int
        readcol,silent=1,files[0],f='f,f,f,f',vel,int,down,up
;        vel = [v1,vel,v2]
;        int = [0,int,0]
        n2 = n_elements(vel)
        
        if (ans eq 's') then begin
            window,0,retain=2
            device,decomposed=0
            loadct,0
            plot,Mvel,Mint,title='LOSVDs for '+bins[0],xrange=[v1,v2],$
              xstyle=1,yrange=[0,0.20],ystyle=1,xtitle='Velocity (km/sec)',$
              charsize=1.2,/nodata
;            oplot,Mvel,Mdown,linestyle=2
;            oplot,Mvel,Mup,linestyle=2
            loadct,27
            ttl = fltarr(29)
            for j=0,n1-1 do begin
;                readcol,silent=1,skipline=1,files[j],format='x,f,f',vel,int
;                vel = [v1,vel,v2]
;                int = [0,int,0]
                readcol,silent=1,files[j],f='f,f,f,f',vel,int,down,up
                scale = total(int)
                int = int/scale
                down = down/scale
                up = up/scale
                print,'Normalization for '+bins[0]+' :'+strn(total(int))
                oplot,vel,int,color=colors[j]
;               oplot,vel,down,linestyle=2,color=colors[j]
;               oplot,vel,up,linestyle=2,color=colors[j]
                xyouts,0.15,(0.90-(j*0.03)),regions[j],color=colors[j],charsize=1.3,/normal
                ttl = [[ttl],[int]]
            endfor
            loadct,0
            oplot,Mvel,Mdown,linestyle=2,thick=1
            oplot,Mvel,Mup,linestyle=2,thick=1
            oplot,Mvel,Mint,thick=1
            if (skip eq 'on') then goto, jump1
            ttl = ttl[*,1:n1]
            std = fltarr(n_elements(ttl[*,0]))
            mn = fltarr(n_elements(ttl[*,0]))
            for k=0,n_elements(ttl[*,0])-1 do begin
                mn[k] = mean(ttl[k,*])
                std[k] = stddev(ttl[k,*])
            endfor
            std = std / sqrt(n_elements(ttl[0,*]))
;           pause
           loadct,4
;          std = std / n_elements(ttl[0,*]) ;TURN ME OFF!
           oplot,vel,mn,thick=2,color=60
           oplot,vel,mn+std,color=60,linestyle=2
           oplot,vel,mn-std,color=60,linestyle=2
           loadct,0
           jump1:
           pause
        endif


        if (ans eq 'f') then begin

            set_plot,'ps'
            device,file=name+'_losvd.ps',/color
            loadct,0
            plot,vel,int,/nodata,xrange=[v1,v2],xstyle=1,$
              yrange=[0,0.15],ystyle=1,xtitle='!3Velocity (km/sec)',charsize=1.2,$
              xthick=5,ythick=5,charthick=3
            loadct,4
            for j=0,n1-1 do begin
                readcol,silent=1,skipline=1,files[j],format='x,f,f',vel,int
                vel = [v1,vel,v2]
                int = [0,int,0]
                if (j eq 0) then begin
                    loadct,0
                    oplot,vel,int,thick=5
                    loadct,4
                endif else begin
                    oplot,vel,int,thick=5,color=colors[j-1]
                endelse
;                xyouts,0.17,(0.88-(j*0.035)),regions[j],color=colors[j],charsize=1.3,/normal,charthick=3
            endfor
            device,/close_file
            set_plot,'x'
        endif
    endfor
endif

if (loop eq 'off') then begin
    ans=''
    print,'Plot to screen or postscript file? (s/f)'
    read,ans
    readcol,silent=1,list,format='a',files
    name = strsplit(list,'_.',/extract)
    name = name[0]
    n1 = n_elements(files)
    bins = strarr(n1)
    regions = bins
        
    colors = floor(255/n1)*indgen(n1)
    for j=0,n1-1 do begin
        temp = strsplit(files[j],'_.',/extract)
        bins[j] = temp[0]
        regions[j] = temp[2]
    endfor
        
    readcol,silent=1,skipline=1,files[0],format='x,f,f',vel,int
    vel = [v1,vel,v2]
    int = [0,int,0]
    n2 = n_elements(vel)
        
    ans2 = ''
    intall = fltarr(n_elements(vel))

    if (ans eq 's') then begin
        print,'Take median? (y/n)'
        read,ans2
        window,0,retain=2
        device,decomposed=0
        loadct,0
        plot,vel,int,title='LOSVDs for '+bins[0],/nodata,xrange=[v1,v2],xstyle=1,$
          yrange=[0,0.15],ystyle=1,xtitle='Velocity (km/sec)',charsize=1.2
        loadct,27
        for j=0,n1-1 do begin
            readcol,silent=1,skipline=1,files[j],format='x,f,f',vel,int
            vel = [v1,vel,v2]
            int = [0,int,0]
            if (ans2 eq 'y') then intall = [[intall],[int]]
            oplot,vel,int,color=colors[j]
            xyouts,0.15,(0.90-(j*0.03)),regions[j],color=colors[j],charsize=1.3,/normal
        endfor
        if (ans2 eq 'y') then begin
            intall = intall[*,1:*]
            mint = median(intall,dimension=2)
            meanint = fltarr(n2)
            for j=0,n2-1 do meanint[j] = mean(intall[j,*])
        endif
        loadct,0
        oplot,vel,meanint,thick=3
    endif

    if (ans eq 'f') then begin
        
        set_plot,'ps'
        device,file=name+'_losvd.ps',/color
        loadct,0
        plot,vel,int,title='LOSVDs for '+bins[0],/nodata,xrange=[v1,v2],xstyle=1,$
          yrange=[0,0.15],ystyle=1,xtitle='Velocity (km/sec)',charsize=1.2,$
          xthick=3,ythick=3,charthick=3
        loadct,4
        for j=0,n1-1 do begin
            readcol,silent=1,skipline=1,files[j],format='x,f,f',vel,int
            vel = [v1,vel,v2]
            int = [0,int,0]
            oplot,vel,int,color=colors[j],thick=3
            xyouts,0.17,(0.88-(j*0.035)),regions[j],color=colors[j],charsize=1.3,/normal,charthick=3
        endfor
        device,/close_file
        set_plot,'x'
    endif
endif
pause
if (ans eq 's') then wdelete

stop
end

PRO Plifetime2, headlist, track=track

; This routine is used to plot the output of lifetime2.pro. It will
; read in the HEADER_1x1_list.txt file and plot various things against
; one another to look for correlations in the fiber output and things
; like temperature, Fracker position, etc.

;********************************************************************
;p2s = 'screen'
p2s = 'ps'
ccdtemp = -32.5
;ccdtemp = -24.5
;track='Dec: 38@2X' ;just gets written onto the plots
kill = 400.0
;********************************************************************

if (n_elements(track) ne 1) then begin
    track = ''
    print,'What track is this?'
    read,track
endif

;readcol,headlist,silent=1,f='a,f,f,f,f,f,f,f,f,f,f',$
;  name,iter,temp,xpos,ypos,rpos,bpos,temp1,temp2,Q1,Q2
readcol,headlist,silent=1,f='a,f,f,f,f,f,f,f,f,f,f,f,f',$
  name,iter,temp,xpos,ypos,rpos,bpos,temp1,temp2,hum1,hum2,Q1,Q2

ihigh = where(temp gt ccdtemp,cnt)
print,'The number of frames with high CCD temperature is '+strn(cnt)
wait,1.0

i = where(Q1 ne -1.0 and Q1 lt kill)
name = name[i]
iter = iter[i]
temp = temp[i]
xpos = xpos[i]
ypos = ypos[i]
rpos = rpos[i]
bpos = bpos[i]
temp1 = temp1[i]
temp2 = temp2[i]
hum1 = hum1[i]
hum2 = hum2[i]
Q1 = Q1[i]
Q2 = Q2[i]

i = bsort(Q1)

colors = Q1[i]
colors = float(colors/max(colors)) * 255.0
;colors = (float(colors/max(colors)) * 245.0) + 10.0

y1t1 = min(temp1)-0.1
y2t1 = max(temp1)+0.1
y1t2 = min(temp2)-0.25
y2t2 = max(temp2)+0.25
y1h1 = min(hum1)-1.0
y2h1 = max(hum1)+1.0
y1h2 = min(hum2)-1.0
y2h2 = max(hum2)+1.0

if (max(xpos) gt 1.0) then osx = 0.25 else osx = 0.02
if (max(ypos) gt 1.0) then osy = 0.25 else osy = 0.05
if (max(rpos) gt 1.0) then osr = 0.25 else osr = 0.05
if (max(bpos) gt 1.0) then osb = 0.25 else osb = 0.05

if (p2s eq 'screen') then begin
    set_plot,'x'
    window,2,retain=2
    device,decomposed=0
    !p.multi = [0,0,0,0,0]

    yl = min(xpos)-osx
    yu = max(xpos)+osx
    loadct,0,/silent
    plot,Q1,xpos,psym=1,/ynozero,title='Fracker X Position',$
      yrange=[yl,yu],/ystyle
    oplot,Q2,xpos,psym=4
    pause

    yl = min(ypos)-osy
    yu = max(ypos)+osy
    plot,Q1,ypos,psym=1,/ynozero,title='Fracker Y Position',$
      yrange=[yl,yu],/ystyle
    oplot,Q2,ypos,psym=4
    pause

    yl = min(rpos)-osr
    yu = max(rpos)+osr
    plot,Q1,rpos,psym=1,/ynozero,title='Fracker Rho Position',$
      yrange=[yl,yu],/ystyle
    oplot,Q2,rpos,psym=4
    pause

    yl = min(bpos)-osb
    yu = max(bpos)+osb
    plot,Q1,bpos,psym=1,/ynozero,title='Lower Stage position',$
      yrange=[yl,yu],/ystyle
    oplot,Q2,bpos,psym=4
    pause

    plot,Q1,temp1,psym=1,/ynozero,title='Temperature 1',$
      yrange=[y1t1,y2t1],/ystyle
    oplot,Q2,temp1,psym=4
    pause

    plot,Q1,temp2,psym=1,/ynozero,title='Temperature 2',$
      yrange=[y1t2,y2t2],/ystyle
    oplot,Q2,temp2,psym=4
    pause

    plot,Q1,hum1,psym=1,/ynozero,title='Humidity 1',$
      yrange=[y1h1,y2h1],/ystyle
    oplot,Q2,hum1,psym=4
    pause

    plot,Q1,hum2,psym=1,/ynozero,title='Humidity 2',$
      yrange=[y1h2,y2h2],/ystyle
    oplot,Q2,hum2,psym=4
    pause

;Now, the coordinates are plotted against themselves.
    loadct,16,/silent
    plot,xpos,ypos,psym=1,/ynozero,title='X vs Y',/nodata,$
      xtitle='X Position',ytitle='Y Position',color=150
    loadct,39,/silent
    for j=0,n_elements(xpos)-1 do plots,xpos[i[j]],ypos[i[j]],psym=1,color=colors[j]
    pause

    loadct,16,/silent
    plot,ypos,rpos,psym=1,/ynozero,title='Y vs Rho',/nodata,$
      xtitle='Y Position',ytitle='Rho Position',color=150
    loadct,39,/silent
    for j=0,n_elements(ypos)-1 do plots,ypos[i[j]],rpos[i[j]],psym=1,color=colors[j]
    pause

    loadct,16,/silent
    plot,xpos,rpos,psym=1,/ynozero,title='X vs Rho',/nodata,$
      xtitle='X Position',ytitle='Rho Position',color=150
    loadct,39,/silent
    for j=0,n_elements(xpos)-1 do plots,xpos[i[j]],rpos[i[j]],psym=1,color=colors[j]
    pause

    loadct,16,/silent
    plot,xpos,bpos,psym=1,/ynozero,title='X vs Bench',/nodata,$
      xtitle='X Position',ytitle='Bench Position',color=150
    loadct,39,/silent
    for j=0,n_elements(xpos)-1 do plots,xpos[i[j]],bpos[i[j]],psym=1,color=colors[j]
    pause

    loadct,16,/silent
    plot,ypos,bpos,psym=1,/ynozero,title='Y vs Bench',/nodata,$
      xtitle='Y Position',ytitle='Bench Position',color=150
    loadct,39,/silent
    for j=0,n_elements(ypos)-1 do plots,ypos[i[j]],bpos[i[j]],psym=1,color=colors[j]
    pause

    loadct,16,/silent
    plot,rpos,bpos,psym=1,/ynozero,title='Rho vs Bench',/nodata,$
      xtitle='Rho Position',ytitle='Bench Position',color=150
    loadct,39,/silent
    for j=0,n_elements(rpos)-1 do plots,rpos[i[j]],bpos[i[j]],psym=1,color=colors[j]
    pause

endif

if (p2s eq 'screen') then begin
    print,'The next ENTER deletes the plots'
    pause
    while !d.window ne -1 do wdelete, !d.window
endif

if (p2s eq 'ps') then begin
    name1 = name[0]
    name2 = name[n_elements(name)-1]
    s = strsplit(headlist,'.',/extract)
    ss1 = s[0]+'.ps'
    ss2 = s[0]+'_tracks.ps'

    set_plot,'ps'
    device,file=ss1,/color
    !p.multi = [0,2,4,0,0]

    yl = min(xpos)-osx
    yu = max(xpos)+osx
    plot,Q1,xpos,psym=1,/ynozero,title='Fracker X Position',$
      yrange=[yl,yu],/ystyle,symsize=0.7,xtitle='Severity of Event',$
      ytitle='X position'
    oplot,Q2,xpos,psym=4,symsize=0.7

    yl = min(ypos)-osy
    yu = max(ypos)+osy
    plot,Q1,ypos,psym=1,/ynozero,title='Fracker Y Position',$
      yrange=[yl,yu],/ystyle,symsize=0.7,xtitle='Severity of Event',$
      ytitle='Y position'
    oplot,Q2,ypos,psym=4,symsize=0.7

    yl = min(rpos)-osr
    yu = max(rpos)+osr
    plot,Q1,rpos,psym=1,/ynozero,title='Fracker Rho Position',$
      yrange=[yl,yu],/ystyle,symsize=0.7,xtitle='Severity of Event',$
      ytitle='Rho position'
    oplot,Q2,rpos,psym=4,symsize=0.7

    yl = min(bpos)-osb
    yu = max(bpos)+osb
    plot,Q1,bpos,psym=1,/ynozero,title='Lower Stage position',$
      yrange=[yl,yu],/ystyle,symsize=0.7,xtitle='Severity of Event',$
      ytitle='Lower Stage position'
    oplot,Q2,bpos,psym=4,symsize=0.7

    plot,Q1,temp1,psym=1,/ynozero,title='Temperature 1',$
      yrange=[y1t1,y2t1],/ystyle,symsize=0.7,xtitle='Severity of Event',$
      ytitle='Temperature 1 (C)'
    oplot,Q2,temp1,psym=4,symsize=0.7

    plot,Q1,temp2,psym=1,/ynozero,title='Temperature 2',$
      yrange=[y1t2,y2t2],/ystyle,symsize=0.7,xtitle='Severity of Event',$
      ytitle='Temperature 2 (C)'
    oplot,Q2,temp2,psym=4,symsize=0.7

    plot,Q1,hum1,psym=1,/ynozero,title='Humidity 1',$
      yrange=[y1h1,y2h1],/ystyle,symsize=0.7,xtitle='Severity of Event',$
      ytitle='Humidity 1'
    oplot,Q2,hum1,psym=4,symsize=0.7

    plot,Q1,hum2,psym=1,/ynozero,title='Humidity 2',$
      yrange=[y1h2,y2h2],/ystyle,symsize=0.7,xtitle='Severity of Event',$
      ytitle='Humidity 2'
    oplot,Q2,hum2,psym=4,symsize=0.7

    loadct,4,/silent
    xyouts,0.50,0.49,'First frame: '+name1,/normal,color=150,charsize=0.6,orient=90
    xyouts,0.515,0.49,'Last frame: '+name2,/normal,color=150,charsize=0.6,orient=90
    xyouts,0.46,0.24,'Number of frames: '+strn(n_elements(Q1)),/normal,$
      charsize=1.0,color=150,charthick=1.5
    xyouts,0.30,0.24,track,/normal,charsize=1.0,color=150,charthick=1.5
    loadct,0,/silent
    
    device,/close_file

    device,file=ss2,/color
    !p.multi = [0,2,3,0,0]

    loadct,0,/silent

    plot,xpos[i],ypos[i],psym=1,/ynozero,title='X vs Y',xtitle='X Position',$
      ytitle='Y Position',symsize=0.7,xrange=[min(xpos[i])-osx,max(xpos[i])+osx],$
      yrange=[min(ypos[i])-osy,max(ypos[i])+osy],/xstyle,/ystyle,/nodata,charsize=1.2
    loadct,39,/silent
    for j=0,n_elements(xpos[i])-1 do plots,xpos[i[j]],ypos[i[j]],psym=1,color=colors[j]

    loadct,0,/silent
    plot,ypos[i],rpos[i],psym=1,/ynozero,title='Y vs Rho',xtitle='Y Position',$
      ytitle='Rho Position',symsize=0.7,xrange=[min(ypos[i])-osy,max(ypos[i])+osy],$
      yrange=[min(rpos[i])-osr,max(rpos[i])+osr],/xstyle,/ystyle,/nodata,charsize=1.2
    loadct,39,/silent
    for j=0,n_elements(xpos[i])-1 do plots,ypos[i[j]],rpos[i[j]],psym=1,color=colors[j]

    loadct,0,/silent
    plot,xpos[i],rpos[i],psym=1,/ynozero,title='X vs Rho',xtitle='X Position',$
      ytitle='Rho Position',symsize=0.7,xrange=[min(xpos[i])-osx,max(xpos[i])+osx],$
      yrange=[min(rpos[i])-osr,max(rpos[i])+osr],/xstyle,/ystyle,/nodata,charsize=1.2
   loadct,39,/silent
    for j=0,n_elements(xpos[i])-1 do plots,xpos[i[j]],rpos[i[j]],psym=1,color=colors[j]

    loadct,0,/silent
    plot,xpos[i],bpos[i],psym=1,/ynozero,title='X vs Bench',xtitle='X Position',$
      ytitle='Bench Position',symsize=0.7,xrange=[min(xpos[i])-0.1,max(xpos[i])+0.1],$
      yrange=[min(bpos[i])-0.1,max(bpos[i])+0.1],/xstyle,/ystyle,/nodata,charsize=1.2
   loadct,39,/silent
    for j=0,n_elements(xpos[i])-1 do plots,xpos[i[j]],bpos[i[j]],psym=1,color=colors[j]

    loadct,0,/silent
    plot,ypos[i],bpos[i],psym=1,/ynozero,title='Y vs Bench',xtitle='Y Position',$
      ytitle='Bench Position',symsize=0.7,xrange=[min(ypos[i])-0.1,max(ypos[i])+0.1],$
      yrange=[min(bpos[i])-0.1,max(bpos[i])+0.1],/xstyle,/ystyle,/nodata,charsize=1.2
   loadct,39,/silent
    for j=0,n_elements(xpos[i])-1 do plots,ypos[i[j]],bpos[i[j]],psym=1,color=colors[j]

    loadct,0,/silent
    plot,rpos[i],bpos[i],psym=1,/ynozero,title='Rho vs Bench',xtitle='Rho Position',$
      ytitle='Bench Position',symsize=0.7,xrange=[min(rpos[i])-0.1,max(rpos[i])+0.1],$
      yrange=[min(bpos[i])-0.1,max(bpos[i])+0.1],/xstyle,/ystyle,/nodata,charsize=1.2
   loadct,39,/silent
    for j=0,n_elements(xpos[i])-1 do plots,rpos[i[j]],bpos[i[j]],psym=1,color=colors[j]

    loadct,4,/silent
    xyouts,0.49,0.49,'First frame: '+name1,/normal,color=150,charsize=0.6,orient=90
    xyouts,0.505,0.49,'Last frame: '+name2,/normal,color=150,charsize=0.6,orient=90
    xyouts,0.46,0.33,'Number of frames: '+strn(n_elements(Q1)),/normal,$
      charsize=1.0,color=150,charthick=1.5
    xyouts,0.30,0.33,track,/normal,charsize=1.0,color=150,charthick=1.5
    loadct,0,/silent

    device,/close_file
    set_plot,'x'

endif
stop
END

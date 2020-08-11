PRO pltfrd3
restore,'vp2data.idl'

eb = fltarr(2,30)
eb337pm = eb
eb365pm = eb
eb400pm = eb
eb600pm = eb
i337=[0,4,8,12,16,20,24,28,32,36]
i365=i337+1
i400=i365+1
i600=i400+1

for j=0,29 do begin
    big = max(data[j,7,i337])
    small = min(data[j,7,i337])
    eb337pm[0,j] = big
    eb337pm[1,j] = small
endfor

for j=0,29 do begin
    big = max(data[j,7,i365])
    small = min(data[j,7,i365])
    eb365pm[0,j] = big
    eb365pm[1,j] = small
endfor

for j=0,29 do begin
    big = max(data[j,7,i400])
    small = min(data[j,7,i400])
    eb400pm[0,j] = big
    eb400pm[1,j] = small
endfor

for j=0,29 do begin
    big = max(data[j,7,i600])
    small = min(data[j,7,i600])
    eb600pm[0,j] = big
    eb600pm[1,j] = small
endfor


ttl='F/#!dout!n for 10 Polymicro VP2 Fibers'
;window,0,retain=2
set_plot,'ps'
device,file='frd3.ps',/color
loadct,0
plot,m337pm,ee337,xrange=[max(m337pm)+0.1,min(m337pm)-0.1],psym=1,$
  yrange=[0.68,1.01],xstyle=1,ystyle=1,thick=2,xthick=2,ythick=2,$
  title=ttl,xtitle='Output F/#',ytitle='Encircled Energy',charsize=1.2,$
  charthick=2,/nodata
loadct,4

;oplot,m337pm,ee337,psym=1,color=60
;for j=0,29 do begin
;    y = [ee337[j],ee337[j]]
;    oplot,eb337pm[*,j],y,color=60
;endfor
oplot,m365pm,ee365,psym=1,color=60
for j=0,29 do begin
    y = [ee365[j],ee365[j]]
    oplot,eb365pm[*,j],y,color=60
endfor
oplot,m400pm,ee400,psym=1,color=110
for j=0,29 do begin
    y = [ee400[j],ee400[j]]
    oplot,eb400pm[*,j],y,color=110
endfor
oplot,m600pm,ee600,psym=1,color=150
for j=0,29 do begin
    y = [ee600[j],ee600[j]]
    oplot,eb600pm[*,j],y,color=150
endfor
lx1=[3.65,3.65]
ly1=[0.6,1.1]
oplot,lx1,ly1,color=60,thick=3
xyouts,3.66,0.73,'INPUT BEAM F/3.65',charthick=2.0,orientation=90,color=60,charsize=1.0

plots,3.5,0.8,psym=1,symsize=2,thick=4,color=60
xyouts,3.47,0.795,': 365 nm',color=60,charthick=2,charsize=1.3
plots,3.5,0.77,psym=1,symsize=2,thick=4,color=110
xyouts,3.47,0.765,': 400 nm',color=110,charthick=2,charsize=1.3
plots,3.5,0.74,psym=1,symsize=2,thick=4,color=150
xyouts,3.47,0.735,': 600 nm',color=150,charthick=2,charsize=1.3

device,/close_file
set_plot,'x'

stop
END

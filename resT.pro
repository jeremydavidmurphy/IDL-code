PRO resT, file

; This function reads in a bin####.res file (the output of
; resolution.pro) and returns a file bin####.resT, which is a text
; file of three numbers, giving the median of the bin####.res file 
; The input file must be in the form
; 4046.55  4077.83  4358.33  4678.15  4799.91  4916.07  5085.82  5154.66  5460.74
;  4.8362   4.8057   4.5354   4.4667   4.4744   4.2674   4.5164   4.8774   4.6218
;  4.8186   4.7709   4.5216   4.4607   4.4585   4.2687   4.4807   5.1480   4.5776
;  4.9689   4.8800   4.4006   4.1535   4.1059   4.0006   4.0438   4.0645   3.9775

yplot = 'on' ;change to 'off' if you want to supress the plotting

readcol,file,b1,b2,b3,g1,g2,g3,r1,r2,r3

pwave = [[b1[0],b2[0],b3[0],g1[0],g2[0],g3[0],r1[0],r2[0],r3[0]]]

wave = fltarr(3)
wave[0] = mean([b1[0],b2[0],b3[0]])
wave[1] = mean([g1[0],g2[0],g3[0]])
wave[2] = mean([r1[0],r2[0],r3[0]])

readcol,file,b1,b2,b3,g1,g2,g3,r1,r2,r3,skipline=1

n0 = n_elements(b1)
alldata = fltarr(9,n0)
for j=0,n0-1 do alldata[*,j] = [b1[j],b2[j],b3[j],g1[j],g2[j],g3[j],r1[j],r2[j],r3[j]]

bmean = mean([[b1],[b2],[b3]])
bmedian = median([[b1],[b2],[b3]],/even)
bvar = variance([[b1],[b2],[b3]])
gmean = mean([[g1],[g2],[g3]])
gmedian = median([[g1],[g2],[g3]],/even)
gvar = variance([[g1],[g2],[g3]])
rmean = mean([[r1],[r2],[r3]])
rmedian = median([[r1],[r2],[r3]],/even)
rvar = variance([[r1],[r2],[r3]])

allmean = [mean(b1),mean(b2),mean(b3),mean(g1),mean(g2),mean(g3),$
           mean(r1),mean(r2),mean(r3)]
allmedian = [median(b1),median(b2),median(b3),median(g1),median(g2),median(g3),$
           median(r1),median(r2),median(r3)]

if (yplot eq 'on') then begin
    colors = intarr(n0)
    for j=0,n0-1 do colors[j] = round((200./n0)*j)
    set_plot,'x'
    window,2,retain=2
    device,decomposed=0
    loadct,0
    plot,pwave,alldata[*,0],title='FWHM for file '+file,xtitle='wavelength (A)',$
      ytitle='FWHM (A) or I.DISP (km/sec)',/nodata,xrange=[4000,5500],xstyle=1,$
      yrange=[min(alldata)-0.1,max(alldata)+0.1],ystyle=1
    oplot,pwave,allmean,psym=2,symsize=2.5,thick=2
    xyouts,0.15,0.80,'MEAN',charsize=1.5,/normal
    loadct,27
    oplot,pwave,allmedian,psym=2,symsize=2.5,thick=2,color=60
    xyouts,0.15,0.85,'MEDIAN',color=60,charsize=1.5,/normal
    for j=0,n0-1 do begin
        oplot,pwave,alldata[*,j],psym=1,color=colors[j]
    endfor
    print,'Next ENTER deletes the plot...'
    pause
    wdelete
endif


stop
END

PRO Pchi, file, SORTON=sorton, OUTNAME=outname
;this routine is used to sort the results of a modeling run. Two
;different plots are generated. Plot 1 is just a X^2 vs. parameter
;plot for the 4 parameters. Plot 2 is a parameter vs. parameter plot
;where the 68% and 95% values are centered upon and plotted.

; After the cres files are catted together, these are sorted in order
; from best to worst chi-squared. This is output as chi2.std.txt

; FILE: The name of the output from the catted-together CRES files.

; SORTON: This is string of either 'stars', 'gc' or 'both' and
; determines which of the three X^2 values gets sorted on. If the
; keyword isn't used then the user is prompted.

;This sets whether the plots are saved as postscript files.
pswitch = 'on'
;pswitch = 'off'

readcol,file,format='f,f,i,i,f,f,f',ml,bhm,vc,rs,x2stars,x2gc,x2tot

if (n_elements(sorton) eq 0) then begin
    sorton = ''
    print,'Enter which chi-squared value to sort on: (stars/gc/both)'
    read,sorton
endif

if (sorton eq 'stars') then begin
    sindex = bsort(x2stars)
    sx2stars = x2stars(sindex)
    temp = sx2stars[1:*]
    x2dn = min(x2stars)-100
    x2up = max(x2stars)+30
    put = sx2stars
endif

if (sorton eq 'gc') then begin
    sindex = bsort(x2gc)
    sx2gc = x2gc(sindex)
    temp = sx2gc[1:*]
    x2dn = min(x2gc)-10
    x2up = max(x2gc)+10
    put = sx2gc
endif

if (sorton eq 'both') then begin
    sindex = bsort(x2tot)
    sx2tot = x2tot(sindex)
    temp = sx2tot[1:*]
    x2dn = min(x2tot)-100
    x2up = max(x2tot)+30
    put = sx2tot
endif

n1 = n_elements(ml)

sml = ml(sindex)
sbhm = bhm(sindex)
svc = vc(sindex)
srs = rs(sindex)
sx2stars = x2stars(sindex)
sx2gc = x2gc(sindex)
sx2tot = x2tot(sindex)

f1 = '(f5.2,2x,e12.5,2x,i4,2x,i3,2x,f9.3,2x,f9.3,2x,f9.3)'

if(n_elements(outname) eq 0) then outname = file+'.'+sorton+'.chi2.std.txt'
openw,5,outname
for j=0,n1-1 do printf,5,sml[j],sbhm[j],svc[j],srs[j],$
  sx2stars[j],sx2gc[j],sx2tot[j],format=f1
free_lun,5

;the breaks between the 68 and 95th percent levels are determined

x268 = where(temp-temp[0] le 1.0)
x295 = where(temp-temp[0] le 2.0)

mlup = max(ml)+0.5
mldn = min(ml)-0.5
bhmup = max(bhm)+1e9
bhmdn = min(bhm)-1e9
vcup = max(vc)+50
vcdn = min(vc)-50
rsup = max(rs)+1
rsdn = min(rs)-1

set_plot,'x'
ans = ''

repeat begin
    window,2,retain=2,xsize=800,ysize=700
    device,decomposed=0
    
    !p.multi = [0,2,2,0,1]
    !y.omargin = [2,4]
    loadct,0
    plot,sml,put,psym=sym(1),symsize=0.5,xtitle='M/L',ytitle='!7v!6!u2!6',$
      yrange=[x2dn,x2up],ystyle=1,xrange=[mldn,mlup],xstyle=1,charsize=1.2
   
    loadct,0
    plot,sbhm,put,psym=sym(1),symsize=0.5,xtitle='Black Hole Mass',ytitle='!7v!6!u2!6',$
      yrange=[x2dn,x2up],ystyle=1,xrange=[bhmdn,bhmup],xstyle=1,charsize=1.2
    
    plot,svc,put,psym=sym(1),symsize=0.5,xtitle='V!dc!n(km/sec)',ytitle='!7v!6!u2!6',$
      yrange=[x2dn,x2up],ystyle=1,xrange=[vcdn,vcup],xstyle=1,charsize=1.2
    
    plot,srs,put,psym=sym(1),symsize=0.5,xtitle='R!dc!n(kpc)',ytitle='!7v!6!u2!6',$
      yrange=[x2dn,x2up],ystyle=1,xrange=[rsdn,rsup],xstyle=1,charsize=1.2
    xyouts,0.5,0.95,/normal,'Chi squared values for '+sorton+': file '+file,alignment=0.5,charsize=2.0
    
    print,'Change the Chi-squared range? ("y" or "n")'
    read,ans
    if (ans eq 'y') then begin
        print,'The current upper limit is '+strn(x2up)
        print,'Enter a new upper limit.'
        read,x2up
        print,'The current lower limit is '+strn(x2dn)
        print,'Enter a new lower limit.'
        read,x2dn
        if (sorton eq 'stars') then i = where(x2stars lt x2up)
        if (sorton eq 'gc') then i = where(x2gc lt x2up)
        if (sorton eq 'both') then i = where(x2tot lt x2up)
        mldn = min(ml[i])-0.5
        mlup = max(ml[i])+0.5
        bhmdn = min(bhm[i])-0.5e9
        bhmup = max(bhm[i])+0.5e9
        vcdn = min(vc[i])-50
        vcup = max(vc[i])+50
        rsdn = min(rs[i])-1
        rsup = max(rs[i])+1
    endif
    
endrep until (ans eq 'n')

if (pswitch eq 'on') then begin
    set_plot,'ps'
    device,file=file+'.'+sorton+'.plotchi1.ps',/color

    !p.multi = [0,2,2,0,1]
    !y.omargin = [2,4]
    loadct,0
    plot,sml,put,psym=sym(1),symsize=0.2,xtitle='M/L',ytitle='!7v!6!u2!6',$
      yrange=[x2dn,x2up],ystyle=1,xthick=3,ythick=3,charthick=2,charsize=0.8,$
      xrange=[mldn,mlup],xstyle=1
    
    plot,sbhm,put,psym=sym(1),symsize=0.2,xtitle='Black Hole Mass',ytitle='!7v!6!u2!6',$
      yrange=[x2dn,x2up],ystyle=1,xthick=3,ythick=3,charthick=2,charsize=0.8,$
     xrange=[bhmdn,bhmup],xstyle=1
    
    plot,svc,put,psym=sym(1),symsize=0.2,xtitle='V!dc!n(km/sec)',ytitle='!7v!6!u2!6',$
      yrange=[x2dn,x2up],ystyle=1,xthick=3,ythick=3,charthick=2,charsize=0.8,$
      xrange=[vcdn,vcup],xstyle=1
    
    plot,srs,put,psym=sym(1),symsize=0.2,xtitle='R!dc!n(kpc)',ytitle='!7v!6!u2!6',$
      yrange=[x2dn,x2up],ystyle=1,xthick=3,ythick=3,charthick=2,charsize=0.8,$
      xrange=[rsdn,rsup],xstyle=1
    xyouts,0.5,0.9,/normal,'Chi squared values for '+sorton+'. File name: '+file,$
      alignment=0.5,charsize=1.3,charthick=2,font=0
    
    device,/close_file
    set_plot,'x'
endif

window,0,retain=2,xsize=800,ysize=700
device,decomposed=0
loadct,0

xyouts,1/6.,5/6.,'M/L',/normal,charsize=4.0,alignment=0.5
xyouts,3/6.,3/6.,'V!dcirc!n',/normal,charsize=4.0,alignment=0.5
xyouts,5/6.,1/6.,'R!dscale!n',/normal,charsize=4.0,alignment=0.5

!p.multi = [8,3,3,0,0]
loadct,0
plot,svc,sml,psym=3,yrange=[mldn,mlup],xrange=[vcdn,vcup],xstyle=1,$
  ystyle=1,charsize=2.0
oplot,svc[x295],sml[x295],psym=sym(1),symsize=0.8
loadct,4
oplot,svc[x268],sml[x268],psym=sym(1),color=150,symsize=1.0

!p.multi = [7,3,3,0,0]
loadct,0
plot,srs,sml,psym=3,yrange=[mldn,mlup],xrange=[rsdn,rsup],xstyle=1,$
  ystyle=1,charsize=2.0
oplot,srs[x295],sml[x295],psym=sym(1),symsize=0.8
loadct,4
oplot,srs[x268],sml[x268],psym=sym(1),color=150,symsize=1.0

!p.multi = [6,3,3,0,0]
loadct,0
plot,sml,svc,psym=3,xrange=[mldn,mlup],yrange=[vcdn,vcup],xstyle=1,$
  ystyle=1,charsize=2.0
oplot,sml[x295],svc[x295],psym=sym(1),symsize=0.8
loadct,4
oplot,sml[x268],svc[x268],psym=sym(1),color=150,symsize=1.0

!p.multi = [4,3,3,0,0]
loadct,0
plot,srs,svc,psym=3,xrange=[rsdn,rsup],yrange=[vcdn,vcup],xstyle=1,$
  ystyle=1,charsize=2.0
oplot,srs[x295],svc[x295],psym=sym(1),symsize=0.8
loadct,4
oplot,srs[x268],svc[x268],psym=sym(1),color=150,symsize=1.0

!p.multi = [3,3,3,0,0]
loadct,0
plot,sml,srs,psym=3,xrange=[mldn,mlup],yrange=[rsdn,rsup],xstyle=1,$
  ystyle=1,charsize=2.0
oplot,sml[x295],srs[x295],psym=sym(1),symsize=0.8
loadct,4
oplot,sml[x268],srs[x268],psym=sym(1),color=150,symsize=1.0

!p.multi = [2,3,3,0,0]
loadct,0
plot,svc,srs,psym=3,xrange=[vcdn,vcup],yrange=[rsdn,rsup],xstyle=1,$
  ystyle=1,charsize=2.0
oplot,svc[x295],srs[x295],psym=sym(1),symsize=0.8
loadct,4
oplot,svc[x268],srs[x268],psym=sym(1),color=150,symsize=1.0

if (pswitch eq 'on') then begin
    set_plot,'ps'
    device,file=file+'.'+sorton+'.plotchi2.ps',/color
    loadct,0
    xyouts,1/6.,5/6.,'M/L',/normal,charsize=3.0,alignment=0.5,charthick=2
    xyouts,3/6.,3/6.,'V!dcirc!n',/normal,charsize=3.0,alignment=0.5,charthick=2
    xyouts,5/6.,1/6.,'R!dscale!n',/normal,charsize=3.0,alignment=0.5,charthick=2
    
    !p.multi = [8,3,3,0,0]
    plot,svc,sml,psym=3,yrange=[mldn,mlup],xrange=[vcdn,vcup],xstyle=1,$
      ystyle=1,charsize=0.8,xthick=3,ythick=3,charthick=2
    oplot,svc[x295],sml[x295],psym=sym(1),symsize=0.4
    loadct,4
    oplot,svc[x268],sml[x268],psym=sym(1),color=150,symsize=0.6
    
    !p.multi = [7,3,3,0,0]
    loadct,0
    plot,srs,sml,psym=3,yrange=[mldn,mlup],xrange=[rsdn,rsup],xstyle=1,$
      ystyle=1,charsize=0.8,xthick=3,ythick=3,charthick=2
    oplot,srs[x295],sml[x295],psym=sym(1),symsize=0.4
    loadct,4
    oplot,srs[x268],sml[x268],psym=sym(1),color=150,symsize=0.6
    
    !p.multi = [6,3,3,0,0]
    loadct,0
    plot,sml,svc,psym=3,xrange=[mldn,mlup],yrange=[vcdn,vcup],xstyle=1,$
      ystyle=1,charsize=0.8,xthick=3,ythick=3,charthick=2
    oplot,sml[x295],svc[x295],psym=sym(1),symsize=0.4
    loadct,4
    oplot,sml[x268],svc[x268],psym=sym(1),color=150,symsize=0.6
    
    !p.multi = [4,3,3,0,0]
    loadct,0
    plot,srs,svc,psym=3,xrange=[rsdn,rsup],yrange=[vcdn,vcup],xstyle=1,$
      ystyle=1,charsize=0.8,xthick=3,ythick=3,charthick=2
    oplot,srs[x295],svc[x295],psym=sym(1),symsize=0.4
    loadct,4
    oplot,srs[x268],svc[x268],psym=sym(1),color=150,symsize=0.6
    
    !p.multi = [3,3,3,0,0]
    loadct,0
    plot,sml,srs,psym=3,xrange=[mldn,mlup],yrange=[rsdn,rsup],xstyle=1,$
      ystyle=1,charsize=0.8,xthick=3,ythick=3,charthick=2
    oplot,sml[x295],srs[x295],psym=sym(1),symsize=0.4
    loadct,4
    oplot,sml[x268],srs[x268],psym=sym(1),color=150,symsize=0.6
    
    !p.multi = [2,3,3,0,0]
    loadct,0
    plot,svc,srs,psym=3,xrange=[vcdn,vcup],yrange=[rsdn,rsup],xstyle=1,$
      ystyle=1,charsize=0.8,xthick=3,ythick=3,charthick=2
    oplot,svc[x295],srs[x295],psym=sym(1),symsize=0.4
    loadct,4
    oplot,svc[x268],srs[x268],psym=sym(1),color=150,symsize=0.6
    
    device,/close_file
    set_plot,'x'

endif

!p.multi = [0]

print,'Next ENTER deletes the plots...'
pause
wdelete,0,2

STOP
END

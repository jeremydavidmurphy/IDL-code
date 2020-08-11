; I have a better looking version of this code (meaning, the output
; figures for chi2.ps are better) in the figures directory of my first
; paper on M87.

PRO Pchi, file, SORTON=sorton, OUTNAME=outname, TYPE=type
;this routine is used to sort the results of a modeling run. Two
;different plots are generated. Plot 1 is just a X^2 vs. parameter
;plot for the 3 or 4 parameters. Plot 2 is a parameter vs. parameter plot
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
;offset = 18.9
offset = 0.0 ;a linear shift applied to the GC values to align with the stellar values.
chistep = 3.0 ;change this if you want to capture a greater deltaX2 range. If set to 2.0, then the code will sort on deltaX2 = 2.0 and 4.0; on 3, then 3.0 and 6.0
alphacorrection = 'on'
STARcutoff = 0.1 ; Set this (and the next) value to where you want to set the cutoff for alpha-correction.
GCcutoff = 0.1   ; Any alpha values BELOW this number are rejected and NO CORRECTION TO THE X2 VALUE IS MADE.
;*********************************************************

if (alphacorrection eq 'on') then begin
readcol,file,format='f,f,f,f,d,d,x,d,d',ml,bhm,vc,rs,x2stars,x2gc,astars,agc
; This new component now makes the corrections to chi2
; based on the convergence of the model (based on the alpha
; parameter). The first column is the alpha parameter from the
; models. The second is a X^2 adjustment that's applied to the
; modeling results. THESE X^2 VALUES GET ADDED TO THE X^2 VALUES
; OUTPUT FROM THE MODELS. 

    readcol,'plotiterSTARS.out',format='d,d',alphaCs,chiCs
    readcol,'plotiterGC.out',format='d,d',alphaCg,chiCg

    x2tot = dblarr(n_elements(astars))
    for j=0,n_elements(astars)-1 do begin
        if (astars[j] gt STARcutoff) then begin
            adjuststars = interpol(chiCs,alphaCs,astars[j],/spline)
            if (adjuststars gt 5.0) then print,adjuststars
            x2stars[j] = x2stars[j] - adjuststars
;            print,'Adjusting frame '+strn(j+1)+' by '+strn(adjuststars)+' in STARS...'
        endif
        if (agc[j] gt GCcutoff) then begin
            adjustgc = interpol(chiCg,alphaCg,agc[j],/spline)
            if (adjustgc gt 5.0) then print,adjustgc
            x2gc[j] = x2gc[j] - adjustgc
;            print,'Adjusting frame '+strn(j+1)+' by '+strn(adjustgc)+' in GC...'
        endif
        x2tot[j] = x2stars[j] + x2gc[j]
    endfor
    x2gc = x2gc + offset
endif else begin
    readcol,file,format='f,f,f,f,d,d,x,d,d',ml,bhm,vc,rs,x2stars,x2gc,astars,agc
    x2tot = x2stars + x2gc
endelse

;*********************************************************
if (n_elements(sorton) eq 0) then begin
    sorton = ''
    print,'Enter which chi-squared value to sort on: (stars/gc/both)'
    read,sorton
endif
if (n_elements(type) eq 0) then begin
    type = ''
    print,'Enter the profile type (nfw or log):'
    read,type
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
sastars = astars[sindex]
sagc = agc[sindex]

f1 = '(f5.2,2x,e12.5,2x,f7.2,2x,f7.2,2x,f11.4,2x,f11.4,2x,f11.4,2x,f9.6,2x,f9.6)'
f2 = '(f5.2,2x,e12.5,2x,f7.2,2x,f7.2,2x,f11.4,2x,f11.4,2x,f11.4)'

if(n_elements(outname) eq 0) then outname = file+'.'+sorton+'.std.txt'
openw,5,outname
for j=0,n1-1 do printf,5,sml[j],sbhm[j],svc[j],srs[j],$
  sx2stars[j],sx2gc[j],sx2tot[j],sastars[j],sagc[j],format=f1
;for j=0,n1-1 do printf,5,sml[j],sbhm[j],svc[j],srs[j],$
;  sx2stars[j],sx2gc[j],sx2tot[j],format=f2
free_lun,5

;the breaks between the 68 and 95th percent levels are determined

x268 = where(temp-temp[0] le chistep * 1.0)
x295 = where(temp-temp[0] le chistep * 2.0)

mlup = max(ml)+0.5
mldn = min(ml)-0.5
bhmup = max(bhm)+1e9
bhmdn = min(bhm)-1e9
if type eq 'log' then vcup = max(vc)+50
if type eq 'log' then vcdn = min(vc)-50
if type eq 'nfw' then vcup = max(vc)+1 
if type eq 'nfw' then vcdn = min(vc)-1
rsup = max(rs)+5
rsdn = min(rs)-5

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
    
    if type eq 'log' then begin
        plot,svc,put,psym=sym(1),symsize=0.5,xtitle='V!dc!n(km/sec)',ytitle='!7v!6!u2!6',$
          yrange=[x2dn,x2up],ystyle=1,xrange=[vcdn,vcup],xstyle=1,charsize=1.2
    endif
    if type eq 'nfw' then begin
        plot,svc,put,psym=sym(1),symsize=0.5,xtitle='C',ytitle='!7v!6!u2!6',$
          yrange=[x2dn,x2up],ystyle=1,xrange=[vcdn,vcup],xstyle=1,charsize=1.2
    endif

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
        if type eq 'log' then vcup = max(vc[i])+50
        if type eq 'log' then vcdn = min(vc[i])-50
        if type eq 'nfw' then vcup = max(vc[i])+1 
        if type eq 'nfw' then vcdn = min(vc[i])-1
        rsdn = min(rs[i])-5
        rsup = max(rs[i])+5
    endif
    
endrep until (ans eq 'n')

if (pswitch eq 'on') then begin
    set_plot,'ps'
    device,file=file+'.'+sorton+'.plotchi1.ps',/color
    !p.multi = [0,2,2,0,1]
    !y.omargin = [2,4]

    loadct,0
    plot,sml,put,psym=sym(1),symsize=0.4,xtitle='M/L',ytitle='!7v!6!u2!6',$
      yrange=[x2dn,x2up],ystyle=1,xthick=3,ythick=3,charthick=2,charsize=0.8,$
      xrange=[mldn,mlup],xstyle=1
    
    plot,sbhm,put,psym=sym(1),symsize=0.4,xtitle='Black Hole Mass',ytitle='!7v!6!u2!6',$
      yrange=[x2dn,x2up],ystyle=1,xthick=3,ythick=3,charthick=2,charsize=0.8,$
      xrange=[bhmdn,bhmup],xstyle=1
    
    if type eq 'log' then begin
        plot,svc,put,psym=sym(1),symsize=0.4,xtitle='V!dc!n(km/sec)',ytitle='!7v!6!u2!6',$
          yrange=[x2dn,x2up],ystyle=1,xthick=3,ythick=3,charthick=2,charsize=0.8,$
          xrange=[vcdn,vcup],xstyle=1
    endif
    if type eq 'nfw' then begin
        plot,svc,put,psym=sym(1),symsize=0.4,xtitle='Concentration',ytitle='!7v!6!u2!6',$
          yrange=[x2dn,x2up],ystyle=1,xthick=3,ythick=3,charthick=2,charsize=0.8,$
          xrange=[vcdn,vcup],xstyle=1
    endif
    plot,srs,put,psym=sym(1),symsize=0.4,xtitle='R!dc!n(kpc)',ytitle='!7v!6!u2!6',$
      yrange=[x2dn,x2up],ystyle=1,xthick=3,ythick=3,charthick=2,charsize=0.8,$
      xrange=[rsdn,rsup],xstyle=1
    xyouts,0.5,0.9,/normal,'Chi squared values for '+sorton+'.'+type+'. File name: '+file,$
      alignment=0.5,charsize=1.3,charthick=2,font=0
    device,/close_file
    set_plot,'x'
endif

window,0,retain=2,xsize=800,ysize=700
device,decomposed=0
loadct,0

xyouts,1/6.,5/6.,'M/L',/normal,charsize=4.0,alignment=0.5
if type eq 'log' then xyouts,3/6.,3/6.,'V!dcirc!n',/normal,charsize=4.0,alignment=0.5
if type eq 'nfw' then xyouts,3/6.,3/6.,'C',/normal,charsize=4.0,alignment=0.5
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
    if type eq 'log' then xyouts,3/6.,3/6.,'V!dcirc!n',/normal,charsize=3.0,alignment=0.5,charthick=2
    if type eq 'nfw' then xyouts,3/6.,3/6.,'C',/normal,charsize=3.0,alignment=0.5,charthick=2
    xyouts,5/6.,1/6.,'R!dscale!n',/normal,charsize=3.0,alignment=0.5,charthick=2
    
    !p.multi = [8,3,3,0,0]
    plot,svc,sml,psym=3,yrange=[mldn,mlup],xrange=[vcdn,vcup],xstyle=1,$
      ystyle=1,charsize=1.0,xthick=3,ythick=3,charthick=2
    oplot,svc[x295],sml[x295],psym=sym(1),symsize=0.4
    loadct,4
    oplot,svc[x268],sml[x268],psym=sym(1),color=150,symsize=0.6
    
    !p.multi = [7,3,3,0,0]
    loadct,0
    plot,srs,sml,psym=3,yrange=[mldn,mlup],xrange=[rsdn,rsup],xstyle=1,$
      ystyle=1,charsize=1.0,xthick=3,ythick=3,charthick=2
    oplot,srs[x295],sml[x295],psym=sym(1),symsize=0.4
    loadct,4
    oplot,srs[x268],sml[x268],psym=sym(1),color=150,symsize=0.6
    
    !p.multi = [6,3,3,0,0]
    loadct,0
    plot,sml,svc,psym=3,xrange=[mldn,mlup],yrange=[vcdn,vcup],xstyle=1,$
      ystyle=1,charsize=1.0,xthick=3,ythick=3,charthick=2
    oplot,sml[x295],svc[x295],psym=sym(1),symsize=0.4
    loadct,4
    oplot,sml[x268],svc[x268],psym=sym(1),color=150,symsize=0.6
    
    !p.multi = [4,3,3,0,0]
    loadct,0
    plot,srs,svc,psym=3,xrange=[rsdn,rsup],yrange=[vcdn,vcup],xstyle=1,$
      ystyle=1,charsize=1.0,xthick=3,ythick=3,charthick=2
    oplot,srs[x295],svc[x295],psym=sym(1),symsize=0.4
    loadct,4
    oplot,srs[x268],svc[x268],psym=sym(1),color=150,symsize=0.6
    
    !p.multi = [3,3,3,0,0]
    loadct,0
    plot,sml,srs,psym=3,xrange=[mldn,mlup],yrange=[rsdn,rsup],xstyle=1,$
      ystyle=1,charsize=1.0,xthick=3,ythick=3,charthick=2
    oplot,sml[x295],srs[x295],psym=sym(1),symsize=0.4
    loadct,4
    oplot,sml[x268],srs[x268],psym=sym(1),color=150,symsize=0.6
    
    !p.multi = [2,3,3,0,0]
    loadct,0
    plot,svc,srs,psym=3,xrange=[vcdn,vcup],yrange=[rsdn,rsup],xstyle=1,$
      ystyle=1,charsize=1.0,xthick=3,ythick=3,charthick=2
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

; Version 2 of this code attempts to add lines showing the X2 minimums
; from the stars and GC. Because of this, the only sorting is on BOTH.

PRO Pchi_v2, file,OUTNAME=outname, TYPE=type

;compile findx2floorf.pro

;this routine is used to sort the results of a modeling run. Two
;different plots are generated. Plot 1 is just a X^2 vs. parameter
;plot for the 3 or 4 parameters. Plot 2 is a parameter vs. parameter plot
;where the 68% and 95% values are centered upon and plotted.

; After the cres files are catted together, these are sorted in order
; from best to worst chi-squared. This is output as chi2.std.txt

; FILE: The name of the output from the catted-together CRES files.

;This sets whether the plots are saved as postscript files.
pswitch = 'on'
;pswitch = 'off'
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
endif else $
readcol,file,format='f,f,f,f,d,d,d,d,d',ml,bhm,vc,rs,x2stars,x2gc,x2tot,astars,agc

;*********************************************************

if (n_elements(type) eq 0) then begin
    type = ''
    print,'Enter the profile type (nfw or log):'
    read,type
endif

sindex = bsort(x2tot)
sx2tot = x2tot(sindex)
temp = sx2tot[1:*]
x2dn = min(x2tot)-100
x2up = max(x2tot)+30
put = sx2tot

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

if(n_elements(outname) eq 0) then outname = file+'.V2.std.txt'
openw,5,outname
for j=0,n1-1 do printf,5,sml[j],sbhm[j],svc[j],srs[j],$
  sx2stars[j],sx2gc[j],sx2tot[j],sastars[j],sagc[j],format=f1
;for j=0,n1-1 do printf,5,sml[j],sbhm[j],svc[j],srs[j],$
;  sx2stars[j],sx2gc[j],sx2tot[j],format=f2
free_lun,5

;the breaks between the 68 and 95th percent levels are determined

x268 = where(temp-temp[0] le chistep * 1.0)
x295 = where(temp-temp[0] le chistep * 2.0)

mlup = max(ml)+0.2
mldn = min(ml)-0.2
bhmup = max(bhm)+1e9
bhmdn = min(bhm)-1e9
if type eq 'log' then vcup = max(vc)+50
if type eq 'log' then vcdn = min(vc)-50
if type eq 'nfw' then vcup = max(vc)+1 
if type eq 'nfw' then vcdn = min(vc)-1
rsup = max(rs)+5
rsdn = min(rs)-5

pmlgc = findx2floorf(ml,x2gc)
pmlst = findx2floorf(ml,x2stars)
pmlbt = findx2floorf(ml,x2tot)
stop
pvcgc = findx2floorf(vc,x2gc)
pvcst = findx2floorf(vc,x2stars)
pvcbt = findx2floorf(vc,x2tot)
prsgc = findx2floorf(rs,x2gc)
prsst = findx2floorf(rs,x2stars)
prsbt = findx2floorf(rs,x2tot)
offset1gc = median(pmlbt[*,1]) - median(pmlgc[*,1])
offset2gc = median(pvcbt[*,1]) - median(pvcgc[*,1])
offset3gc = median(prsbt[*,1]) - median(prsgc[*,1])
offset1st = median(pmlbt[*,1]) - median(pmlst[*,1])
offset2st = median(pvcbt[*,1]) - median(pvcst[*,1])
offset3st = median(prsbt[*,1]) - median(prsst[*,1])

set_plot,'x'
ans = ''

repeat begin
    window,2,retain=2,xsize=800,ysize=700
    device,decomposed=0
    
    !p.multi = [0,2,2,0,1]
    !y.omargin = [2,4]
    loadct,0,/silent
    plot,sml,put,psym=sym(1),symsize=0.5,xtitle='M/L',ytitle='!7v!6!u2!6',$
      yrange=[x2dn,x2up],ystyle=1,xrange=[mldn,mlup],xstyle=1,charsize=1.2
    loadct,4,/silent
    oplot,pmlgc[*,0],pmlgc[*,1]+offset1gc,color=110,thick=2
    oplot,pmlst[*,0],pmlst[*,1]+offset1st,color=150,thick=2
    loadct,0,/silent

    plot,sbhm,put,psym=sym(1),symsize=0.5,xtitle='Black Hole Mass',ytitle='!7v!6!u2!6',$
      yrange=[x2dn,x2up],ystyle=1,xrange=[bhmdn,bhmup],xstyle=1,charsize=1.2
    
    if type eq 'log' then begin
        plot,svc,put,psym=sym(1),symsize=0.5,xtitle='V!dc!n(km/sec)',ytitle='!7v!6!u2!6',$
          yrange=[x2dn,x2up],ystyle=1,xrange=[vcdn,vcup],xstyle=1,charsize=1.2
        loadct,4,/silent
        oplot,pvcgc[*,0],pvcgc[*,1]+offset2gc,color=110,thick=2
        oplot,pvcst[*,0],pvcst[*,1]+offset2st,color=150,thick=2
        loadct,0,/silent
    endif
    if type eq 'nfw' then begin
        plot,svc,put,psym=sym(1),symsize=0.5,xtitle='C',ytitle='!7v!6!u2!6',$
          yrange=[x2dn,x2up],ystyle=1,xrange=[vcdn,vcup],xstyle=1,charsize=1.2
        loadct,4,/silent
        oplot,pvcgc[*,0],pvcgc[*,1]+offset2gc,color=110,thick=2
        oplot,pvcst[*,0],pvcst[*,1]+offset2st,color=150,thick=2
        loadct,0,/silent
    endif

    plot,srs,put,psym=sym(1),symsize=0.5,xtitle='R!dc!n(kpc)',ytitle='!7v!6!u2!6',$
      yrange=[x2dn,x2up],ystyle=1,xrange=[rsdn,rsup],xstyle=1,charsize=1.2
    loadct,4,/silent
    oplot,prsgc[*,0],prsgc[*,1]+offset3gc,color=110,thick=2
    oplot,prsst[*,0],prsst[*,1]+offset3st,color=150,thick=2
    loadct,0,/silent

    xyouts,0.5,0.95,/normal,'Chi squared values for the V2: file '+file,alignment=0.5,charsize=2.0
    
    print,'Change the Chi-squared range? ("y" or "n")'
    read,ans
    if (ans eq 'y') then begin
        print,'The current upper limit is '+strn(x2up)
        print,'Enter a new upper limit.'
        read,x2up
        print,'The current lower limit is '+strn(x2dn)
        print,'Enter a new lower limit.'
        read,x2dn
        i = where(x2tot lt x2up)
        mldn = min(ml[i])-0.2
        mlup = max(ml[i])+0.2
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
    device,file=file+'.V2.plotchi1.ps',/color
    !p.multi = [0,2,2,0,1]
    !y.omargin = [2,4]

    loadct,0,/silent
    plot,sml,put,psym=sym(1),symsize=0.4,xtitle='M/L',ytitle='!7v!6!u2!6',$
      yrange=[x2dn,x2up],ystyle=1,xthick=3,ythick=3,charthick=2,charsize=0.8,$
      xrange=[mldn,mlup],xstyle=1
    loadct,4,/silent
    oplot,pmlgc[*,0],pmlgc[*,1]+offset1gc,color=110,thick=2
    oplot,pmlst[*,0],pmlst[*,1]+offset1st,color=150,thick=2
    loadct,0,/silent

    plot,sbhm,put,psym=sym(1),symsize=0.4,xtitle='Black Hole Mass',ytitle='!7v!6!u2!6',$
      yrange=[x2dn,x2up],ystyle=1,xthick=3,ythick=3,charthick=2,charsize=0.8,$
      xrange=[bhmdn,bhmup],xstyle=1
    
    if type eq 'log' then begin
        plot,svc,put,psym=sym(1),symsize=0.4,xtitle='V!dc!n(km/sec)',ytitle='!7v!6!u2!6',$
          yrange=[x2dn,x2up],ystyle=1,xthick=3,ythick=3,charthick=2,charsize=0.8,$
          xrange=[vcdn,vcup],xstyle=1
        loadct,4,/silent
        oplot,pvcgc[*,0],pvcgc[*,1]+offset2gc,color=110,thick=2
        oplot,pvcst[*,0],pvcst[*,1]+offset2st,color=150,thick=2
        loadct,0,/silent
    endif
    if type eq 'nfw' then begin
        plot,svc,put,psym=sym(1),symsize=0.4,xtitle='Concentration',ytitle='!7v!6!u2!6',$
          yrange=[x2dn,x2up],ystyle=1,xthick=3,ythick=3,charthick=2,charsize=0.8,$
          xrange=[vcdn,vcup],xstyle=1
        loadct,4,/silent
        oplot,pvcgc[*,0],pvcgc[*,1]+offset2gc,color=110,thick=2
        oplot,pvcst[*,0],pvcst[*,1]+offset2st,color=150,thick=2
        loadct,0,/silent
    endif

    plot,srs,put,psym=sym(1),symsize=0.4,xtitle='R!dc!n(kpc)',ytitle='!7v!6!u2!6',$
      yrange=[x2dn,x2up],ystyle=1,xthick=3,ythick=3,charthick=2,charsize=0.8,$
      xrange=[rsdn,rsup],xstyle=1
    loadct,4,/silent
    oplot,prsgc[*,0],prsgc[*,1]+offset3gc,color=110,thick=2
    oplot,prsst[*,0],prsst[*,1]+offset3st,color=150,thick=2
    loadct,0,/silent

    xyouts,0.5,0.9,/normal,'Chi squared values for V2.'+type+'. File name: '+file,$
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
    device,file=file+'.V2.plotchi2.ps',/color
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

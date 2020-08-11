pro fixweakfib, frame1, frame2, frame3

;window, 0, retain=2, xsize=500, ysize=500
;loadct, 39
;device, decompose=0

print, frame1

f1=mrdfits(frame1+'pe.fits', 0, h)
f2=mrdfits(frame2+'pe.fits', 0, h)
f3=mrdfits(frame3+'pe.fits', 0, h)
flat=mrdfits(frame1+'pef.fits', 0, h)

f=dblarr(1024,2048)
for i=0, 1024-1 do begin
    for j=0, 2048-1 do begin
        f[i,j]=median([f1[i,j],f2[i,j],f3[i,j]])
    endfor
endfor

f=f/flat
f=f/median(f)

readcol, 'ptow_Mar09_n1.dat', a0, a1, a2, a3, a4
lambda=dblarr(1024, 246)
x=dindgen(1024)
for i=0, 246-1 do lambda[*,i]=a0[i]+a1[i]*x+a2[i]*x^2+a3[i]*x^3+a4[i]*x^4

fc=dblarr(1024, 246)
for i=0, 246-1 do begin
    ymin=(i+1)*8-2-1
    ymax=(i+1)*8+2-1
    for j=0, 1024-1 do fc[j,i]=median(f[j,ymin:ymax])
endfor

A5577=dblarr(246)

for i=0, 246-1 do begin
    auxl0=lambda[*,i]
    auxf0=fc[*,i]
    sel=where(auxl0 ge 5555 and auxl0 le 5600)
    auxl=auxl0[sel]
    auxf=auxf0[sel]
    a=[max(auxf), 5577.0, 1.0, auxf[0]]
    fit=gaussfit(auxl, auxf, a, nterms=4)
;   print, a
    A5577[i]=a[0]*sqrt(2*!pi*a[2]^2)
 ;   plot, auxl, auxf
 ;   oplot, auxl, fit, color=150
;if (i eq 42) then stop
endfor



flatc=dblarr(1024, 246)
flatcc=dblarr(246)
for i=0, 246-1 do begin
    ymin=(i+1)*8-2-1
    ymax=(i+1)*8+2-1
    for j=0, 1024-1 do begin
        aux0=flat[j,ymin:ymax]
        sel=where(aux0 ne -666)
        if (sel eq [-1]) then begin
            flatc[j,i]=-666
        endif else begin
        aux1=aux0[sel]
        flatc[j,i]=mean(aux1)
        endelse
    endfor
    aux2=flatc[*,i]
    flatcc[i]=mean(aux2[where(aux2 ne -666)])
endfor


;window, 1, retain=2, xsize=1800, ysize=500
;plot, indgen(246), a5577/median(a5577), psym=-1
;oplot, indgen(246), flatcc/median(flatcc), color=250, psym=-2

corr=a5577/median(a5577)

fcorr=dblarr(1024,2048)
for i=0, 246-1 do begin 
    ymin=(i+1)*8-2-1
    ymax=(i+1)*8+2-1
    fcorr[*,ymin:ymax]=corr[i]
endfor

sky=mrdfits(frame1+'pefy.fits', 0, h)

badframe=mrdfits(frame1+'pefs.fits', 0, h) ; for copying header

fixed=dblarr(1024,2048)
sel=where(flat ne 0)
fixed[sel]=f1[sel]/flat[sel]/fcorr[sel]-sky[sel]
mwrfits, fixed , frame1+'pefsf.fits', h, /create



;stop
end

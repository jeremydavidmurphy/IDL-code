; This code is designed to read in two lists, containing the dM.dat
; and dMhalo.dat output of the dynamical models, over the uncertainty
; range of the modeling. It then output a text file of the form

; R, light, tot_mass, tot_U, tot_L, DM_mass, DM_U, DM_L

; the "light" will need to get multipled by M/L to get a mass out.

; THE FIRST MODELS IN THE LIST ARE ASSUMED TO BE THE BEST FIT. IF THEY
; ARE NOT, THEN THE COLUMNS WITH "TOTALS" DON'T MEAN ANYTHING, AND
; ONLY THE UNCERTAINTY SHOULD BE USED.

; This needs the dL.dat, and both bindemo_r.out and bindemo_v.out in
; the calling directory

; NOTE: You need to put in the proper black hole mass into this before
; running it...

PRO error_calc

;*********************************************************
angleRmax = 2000.0
BHlower = 5.9e9
BHmass = 6.4e9
BHupper = 6.9e9
;*********************************************************

readcol,'bindemo_r.out',silent=1,format='i,x,d,d,d,x,x',$
  skipline=4,rad,radl,radm,radh
;the initial zeros are dumped as they don't make it into dL.dat,
;dM.dat, etc
rad = rad[1:*]
radl = radl[1:*]
radm = radm[1:*]
radh = radh[1:*]
readcol,'bindemo_v.out',silent=1,format='i,x,d,d,d',$
  skipline=4,angle,angl,angm,angh


readcol,'dMhalo.list',silent=1,format='a',dmhalolist
readcol,'dM.list',silent=1,format='a',dmlist
readcol,'dL.dat',silent=1,format='i,i,d',l1,l2,l3

nlist = n_elements(dmlist)
readcol,dmlist[0],f='i,x,x',a    
n0 = n_elements(a)

if (n0 ne n_elements(l1)) then stop

dHarr = dblarr(nlist,n0)
dMarr = dblarr(nlist,n0)
light = dblarr(n0)

for j=0,nlist-1 do begin
    if (j eq 0) then begin
        readcol,dmlist[j],f='i,i,d',t1,t2,t3,silent=1
        Mr = t1
        Ma = t2
        dMarr[j,*] = t3
        readcol,dmhalolist[j],f='i,i,d',t1,t2,t3,silent=1
        dHarr[j,*] = t3
    endif else begin
        readcol,dmlist[j],f='x,x,d',t3,silent=1
        dMarr[j,*] = t3
        readcol,dmhalolist[j],f='x,x,d',t3,silent=1
        dHarr[j,*] = t3
    endelse
endfor

radl = radl * angleRmax
radm = radm * angleRmax
radh = radh * angleRmax
angl = asin(angl) * 360.0/(2*!pi)
angm = asin(angm) * 360.0/(2*!pi)
angh = asin(angh) * 360.0/(2*!pi)

na = n_elements(angle)
nr = n_elements(rad)

Mtot = dblarr(nr,nlist)
Mdark = Mtot
EtotU = dblarr(nr)
EtotL = dblarr(nr)
EdarkU = EtotU
EdarkL = EtotL

for j=0,nr-1 do begin
    ir = where(Mr eq rad[j])
    if (j eq 0) then light[j] = total(l3[ir]) else light[j] = light[j-1] + total(l3[ir])
    for k=0,nlist-1 do begin
        Mpm = total(dMarr[k,ir])
        Mph = total(dHarr[k,ir])
        if (j eq 0) then begin
            Mtot[j,k] = Mpm
            Mdark[j,k] = Mph
        endif else begin
            Mtot[j,k] = Mtot[j-1,k] + Mpm
            Mdark[j,k] = Mdark[j-1,k] + Mph
        endelse
    endfor
endfor

for j=0,nr-1 do begin
    EtotU[j] = max(Mtot[j,*]) + BHupper
    EtotL[j] = min(Mtot[j,*]) + BHlower
    EdarkU[j] = max(Mdark[j,*])
    EdarkL[j] = min(Mdark[j,*])
endfor

form1 = '(f7.2,1x,e14.8,1x,e14.8,1x,e14.8,1x,e14.8,1x,e14.8,1x,e14.8,1x,e14.8)'

free_lun,5
openw,5,'model_uncertainty.txt'
printf,5,'   R,   lumosity,      tot_mass,      tot_U,         tot_L,         DM_mass,       DM_U,          DM_L'
for j=0,nr-1 do begin
    printf,5,radm[j],light[j],Mtot[j,0]+BHmass,EtotU[j],$
      EtotL[j],Mdark[j,0],EdarkU[j],EdarkL[j],format=form1
endfor
free_lun,5

window,0,retain=2

plot,radm,Mtot[*,0]+BHmass,/xlog,/ylog,xrange=[1,2000],xstyle=1,$
  xtitle='Radius (arcsec)',ytitle='Enclosed Mass'
oplot,radm,EtotU,linestyle=2
oplot,radm,EtotL,linestyle=2

oplot,radm,Mdark[*,0]
oplot,radm,EdarkU,linestyle=2
oplot,radm,EdarkL,linestyle=2
stop
END

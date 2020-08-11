PRO findx2floor

; This routine is used to find the minimum X2 values at each
; parameter. It returns text files

;findfloor = 'yes' ;set to 'yes' if you haven't found the X^2 min for fitting. Set to 'no' if you have a want the code to read in a file.
findfloor = 'no'
;*****************************************************************

;readcol,cresfile,format='f,f,i,i,f,f,f,f,f',ml,bhm,vc,rs,x2stars,x2gc,x2tot,astars,agc
readcol,cresfile,format='f,f,f,f,f,f,f,x,x',ml,bhm,vc,rs,x2stars,x2gc,x2tot

n0 = n_elements(ml)

;a test to see if you've fit for the BH.
bhm1 = bhm[0]
ibhm = where(bhm eq bhm1,cbhm)
if (cbhm eq n0) then begin
    bhfixed = 'yes'
    print,'The BH is fixed! You are making 3 figures!'
endif else begin
    bhfixed = 'no'
    print,'The BH is NOT fixed! You are making 4 figures!'
endelse

; the lowest X2 values are fit with a polynomial...
if (findfloor eq 'yes') then begin
    ans = ''

    ;ITERATION ON M/L
    iml = bsort(ml)
    p = ml[iml]
    x2ml0 = x2tot[iml]
    x2ml1 = x2stars[iml]
    x2ml2 = x2gc[iml]
    mlx0 = [0.0]
    mlx1 = [0.0]
    mlx2 = [0.0]
    mly0 = [0.0]
    mly1 = [0.0]
    mly2 = [0.0]
    i = 0
    repeat begin
        p1 = p[i]
        p2 = p[i+1]
        if (p2 gt p1) then begin
            ii = where(ml eq p1)
            mlx0 = [mlx0,p1]
            mlx1 = [mlx1,p1]
            mlx2 = [mlx2,p1]
            mly0 = [mly0,min(x2tot[ii])]
            mly1 = [mly1,min(x2stars[ii])]
            mly2 = [mly2,min(x2gc[ii])]
        endif
        i = i + 1
    endrep until (p2 ge max(ml))
    ii = where(ml eq p2)
    mlx0 = [mlx0,p2]
    mlx1 = [mlx1,p2]
    mlx2 = [mlx2,p2]
    mly0 = [mly0,min(x2tot[ii])]
    mly1 = [mly1,min(x2stars[ii])]
    mly2 = [mly2,min(x2gc[ii])]
    mlx0 = mlx0[1:*]
    mlx1 = mlx1[1:*]
    mlx2 = mlx2[1:*]
    mly0 = mly0[1:*]
    mly1 = mly1[1:*]
    mly2 = mly2[1:*]
    
    free_lun,5
    openw,5,namechi+'.mltA'
    for k=0,n_elements(mlx0)-1 do printf,5,mlx0[k],mly0[k]
    free_lun,5
    openw,5,namechi+'.mlsA'
    for k=0,n_elements(mlx1)-1 do printf,5,mlx1[k],mly1[k]
    free_lun,5
    openw,5,namechi+'.mlgA'
    for k=0,n_elements(mlx2)-1 do printf,5,mlx2[k],mly2[k]
    free_lun,5

    ;ITERATION ON V_CIRC OR CONCENTRATION
    ivc = bsort(vc)
    p = vc[ivc]
    x2vc0 = x2tot[ivc]
    x2vc1 = x2stars[ivc]
    x2vc2 = x2gc[ivc]
    vcx0 = [0.0]
    vcx1 = [0.0]
    vcx2 = [0.0]
    vcy0 = [0.0]
    vcy1 = [0.0]
    vcy2 = [0.0]
    i = 0
    repeat begin
        p1 = p[i]
        p2 = p[i+1]
        if (p2 gt p1) then begin
            ii = where(vc eq p1)
            vcx0 = [vcx0,p1]
            vcx1 = [vcx1,p1]
            vcx2 = [vcx2,p1]
            vcy0 = [vcy0,min(x2tot[ii])]
            vcy1 = [vcy1,min(x2stars[ii])]
            vcy2 = [vcy2,min(x2gc[ii])]
        endif
        i = i + 1
    endrep until (p2 ge max(vc))
    ii = where(vc eq p2)
    vcx0 = [vcx0,p2]
    vcx1 = [vcx1,p2]
    vcx2 = [vcx2,p2]
    vcy0 = [vcy0,min(x2tot[ii])]
    vcy1 = [vcy1,min(x2stars[ii])]
    vcy2 = [vcy2,min(x2gc[ii])]
    vcx0 = vcx0[1:*]
    vcx1 = vcx1[1:*]
    vcx2 = vcx2[1:*]
    vcy0 = vcy0[1:*]
    vcy1 = vcy1[1:*]
    vcy2 = vcy2[1:*]

    free_lun,5
    openw,5,namechi+'.vctA'
    for k=0,n_elements(vcx0)-1 do printf,5,vcx0[k],vcy0[k]
    free_lun,5
    openw,5,namechi+'.vcsA'
    for k=0,n_elements(vcx1)-1 do printf,5,vcx1[k],vcy1[k]
    free_lun,5
    openw,5,namechi+'.vcgA'
    for k=0,n_elements(vcx2)-1 do printf,5,vcx2[k],vcy2[k]
    free_lun,5

    ;ITERATION OVER THE SCALE RADIUS
    irs = bsort(rs)
    p = rs[irs]
    x2rs0 = x2tot[irs]
    x2rs1 = x2stars[irs]
    x2rs2 = x2gc[irs]
    rsx0 = [0.0]
    rsx1 = [0.0]
    rsx2 = [0.0]
    rsy0 = [0.0]
    rsy1 = [0.0]
    rsy2 = [0.0]
    i = 0
    repeat begin
        p1 = p[i]
        p2 = p[i+1]
        if (p2 gt p1) then begin
            ii = where(rs eq p1)
            rsx0 = [rsx0,p1]
            rsx1 = [rsx1,p1]
            rsx2 = [rsx2,p1]
            rsy0 = [rsy0,min(x2tot[ii])]
            rsy1 = [rsy1,min(x2stars[ii])]
            rsy2 = [rsy2,min(x2gc[ii])]
        endif
        i = i + 1
    endrep until (p2 ge max(rs))
    ii = where(rs eq p2)
    rsx0 = [rsx0,p2]
    rsx1 = [rsx1,p2]
    rsx2 = [rsx2,p2]
    rsy0 = [rsy0,min(x2tot[ii])]
    rsy1 = [rsy1,min(x2stars[ii])]
    rsy2 = [rsy2,min(x2gc[ii])]
    rsx0 = rsx0[1:*]
    rsx1 = rsx1[1:*]
    rsx2 = rsx2[1:*]
    rsy0 = rsy0[1:*]
    rsy1 = rsy1[1:*]
    rsy2 = rsy2[1:*]

    free_lun,5
    openw,5,namechi+'.rstA'
    for k=0,n_elements(rsx0)-1 do printf,5,rsx0[k],rsy0[k]
    free_lun,5
    openw,5,namechi+'.rssA'
    for k=0,n_elements(rsx1)-1 do printf,5,rsx1[k],rsy1[k]
    free_lun,5
    openw,5,namechi+'.rsgA'
    for k=0,n_elements(rsx2)-1 do printf,5,rsx2[k],rsy2[k]
    free_lun,5

endif else begin
;total
    readcol,namechi+'.mltT',silent=1,f='f,f',mlx0,mly0
    readcol,namechi+'.vctT',silent=1,f='f,f',vcx0,vcy0
    readcol,namechi+'.rstT',silent=1,f='f,f',rsx0,rsy0
;stars
    readcol,namechi+'.mlsT',silent=1,f='f,f',mlx1,mly1
    readcol,namechi+'.vcsT',silent=1,f='f,f',vcx1,vcy1
    readcol,namechi+'.rssT',silent=1,f='f,f',rsx1,rsy1
;gc
    readcol,namechi+'.mlgT',silent=1,f='f,f',mlx2,mly2
    readcol,namechi+'.vcgT',silent=1,f='f,f',vcx2,vcy2
    readcol,namechi+'.rsgT',silent=1,f='f,f',rsx2,rsy2
endelse

;SMOOTHING AND FITTING FOR M/L
d = (max(mlx0) - min(mlx0)) * 100.0 ;total
xml0 = fltarr(d)
for j=0,d-1 do xml0[j] = min(mlx0) + (j*0.01)
smly = smooth(mly0,1)
yml0 = spline(mlx0,smly,xml0)
d = (max(mlx1) - min(mlx1)) * 100.0 ;stars
xml1 = fltarr(d)
for j=0,d-1 do xml1[j] = min(mlx1) + (j*0.01)
smly = smooth(mly1,1)
yml1 = spline(mlx1,smly,xml1)
d = (max(mlx2) - min(mlx2)) * 100.0 ;GC
xml2 = fltarr(d)
for j=0,d-1 do xml2[j] = min(mlx2) + (j*0.01)
smly = smooth(mly2,1)
yml2 = spline(mlx2,smly,xml2)

;SMOOTHING AND FITTING FOR V_CIRC OR CONCENTRATION
if (halotype eq 'log') then begin
    d = (max(vcx0) - min(vcx0)) ;total
    xvc0 = fltarr(d)
    for j=0,d-1 do xvc0[j] = min(vcx0) + j
    svcy = smooth(vcy0,1)
endif
if (halotype eq 'nfw') then begin
    d = (max(vcx0) - min(vcx0))*100
    xvc0 = fltarr(d)
    for j=0,d-1 do xvc0[j] = min(vcx0) + (j*0.01)
    svcy = smooth(vcy0,1)
endif
yvc0 = spline(vcx0,svcy,xvc0)
d = (max(vcx1) - min(vcx1)) ;stars
xvc1 = fltarr(d)
for j=0,d-1 do xvc1[j] = min(vcx1) + j
svcy = smooth(vcy1,1)
yvc1 = spline(vcx1,svcy,xvc1)
d = (max(vcx2) - min(vcx2)) ;GC
xvc2 = fltarr(d)
for j=0,d-1 do xvc2[j] = min(vcx2) + j
svcy = smooth(vcy2,1)
yvc2 = spline(vcx2,svcy,xvc2)

;SMOOTHING AND FITTING FOR THE SCALE RADIUS
d = (max(rsx0) - min(rsx0)) * 10.0 ;total
xrs0 = fltarr(d)
for j=0,d-1 do xrs0[j] = min(rsx0) + (j*0.1)
srsy = smooth(rsy0,1)                      ;
yrs0 = spline(rsx0,srsy,xrs0)
d = (max(rsx1) - min(rsx1)) * 10.0 ;stars
xrs1 = fltarr(d)
for j=0,d-1 do xrs1[j] = min(rsx1) + (j*0.1)
srsy = smooth(rsy1,1)
yrs1 = spline(rsx1,srsy,xrs1)
d = (max(rsx2) - min(rsx2)) * 10.0 ;GC
xrs2 = fltarr(d)
for j=0,d-1 do xrs2[j] = min(rsx2) + (j*0.1)
srsy = smooth(rsy2,1)                      ;
yrs2 = spline(rsx2,srsy,xrs2)

;and the linear shifts are applied...
yml1 = yml1 + os0s
yml2 = yml2 + os0g
yvc1 = yvc1 + os1s
yvc2 = yvc2 + os1g
yrs1 = yrs1 + os2s
yrs2 = yrs2 + os2g

set_plot,'ps'
device,file=nameout1,/color

!p.multi = [0,1,3,0,1]
!y.omargin = [2,4]
loadct,0
plot,ml,x2tot,psym=sym(1),symsize=0.2,xtitle='!3M/L',ytitle='!7v!6!u2!6',$
  yrange=[x2dn,x2up],ystyle=1,xthick=2.5,ythick=2.5,charthick=2.5,charsize=1.2,$
  xrange=[mldn,mlup],xstyle=1,position=[pw0,0.71,pw1,0.96],yminor=2,xticklen=0.05
oplot,xml0,yml0,thick=3
oplot,mlx0,mly0,psym=sym(6),symsize=0.5
oplot,mlx0,mly0,psym=sym(1),symsize=0.3
loadct,4,/silent
oplot,xml1,yml1,color=150
oplot,xml2,yml2,color=60 
loadct,0
plot,vc,x2tot,psym=sym(1),symsize=0.2,xtitle=xtitle1,ytitle='!7v!6!u2!6',$
  yrange=[x2dn,x2up],ystyle=1,xthick=2.5,ythick=2.5,charthick=2.5,charsize=1.2,$
  xrange=[vcdn,vcup],xstyle=1,position=[pw0,0.39,pw1,0.64],yminor=2,xticklen=0.05
oplot,xvc0,yvc0,thick=3
oplot,vcx0,vcy0,psym=sym(6),symsize=0.5
oplot,vcx0,vcy0,psym=sym(1),symsize=0.3
loadct,4,/silent
oplot,xvc1,yvc1,color=150
oplot,xvc2,yvc2,color=60 
loadct,0
plot,rs,x2tot,psym=sym(1),symsize=0.2,xtitle=xtitle2,ytitle='!7v!6!u2!6',$
  yrange=[x2dn,x2up],ystyle=1,xthick=2.5,ythick=2.5,charthick=2.5,charsize=1.2,$
  xrange=[rsdn,rsup],xstyle=1,position=[pw0,0.07,pw1,0.32],xticklen=0.05,$
  yminor=2
oplot,xrs0,yrs0,thick=3
oplot,rsx0,rsy0,psym=sym(6),symsize=0.5
oplot,rsx0,rsy0,psym=sym(1),symsize=0.3
loadct,4,/silent
oplot,xrs1,yrs1,color=150
oplot,xrs2,yrs2,color=60 

device,/close_file
jump1:

set_plot,'ps'
!p.multi = [0,1,1,0,0]

G = 4.516e-30
rho = fltarr(n0)
for j=0,n0-1 do begin
    top = 3 * vc[j]^2
    r = rs[j]*3.086e13
    bot = 4 * 3.14159 * G * r^2
    rho[j] = top / bot
endfor

device,file='chidensity.ps'
plot,rho,x2tot,psym=sym(1),symsize=0.2,xtitle='!3Central Density (M_sun/pc)',$
  xthick=4,ythick=4,charthick=4,pos=[0.13,0.13,0.93,0.89],$
  ytitle='!7v!6!u2!6'
device,/close


STOP
END

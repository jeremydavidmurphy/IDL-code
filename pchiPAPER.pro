PRO pchiPAPER, cresfile

; This is a generic template used to plot the resulting X^2 values
; from your dynamical modeling results. After running Pchi.pro, which
; does the initial sorting and adjustment for the alpha parameter, 

; CRESFILE: The SORTED cres results.

pw0 = 0.23
pw1 = 0.48
pw2 = 0.56
pw3 = 0.81

x2up = 1950
x2dn = 1850

halotype = 'log' ;change to 'log' or 'nfw'
findfloor = 'yes' ;set to 'yes' if you haven't found the X^2 min for fitting. Set to 'no' if you have a want the code to read in a file.
;*****************************************************************

if (halotype eq 'log') then begin
    xtitle1 = '!3V!dc!n (km/sec)'
    xtitle2 = '!3R!ds!n (kpc)'
    nameout1 = 'chi1LOG.ps'
    nameout2 = 'chi2LOG.ps'
    namechi = 'x2minLOG'
endif
if (halotype eq 'nfw') then begin
    xtitle1 = '!3Concentration'
    xtitle2 = '!3R!ds!n (kpc)'
    nameout1 = 'chi1NFW.ps'
    nameout2 = 'chi2NFW.ps'
    namechi = 'x2minNFW'
endif

;readcol,cresfile,format='f,f,i,i,f,f,f,f,f',ml,bhm,vc,rs,x2stars,x2gc,x2tot,astars,agc
readcol,cresfile,format='f,f,i,i,f,f,f',ml,bhm,vc,rs,x2stars,x2gc,x2tot

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

mlup = max(ml)+0.5
mldn = min(ml)-0.5
bhmup = max(bhm)+1e9
bhmdn = min(bhm)-1e9
vcup = max(vc)+50
vcdn = min(vc)-50
rsup = max(rs)+1
rsdn = min(rs)-1

i = where(x2tot lt x2up,count)
if (count ne 0) then mldn = min(ml[i]);-0.5
if (count ne 0) then mlup = max(ml[i])+0.15
if (count ne 0) then bhmdn = min(bhm[i]);-0.5e9
if (count ne 0) then bhmup = max(bhm[i]);+0.5e9
if (count ne 0) then vcdn = min(vc[i])-25
if (count ne 0) then vcup = max(vc[i]);+50
if (count ne 0) then rsdn = min(rs[i]);-1
if (count ne 0) then rsup = max(rs[i])+0.5

; the lowest X2 values are fit with a polynomial...
if (findfloor eq 'yes') then begin
    ans = ''

    ;ITERATION ON M/L
    iml = bsort(ml)
    p = ml[iml]
    x2ml = x2tot[iml]
    mlx = [0.0]
    mly = [0.0]
    i = 0
    repeat begin
        p1 = p[i]
        p2 = p[i+1]
        if (p2 gt p1) then begin
            ii = where(ml eq p1)
            mlx = [mlx,p1]
            mly = [mly,min(x2tot[ii])]
        endif
        i = i + 1
    endrep until (p2 ge max(ml))
    ii = where(ml eq p2)
    print,p2,min(x2tot[ii])
    mlx = [mlx,p2]
    mly = [mly,min(x2tot[ii])]
    mlx = mlx[1:*]
    mly = mly[1:*]

    free_lun,5
    openw,5,namechi+'.ml'
    for k=0,n_elements(mlx)-1 do printf,5,mlx[k],mly[k]
    free_lun,5

    ;ITERATION ON V_CIRC OR CONCENTRATION
    ivc = bsort(vc)
    p = vc[ivc]
    x2vc = x2tot[ivc]
    vcx = [0.0]
    vcy = [0.0]
    i = 0
    repeat begin
        p1 = p[i]
        p2 = p[i+1]
        if (p2 gt p1) then begin
            ii = where(vc eq p1)
            vcx = [vcx,p1]
            vcy = [vcy,min(x2tot[ii])]
        endif
        i = i + 1
    endrep until (p2 ge max(vc))
    ii = where(vc eq p2)
    print,p2,min(x2tot[ii])
    vcx = [vcx,p2]
    vcy = [vcy,min(x2tot[ii])]
    vcx = vcx[1:*]
    vcy = vcy[1:*]
    
    free_lun,5
    openw,5,namechi+'.vc'
    for k=0,n_elements(mlx)-1 do printf,5,vcx[k],vcy[k]
    free_lun,5

    ;ITERATION OVER THE SCALE RADIUS
    irs = bsort(rs)
    p = rs[irs]
    x2rs = x2tot[irs]
    rsx = [0.0]
    rsy = [0.0]
    i = 0
    repeat begin
        p1 = p[i]
        p2 = p[i+1]
        if (p2 gt p1) then begin
            ii = where(rs eq p1)
            rsx = [rsx,p1]
            rsy = [rsy,min(x2tot[ii])]
        endif
        i = i + 1
    endrep until (p2 ge max(rs))
    ii = where(rs eq p2)
    print,p2,min(x2tot[ii])
    rsx = [rsx,p2]
    rsy = [rsy,min(x2tot[ii])]
    rsx = rsx[1:*]
    rsy = rsy[1:*]

    free_lun,5
    openw,5,namechi+'.rs'
    for k=0,n_elements(mlx)-1 do printf,5,rsx[k],rsy[k]
    free_lun,5

endif

readcol,namechi+'.ml',silent=1,f='f,f',mlx,mly
readcol,namechi+'.vc',silent=1,f='f,f',vcx,vcy
readcol,namechi+'.rs',silent=1,f='f,f',rsx,rsy


;SMOOTHING AND FITTING FOR M/L
d = (max(mlx) - min(mlx)) * 100.0
xml = fltarr(d)
for j=0,d-1 do xml[j] = min(mlx) + (j*0.01)
smly = smooth(mly,3)
yml = spline(mlx,smly,xml)

;SMOOTHING AND FITTING FOR V_CIRC OR CONCENTRATION
d = (max(vcx) - min(vcx))
xvc = fltarr(d)
for j=0,d-1 do xvc[j] = min(vcx) + j
svcy = smooth(vcy,2)
yvc = spline(vcx,svcy,xvc)

;SMOOTHING AND FITTING FOR THE SCALE RADIUS
d = (max(rsx) - min(rsx)) * 10.0
xrs = fltarr(d)
for j=0,d-1 do xrs[j] = min(rsx) + (j*0.1)
srsy = rsy                      ;
yrs = spline(rsx,srsy,xrs)

set_plot,'ps'
device,file=nameout1,/color

!p.multi = [0,1,3,0,1]
!y.omargin = [2,4]
loadct,0
plot,ml,x2tot,psym=sym(1),symsize=0.2,xtitle='!3M/L',ytitle='!7v!6!u2!6',$
  yrange=[x2dn,x2up],ystyle=1,xthick=1.5,ythick=1.5,charthick=1.5,charsize=1.2,$
  xrange=[mldn,mlup],xstyle=1,position=[pw0,0.71,pw1,0.96],yminor=2,xticklen=0.05
oplot,xml,yml

plot,vc,x2tot,psym=sym(1),symsize=0.2,xtitle=xtitle1,ytitle='!7v!6!u2!6',$
  yrange=[x2dn,x2up],ystyle=1,xthick=1.5,ythick=1.5,charthick=1.5,charsize=1.2,$
  xrange=[vcdn,vcup],xstyle=1,position=[pw0,0.39,pw1,0.64],yminor=2,xticklen=0.05
oplot,xvc,yvc

plot,rs,x2tot,psym=sym(1),symsize=0.2,xtitle=xtitle2,ytitle='!7v!6!u2!6',$
  yrange=[x2dn,x2up],ystyle=1,xthick=1.5,ythick=1.5,charthick=1.5,charsize=1.2,$
  xrange=[rsdn,rsup],xstyle=1,position=[pw0,0.07,pw1,0.32],xticklen=0.05,$
  yminor=2
oplot,xrs,yrs

device,/close_file
set_plot,'x'

!p.multi = [0,1,1,0,0]

STOP
END

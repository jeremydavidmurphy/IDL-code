; This function is used in BINGEN.PRO to generate the bin#### and
; bnn#### lists. It's a modified version of ROTATE.PRO

; NOTES ON THE DIRECTION OF ROTATION:
; Currently this transformation is a bit screwy. The correct PA is to
; enter the degree of rotation counterclockwise from the North
; direction in positve degrees. Example: The PA of N4472 is given as
; +155 or -25. This would be entered as 25 ('bin' ends up South) or
; 180+25=205 ('bin' ends up North).

FUNCTION rotateF, offsetfile, pa, letter

plotting = 'on'
readcol,offsetfile,silent=1,format='i,f,f,f',fibF,radF,dRAF,dDecF

theta = -(90.0 - pa) ;converted to Cartesian coordinates
;theta = (90.0 + pa)

print,theta
theta = theta/(360.0/(2*!pi)) ;converted to radians)
print,theta

n00 = n_elements(fibF)
;the next 4 arrays give the values of their location IN THE REFERENCE
;FRAME OF THE GALAXY! (this is both RA and Dec offsets and also in
;polar coordinates (Rp and Ap). Once this coordinate transformation is
;made then I can pick out which fibers are where on the galaxy.
dRAp = dblarr(n00)
dDecp = dblarr(n00)
Rp = dblarr(n00)
Ap = dblarr(n00)

for jj=0,n00-1 do begin
    dRAp[jj] = dRAF[jj]*cos(theta) - dDecF[jj]*sin(theta)
    dDecp[jj] = dRAF[jj]*sin(theta) + dDecF[jj]*cos(theta)
    Ap[jj] = atan(dDecp[jj] / dRAp[jj])
    Rp[jj] = dRAp[jj] / cos(Ap[jj])
    print,jj+1,Rp[jj],Ap[jj]
endfor

; dRAp and dDecp are the rotated coordinates in RA and Dec offsets
; Rp and Ap are the rotated polar coordinates
; This new coordinate system, dRAp, dDecp, Rp and Ap are now in the
; reference frame of the galaxy.

if (plotting eq 'on') then begin
    ramx = max([max(dRAF),max(dRAp)])
    ramn = min([min(dRAF),min(dRAp)])
    decmx = max([max(dDecF),max(dDecp)])
    decmn = min([min(dDecF),min(dDecp)])
    window,2,retain=2,xsize=500,ysize=500
    device,decomposed=0
    
    loadct,0
    plot,dRAF,dDecF,psym=sym(0),xrange=[ramn-3.0,ramx+3.0],$
      yrange=[decmn-3.0,decmx+3.0],/xstyle,/ystyle,/isotropic,$
      title='Fibers in Galaxy Coordinate System'
    print,'Now transforming into galaxy reference frame...'
    oplot,dRAp,dDecp,psym=sym(6)
    loadct,4
endif

readcol,'bin_r.out',format='i,f,f,f,f',skipline=2,rbin,rlow,rmid,rup,rsize,silent=1
readcol,'bin_v.out',format='i,f,f,f,f',skipline=2,abin,alow,amid,aup,asize,silent=1

n11 = n_elements(rbin)
n22 = n_elements(abin)
alow = alow/(360.0/(2*!pi))
aup = aup/(360.0/(2*!pi))

cntr = 0
color = [60,110,180,255,150]
for jj=0,n22-1 do begin ;a loop through the angular bins
    a = abin[jj]
    al = alow[jj]
    ah = aup[jj]
    for k=0,n11-1 do begin ;a loop through the radial bins
        r = rbin[k]
        rl = rlow[k]
        rh = rup[k]
        index = where(abs(Rp) ge rl and abs(Rp) lt rh and $
                      abs(Ap) ge al and abs(Ap) lt ah,count)
        if (count ne 0) then begin

; At this point the fibers that fall into a given radial and angular
; bin are selected. however, this doesn't account for where on the
; galaxy it is. For axisymmetric models, just the bin and bnn
; distinctions can be made. But when you want to break this up between
; the four quadrants you'll have to work on the following section by
; using the A + and - values to subdivide the bin and bnn outputs. The
; coordinate convention is as follows:
; Q1(0 to 90): R is +, A is +
; Q2(90 to 180): R is +, A is -
; Q3(180 to 270): R is -, A is +
; Q4(270 to 360): R is -, A is - 

            ibin = where(Rp[index] ge 0.0,cplus)
            ibnn = where(Rp[index] lt 0.0,cminus)

            if (cplus gt 0) then begin
                if (r le 9) then begin
                    tempnameF = 'bin'+strn(0)+strn(r)+strn(0)+strn(a)
                    readcol,tempnameF,f='a',preexisting,silent=1
                    openw,55,tempnameF
                    for mm=0,n_elements(preexisting)-1 do $
                      printf,55,preexisting[mm]
                    for mm=0,cplus-1 do $
                      printf,55,strn(index[ibin[mm]]+1)+'_'+letter
                endif else begin
                    tempnameF = 'bin'+strn(r)+strn(0)+strn(a)
                    readcol,tempnameF,f='a',preexisting,silent=1
                    openw,55,tempnameF
                    for mm=0,n_elements(preexisting)-1 do $
                      printf,55,preexisting[mm]
                    for mm=0,cplus-1 do $
                      printf,55,strn(index[ibin[mm]]+1)+'_'+letter
                endelse
                free_lun,55
            endif

            if (cminus gt 0) then begin
                if (r le 9) then begin
                    tempnameF = 'bnn'+strn(0)+strn(r)+strn(0)+strn(a)
                    readcol,tempnameF,f='a',preexisting,silent=1
                    openw,55,tempnameF
                    for mm=0,n_elements(preexisting)-1 do $
                      printf,55,preexisting[mm]
                    for mm=0,cminus-1 do $
                      printf,55,strn(index[ibnn[mm]]+1)+'_'+letter
                endif else begin
                    tempnameF = 'bnn'+strn(r)+strn(0)+strn(a)
                    readcol,tempnameF,f='a',preexisting,silent=1
                    openw,55,tempnameF
                    for mm=0,n_elements(preexisting)-1 do $
                      printf,55,preexisting[mm]
                    for mm=0,cminus-1 do $
                      printf,55,strn(index[ibnn[mm]]+1)+'_'+letter
                endelse
                free_lun,55
            endif
            temp = where(Rp[index] eq -1,count)
            if (count eq 0 and plotting eq 'on') then begin
                print,'Plotting axis # '+strn(cntr+1)
                wait,0.2
                oplot,Rp[index],Ap[index],psym=sym(1),/polar,$
                  color=color[cntr],symsize=1.5
            endif
        endif
    endfor
    cntr = cntr + 1
endfor
                
if (plotting eq 'on') then begin
    pause
    wdelete
endif

return,'test'
END

; This routine is used to generate the bins that go into PIPE2. It
; uses the dotR files to get the dRA and dDec of each fiber, then
; makes a rotation (based on the galaxy PA). The PA (position angle)
; is read from the north.

; This code will look for the bin_r.out and bin_v.out files.

; APRIL 25TH, 2010: bingen.pro now accepts the data.list file (see
; below for format of this list), galRA, galDec, and PA, and returns
; bin bin#### and bnn#### lists needed for PIPE2.pro. IT USES A
; FUNCTIONAL FORM OF ROTATE named rotateF
PRO rotate, dotR, pa

plotting = 'on' ;change this to supress the plotting to screen
letter = 'a' ; gets subscripted to the end of the fiber
slow = 'on' ;this turns off or on some pauses in the code

readcol,dotR,silent=1,format='i,f,f,f,x,x',skipline=1,$
  fib,rad,dRA,dDec

ibad = where(fib eq -1,count)
igood = where(fib ne -1)

theta = (90.0 - pa) ;converted to Cartesian coordinates
theta = theta/(360.0/(2*!pi)) ;converted to radians)

n0 = n_elements(fib)
;the next 4 arrays give the values of their location IN THE REFERENCE
;FRAME OF THE GALAXY! (this is both RA and Dec offsets and also in
;polar coordinates (Rp and Ap). Once this coordinate transformation is
;made then I can pick out which fibers are where on the galaxy.
dRAp = dblarr(n0)
dDecp = dblarr(n0)
Rp = dblarr(n0)
Ap = dblarr(n0)

for j=0,n0-1 do begin
    if (fib[j] eq -1) then begin
        dRAp[j] = -1.0
        dDecp[j] = -1.0
        Ap[j] = -1.0
        Rp[j] = -1.0
        goto, jump1
    endif
    dRAp[j] = dRA[j]*cos(theta) - dDec[j]*sin(theta)
    dDecp[j] = dRA[j]*sin(theta) + dDec[j]*cos(theta)
    Ap[j] = atan(dDecp[j] / dRAp[j])
    Rp[j] = dRAp[j] / cos(Ap[j])
    print,j+1,Rp[j],Ap[j]
jump1:
endfor

; dRAp and dDecp are the rotated coordinates in RA and Dec offsets
; Rp and Ap are the rotated polar coordinates
; This new coordinate system, dRAp, dDecp, Rp and Ap are now in the
; reference frame of the galaxy.

if (plotting eq 'on') then begin
    ramx = max([max(dRA),max(dRAp)])
    ramn = min([min(dRA),min(dRAp)])
    decmx = max([max(dDec),max(dDecp)])
    decmn = min([min(dDec),min(dDecp)])
    window,2,retain=2,xsize=500,ysize=500
    device,decomposed=0
    
    loadct,0
    plot,dRA[igood],dDec[igood],psym=sym(0),xrange=[ramn-3.0,ramx+3.0],$
      yrange=[decmn-3.0,decmx+3.0],/xstyle,/ystyle,/isotropic
;    oplot,dRA[igood],dDec[igood],psym=sym(6)
    pause
;    loadct,0
    oplot,dRAp[igood],dDecp[igood],psym=sym(6)
    loadct,4
;    oplot,Rp[igood],Ap[igood],psym=1,color=60,/polar
endif

readcol,'bin_r.out',format='i,f,f,f,f',skipline=2,rbin,rlow,rmid,rup,rsize
readcol,'bin_v.out',format='i,f,f,f,f',skipline=2,abin,alow,amid,aup,asize

n1 = n_elements(rbin)
n2 = n_elements(abin)
alow = alow/(360.0/(2*!pi))
aup = aup/(360.0/(2*!pi))

cntr = 0
color = [60,110,180,255,150]
for j=0,n2-1 do begin
    a = abin[j]
    al = alow[j]
    ah = aup[j]
    for k=0,n1-1 do begin
        r = rbin[k]
        rl = rlow[k]
        rh = rup[k]
        index = where(abs(Rp) ge rl and abs(Rp) lt rh and $
                      abs(Ap) ge al and abs(Ap) lt ah,count)
        if (count ne 0) then begin

; at this point the fibers that fall into a given radial and angular
; bin are selected. however, this doesn't account for where on the
; galaxy it is. for axisymmetric models the just the bin and bnn
; distinctions can be made. but when you want to break this up between
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
                openw,5,'bin'+strn(0)+strn(r)+strn(0)+strn(a)
                print,'bin'+strn(0)+strn(r)+strn(0)+strn(a)
                for m=0,cplus-1 do begin
                    printf,5,strn(index[ibin[m]]+1)+'_'+letter
                    print,strn(index[ibin[m]]+1)+'_'+letter
                endfor
                free_lun,5
            endif

            if (cminus gt 0) then begin
                openw,5,'bnn'+strn(0)+strn(r)+strn(0)+strn(a)
                print,'bnn'+strn(0)+strn(r)+strn(0)+strn(a)
                for m=0,cminus-1 do begin
                    printf,5,strn(index[ibnn[m]]+1)+'_'+letter
                    print,strn(index[ibnn[m]]+1)+'_'+letter
                endfor
                free_lun,5
            endif
            print,rl,rh
            print,Rp[index]
            print,''
            print,al,ah
            print,Ap[index]
            print,''
            temp = where(Rp[index] eq -1,count)
            if (count eq 0) then begin ;all the values are good
                oplot,Rp[index],Ap[index],psym=sym(1),/polar,$
                  color=color[cntr],symsize=1.5
;                openw,5,
            endif
                
            if (slow eq 'on') then pause
        endif
    endfor
    cntr = cntr + 1
endfor
;        i1 = where(abs(Rp) ge rl and abs(Rp) lt rh,count1)
;        if (count1 ne -1) then begin ;restricted to annulus
;            i2 = where(abs(Ap) ge al and abs(Ap) lt ah,count2)
;            if (count2 ne -1) then begin ;restricted to slice
;                for m=0,count1-1 do begin
;                    rs = Rp[i1[m]]
;                    for n=0,count2-1 do as = Ap[i2[n]]
;                    if (rs ge 0.0) then begin
;                        if (as ge 0.0) then b1
                
wdelete
stop          
for j=0,n1-1 do begin ;a loop through each radial bin
    print,'radial bin #: '+strn(j+1)
    r = rbin(j)
    rl = rlow[j]
    rh = rup[j]
    for k=0,n2-1 do begin ;a loop through each angular bin
        print,'angular bin #: '+strn(k)
        a = abin[k]
        al = alow[k]
        ah = aup[k]
        iralow = where(dRAp gt rl*cos(al))
        irahigh = where(dRAp lt rh*cos(al))
        i = where(iralow eq irahigh,count)
        if (count ne 0) then begin
            ira = iralow[where(iralow eq irahigh)]
            print,dRAp[ira]
        endif
        ideclow = where(dDecp gt rl*sin(al))
        idechigh = where(dDecp lt rh*sin(al))
        i = where(ideclow eq idechigh,count)
        if (count ne 0) then begin
            idec = ideclow[where(ideclow eq idechigh)]
            print,dRAp[idec]
        endif

    endfor
endfor


wdelete

stop
END

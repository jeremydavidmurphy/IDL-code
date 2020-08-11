PRO bingen, galRA, galDec, galPA, coordlist

; This routine is used to create the lists that form the bins. It
; reads in both the bin_r.out and bin_v.out for the model/galaxy, the
; IFUcen.txt file, and the coordinate files for each of the pointings
; that goes into a given bin.

; galRA: the right ascension of the galaxy. hh:mm:ss
; galDec: the declination of the galaxy. dd:mm:ss
; galPA: the position angle of the galaxy. the value from NED, that is
; a degree rotation of the major axis from N and towards the E is proper.

;COORDLIST: A list of the coordinate files. These are the fiber RA and
;Dec values from the finder_code output. (i.e. the galaxy_D1_coords.txt file)

;COMPILE: rotateF.pro

;******************************************************************

;the galaxy RA and Dec are read in, then converted to decimal coordinates.
galRA = strn(galRA) & galDec = strn(galDec)
galRA = strsplit(galRA,':',/EXTRACT)
galDec = strsplit(galDec,':',/EXTRACT)

galRA = double((galRA[0]+galRA[1]/60.+galRA[2]/3600.)*15)
;The + or - sign of the declination is accounted for
IF (galDec[0] LT 0) THEN galDec = double(galDec[0] - galDec[1]/60. - galDec[2]/3600.)
IF (galDec[0] GE 0) THEN galDec = double(galDec[0] + galDec[1]/60. + galDec[2]/3600.)

readcol,'bin_r.out',format='i,f,x,f,x',skipline=2,binR,rlow,rhigh,silent=1
readcol,'bin_v.out',format='i,f,x,f,x',skipline=2,binV,vlow,vhigh,silent=1

;the fake bins are created. these are bins for all possible values
;with a -1 value in them. 
for j=0,n_elements(binR)-1 do begin
    for k=0,n_elements(binV)-1 do begin
        free_lun,10
        if (binR[j] le 9) then begin
            openw,10,'bin'+strn(0)+strn(binR[j])+strn(0)+strn(binV[k])
            printf,10,-1
            free_lun,10
            openw,10,'bnn'+strn(0)+strn(binR[j])+strn(0)+strn(binV[k])
            printf,10,-1
            free_lun,10
        endif else begin
            openw,10,'bin'+strn(binR[j])+strn(0)+strn(binV[k])
            printf,10,-1
            free_lun,10
            openw,10,'bnn'+strn(binR[j])+strn(0)+strn(binV[k])
            printf,10,-1
            free_lun,10
        endelse
    endfor
endfor

readcol,coordlist,format='a,a',list,pointing,silent=1
n0 = n_elements(list)

form1 = '(i3,2x,f10.3,f10.3,f10.3)'
for j=0,n0-1 do begin
    Cname = list[j]
    ptg = pointing[j]
    readcol,Cname,format='i,a,a',fnum,fRA,fDec,silent=1
    n1 = n_elements(fnum)
    free_lun,5
    tname = strsplit(Cname,'.',/extract)
    tname = tname[0]
    tname = tname+'.std'
    openw,5,tname
    for k=0,n1-1 do begin
        oneRA = fRA[k]
        oneDec = fDec[k]
        oneRA = strsplit(oneRA,':',/EXTRACT)
        oneDec = strsplit(oneDec,':',/EXTRACT)
        RA = double((oneRA[0]+oneRA[1]/60.+oneRA[2]/3600.)*15.)
        if (oneDec[0] lt 0) then Dec = double(oneDec[0]-oneDec[1]/60.-oneDec[2]/3600.)
        if (oneDec[0] gt 0) then Dec = double(oneDec[0]+oneDec[1]/60.+oneDec[2]/3600.)
        dRA = double((RA - galRA) * 3600.0)
        dDec = double((Dec - galDec) * 3600.0)
        rad = float(sqrt(dRA^2 + dDec^2))
        out = [fnum[k],rad,dRA,dDec]
;        print,out
        printf,5,out,format=form1
    endfor
    free_lun,5
    what = rotateF(tname,galPA,ptg)
endfor

stop
END

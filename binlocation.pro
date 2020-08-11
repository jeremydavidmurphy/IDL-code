PRO binlocation, galra, galdec

; This routine uses the output of pipe1.pro (the bin####.coord files)
; to determine a light-weighted center for a given spatial bin.
;
;M49: '12:29:46.7','08:00:02'

;*******************************************************************
RA = strn(galRA) & Dec = strn(galDec)
RA = strsplit(RA,':',/EXTRACT)
Dec = strsplit(Dec,':',/EXTRACT)
RA = double((RA[0]+RA[1]/60.+RA[2]/3600.)*15)
;The + or - sign of the declination is accounted for
IF (Dec[0] LT 0) THEN Dec = double(Dec[0]-Dec[1]/60.-Dec[2]/3600.)
IF (Dec[0] GE 0) THEN Dec = double(Dec[0]+Dec[1]/60.+Dec[2]/3600.)

readcol,'coord.list',f='a',coordfiles
n0 = n_elements(coordfiles)

free_lun,5
openw,5,'binlocation.txt'
for j=0,n0-1 do begin
    bin = coordfiles[j]
    readcol,bin,f='x,d,d,x,d',skip=1,dra,ddec,flux
    dra = dra/3600.0
    ddec = ddec/3600.0
    ara = ra + dra
    adec = dec + ddec
    n1 = n_elements(dra)
    top = fltarr(n1)
    bot = flux
    for k=0,n1-1 do top[k] = ara[k] * flux[k]
    wra = total(top) / total(bot)
    for k=0,n1-1 do top[k] = adec[k] * flux[k]
    wdec = total(top) / total(bot)
    wra = string(wra)
    wdec = string(wdec)
    printf,5,[bin,wra,wdec]
endfor
free_lun,5

stop
END


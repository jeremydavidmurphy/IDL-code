;THIS CODE HAS BEEN SUPERSEDED BY THE NEW PIPE1 THAT INCLUDES THIS IN
;THE REDUCTION. THE DELTA-RA AND DELTA-DEC VALUES ARE NOW OUTPUT INTO
;THE dotR FILE.

PRO delfib, name, RA, Dec

;This code is a slight modification of the regular radius.pro
;routine. It includes the regular radial position for each fiber, but
;now it also outputs a file named "name.RAnDEC" which is the delta X
;and Y offsets of each fiber in terms of arcseconds.  This is just the
;legs of the triangle, rather than the hypotenuse...

;RA and Dec values (which are the center of the galaxy) must be seperated by colons, i.e. 04:34:45.6
;the 'name' is the name you orignially gave your fields when
;generating the finder charts.
;THE COORDINATE TEXT FILE MUST BE IN THE CALLING DIRECTORY

;cd,'radius'
file = name+'_D1_coords.txt'
RA = strn(RA) & Dec = strn(Dec)
RA = strsplit(RA,':',/EXTRACT)
Dec = strsplit(Dec,':',/EXTRACT)
rad = dblarr(2,247)
rads = rad
darc = dblarr(3,247)

RA = double((RA[0]+RA[1]/60.+RA[2]/3600.)*15)
IF (Dec[0] LT 0) THEN Dec = double(Dec[0]-Dec[1]/60.-Dec[2]/3600.)
IF (Dec[0] GE 0) THEN Dec = double(Dec[0]+Dec[1]/60.+Dec[2]/3600.)


openr,5,file
coord = strarr(1,247)
readf,5,coord
free_lun,5
cor = dblarr(3,247)
FOR j=0,246 DO BEGIN
    fc = coord[0,j]
    fc = strsplit(fc,' :',/EXTRACT)
    cor[0,j] = fc[0]
    cor[1,j] = (fc[1]+fc[2]/60.+fc[3]/3600.)*15.
    IF (fc[4] LT 0) THEN BEGIN
        cor[2,j] = fc[4]-fc[5]/60.-fc[6]/3600.
    ENDIF ELSE BEGIN
        cor[2,j] = fc[4]+fc[5]/60.+fc[6]/3600.
    ENDELSE
;    print,coord[0,j]
;    print,transpose(fc)
;    print,'Fiber '+strn(j+1)+' has an RA and Dec of: ',cor[1,j],cor[2,j]
ENDFOR

FOR j=0,246 DO BEGIN
    dx = double((cor[1,j]-RA))
    dxs = double((cor[1,j]-RA)^2)
    dy = double((cor[2,j]-Dec))
    dys = double((cor[2,j]-Dec)^2)
print,dx
print,dy
    darc[0,j] = j+1
    darc[1,j] = dx*3600
    darc[2,j] = dy*3600
    rad[0,j] = j+1
    rad[1,j] = sqrt(dxs+dys)*3600.
    Print, 'The radius for fiber '+strn(j+1)+ ' is: ',strn(rad[1,j])
ENDFOR

indy = bsort(rad[1,*])
rads[0,*] = rad[0,indy]
rads[1,*] = rad[1,indy]

OpenW,6,name+'.rad'
PrintF,6,rad
Free_Lun,6

OpenW,7,name+'.std.rad'
PrintF,7,rads
Free_Lun,7

OpenW,8,name+'.RAnDEC'
PrintF,8,darc
Free_Lun,8
;cd,'../'
END


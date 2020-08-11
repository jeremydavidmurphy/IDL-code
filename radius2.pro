;This code was re-written to output the delta RA and delta Dec values
;along with the radial positions of each fiber relative to the center
;of the galaxy. The coordinate text (NAME_D1_COORDS.TXT) must be in
;the calling directory
;RA and Dec values (which are the center of the galaxy) must be
;seperated by colons, i.e. 04:34:45.6 the 'name' is the name you
;orignially gave your fields when generating the finder charts.

;The output is named NAME.RAD and NAME.STD.RAD Each have 4 colummns
;with the only difference is that std.rad has been sorted by radial
;position.
;  fiber #      radial position      RA offset      Dec offset

PRO radius, name, RA, Dec

; EX: radius,'CGCG137-019_D1','16:02:30.5','21:07:15'

file = name+'_coords.txt'
RA = strn(RA) & Dec = strn(Dec)
RA = strsplit(RA,':',/EXTRACT)
Dec = strsplit(Dec,':',/EXTRACT)
RA = double((RA[0]+RA[1]/60.+RA[2]/3600.)*15)
;The + or - sign of the declination is accounted for
IF (Dec[0] LT 0) THEN Dec = double(Dec[0]-Dec[1]/60.-Dec[2]/3600.)
IF (Dec[0] GE 0) THEN Dec = double(Dec[0]+Dec[1]/60.+Dec[2]/3600.)

;The RA and Dec for each fiber is read in.
readcol,file,format='I,A,A',fib,fra,fdec
n1 = n_elements(fib)
rad = dblarr(4,n1) ;The radius array (fiber#, rad, deltaRA, deltaDec)
cor = dblarr(3,n1) ;The coordinate array (fiber#, RA, Dec)
rads = rad

FOR j=0,n1-1 DO BEGIN
    cor[0,j] = fib[j]
    fc = fra[j]
    fd = fdec[j]
    fc = strsplit(fc,':',/EXTRACT)
    fd = strsplit(fd,':',/EXTRACT)
;The coordinates are converted to decimal degrees
    cor[1,j] = (fc[0]+fc[1]/60.+fc[2]/3600.)*15.
    IF (fd[0] LT 0) THEN cor[2,j] = double(fd[0]-fd[1]/60.-fd[2]/3600.)
    IF (fd[0] GT 0) THEN cor[2,j] = double(fd[0]+fd[1]/60.+fd[2]/3600.)
ENDFOR

;The radial distance from the galaxy center, deltaRA and deltaDec are
;calculated, then converted back to arcseconds.
FOR j=0,n1-1 DO BEGIN
    dx = double((cor[1,j]-RA))
    dy = double((cor[2,j]-Dec))
    rad[0,j] = j+1
    rad[1,j] = sqrt(dx^2+dy^2)*3600.
    rad[2,j] = dx*3600.
    rad[3,j] = dy*3600.
ENDFOR

indy = bsort(rad[1,*])
rads[0,*] = rad[0,indy]
rads[1,*] = rad[1,indy]
rads[2,*] = rad[2,indy]
rads[3,*] = rad[3,indy]

window,0,retain=2
device,decomposed=0
loadct,0
plot,cor[1,*],cor[2,*],psym=3,xtitle='RA',ytitle='Dec',/ynozero,$
  title='Radial positions of each fiber relative to galaxy center',$
  charsize=1.5
loadct,4
for j=0,n1-1 do xyouts,cor[1,j],cor[2,j],strn(round(rad[1,j])),color=150,charsize=1.2

openw,6,name+'.rad'
printf,6,['Fib','  Rad      ','dRA      ','dDec']
for j=0,n1-1 do printf,6,uint(rad[0,j]),rad[1,j],rad[2,j],rad[3,j],format='(I3,2x,f8.4,2x,f8.4,2x,f8.4)'
free_lun,6

openw,6,name+'.rad.std'
printf,6,['Fib','  Rad      ','dRA      ','dDec']
for j=0,n1-1 do printf,6,uint(rads[0,j]),rads[1,j],rads[2,j],rads[3,j],format='(I3,2x,f8.4,2x,f8.4,2x,f8.4)'
free_lun,6

wait,2
wdelete,0
END


; This code is used to plot the major-axis flux for an elliptical
; galaxy. It reads in a processed data frame (up to cosmic ray
; rejection), and several other files, and returns a plot of the
; galaxy as a function of major axis radius. This essentially corrects
; for ellipticity in the plotting of the flux as a function of radius.

; The first part of this routine has been stripped from
; sumtemp.pro. If you ever want to just get the counts without
; calculating the R_major values then you can run that code.

; DATAFILE: The extracted data file
; MONTH+night: The month and year the data was taken. This is used to
; locate the ptow and mask files. EX: jan08_n3
; GALRA & GALDEC: The galaxy center ex:12:34:45.6
; PA: Position angle. This must be converted to the angle in the first
; quadrant. (i.e PA: 235.0 should be entered as 235-180=55)
; AXISRATIO: The axis ratio of the galaxy = minor/major (in diameter)
; WAVE1 AND 2: This is the wavelength range over which to sum the
; fibers. If no values are provided then the default of W1=4400,
; W2=4700 is applied
; OUTNAME= The name of the output text file. If no name is given the
; datafile is parsed on the first '_' and the first piece is
; subscripted with '.txt' and used. The format of the output is as
; follows:

PRO Prad_major, datafile, month_night, galRA, galDec, PA, coordfile,$
                axisratio, W1=wave1, W2=wave2,OUTNAME=outname

;COMPILE: realignF.pro

;***********************************************************************
;raoffset = 0.00331116
;decoffset = 0.0017014
raoffset = 0.001633 ;SUBTRACT FROM HET RA AND DEC
decoffset = 0.001084 ;SUBTRACT FROM HET RA AND DEC
;offsetcorrect = 'on'
offsetcorrect = 'off'
;listswitch = 'off'
listswitch = 'on'
wi = 3550.0
disp = 1.125
;disp = 1.125 * 2.0
aperture = 5
plotting = 'on' ;change to 'on' to plot various things to screen
skysub = 'no';change to 'no' if you don't want to attempt a
                   ;sky subtraction using a list of fibers (skyfibers.txt)
;***********************************************************************

ptow = 'ptow_'+month_night+'.dat'
mask = 'mask_'+month_night+'.dat'

if(n_elements(wave1) eq 0) then wave1 = 4400.00
if(n_elements(wave2) eq 0) then wave2 = 4700.00

step = floor((aperture -1) / 2.0)

if (listswitch eq 'on') then begin
    readcol,datafile,f='a',list
    nnn = n_elements(list)
    for j=0,n_elements(list)-1 do begin
        data = readfits(list[j],header,/silent)
        n0 = n_elements(data[*,0])
        n1 = n_elements(data[0,*])
        if (j eq 0) then dataarr = fltarr(n0,n1,nnn)
        wavezp = 'on'
        cntr3 = 0
        repeat begin
            temp = strsplit(header[cntr3],' ',/extract)
            ans = temp[0]
            cntr3 = cntr3 + 1
        endrep until (ans eq 'WAVEZP')
        wzp = double(temp[2])
        print,wzp
        dataarr[*,*,j] = realignF(data,ptow,mask,wi,disp,aperture,wzp)
    endfor
    dataA = median(dataarr,dimension=3,/even)
;writefits,'het_d1.fits',dataA
;writefits,'hjs_d1.fits',dataA
;delete once this runs...
endif

if (listswitch ne 'on') then begin
    data = readfits(datafile,header,/silent)
    n0 = n_elements(data[*,0])
    n1 = n_elements(data[0,*])
    
    wavezp = 'on'
    cntr3 = 0
    repeat begin
        temp = strsplit(header[cntr3],' ',/extract)
        ans = temp[0]
        cntr3 = cntr3 + 1
    endrep until (ans eq 'WAVEZP')
    wzp = double(temp[2])
    
    dataA = realignF(data,ptow,mask,wi,disp,aperture,wzp)
endif

readcol,coordfile,format='i,a,a',fibnum,RA,Dec
n2 = n_elements(RA)

wave = dblarr(n0)
for j=0,n0-1 do wave[j] = wi + j*disp 

readcol,mask,format='i',yindex
yindex = yindex - 1

;a sloppy work-a-round is introduced to handle sky-subtraction
if (skysub eq 'yes') then begin
    readcol,'skyfibers.txt',format='i',skyfibers
    isky = skyfibers - 1
    skyarr = dblarr(n0,5,n_elements(isky))
    for j=0,n_elements(isky)-1 do begin
        skyarr[*,*,j] = dataA[*,yindex[isky[j]]-step:yindex[isky[j]]+step]
    endfor
    skyestimate = median(skyarr,dimension=3,/even)
endif

;some numbers you may not even care about
i1 = where(wave ge wave1)
i1 = i1[0]
wout1 = wave[i1]
pout1 = i1 + 1
i2 = where(wave ge wave2)
i2 = i2[0]
wout2 = wave[i2]
pout2 = i2 + 1
deltaw = wout2 - wout1
deltap = pout2 - pout1

waveP = wave[i1:i2]
n3 = n_elements(waveP)
get_lun,lun
openw,lun,'temp.txt'
form1 = '(i3,a13,a13,2x,f6.1,2x,f6.1,2x,f6.1,i5,i5,i5,i5,f10.2,f10.2)'
form2 = '(i3,a13,a13,i7,f11.2,1x,f9.2)'

for j=0,n2-1 do begin
    fn = fibnum[j]
    raP = RA[j]
    decP = Dec[j]
    dataP = dataA[i1:i2,yindex[j]-step:yindex[j]+step]
    ibad = where(dataP le -500,nbad)
    if (nbad eq n_elements(dataP)) then begin
        ttl = 0.0
        ttlN = 0.0
        goto,jump1
    endif
    if (skysub eq 'yes') then dataP = dataP - skyestimate
    if (nbad ne 0) then dataP[ibad] = 0.0
    ttl = total(dataP)
    ttlN = ttl/(n_elements(dataP)-nbad)
    jump1:
    printf,lun,fn,raP,decP,nbad,ttl,ttlN,format=form2
endfor
free_lun,lun

if (plotting eq 'on') then begin
    window,0,retain=2
    device,decomposed=0
    readcol,'temp.txt',format='x,x,x,f,f,f',nbad,ttl,ttlN
    plot,ttlN
    pause
endif


; THE START OF THE SECOND BIT THAT DOES THE R_MAJOR CALCULATION
PA = PA/(360.0/(2*!pi))
readcol,'temp.txt',format='i,a,a,i,f,f',fib,fRA,fDec,nbad,ttl,ttln

if(n_elements(fib) ne n2) then stop

RA = strn(galRA) & Dec = strn(galDec)
RA = strsplit(RA,':',/EXTRACT)
Dec = strsplit(Dec,':',/EXTRACT)
galRA = double((RA[0]+RA[1]/60.+RA[2]/3600.)*15)
;The + or - sign of the declination is accounted for
IF (Dec[0] LT 0) THEN galDec = double(Dec[0]-Dec[1]/60.-Dec[2]/3600.)
IF (Dec[0] GE 0) THEN galDec = double(Dec[0]+Dec[1]/60.+Dec[2]/3600.)

fibcoords = dblarr(3,n2) ;fiber #, fiberRA, fiberDec
fiboffsets = dblarr(4,n2) ;fiber #, radius, deltaRA, deltaDec
galoffsets = dblarr(4,n2)
Rmajor = dblarr(n2)

FOR j=0,n2-1 DO BEGIN ;A loop through each fiber
    fibcoords[0,j] = fib[j]
    fc = fra[j]
    fd = fdec[j]
    fc = strsplit(fc,':',/EXTRACT)
    fd = strsplit(fd,':',/EXTRACT)
;The coordinates are converted to decimal degrees
    fibcoords[1,j] = double((fc[0]+fc[1]/60.+fc[2]/3600.)*15.)
    if (offsetcorrect eq 'on') then fibcoords[1,j] = fibcoords[1,j] - raoffset
    IF (fd[0] LT 0) THEN fibcoords[2,j] = double(fd[0]-fd[1]/60.-fd[2]/3600.)
    IF (fd[0] GT 0) THEN fibcoords[2,j] = double(fd[0]+fd[1]/60.+fd[2]/3600.)
    if (offsetcorrect eq 'on') then fibcoords[2,j] = fibcoords[2,j] - decoffset
ENDFOR

FOR j=0,n2-1 DO BEGIN
;    dx = abs(double((fibcoords[1,j]-galRA)))
;    dy = abs(double((fibcoords[2,j]-galDec)))
    dx = double((fibcoords[1,j]-galRA))
    dy = double((fibcoords[2,j]-galDec))
    fiboffsets[0,j] = j+1
    fiboffsets[1,j] = sqrt(dx^2 + dy^2)*3600.0
    fiboffsets[2,j] = dx*3600.0
    fiboffsets[3,j] = dy*3600.0
;    Print, 'The radius for fiber '+strn(j+1)+ ' is: ',strn(rad[1,j])
ENDFOR

; convert to the galaxy coordinate system
for j=0,n2-1 do begin
    r = sqrt(fiboffsets[2,j]^2 + fiboffsets[3,j]^2)
    galoffsets[0,j] = fibcoords[0,j]
    theta = acos(fiboffsets[2,j]/r)
    phi = theta + PA + (90.0/(360.0/(2*!pi)))
    galX = r * cos(phi)
    galY = r * sin(phi)
    galoffsets[2,j] = galX
    galoffsets[3,j] = galY
    galR = sqrt(galX^2 + galY^2)
    galoffsets[1,j] = galR
    Rmajor[j] = sqrt(galX^2 + (galY^2/axisratio^2))
endfor

if (plotting eq 'on') then begin
    window,0,retain=2
    device,decomposed=0
    loadct,0
    plot,Rmajor,ttln,psym=1
    ;pause
    plot,galoffsets[1,*]/rmajor,/ynozero,psym=1,$
      title='R_galaxy / R_major'
    ;pause
    wdelete
endif

if (n_elements(outname) eq 0) then begin
    outname = strsplit(datafile,'_',/extract)
    outname = outname[0]+'.txt'
endif

name=['Fib#','RA','Dec','RA.0','Dec.0','DeltaRA','DeltaDec','Radius','GalRA','GalDec','Nbad','Totals','N_ttl','R_Maj']
form4 = '(a4,a7,a14,a11,a12,a12,a13,a8,a6,a7,a7,a11,a11,a7)'
;form3 = '(i3,a13,a13,2x,f10.5,f10.5,f11.5,f11.5,f7.2,f7.2,f7.2,i6,f14.2,f9.2,f7.2)'
form5 = '(f10.5,f10.5,f9.2)'
get_lun,lun
openw,lun,outname
;printf,lun,name,format=form4

;for j=0,n2-1 do begin
;    printf,lun,j+1,fra[j],fdec[j],fibcoords[1,j],fibcoords[2,j],fiboffsets[2,j],$
;      fiboffsets[3,j],fiboffsets[1,j],galoffsets[2,j],galoffsets[3,j],$
;      nbad[j],ttl[j],ttln[j],rmajor[j],format=form3
;endfor
for j=0,n2-1 do begin
    printf,lun,fibcoords[1,j],fibcoords[2,j],ttln[j],format=form5
endfor
free_lun,lun

stop

END


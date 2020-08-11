;This purpose of this routine is to create the final output values for
;each bin. The output of this routine is to return a text output of
;the following form:

;BIN NAME   RA    DEC    VELOCITY    DISPERSION    H3    H4

;Where all six of the values are WEIGHTED BY EITHER THE FIBER
;INTENSITY (in the case of the RA and DEC) OR THE RMS (as output from
;rfitlov). In this way, what's returned is a type of average for each
;radial bin.

;**********************************************************************
pro binoffset_vp1,data
;**********************************************************************

;This routine runs after you've created the ####.bin files. It
;requires these files (AND NO OTHERS WITH A .BIN SUFFIX!), the
;IFUcen.txt, the coords.txt and the data you're binning to be in the
;calling directory. It also needs a file called "rms.txt"- this is the
;screen output from running rfitlov, piped into a text file and edited to
;contain only the rms values... 

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;THIS CODE IS FOR VP1!!!!
;When changing this, the central fiber that makes the alignment with
;the galaxy center will change...
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;the offsets from the central fiber are read in (VP1)
readcol,'IFUcen.txt',fibzero,ra,dec,FORMAT='I,F,F'
;the actual positions on sky in RA and Dec are read in
readcol,'coords.txt',fibloc,fra,fdec,FORMAT='I,A,A'

files = findfile('*.bin',count=n1)
n3 = n_elements(ra)

data = readfits(data)

;n1 = number of bins
;n3 = number of fibers (247 or 246)

temp = data[700:720,*]
ttl = total(temp,1)
zeroindx = where(ttl eq 0)

frame = data[100:n_elements(data[*,0])-100,*]
ttls = dblarr(n3)
for j=0,n3-1 do ttls[j] = total(frame[*,j])

raoffset = fltarr(n3)
decoffset = fltarr(n3)
ra4bin = fltarr(n1)
dec4bin = fltarr(n1)

galra = strarr(1)
galdec = strarr(1)

;------------------------------------------------------
;UNCOMMENT THESE TO BE ABLE TO INPUT THE GALAXY VALUES FROM THE
;COMMAND LINE

;print,'Enter the galaxy RA (HH:MM:SS.S):'
;read,galra
;print,'Enter the galaxy Dec (DD:MM:SS.S):'
;read,galdec
;print,'Enter the galaxy PA:'
;read,pa

galra = '04:31:39.8'
galdec='-05:05:10'
pa = 5
;---------------------------------------------------------

galra = strsplit(galra,':',/extract)
galdec = strsplit(galdec,':',/extract)
f124ra = strsplit(fra[123],':',/extract)
f124dec = strsplit(fdec[123],':',/extract)
f8ra = strsplit(fra[7],':',/extract)
f8dec = strsplit(fdec[7],':',/extract)

;galaxy RA and Dec are converted into decimal
galra = double((galra[0]+galra[1]/60.+galra[2]/3600.)*15)
if (galdec[0] lt 0) then galdec = double(galdec[0]-galdec[1]/60.-galdec[2]/3600.)
if (galdec[0] ge 0) then galdec = double(galdec[0]+galdec[1]/60.+galdec[2]/3600.)

;fiber positions are converted into decimal
for j=0,n3-1 do begin
    temp = strsplit(fra[j],':',/extract)
    fra[j] = double((temp[0]+temp[1]/60.+temp[2]/3600.)*15)
    temp = strsplit(fdec[j],':',/extract)
    if (temp[0] lt 0) then fdec[j] = double(temp[0]-temp[1]/60.-temp[2]/3600.)
    if (temp[0] ge 0) then fdec[j] = double(temp[0]+temp[1]/60.+temp[2]/3600.)
endfor

;the central fiber is converted to decimal
f124ra = double((f124ra[0]+f124ra[1]/60.+f124ra[2]/3600.)*15)
if (f124dec[0] lt 0) then f124dec = double(f124dec[0]-f124dec[1]/60.-f124dec[2]/3600.)
if (f124dec[0] ge 0) then f124dec = double(f124dec[0]+f124dec[1]/60.+f124dec[2]/3600.)

f8ra = double((f8ra[0]+f8ra[1]/60.+f8ra[2]/3600.)*15)
if (f8dec[0] lt 0) then f8dec = double(f8dec[0]-f8dec[1]/60.-f8dec[2]/3600.)
if (f8dec[0] ge 0) then f8dec = double(f8dec[0]+f8dec[1]/60.+f8dec[2]/3600.)

for j=0,n3-1 do raoffset[j] = fra[j] - galra
for j=0,n3-1 do decoffset[j] = fdec[j] - galdec    

galcentra = galra-f124ra
galcentdec = galdec-f124dec

delra = (abs(ra[7]-ra[123]))^2
deldec = (abs(dec[7]-dec[123]))^2
l = sqrt(delra+deldec)

galpa = pa/(360/(2*!pi))
s1 = l*tan(galpa)
s2 = -s1

mjraxis = [[s2,-l],[s1,l]]
mnraxis =  [[-l,s1],[l,s2]]


cntr3 = 0
for k=0,n1-1 do begin
    file = files[k]
    readcol,file,bin,format='I'
    print,'This is the bin'
    print,file
    print,bin
    n2 = n_elements(bin)
    indy = bin - 1
    wgtarray = fltarr(n2)
    ran = fltarr(n2)
    decn = fltarr(n2)
    all = total(ttls[indy])
    for j=0,n2-1 do wgtarray[j] = ttls[indy[j]]/all
    for j=0,n2-1 do begin
        ran[j] = raoffset[indy[j]] * wgtarray[j]
        decn[j] = decoffset[indy[j]] * wgtarray[j]
    endfor
    ra4bin[cntr3] = total(ran)
    dec4bin[cntr3] = total(decn)
    print,ra4bin[cntr3]
    print,dec4bin[cntr3]
    cntr3 = cntr3 + 1
endfor

raoffset = raoffset*3600
decoffset = decoffset*3600
ra4bin = ra4bin*3600
dec4bin = dec4bin*3600

openw,5,'raoffset.txt'
printf,5,transpose(raoffset)
openw,6,'decoffset.txt'
printf,6,transpose(decoffset)
free_lun,5,6


readcol,'pfitlov.out',temp,vel,disp,junk,h3,h4,format='a,f,f,f,f,f'
n4 = n_elements(temp)
readcol,'rms.txt',rms,junk,format='f,f'

openw,5,'rms.out'
for j=0,n4-1 do begin
    printf,5,temp[j],rms[j]
endfor
free_lun,5

names = strarr(n4)
for j=0,n4-1 do begin 
    jnk = strsplit(temp[j],'f',/EXTRACT)
    names[j] = jnk[0]
endfor


if (n_elements(rms) ne n4) then begin
    print,'You do not have the correct number of RMS files!'
    stop
endif

n5 = n4/n1
s1=fltarr(n5)
s2 = s1 & s3 = s1 & s4 = s1
rmswgt = fltarr(n5)
weights = fltarr(n5,n1)
out = fltarr(7,n1)

;the weight array is calculated and parts of the final array are set
;into place
cntr = 0
for k=0,n4-1,n5 do begin
    chunk = rms[k:k+4]
    chunknames = names[k:k+4]
    chunkvel = vel[k:k+4]
    chunkdisp = disp[k:k+4]
    chunkh3 = h3[k:k+4]
    chunkh4 = h4[k:k+4]
    for j=0,n5-1 do rmswgt[j] = 1/chunk[j]
    for j=0,n5-1 do begin
        t = total(rmswgt)
        weights[j,cntr] = rmswgt[j]/t
    endfor
    for j=0,n5-1 do begin
        s1[j] = weights[j,cntr]*chunkvel[j]
        s2[j] = weights[j,cntr]*chunkdisp[j]
        s3[j] = weights[j,cntr]*chunkh3[j]
        s4[j] = weights[j,cntr]*chunkh4[j]
    endfor
    out[0,cntr] = chunknames[0]
    out[1,cntr] = ra4bin[cntr]
    out[2,cntr] = dec4bin[cntr]
    out[3,cntr] = total(s1)
    out[4,cntr] = total(s2)
    out[5,cntr] = total(s3)
    out[6,cntr] = total(s4)
    cntr = cntr + 1
endfor

nm = ['Bin','RA','Dec','Velocity','Dispersion','H_3','H_4']
openw,5,'binvalues.out'
for j=0,n1-1 do begin
    bn = fix(out[0,j])
printf,5,bn,out[1,j],out[2,j],out[3,j],out[4,j],out[5,j]
endfor

free_lun,5

set_plot,'ps'
device,filename='bin_location.ps',/color
loadct,0
plot,out[1,*],out[2,*],psym=4,xtitle='RA (arcsec)',ytitle='Dec (arcsec)',$
  title='Radial Bins and Galaxy Axes',xthick=2,ythick=2,thick=2,charthick=2
loadct,4
oplot,mjraxis[0,*],mjraxis[1,*],thick=2,color=150
oplot,mnraxis[0,*],mnraxis[1,*],linestyle=2,thick=2,color=150
for j=0,n1-1 do begin
    s = fix(out[0,j])
    s = ' '+strn(s)
    xyouts,out[1,j],out[2,j],s,charthick=2,charsize=0.7
endfor
device,/close_file
set_plot,'x'

print,'SUCCESS!'
stop
end

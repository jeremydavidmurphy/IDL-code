pro makebins, name, maskname 

nprofile = 5

; NAME should not include the pefsm.fits subscript. It should be
; everything up to the _pefsm.fits and _pefw.fits subscript.

;THIS CODE IS USED AS A VISUAL AID IN DETERMINING THE BINNING REGIONS
;WHEN THE POINTING IS CENTERED IN THE GALAXY. IT SHOWS THE MAJOR AND
;MINOR AXIS IF THE GALAXY, WHICH FIBERS ARE TAKEN UP IN OTHER BINS,
;AND WHICH ARE EITHER DEAD OR REJECTED. AFTER ENTERING THE FIBERS FOR
;A GIVEN REGION, THEY ARE WRITTEN OUT TO A TEXT FILE CALLED
;'NAME.BIN'. A FITS FILE IS ALSO CREATED BY A SIMPLE SUM OF THE
;SPECTRA FROM EACH FIBER INCLUDED IN THE BIN (there's an assumption
;here that the spectra fed in is fully reduced)

;The "name" file above is your fully reduced data frame. It's used to
;determine which fibers are not in use.

;It requires that IFUcen.txt and the corresponding coordinate file
;exist in the calling directory, renamed as coord.txt.

;It also will open your final, combined data file and place "X"s over
;all fibers for which there is no data (either dead or contaminated
;with foreground/background objects.)

; Modified on 07/28/09: The code now reads in the uncollapsed data and
; the weight file. The code now expects the NAME to be the name
; of the file, WITHOUT any subscript. The collapse is now done within
; the routine wgtcollapse.pro

; collapsed data = wgtcollapse(data,weight,maskname,nprofile)
; NAME: the actual uncollapsed data (pefsm)
; WEIGHT: the actual weight (pefw) file
; MASKNAME: the name of the mask file.
; NPROFILE: the number of fibers (5 or 7) in the extracted data.

readcol,'IFUcen.txt',fibzero,ra,dec,FORMAT='I,F,F'
n1 = n_elements(fibzero)
readcol,'coords.txt',fibloc,fra,fdec,FORMAT='I,A,A'

dname = name+'_pefsm.fits'
wname = name+'_pefw.fits'
data = readfits(dname,/silent)
weight = readfits(wname,/silent)

coll = wgtcollapse(data,weight,maskname,nprofile)
writefits,'testC.fits',coll

coll[where(coll eq -666)] = 0.0
n2 = n_elements(coll[*,0])
frame = coll[100:n2-100,*]
ttl = total(frame,1)
zeroindx = where(ttl eq 0)

ttls = dblarr(n1) ; The approximate totals for each fiber...
for j=0,n1-1 do ttls[j] = total(frame[*,j])

raoffset = fltarr(n1)
decoffset = fltarr(n1)
ra4bin = fltarr(n1);these will be cut down later to match the number of bins
dec4bin = fltarr(n1)

galra = strarr(1)
galdec = strarr(1)
pa = fltarr(1)
nangbins = 1

;------------------------------------------------------
;UNCOMMENT THESE TO BE ABLE TO INPUT THE GALAXY VALUES 

;print,'Enter the galaxy RA (HH:MM:SS.S):'
;read,galra
;print,'Enter the galaxy Dec (DD:MM:SS.S):'
;read,galdec
;print,'Enter the galaxy PA:'
;read,pa

;galra = '04:31:39.8' (N1600)
;galdec='-05:05:10'   (N1600)
;pa = 5               (N1600)
pa = 165.0
galra = '12:29:46.7'
galdec = '08:00:02'
;---------------------------------------------------------

galra = strsplit(galra,':',/extract)
galdec = strsplit(galdec,':',/extract)
f124ra = strsplit(fra[123],':',/extract)
f124dec = strsplit(fdec[123],':',/extract)
f8ra = strsplit(fra[7],':',/extract)
f8dec = strsplit(fdec[7],':',/extract)

galra = double((galra[0]+galra[1]/60.+galra[2]/3600.)*15)
if (galdec[0] lt 0) then galdec = double(galdec[0]-galdec[1]/60.-galdec[2]/3600.)
if (galdec[0] ge 0) then galdec = double(galdec[0]+galdec[1]/60.+galdec[2]/3600.)

for j=0,n1-1 do begin
    temp = strsplit(fra[j],':',/extract)
    fra[j] = double((temp[0]+temp[1]/60.+temp[2]/3600.)*15)
    temp = strsplit(fdec[j],':',/extract)
    if (temp[0] lt 0) then fdec[j] = double(temp[0]-temp[1]/60.-temp[2]/3600.)
    if (temp[0] ge 0) then fdec[j] = double(temp[0]+temp[1]/60.+temp[2]/3600.)
endfor

f124ra = double((f124ra[0]+f124ra[1]/60.+f124ra[2]/3600.)*15)
if (f124dec[0] lt 0) then f124dec = double(f124dec[0]-f124dec[1]/60.-f124dec[2]/3600.)
if (f124dec[0] ge 0) then f124dec = double(f124dec[0]+f124dec[1]/60.+f124dec[2]/3600.)

f8ra = double((f8ra[0]+f8ra[1]/60.+f8ra[2]/3600.)*15)
if (f8dec[0] lt 0) then f8dec = double(f8dec[0]-f8dec[1]/60.-f8dec[2]/3600.)
if (f8dec[0] ge 0) then f8dec = double(f8dec[0]+f8dec[1]/60.+f8dec[2]/3600.)

for j=0,n1-1 do raoffset[j] = fra[j] - galra
for j=0,n1-1 do decoffset[j] = fdec[j] - galdec    

galcentra = galra-f124ra
galcentdec = galdec-f124dec
print,galcentra,galcentdec

delra = (abs(ra[7]-ra[123]))^2
deldec = (abs(dec[7]-dec[123]))^2
l = sqrt(delra+deldec)

galpa = pa/(360/(2*!pi))
s1 = l*tan(galpa)
s2 = -s1

mjraxis = [[s2,-l],[s1,l]]
mnraxis =  [[-l,s1],[l,s2]]

;print,'How many angular bins?'
;read,nangbins
nangbins=12
ang = 360/nangbins
angbins = fltarr(nangbins)
radeus = angbins
npa = 90-pa
for j=0,nangbins-1 do begin
    angbins[j] = (npa+(j*ang))/(360/(2*!pi))
    radeus[j] = l+l
endfor

set_plot,'x'

window,0,retain=2,xsize=600,ysize=550
plot,ra,dec,/nodata,title='Bin Selection ("*"s are bad fibers)',$
  xtitle='RA (arcsec)',ytitle='Dec (arcsec)',charsize=1.5

cntr = 0
for j=0,n1-1 do begin
    nm = strn(fibzero[j])
    xyouts,RA[j]-2.5,Dec[j],nm,charsize=1.2
    if (j eq zeroindx[cntr]) then begin
        plots,ra[j],dec[j],psym=2,color=255,symsize=2.5,thick=2
        if (n_elements(zeroindx)-1 gt cntr) then cntr = cntr + 1
    endif
endfor

oplot,mjraxis[0,*],mjraxis[1,*],thick=2,color=255
oplot,mnraxis[0,*],mnraxis[1,*],linestyle=2,thick=2,color=255
for j=0,19 do oplot,radeus-(j*5),angbins,psym=1,symsize=0.5,color=255,/polar
for j=0,7 do oplot,radeus-(j*15),angbins,psym=-1,symsize=0.5,color=255,/polar
ans1 = ''
biname=''
new = 1

cntr3 = -1
repeat begin
array = 1
print,'Enter "n" for a new bin, "c" to clear, "a" for another bin and "d" when done:'
read,ans1
;------------------------------------------------------------------
;CLEAR
;------------------------------------------------------------------
if (ans1 eq 'c') then begin
    cntr = 0
    wset,0
    plot,ra,dec,/nodata,title='Bin Selection ("*"s are bad fibers)',$
      xtitle='RA (arcsec)',ytitle='Dec (arcsec)',charsize=1.5
    for j=0,n1-1 do begin
        nm = strn(fibzero[j])
        xyouts,RA[j]-2.5,Dec[j]-1,nm,charsize=1.2
        if (j eq zeroindx[cntr]) then begin
            plots,ra[j],dec[j],psym=2,color=255,symsize=2.5,thick=2
            if (n_elements(zeroindx)-1 gt cntr) then cntr = cntr + 1
        endif
    endfor
    oplot,mjraxis[0,*],mjraxis[1,*],thick=2,color=255
    oplot,mnraxis[0,*],mnraxis[1,*],linestyle=2,thick=2,color=255
    for j=0,19 do oplot,radeus-(j*5),angbins,psym=1,symsize=0.5,color=255,/polar
    for j=0,7 do oplot,radeus-(j*15),angbins,psym=-1,symsize=0.5,color=255,/polar
endif
;------------------------------------------------------------------
;NEW
;------------------------------------------------------------------
if (ans1 eq 'n') then begin
    cntr = 0
    array = intarr(100)
    ans2 = 1
    repeat begin
        cntr2 = 0
        print,'Enter a fiber number, or "0" when done:'
        read,ans2
        if (ans2 eq 0) then goto,jump1
        array[cntr] = ans2
        cntr = cntr + 1
;quick plot occurs
        ind = array[where(array ne 0)]
        ind = ind-1
        piece = dblarr(n2,n_elements(ind))
        for j=0,n_elements(ind)-1 do piece[*,j] = coll[*,ind[j]]
        if n_elements(ind) eq 1 then tempspec = piece else begin
            print,'BIWEIGHT!'            
            wgts = dblarr(n_elements(ind))
            for l=0,n_elements(ind)-1 do wgts[l] = median(piece[*,l])
            wgts = wgts/max(wgts)
            print,wgts & pause
            tempspec = biweight(piece,N=wgts)
;            tempspec = total(piece,2)
        endelse
        junk = qkplt(tempspec)
;------------------------------------------------------------------
;replot
;------------------------------------------------------------------
        wset,0
        plot,ra,dec,/nodata,title='Bin Selection ("*"s are bad fibers)',$
          xtitle='RA (arcsec)',ytitle='Dec (arcsec)',charsize=1.5
        for j=0,n1-1 do begin
            nm = strn(fibzero[j])
            xyouts,RA[j]-2.5,Dec[j]-1,nm,charsize=1.2
            if (j eq zeroindx[cntr2]) then begin
                plots,ra[j],dec[j],psym=2,color=255,symsize=2.5,thick=2
                if (n_elements(zeroindx)-1 gt cntr2) then cntr2= cntr2 + 1
            endif
        endfor
        oplot,mjraxis[0,*],mjraxis[1,*],thick=2,color=255
        oplot,mnraxis[0,*],mnraxis[1,*],linestyle=2,thick=2,color=255
        for j=0,19 do oplot,radeus-(j*5),angbins,psym=1,symsize=0.5,color=255,/polar
        for j=0,7 do oplot,radeus-(j*15),angbins,psym=-1,symsize=0.5,color=255,/polar
        for j=0,n_elements(ind)-1 do begin
            pra = ra[array[j]-1]
            pdec = dec[array[j]-1]
            plots,pra,pdec,psym=6,color=255,symsize=2.5,thick=2
        endfor

        jump1:
    endrep until (ans2 eq 0)

    if (n_elements(array[where(array eq 0)]) eq 100) then begin
        ans1 = 'd'
        goto,jump3
    endif

    howmany = n_elements(array[where(array ne 0)])
    array = array[0:howmany-1]
    print,'Name the bin:'
    print,'(The fits file will be automatically subscripted...)'
    read,biname
    txt = biname
    fts = biname+'.fits'
    openw,5,txt
    printf,5,transpose(array)
    free_lun,5
    writefits,fts,tempspec
    if (new eq 1) then oldarr = array else oldarr = [oldarr,array]
    new = 2

;the weighted ra and dec positions for the bin are determined.
    indy = array - 1
    wgtarray = fltarr(howmany)
    ran = fltarr(howmany)
    decn = fltarr(howmany)
    all = total(ttls[indy])
    for k=0,howmany-1 do wgtarray[k] = ttls[indy[k]]/all
    for k=0,howmany-1 do begin
        ran[k] = raoffset[indy[k]] * wgtarray[k]
        decn[k] = decoffset[indy[k]] * wgtarray[k]
    endfor
    cntr3 = cntr3 + 1
    ra4bin[cntr3] = total(ran)
    dec4bin[cntr3] = total(decn)
    print,ra4bin[cntr3]
    print,dec4bin[cntr3]
endif

;ANOTHER----------------------------------------------------------------------

if (ans1 eq 'a') then begin
;the existing bins are marked...
    cntr = 0
    wset,0
    plot,ra,dec,/nodata,title='The "X"s are already within another bin',$
      xtitle='RA (arcsec)',ytitle='Dec (arcsec)',charsize=1.5
    for j=0,n1-1 do begin
        nm = strn(fibzero[j])
        xyouts,RA[j]-2.5,Dec[j]-1,nm,charsize=1.2
        if (j eq zeroindx[cntr]) then begin
            plots,ra[j],dec[j],psym=2,color=255,symsize=2.5,thick=2
            if (n_elements(zeroindx)-1 gt cntr) then cntr = cntr + 1
        endif
    endfor
    oplot,mjraxis[0,*],mjraxis[1,*],thick=2,color=255
    oplot,mnraxis[0,*],mnraxis[1,*],linestyle=2,thick=2,color=255
    for j=0,19 do oplot,radeus-(j*5),angbins,psym=1,symsize=0.5,color=255,/polar
    for j=0,7 do oplot,radeus-(j*15),angbins,psym=-1,symsize=0.5,color=255,/polar
    for j=0,n_elements(oldarr)-1 do begin
        pra = ra[oldarr[j]-1]
        pdec = dec[oldarr[j]-1]
        plots,pra,pdec,psym=7,color=255,symsize=2.5,thick=2
    endfor

    array = intarr(100)
    ans2 = intarr(1)
    cntr = 0
    repeat begin
        cntr2 = 0
        print,'Enter a fiber number, or "0" when done:'
        read,ans2
        if (ans2 eq 0) then goto,jump2
        array[cntr] = ans2
        cntr = cntr + 1
        pra = ra[ans2-1]
        pdec = dec[ans2-1]
        wset,0
        plots,pra,pdec,psym=6,color=255,symsize=2.5,thick=2
;a quick plot of the spectra occurs
        ind = array[where(array ne 0)]
        ind = ind-1
        piece = dblarr(n2,n_elements(ind))
        for j=0,n_elements(ind)-1 do piece[*,j] = coll[*,ind[j]]
        if n_elements(ind) eq 1 then tempspec = piece else begin
            print,'BIWEIGHT!'
            wgts = dblarr(n_elements(ind))
            for l=0,n_elements(ind)-1 do wgts[l] = median(piece[*,l])
            wgts = wgts/max(wgts)
            print,wgts & pause
            tempspec = biweight(piece,N=wgts)
        endelse
;            tempspec = total(piece,2)
        junk = qkplt(tempspec)
       
;replot
        wset,0
        plot,ra,dec,/nodata,title='Bin Selection ("*"s are bad fibers)',$
          xtitle='RA (arcsec)',ytitle='Dec (arcsec)',charsize=1.5
        for j=0,n1-1 do begin
            nm = strn(fibzero[j])
            xyouts,RA[j]-2.5,Dec[j]-1,nm,charsize=1.2
            if (j eq zeroindx[cntr2]) then begin
                plots,ra[j],dec[j],psym=2,color=255,symsize=2.5,thick=2
                if (n_elements(zeroindx)-1 gt cntr2) then cntr2= cntr2 + 1
            endif
        endfor
        oplot,mjraxis[0,*],mjraxis[1,*],thick=2,color=255
        oplot,mnraxis[0,*],mnraxis[1,*],linestyle=2,thick=2,color=255
        for j=0,19 do oplot,radeus-(j*5),angbins,psym=1,symsize=0.5,color=255,/polar
        for j=0,7 do oplot,radeus-(j*15),angbins,psym=-1,symsize=0.5,color=255,/polar
        for j=0,n_elements(ind)-1 do begin
            pra = ra[array[j]-1]
            pdec = dec[array[j]-1]
            plots,pra,pdec,psym=6,color=255,symsize=2.5,thick=2
        endfor
        for j=0,n_elements(oldarr)-1 do begin
            pra = ra[oldarr[j]-1]
            pdec = dec[oldarr[j]-1]
            plots,pra,pdec,psym=7,color=255,symsize=2.5,thick=2
        endfor

        jump2:
    endrep until (ans2 eq 0)

    if (n_elements(array[where(array eq 0)]) eq 100) then begin
        ans1 = 'd'
        goto,jump3
    endif

    howmany = n_elements(array[where(array ne 0)])
    array = array[0:howmany-1]
    print,'Name the bin- it will automatically be subscripted by ".bin":'
    read,biname
    txt = biname+'.bin'
    fts = biname+'.fits'
    openw,5,txt
    printf,5,transpose(array)
    free_lun,5
    writefits,fts,tempspec
    oldarr = [oldarr,array]

;the weighted ra and dec positions for the bin are determined.
    indy = array - 1
    wgtarray = fltarr(howmany)
    ran = fltarr(howmany)
    decn = fltarr(howmany)
    all = total(ttls[indy])
    for k=0,howmany-1 do wgtarray[k] = ttls[indy[k]]/all
    for k=0,howmany-1 do begin
        ran[k] = raoffset[indy[k]] * wgtarray[k]
        decn[k] = decoffset[indy[k]] * wgtarray[k]
    endfor
    cntr3 = cntr3 + 1
    ra4bin[cntr3] = total(ran)
    dec4bin[cntr3] = total(decn)
    print,ra4bin[0:cntr3]
    print,dec4bin[0:cntr3]

endif

jump3:
endrep until (ans1 eq 'd')

ra4bin = ra4bin[0:cntr3]
dec4bin = dec4bin[0:cntr3]

openw,5,'ra_offsets.txt'
printf,5,transpose(ra4bin)
free_lun,5
openw,5,'dec_offsets.txt'
printf,5,transpose(dec4bin)
free_lun,5
print,'The RA and Dec offsets have been written.'
print,''

print,'The last version of the display will be saved as "binselect.ps".'
set_plot,'ps'
device,filename='binselect.ps',/color

plot,ra,dec,/nodata,title='Bin Selection ("*"s are bad fibers)',$
  xtitle='RA (arcsec)',ytitle='Dec (arcsec)',charsize=1.2,$
  xthick=2,ythick=2,charthick=2
;the numbers
for j=0,n1-1 do begin
    nm = strn(fibzero[j])
    xyouts,ra[j]-2.5,dec[j]-1,nm,charsize=0.8,charthick=2
endfor
;the bad fibers
for j=0,n_elements(zeroindx)-1 do begin
    ind = zeroindx[j]
    plots,ra[ind],dec[ind],psym=2,symsize=1.5,thick=2
endfor
;the galaxy axes
oplot,mjraxis[0,*],mjraxis[1,*],thick=2
oplot,mnraxis[0,*],mnraxis[1,*],linestyle=2,thick=2
;the angular bins
for j=0,19 do oplot,radeus-(j*5),angbins,psym=1,symsize=0.5,color=255,/polar
for j=0,7 do oplot,radeus-(j*15),angbins,psym=-1,symsize=0.5,color=255,/polar
;the binned fibers
for j=0,n_elements(oldarr)-1 do begin
    pra = ra[oldarr[j]-1]
    pdec = dec[oldarr[j]-1]
    plots,pra,pdec,psym=7,symsize=1.5,thick=2
endfor
device,/close_file
set_plot,'x'

print,'Next ENTER deletes the plot:'
pause
wdelete,2,0

stop
end

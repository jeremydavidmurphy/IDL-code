;
;+
;
;PURPOSE:
;
;      This procedure produces collapsed multi-extension data 'cubes'
;      for VIRUS-P IFU data for a single frame or interpolates over a
;      common wavelength grid for multiple data frames.  This program
;      assumes the data is flux calibrated in units of ergs/s/cm^2/A
;      (output from fluxCal_refImg.pro)
;
;
;CALLING SEQUENCE:
;      mkvpcube, inputfile, output
;INPUTS:
;       inputfile - file that lists where the input files are located
;                   (see example below)
;
;
;OUTPUTS:
;       output - multi-extension data file that includes:
;
;       1.  Flux in ergs/s/cm^2/A
;       2.  Error in the Flux in ergs/s/cm^2/Ang
;       3.  Wavelength of each pixel in Ang
;       4.  Coordinates (Deg)
;       5.  Instrumental dispersion/FWHM in Ang
;
;KEYWORDS:
;        
;      BIN2X1 - Set this keyword if data was binned 2x1, otherwise
;               code assumes it data has 1x1 binning.
;
;      LOGBIN - Set this keyword if you would like to output a
;               datacube with natural log (wavelength)
;
;      MULTI  - Set this keyword if you have multiple frames to
;               collapse and interpolate onto a single wavelength
;               grid. If /MULTI is not set, program will assume you
;               want to collapse a single frame.
;
;
;EXAMPLE:
;
;    inputfile should have the following:
;
;       wavelength solution (ptowfile) for data frame 1
;               "       "        "         "    frame 2
;               "       "        "         "    frame 3  
;       wavelength solution (ptowfile) for data frame N
;       coordinate file (coordfile) for object dither frame
;       instrument fwhm (instresfile) file
;       location of data frame 1
;         "          "   frame 2
;         "          "   frame 3
;       location of data frame N
;
;       an example input file for 3 data frames taken on galaxy ARP299
;        called: '~/mkcubes/mkvpcube_arp299p1d1r.param' will look like
;        the following:
;
;       ~/vaccine/proc_data/ptow_AHMay09_n2.dat
;       ~/vaccine/proc_data/ptow_AHMay09_n2.dat
;       ~/vaccine/proc_data/ptow_AHMay09_n4.dat
;       ~/mkcubes/coords/ARP299_P1_D1_coords.txt
;       ~/mkcubes/vpinstfwhm_Jan09.fits
;       ~/vaccine/proc_data/AH0136
;       ~/vaccine/proc_data/AH0138
;       ~/vaccine/proc_data/AH0360
;
;
;      Run mkvpcube on ARP 299 that was observed using the 4580-6800
;      Ang red set up:
;
;      IDL>   mkvpcube, '~/mkcubes/mkvpcube_arp299p1d1r.param',$
;             '~/mkcubes/data/ARP299.p1d1'
;
;      the output fits file ARP299.p1d1.fits (note that '.fits' will be
;      appended) will be the multi-extension fits file for dither 1 of
;      pointing 1.
;
;
;MODIFICATION HISTORY: Written by G. Blanc 2009. Updated by
;       A. Heiderman November 2010 to use on data observed on
;       different dates and calibration files (wavelength solutions,
;       flux calibration), scaled frames and errors before averaging,
;       corrected for heliocentric velocity, updated header parameters
;       to include units and wavelength array.  AH also added
;       documentation. Mostly rewritten by A. Heiderman February 2011,
;       took out airmass correction, flux calibration, tracked bad
;       pixels throughout with -666 values, updated for general use on
;       a single or multiple frames.
;
;-
;
;


pro mkvpcube, file, output, BIN2X1=BIN2X1, LOGBIN=LOGBIN, MULTI=MULTI



;READ PARAMETER FILE FOR ONE DITHER
; if the parameter file is for an object observed on the same night
nfiles=file_lines(file)
nframes=(nfiles-2)/2 ; number of lines in param file minus fwhm + coordfile /2 = number of frames

;openr, lun, file, /get_lun
readcol, file, paramlist,format=('a')
ptowfile=strarr(nframes)
for beer=0, nframes-1 do begin
    ptowname=paramlist[beer]
    ptowfile[beer]=ptowname
endfor
coordfile=paramlist[nframes]
instresfile=paramlist[nframes+1]

frames=strarr(nframes)
ind_frame=fix(nframes+2)

for beer=0, nframes-1 do begin
    beer2=beer+ind_frame
    framename=paramlist[beer2]
    frames[beer]=framename
    beer2=beer+1
endfor

;READ INSTRUMENTAL RESOLUTION FILE

fwhm=mrdfits(instresfile, 0, h)

;READ PTOW FILE, CHECKING FOR ORDER OF WAVELENGTH SOLUTION

ncols=intarr(nframes)
nfib=intarr(nframes)

for moon=0, nframes-1 do begin
    openr, lun, ptowfile[moon], /get_lun
    line=''
    readf, lun, line
    free_lun, lun
    ncols1=fix(n_elements(strsplit(line, /regex, /extract))) ;order of polynomial for wavelength solution
    ncols[moon]=ncols1
    nfib1=fix(file_lines(ptowfile[moon]))
    nfib[moon]=nfib1
endfor
if where(uniq(ncols)) eq 0l then begin 
    nfib=nfib[0]
endif else begin
    nfib=max(nfib)
endelse


if where(uniq(ncols)) eq 0l then begin 
    ncols=ncols[0]
    acoef=dblarr(ncols,nfib,nframes) ; array with the wavelength solution coefficients for each fiber
endif else begin
    ncols=max(ncols)
    acoef=dblarr(ncols,nfib,nframes)
endelse
acoef1=dblarr(ncols,nfib) 
for moon=0, nframes-1 do begin
    openr, lun1, ptowfile[moon], /get_lun
    readf, lun1, acoef1
    free_lun, lun1
    acoef[*,*,moon]=acoef1
endfor


;READ COORDINATES FILE

bb=strarr(nfib)
openr, lun, coordfile, /get_lun
readf, lun, bb
free_lun, lun
aux=strarr(3,nfib)
for i=0, nfib-1 do aux[*,i]=strsplit(bb[i], /regex, /extract)
coords=dblarr(2,nfib)
for i=0, nfib-1 do begin get_coords, junk, instring=aux[1,i]+' '+aux[2,i] & coords[*,i]=junk & endfor
coords[0,*]=15d*coords[0,*] ;so both ra and dec are in decimal degrees


;READ DATA CHECKING FOR BINNING, GET DATE, UT from header for velocity correction

month=intarr(nframes)
year=intarr(nframes)
day=intarr(nframes)
ut=fltarr(nframes)
for moon=0, nframes-1 do begin
    iframe=mrdfits(frames[moon]+'pefsma.fits', 0, hdrf)
    date1=sxpar(hdrf, 'DATE-OBS')
    month[moon]=strmid(date1,5,2)
    year[moon]=strmid(date1, 0,4)
    day[moon]=strmid(date1,8,2)
    ut1=sxpar(hdrf, 'UT')
    ut[moon]=strmid(ut1,0,2)+strmid(ut1,3,2)/60d +  strmid(ut1,6,5)/3600d
endfor

aux=mrdfits(frames[0]+'pefsma.fits', 0, hdrf1)
nx=(size(aux))[1]
ny=(size(aux))[2]

data=dblarr(nx,ny,nframes)
weight=dblarr(nx,ny,nframes)
error=dblarr(nx,ny,nframes)
airmass=dblarr(nx,ny,nframes)
exptime=dblarr(nx,ny,nframes)
for i=0, nframes-1 do begin 
    data[*,*,i]=mrdfits(frames[i]+'pefsma.fits', 0, hdr) ;images for 3 dithers
    airmass[*,*,i]=sxpar(hdr, 'airmass')
    exptime[*,*,i]=sxpar(hdr, 'EXPTIME')
    weight[*,*,i]=mrdfits(frames[i]+'pefwea.fits', 0, hdr2)
endfor


nogap=where(weight ne 0)

weight[nogap]=1.0/weight[nogap]^2 ; Mimi's code now outputs error, not weights!!!!

error[nogap]=1.0/sqrt(weight[nogap]) ;error frames, where weight is 1/noise^2 or 1/error^2

;CREATE WAVELENGTH MAP


lambdaarr=dblarr(nx, nfib,nframes)

   
for moon=0, nframes-1 do begin
    xarr=dindgen(nx+1)
    lambda0=dblarr(nx+1, nfib)
    lambda=dblarr(nx, nfib)
for i=0, nfib-1 do begin & for j=0, ncols-1 do lambda0[*,i]=lambda0[*,i]+acoef[j,i,moon]*xarr^j & endfor
    
for i=0, nx-1 do lambda[i,*]=(lambda0[i,*]+lambda0[i+1,*])/2. ;place lambda at center of pixel

    lambdaarr[*,*,moon]=lambda
endfor



;----------- SCALE FRAMES---------- (prints scaling factor, not applied)


mframe=dblarr(nframes)
for beer=0, nframes-1 do begin
    dat=data[*,*,beer]
    errdat=error[*,*,beer]
    indat=where(dat ne -666 and dat ne 0.0)
    med=median(dat[indat], /double)
    mframe[beer]=med
endfor
bestframe=where(mframe eq max(mframe))
medbest=mframe[bestframe]
medbest=medbest[0]
;scales=medbest/mframe
scales=mframe/median(mframe)
print, 'scaling factors between frames are:  ',scales, '  for frames 1 to', nframes


;make masks to preserve the -666 flag for bad pixels (1.0 for good data):

msk=fltarr(nx,ny,nframes)+1.
;preserve flags:
for moon=0, nframes-1 do begin
    dat=data[*,*,moon]
    datindfa=where(dat eq -666)
    datind2dfa=array_indices(dat,datindfa) ;2d indices
    msk[datind2dfa[0,*],datind2dfa[1,*],moon]=-666
endfor




;---------- CORRECT WAVELENGTH ARRAY FOR BARYCENTRIC VELOCITY ----------

lambdac=dblarr(nx,nfib,nframes)
for moon=0, nframes-1 do begin

jdcnv, year[moon], month[moon], day[moon], UT[moon], jd  ;get julian date at time object was observed(roughly for UT in middle of the night, UT=6 hr)   
baryvel, jd, 2000, vhel, vbary  ;compute heliocentric & barycentric velocity (in km/s) at jd
;extract RA & DEC of object from image frame header:
ra1=sxpar(hdr,'RA')
dec1=sxpar(hdr,'DEC')

ra2 = ten(float(strmid(ra1, 0,2)), float(strmid(ra1, 3,2)),$
          float(strmid(ra1, 6,5)))*15./!RADEG ;RA  of object in radians
if strmid(dec1,0,1) eq '-' then begin
    dec2 = (-1.)*ten(float(strmid(dec1, 1,2)),float(strmid(dec1, 4,2)),$
                     float(strmid(dec1, 7,4)))/!RADEG  
endif else begin
    dec2 = ten(float(strmid(dec1, 1,2)),float(strmid(dec1, 4,2)),$
               float(strmid(dec1, 7,4)))/!RADEG  
endelse
;velocity correction in km/s: 
vcor = vbary[0]*cos(dec2)*cos(ra2) + $ ;Project velocity toward star
  vbary[1]*cos(dec2)*sin(ra2) + vbary[2]*sin(dec2) 
cc=2.99792458d5                 ; speed of light in km/s

lambdac[*,*,moon]= (vcor/cc)*lambdaarr[*,*,moon] ;barycentric wavelength correction

lambdaarr[*,*,moon]=lambdac[*,*,moon]+ lambdaarr[*,*,moon] ; corrected wavelength array

endfor



;REBIN WAVELENGTH ARRAYS TO A DELTA_LAMBDA=1.1 Ang (1X1 BINNING); 2.2
;Ang (2X1 BINNING) FOR A RANGE OF 4580-6830 (VENGA DATA) 4630-6850 (VIXENS DATA)

if keyword_set(bin2x1) then begin
    delta_lambda=2.2            ;Ang
endif else begin
    delta_lambda=1.1            ;Angstroms 1x1 binning
endelse


if keyword_set(multi) then begin
;lambda_all=dblarr[
    lambda_all=[lambdaarr[*,*,0],lambdaarr[*,*,1],lambdaarr[*,*,2]] 
endif else begin 
    lambda_all=lambdaarr[*,*,0] 
endelse

if keyword_set(multi) then begin
    minn=dblarr(nframes)
    maxx=dblarr(nframes)
    for moon=0, nframes-1 do begin ; find minimum (bluest) wavelength to start wavelength array
        ll=lambdaarr[*,*,moon]
        tnozero=where(ll ne 0.0)
        mn=min(ll[tnozero])
        mx=max(ll[tnozero])
        minn[moon]=mn
        maxx[moon]=mx
    endfor
    lambdastart=float(floor(min(minn)))
    lambdaend=float(ceil(max(maxx)))


;define number of pixels by range of wavelength in all frames:
    npix=round((lambdaend-lambdastart)/delta_lambda)

    print, 'wavelength range:  ', lambdastart, ' to ',lambdastart+(npix-1)*delta_lambda, '  Ang with N_pixels =', npix
;make new wavelength array for each fiber starting from bluest
;wavelength in data set
    lambdanew=fltarr(npix,nfib)
    for beer=0, npix-1 do begin
        wavel=lambdastart+delta_lambda*beer
        lambdanew[beer,*]=wavel
    endfor

    if keyword_set(logbin) then begin
        for i=0, nfib-1 do  lambdanew[*,i]=exp(range(alog(lambdanew[0,i]), alog(lambdanew[npix-1,i]), npix))  
        c=299792.458 ; speed of light in km/s
        velscale=(alog(lambdanew[npix-1,0])-alog(lambdanew[0,0]))/npix*c
    endif
        
;COMBINE DATA AND ERRORS 



    naper=5                     ; extraction aperture height in pixels
    nval=naper


;INTERPOLATE EACH FRAME TO NEW WAVELENGTH SCALE:
;data values are assumed to be flux calibrated
;stop
    nynew=naper*nfib ;-> make a new data file interpolated to new wavelength with just extracted fiber values (naper=number of pixels to collapose for each fiber pixel)
    fluxnew=dblarr(npix,nynew,nframes)
    fluxerrornew=dblarr(npix,nynew,nframes)
    msknew=intarr(npix,nynew,nframes)
    j=0
    print, 'Interpolating flux, flux errors, and bad pixel masks...'
    wait, 1
    for beer=0, nframes-1 do begin
        print, 'Frame number  ', (beer+1)
        wait, 1
        j=0
        for i=0, nfib-1 do begin
            print, "FIBER: ", (i+1)
            ymin=(i+1)*8-2-1    ; fibers are 8 pixels apart
                                ;ymax=(i+1)*8+2-1
;INTERPOLATE FLUX TO NEW WAVELENGTH SCALE:
            for moon=0, naper-1 do begin
                datanow=data[*,ymin+moon,beer]
                bad1=where(datanow eq -666, complement=good00)
                if bad1[0] ne -1l then datanow[bad1]=!Values.F_NAN ;flag bad pixels with NaN so to preserve them in interpolation
;                fluxnew1=rebin_spectrum(datanow,lambdaarr[*,i,beer],lambdanew[*,i])
                auxl0=lambdaarr[*,i,beer]
                auxf0=datanow
                auxl1=lambdanew[*,i]
                outlam=where((auxl1 gt max(auxl0)) or (auxl1 lt min(auxl0)))
;                linterp, auxl0, auxf0, auxl1, fluxnew1
                if (mean(auxl0) ne 0.0) then fluxnew1=INTERPOL(auxf0, auxl0, auxl1, /SPLINE) else fluxnew1=replicate(!Values.F_NAN, npix)     
                if (outlam ne [-1]) then fluxnew1[outlam]=!Values.F_NAN
                                ;fluxnew[*,j,beer]=-666
;INTERPOLATE FLUX ERROR TO NEW WAVELENGTH SCALE (SQUARE OF ERROR MAPS, TAKE SQUARE ROOT):
                errornowsq=(error[*,ymin+moon,beer])^2
                if bad1[0] ne -1l then errornowsq[bad1]=0 ;flag bad pixels with 0
;                ferrnewsq1=rebin_spectrum(ferrnewsq1,lambdaarr[*,i,beer],lambdanew[*,i])
                auxe0=errornowsq
;                linterp, auxl0, auxe0, auxl1, ferrnewsq1
                if (mean(auxl0) ne 0.0) then ferrnewsq1=INTERPOL(auxe0[good00], auxl0[good00], auxl1, /SPLINE) else ferrnewsq1=replicate(0., npix)
                if (outlam ne [-1]) then ferrnewsq1[outlam]=0
                if (where(ferrnewsq1 lt 0) ne [-1]) then begin; sometimes spline goes negative at the edge of bad pixels which have err=0, so I'm flagging those too
                    fluxnew1[where(ferrnewsq1 lt 0)]=!Values.F_NAN 
                    ferrnewsq1[where(ferrnewsq1 lt 0)]=0
                endif                                             
                ferrnew1=sqrt(ferrnewsq1)
                fluxnew[*,j,beer]=fluxnew1
                fluxerrornew[*,j,beer]=ferrnew1
                j=j+1
            endfor
        endfor
    endfor
    
    
    
;weight=1/error^2, make new weights:
    weightnew=1d/(fluxerrornew)^2.
    weightnew[where(fluxerrornew eq 0)]=0


;FLAG FLUX, ERROR, AND WEIGHT FRAMES WITH -666 (NaN -> -666)
    bad2=where(finite(fluxnew, /nan)) ;find where values are NaN
;    in2dbad2=array_indices(fluxnew,bad2)
;    fluxnew[in2dbad2[0,*],in2dbad2[1,*]]=-666
;    fluxerrornew[in2dbad2[0,*],in2dbad2[1,*]]=-666
;    weightnew[in2dbad2[0,*],in2dbad2[1,*]]=-666
    fluxnew[bad2]=-666
    fluxerrornew[bad2]=-666
    weightnew[bad2]=-666
 
    
    flux=dblarr(npix,nfib)
    ferr=dblarr(npix,nfib)
    print, 'Collapsing data by ',naper, ' rows per fiber and taking weighted mean'
    wait, 1
    
;    for beer=0, nframes-1 do begin
 ;       print, 'Frame number  ', (beer+1)
        wait, 1
        for i=0, nfib-1 do begin
            print, "Fiber: ", (i+1)
            ymin=i*naper ; fibers are naper pixels apart on interpolated data (fluxnew)
            ymax=(i+1)*naper-1
            for j=0, npix-1 do begin
                val=reform(fluxnew[j,ymin:ymax,*], naper*nframes) ; for each row, take pixels for each fiber
                wei=reform(weightnew[j,ymin:ymax,*], naper*nframes)
                err=reform(fluxerrornew[j,ymin:ymax,*], naper*nframes)
                good0=where(val ne -666)
                                ; SIGMA CLIP VALUES TO REJECT COSMIC RAYS2
                if (good0[0] ne -1l) then med1=median(val[good0]) else med1=-666
                                ;mede=median(err[good0])
                if (good0[0] ne -1l) then mede1=robust_sigma(val[good0]) else mede1=-666
;                if (n_elements(good0) gt 1) then mede1=stddev(val[good0]) else mede1=-666
                if (n_elements(good0) eq 1 and good0[0] ne -1l) then mede1=0
                                ; IF ALL DATA POINTS ARE BAD WRITE -666 AND EXIT
                good1=where(val ne -666 and val le med1+3.*mede1 and val ge med1-3.*mede1 ) ; bad pixels + sigma clipped pixels
                                ; IF ALL DATA POINTS ARE BAD WRITE -666 AND EXIT
                if (good1[0] ne -1l or med1 eq 0) then flux[j,i]=total(val[good1]*wei[good1])/total(wei[good1]) else flux[j,i]=-666 ;take weighted mean
;                if (good1[0] ne -1l or med1 eq 0) then ferr[j,i]=sqrt((float(n_elements(good1))*total(wei[good1]))^(-2.)*total(wei[good1]^2.*err[good1]^2)) else ferr[j,i]=-666
                if (good1[0] ne -1l or med1 eq 0) then ferr[j,i]=sqrt(total((wei[good1]*err[good1]/total(wei[good1]))^2)) else ferr[j,i]=-666

;if  (flux[j,i] eq -666) then  print, n_elements(val), n_elements(good0), n_elements(good1), med1, mede1, flux[j,i], ferr[j,i]

            endfor
            
        endfor
        
        
;    endfor
    
    

;FLAG FLUX AND ERROR FRAMES WITH -666 (NaN -> -666)
    bad33=where(finite(flux, /nan)) ;find where values are NaN
    if (bad33 ne [-1]) then begin
;        in2dbad33=array_indices(flux,bad33)
;        flux[in2dbad33[0,*],in2dbad33[1,*]]=-666
;        ferr[in2dbad33[0,*],in2dbad33[1,*]]=-666
        flux[bad33]=-666
        ferr[bad33]=-666
    endif
    

    
    
endif else begin
;for a single frame, collapse, without interpolating:

    flux=dblarr(nx, nfib)
    ferr=dblarr(nx, nfib)
    naper=5
    nval=naper



    for i=0, nfib-1 do begin
        print, "Fiber: ", (i+1)
        ymin=(i+1)*8-2-1        ; fibers are 8 pixels apart
        ymax=(i+1)*8+2-1
        for j=0, nx-1 do begin
            val=reform(data[j,ymin:ymax], nval)
            wei=reform(weight[j,ymin:ymax], nval)
            err=reform(error[j,ymin:ymax], nval)
            good0=where(val ne -666) ; bad pixels
        if (good0 eq [-1]) then begin flux[j,i]=-666 & flux[j,i]=-666 & continue & endif
                                ; SIGMA CLIP VALUES TO REJECT COSMIC RAYS
            medf=median(val[good0])
                                ;mede=median(err[good0])
            mede=robust_sigma(val[good0])
                                ; IF ALL DATA POINTS ARE BAD WRITE -666 AND EXIT
            good1=where(val ne -666 and val le medf+3.*mede and val ge medf-3.*mede ) ; bad pixels + sigma clipped pixels
                                ; IF ALL DATA POINTS ARE BAD WRITE -666 AND EXIT
            if (good1 eq [-1] or medf eq 0) then flux[j,i]=-666 
            if (good1[0] ne -1l or medf eq 0) then flux[j,i]=total(val[good1]*wei[good1])/total(wei[good1]) else flux[j,i]=-666 ;take weighted mean       
;            if (good1[0] ne -1l or medf eq 0) then ferr[j,i]=sqrt((float(n_elements(good1))*total(wei[good1]))^(-2.)*total(wei[good1]^2.*err[good1]^2)) else ferr[j,i]=-666
            if (good1[0] ne -1l or medf eq 0) then ferr[j,i]=sqrt(total((wei[good1]*err[good1]/total(wei[good1]))^2)) else ferr[j,i]=-666
            endfor
        
    endfor
    
    
    lambdanew=lambdaarr 
    
endelse





;WRITE MULTI EXTENSION FITS FILE
hnew=hdr
if keyword_set(logbin) then output=output+'_log.fits' else output=output+'.fits' 
sxaddpar,hnew,'HISTORY','------ Data cube made with IDL program mkvpcube.pro -------'

objname=sxpar(hdr, 'OBJECT') ; get object name from header
;STORE WAVELENGTH INFORMATION IN THE HEADER:
if keyword_set(multi) then begin
    ll=where(lambdanew ne 0.)
    lam=lambdanew[ll]
    sxaddpar,hnew,'CRVAL1', min(lam), "the value at the reference pixel" ;starting wavelength (Ang)
    sxaddpar,hnew,'CDELT1', delta_lambda , "the per pixel increment " ; increment (Ang) or bin size
    sxaddpar,hnew,'CTYPE1', 'LINEAR',   "the type of transform" 
    sxaddpar,hnew,'CRPIX1', 1.0, "use first pixel as reference pixel"
    sxaddpar,hnew,'CD1_1', delta_lambda
    sxaddpar,hnew,'CUNIT1','Angstrom'," Units of Wavelength"
endif
;ADD HEADER PARAMS FOR EACH EXTENSION, INCLUDE UNITS, CHANGE OBJECT NAME
sxdelpar,hnew,'OBJECT'      ;delete old object name
sxaddpar,hnew,'OBJECT',objname + " Flux","Object name" ; add object name + mask to header
sxaddpar,hnew,'BUNIT','ergs/s/cm^2/A'," Units of Flux"
mwrfits, flux, output, hnew,/create
sxdelpar,hnew,'OBJECT'      ;
sxaddpar,hnew,'OBJECT',objname + " Error in Flux","Object name" 
sxdelpar,hnew,'BUNIT'
sxaddpar,hnew,'BUNIT','ergs/s/cm^2/A'," Units of Flux"
mwrfits, ferr, output, hnew
sxdelpar,hnew,'BUNIT'
sxaddpar,hnew,'BUNIT','Angstrom'," Units of Wavelength"
sxdelpar,hnew,'OBJECT'      ;
sxaddpar,hnew,'OBJECT',objname + " Wavelength at pixel center","Object name" 
if keyword_set(logbin) then sxaddpar,hnew,'VELSCALE', velscale, "Pixel Velocity Scale in km/s"
mwrfits, lambdanew, output, hnew
sxdelpar,hnew,'BUNIT'
sxaddpar,hnew,'BUNIT','Deg'," RA,DEC, decimal degrees"
sxdelpar,hnew,'OBJECT'      ;
sxaddpar,hnew,'OBJECT',objname + " Coordinates","Object name" 
mwrfits, coords, output, hnew
sxdelpar,hnew,'BUNIT'
sxaddpar,hnew,'BUNIT','Angstrom'," Units of Instrumental Resolution"
sxdelpar,hnew,'OBJECT'      ;
sxaddpar,hnew,'OBJECT',objname + " Instrumental FWHM","Object name" 
mwrfits, fwhm, output, hnew
print, 'Data file named  ' + output + '  is made.'



;stop
end




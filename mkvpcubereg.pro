;
;+
;
;PURPOSE
;
;     This procedure combines VIRUS-P science frames into a final data
;     cube. It interpolates the spectra to a common grid in
;     wavelength, and then applies a Gaussian filter weighting scheme
;     to merge different exposures into a common regularly spaced
;     spatial grid.
;
;CALLING SEQUENCE
;
;     mkvpcube, file, output, npointings [, /LOGBIN]
;
;INPUTS:
;
;     file - file listing paths to all input files, including
;                 frames, ptow files, coordinate files, and arc lamp
;                 frames.
;
;OUTPUTS:
;
;     
;  1)   output+'_cube.fits'
;       Datacube with flux in ergs/s/cm^2 for visualiztion purposes
;
;  2)   output+'.fits'
;       RSSmulti-extension fits file that includes:
;
;       1.  Flux in ergs/s/cm^2/A
;       2.  Error in the Flux in ergs/s/cm^2/Ang
;       3.  Wavelength of each pixel in Ang
;       4.  Coordinates (Deg)
;       5.  Instrumental dispersion/FWHM in Ang
;
;KEYWORDS:
;
;     LOGBIN - Set this keyword if you would like to output a
;               datacube regularly spaced in natural log (wavelength)
;               instead of in wavelength. 
;
;
;

pro mkvpcubereg, file, output, LOGBIN=LOGBIN

; ------------- MODIFY THESE PARAMETERS ---------------
npointings=3.
delta_space=2.0             ; spatial size of output pixels [arcsec]
fwhm_gauss=4.235            ; FWHM of spatial filter [arcsec]
delta_lambda=1.1            ; spectral size of output pixels [Angstroms] 
delta_vel=60d               ; spectral size of output pixels [km/s] (if /LOGBIN is used)
; -----------------------------------------------------

print, '---------- OUTPUT SPATIAL RESOLUTION FWHM: ', sqrt(4.235^2+fwhm_gauss^2)

dmax=sqrt(-2d*(fwhm_gauss/2.355)^2*alog(1d-2)) ; maximum distance at which a fiber is conidered for a given pixel [arcsec], i.e. the fiber contributes to less than a 1% after gaussian filter
c = 299792.458d ; Speed of light in km/s
delta_loglam=delta_vel/c 



;------------------- READ PARAMETER FILE ----------------------------------------------------
nfiles=file_lines(file)         ; number of lines in param file
nframes=(nfiles)/4              ; number of frames

readcol, file, paramlist,format=('a')

ptowfile=strarr(nframes)
coordfile=strarr(nframes)
frames=strarr(nframes)
arcs=strarr(nframes)

for beer=0, nframes-1 do begin
    ptowfile[beer]=paramlist[beer]
    coordfile[beer]=paramlist[beer+fix(nframes)]
    frames[beer]=paramlist[beer+fix(2*nframes)]
    arcs[beer]=paramlist[beer+fix(3*nframes)]
endfor
; ---------------------------------------------------------------------------------------------



;------------------- READ PTOW FILE, CHECKING FOR ORDER OF WAVELENGTH SOLUTION
print, "---------- READING PTOW FILES ---------"

ncols=intarr(nframes)
nfib=intarr(nframes)

for moon=0, nframes-1 do begin
    openr, lun, ptowfile[moon], /get_lun
    line=''
    readf, lun, line
    free_lun, lun
    ncols[moon]=fix(n_elements(strsplit(line, /regex, /extract))) ;order of polynomial for wavelength solution
    nfib[moon]=fix(file_lines(ptowfile[moon]))
endfor

maxnfib=max(nfib)
maxncols=max(ncols)
acoef=dblarr(maxncols,maxnfib,nframes); array with the wavelength solution coefficients for each fiber

for moon=0, nframes-1 do begin
    acoef1=dblarr(ncols[moon],nfib[moon]) 
    openr, lun1, ptowfile[moon], /get_lun
    readf, lun1, acoef1
    free_lun, lun1
    acoef[0:ncols[moon]-1,0:nfib[moon]-1,moon]=acoef1
endfor

; ----------------------------      READ COORDINATES FILE --------------------------------------------
print, '---------- READING COORDINATE FILES ---------'

bb=strarr(maxnfib, nframes)
coords=dblarr(2, maxnfib, nframes)
for moon=0, nframes-1 do begin
    bb1=strarr(nfib[moon])
    openr, lun, coordfile[moon], /get_lun
    readf, lun, bb1
    free_lun, lun
    bb[0:nfib[moon]-1,moon]=bb1
    aux=strarr(3,nfib[moon])
    for i=0, nfib[moon]-1 do aux[*,i]=strsplit(bb[i,moon], /regex, /extract)
    coords1=dblarr(2,nfib[moon])
    for i=0, nfib[moon]-1 do begin get_coords, junk, instring=aux[1,i]+' '+aux[2,i] & coords1[*,i]=junk & endfor
    coords1[0,*]=15d*coords1[0,*] ;so both ra and dec are in decimal degrees
    coords[*,0:nfib[moon]-1,moon]=coords1
endfor
; ----------------------------------------------------------------------------------------------------


;--------------- GET DATE, UT, WAVEZP, AIRMASS, EXPTIME from header ---------------

print, '---------- READING UT ---------'

month=intarr(nframes)
year=intarr(nframes)
day=intarr(nframes)
ut=fltarr(nframes)
wavezp=dblarr(nframes)
wavezperr=dblarr(nframes)
airmass=dblarr(nframes)
exptime=dblarr(nframes)

for moon=0, nframes-1 do begin
    print, 'READING UT FOR FRAME: ', moon+1, ' / ', strtrim(string(nframes),2)
;    iframe=mrdfits(frames[moon]+'pefsma.fits', 0, hdrf)
    hdrf=HEADFITS(frames[moon]+'pefsma.fits')
    date1=sxpar(hdrf, 'DATE-OBS')
    month[moon]=strmid(date1,5,2)
    year[moon]=strmid(date1, 0,4)
    day[moon]=strmid(date1,8,2)
    ut1=sxpar(hdrf, 'UT')
    ut[moon]=strmid(ut1,0,2)+strmid(ut1,3,2)/60d +  strmid(ut1,6,5)/3600d
    wavezp[moon]=sxpar(hdrf, 'WAVEZP')
    wavezperr[moon]=sxpar(hdrf, 'WAVEERR')
    airmass[moon]=sxpar(hdrf, 'AIRMASS')
    exptime[moon]=sxpar(hdrf, 'EXPTIME')
endfor



; --------------- CREATE WAVELENGTH MAP --------------------------

aux=mrdfits(frames[0]+'pefsma.fits', 0, hdrf1)
nx=(size(aux))[1]  
ny=(size(aux))[2]

lambdaarr=dblarr(nx,maxnfib,nframes)

for moon=0, nframes-1 do begin
    print, 'MAKING WAVELENGTH ARRAYS FOR FRAME', moon+1, ' / ', strtrim(string(nframes),2)
    xarr=dindgen(nx+1)
    lambda0=dblarr(nx+1, nfib[moon])
    lambda=dblarr(nx, nfib[moon])
    for i=0, nfib[moon]-1 do begin & for j=0, ncols[moon]-1 do lambda0[*,i]=lambda0[*,i]+acoef[j,i,moon]*xarr^j & endfor
    for i=0, nx-1 do lambda[i,*]=(lambda0[i,*]+lambda0[i+1,*])/2. ;place lambda at center of pixel
    lambdaarr[*,*,moon]=lambda
endfor

lambdaarr0=lambdaarr


;---------- CORRECT WAVELENGTH ARRAY FOR BARYCENTRIC VELOCITY ----------

lambdac=dblarr(nx,maxnfib,nframes)
for moon=0, nframes-1 do begin
    print, 'BARYCENTRIC CORRECTION FOR FRAME', moon+1, ' / ', strtrim(string(nframes),2)
    jdcnv, year[moon], month[moon], day[moon], UT[moon], jd ;get julian date at time object was observed(roughly for UT in middle of the night, UT=6 hr)   
    baryvel, jd, 2000, vhel, vbary ;compute heliocentric & barycentric velocity (in km/s) at jd
;extract RA & DEC of object from image frame header:


ra2=mean(coords[0,*,moon])/!RADEG
dec2=mean(coords[1,*,moon])/!RADEG

;velocity correction in km/s: 
    vcor = vbary[0]*cos(dec2)*cos(ra2) + $ ;Project velocity toward star
      vbary[1]*cos(dec2)*sin(ra2) + vbary[2]*sin(dec2) 
    cc=2.99792458d5             ; speed of light in km/s
    lambdac[*,*,moon]= (vcor/cc)*lambdaarr[*,*,moon] ;barycentric wavelength correction
    lambdaarr[*,*,moon]=lambdac[*,*,moon]+ lambdaarr[*,*,moon] ; corrected wavelength array
    for j=0, nfib[moon]-1 do if (wavezp[moon] ne 0.0 and mean(lambdaarr[*,j,moon]) ne 0) then lambdaarr[*,j,moon]=lambdaarr[*,j,moon]+wavezp[moon] ; if lambda-ZP from skyline is in header use it
endfor



;------------- CREATE OUTPUT WAVELENGTH ARRAY ----------------------

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
    lambdastart=double(floor(min(minn)))
    lambdaend=double(ceil(max(maxx)))
;define number of pixels by range of wavelength in all frames:
    npix=round((lambdaend-lambdastart)/delta_lambda)
    print, '---------- OUTPUT WAVELENGTH GRID ----------'
    print, 'WAVELENGTH RANGE:  ', lambdastart, ' to ',lambdastart+(npix-1)*delta_lambda, '  Ang with N_pixels =', npix
;make new wavelength array for each fiber starting from bluest
;wavelength in data set
    lambdanew=dblarr(npix)
    for beer=0, npix-1 do begin
        wavel=lambdastart+delta_lambda*beer
        lambdanew[beer]=wavel
    endfor


   if keyword_set(logbin) then begin
        npix=round((alog(lambdaend)-alog(lambdastart))/delta_loglam)
        lambdanew=dblarr(npix)
        lambdanew=exp(alog(lambdastart)+delta_loglam*dindgen(npix))  
        velscale=(alog(lambdanew[npix-1])-alog(lambdanew[0]))/(npix-1)*c
    endif
      
; ------------ CREATE OUTPUT DATA CUBE ----------------


maxra=max(coords[0,*,*])+dmax/3600d
minra=min(coords[0,*,*])-dmax/3600d
maxdec=max(coords[1,*,*])+dmax/3600d
mindec=min(coords[1,*,*])-dmax/3600d



npixdec=round((maxdec-mindec)*3600d/delta_space)
meddec=median(coords[1,*,*])
npixra=round((maxra-minra)*3600d*cos(meddec*!pi/180d)/delta_space)

decarr=dindgen(npixra, npixdec)
raarr=dindgen(npixra, npixdec)

for i=0, npixdec-1 do  decarr[*,i]=mindec+i*delta_space/3600d
for i=0, npixdec-1 do begin &  for j=0, npixra-1 do raarr[j,i]=maxra-j*delta_space/3600d/cos(decarr[j,i]*!pi/180d) & endfor


fluxout=dblarr(npixra, npixdec, npix)
ferrout=dblarr(npixra, npixdec, npix)
arcout=dblarr(npixra, npixdec, npix)
Npixout=dblarr(npixra, npixdec, npix)
fwhmout=dblarr(npixra, npixdec, npix)
Nneararr=dblarr(npixra, npixdec)


; ------------- INTERPOLATE AND COMBINE DATA BY BLOCKS IN WAVELENGTH
;               ---------------------------


Nblocks=9
dpix=round(Npix/Nblocks)
pix1=indgen(Nblocks)*dpix
pix2=indgen(Nblocks)
for i=0, Nblocks-2 do pix2[i]=pix1[i+1]-1
pix2[Nblocks-1]=npix-1

naper=5                         ; extraction aperture height in pixels
nval=naper
nynew=naper*maxnfib ;-> data file interpolated to new wavelength with just extracted fiber values (naper=number of pixels to collapose for each fiber pixel)



for l=0, Nblocks-1 do begin
    print, 'DOING BLOCK ', l+1, ' / ', strtrim(string(Nblocks),2)
    
    lam1=lambdanew[pix1[l]]
    lam2=lambdanew[pix2[l]]
    auxnpix=pix2[l]-pix1[l]+1
    
    data=dblarr(nx,ny)
    comp=dblarr(nx,ny)
    weight=dblarr(nx,ny)
    error=dblarr(nx,ny)
    fluxnew=dblarr(auxnpix,nynew,nframes)
    arcnew=dblarr(auxnpix,nynew,nframes)
    fluxerrornew=dblarr(auxnpix,nynew,nframes)
    
    for beer=0, nframes-1 do begin 

        print, 'READING FRAME: ', beer+1, ' / ', nframes 
        data=mrdfits(frames[beer]+'pefsma.fits', 0, hdr)



        comp=mrdfits(arcs[beer], 0, hdrarc)
       ; NORMALIZE ARCS
        comp=comp/mean(comp)



        weight=mrdfits(frames[beer]+'pefwea.fits', 0, hdr2)
        ; FIX WEIGHT AND ERROR MAPS
        weight=1.0/(weight)^2 ; Mimi's code now outputs error, not weights!!!!
        error=1.0/sqrt(weight) ;error frames, where weight is 1/noise^2 or 1/error^2
        
  

    
        print, 'INTERPOLATING FRAME: ', (beer+1), ' / ', nframes
        j=0
        for i=0, nfib[beer]-1 do begin
;           print, "FIBER: ", (i+1)
            ymin=(i+1)*8-2-1    ; fibers are 8 pixels apart
            ; --INTERPOLATE FLUX TO NEW WAVELENGTH SCALE:
            for moon=0, naper-1 do begin
                datanow=data[*,ymin+moon]
                bad1=where(datanow eq -666, complement=good00)
                if bad1[0] ne -1l then datanow[bad1]=!Values.F_NAN ;flag bad pixels with NaN so to preserve them in interpolation
                auxl0=lambdaarr[*,i,beer]
                auxl0arc=lambdaarr0[*,i,beer] ; arcs don't have barycentric correction
                auxf0=datanow
                auxa0=comp[*,ymin+moon]
                auxl1=lambdanew[pix1[l]:pix2[l]]
                
                outlam=where((auxl1 gt max(auxl0)) or (auxl1 lt min(auxl0)))
               ;INTERPOLATE ARC
;               if (mean(auxl0) ne 0.0) then arcnew1=INTERPOL(auxa0, auxl0arc, auxl1, /SPLINE) else arcnew1=replicate(!Values.F_NAN, auxnpix)     
                if (mean(auxl0) ne 0.0) then arcnew1=INTERPOL(auxa0, auxl0arc, auxl1) else arcnew1=replicate(!Values.F_NAN, auxnpix)     
               ;INTERPOLATE FLUX
;               if (mean(auxl0) ne 0.0) then fluxnew1=INTERPOL(auxf0, auxl0, auxl1, /SPLINE) else fluxnew1=replicate(!Values.F_NAN, auxnpix)     
                if (mean(auxl0) ne 0.0) then fluxnew1=INTERPOL(auxf0, auxl0, auxl1) else fluxnew1=replicate(!Values.F_NAN, auxnpix)     
                if (outlam ne [-1]) then fluxnew1[outlam]=!Values.F_NAN
               ;INTERPOLATE FLUX ERROR TO NEW WAVELENGTH SCALE (SQUARE OF ERROR MAPS, TAKE SQUARE ROOT):
                errornowsq=(error[*,ymin+moon])^2
                if bad1[0] ne -1l then errornowsq[bad1]=0 ;flag bad pixels with 0
                auxe0=errornowsq
;               if (mean(auxl0) ne 0.0) then ferrnewsq1=INTERPOL(auxe0[good00], auxl0[good00], auxl1, /SPLINE) else ferrnewsq1=replicate(0., auxnpix)
                if (mean(auxl0) ne 0.0) then ferrnewsq1=INTERPOL(auxe0[good00], auxl0[good00], auxl1) else ferrnewsq1=replicate(0., auxnpix)
                if (outlam ne [-1]) then ferrnewsq1[outlam]=0
                if (where(ferrnewsq1 lt 0) ne [-1]) then begin; sometimes spline goes negative at the edge of bad pixels which have err=0, so I'm flagging those too
                    fluxnew1[where(ferrnewsq1 lt 0)]=!Values.F_NAN 
                    ferrnewsq1[where(ferrnewsq1 lt 0)]=0
                endif                                             
                ferrnew1=sqrt(ferrnewsq1)
                fluxnew[*,j,beer]=fluxnew1
                arcnew[*,j,beer]=arcnew1
                fluxerrornew[*,j,beer]=ferrnew1
                j=j+1
            endfor
        endfor





    endfor

; Release memory
    data=''
    comp=''
    error=''
    weight=''
  
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
    arcnew[bad2]=-666
 


; ------------ COMBINE THE INTERPOLATED SPECTRA


    for i=0, npixra-1 do begin
        print, 'DOING COLUMN: ',i+1 , ' / ', npixra , ' of Section ', l
        for j=0, npixdec-1 do begin
                                ; DEFINE WHICH FIBERS WILL BE CONSIDERED
            d=sqrt(((coords[0,*,*]-raarr[i,j])*3600d*cos(0.5*(decarr[i,j]+coords[1,*,*])*!pi/180d))^2+((coords[1,*,*]-decarr[i,j])*3600d)^2)
            
            aux=where(d le dmax)
            
            if (aux[0] ne -1l) then begin
                near=array_indices(d, aux) 
                near=near[1:2,*] ; near[0,*] has fiber numbers and near[1,*] has frame numbers
                Nnear=n_elements(near[0,*])
                Nneararr[i,j]=Nnear
        ;        print,'Number of Fibers to Combine: ', Nnear
                dnear=d[where(d le dmax)]
                fluxnew2=dblarr(auxnpix,Nnear*naper)
                fluxerrornew2=dblarr(auxnpix,Nnear*naper)
                weightnew2=dblarr(auxnpix,Nnear*naper)
                arcnew2=dblarr(auxnpix,Nnear*naper)
                
                for l2=0, Nnear-1 do begin
                    ymin0=l2*naper
                    ymax0=l2*naper+naper-1
                    ymin1=near[0,l2]*naper
                    ymax1=near[0,l2]*naper+naper-1
                    fluxnew2[*,ymin0:ymax0]=fluxnew[*,ymin1:ymax1,near[1,l2]]
                    fluxerrornew2[*,ymin0:ymax0]=fluxerrornew[*,ymin1:ymax1,near[1,l2]]
                    weightnew2[*,ymin0:ymax0]=weightnew[*,ymin1:ymax1,near[1,l2]]*exp(-0.5*dnear[l2]^2/(fwhm_gauss/2.355)^2) ; introduce gaussian filter
                    arcnew2[*,ymin0:ymax0]=arcnew[*,ymin1:ymax1,near[1,l2]]           
                endfor
                
                
                for k=0, auxnpix-1 do begin
                    val=reform(fluxnew2[k,*], naper*Nnear) ; for each wavelength take all values
                    wei=reform(weightnew2[k,*], naper*Nnear)
                    err=reform(fluxerrornew2[k,*], naper*Nnear)
                    good0=where(val ne -666)
                    
; SIGMA CLIP VALUES TO REJECT COSMIC RAYS2
                    if (good0[0] ne -1l) then med1=median(val[good0]) else med1=-666
                    if (good0[0] ne -1l) then mede1=robust_sigma(val[good0]) else mede1=-666
                    if (n_elements(good0) eq 1 and good0[0] ne -1l) then mede1=0
                    good1=where(val ne -666 and val le med1+3.*mede1 and val ge med1-3.*mede1, Ngood) ; bad pixels + sigma clipped pixels
                    if (good1[0] ne -1l or med1 eq 0) then fluxout[i,j,k+pix1[l]]=total(val[good1]*wei[good1])/total(wei[good1]) else fluxout[i,j,k+pix1[l]]=-666 ;take weighted mean
                    if (good1[0] ne -1l or med1 eq 0) then ferrout[i,j,k+pix1[l]]=sqrt(total((wei[good1]*err[good1]/total(wei[good1]))^2)) else ferrout[i,j,k+pix1[l]]=-666
                    Npixout[i,j,k+pix1[l]]=Ngood
                    
                    arcval=reform(arcnew2[k,*], naper*Nnear) ; for each row, take pixels for each fiber
                    agood0=where(arcval ne -666)
                    if (agood0[0] ne -1l) then amed1=median(arcval[agood0]) else amed1=-666
                    if (agood0[0] ne -1l) then amede1=robust_sigma(arcval[agood0]) else amede1=-666
                    agood1=where(arcval ne -666 and arcval le amed1+3.*amede1 and arcval ge amed1-3.*amede1 ) ; bad pixels + sigma clipped pixels
                    if (agood1[0] ne -1l or amed1 eq 0) then arcout[i,j,k+pix1[l]]=total(arcval[agood1]*wei[agood1])/total(wei[agood1]) else arcout[i,j,k+pix1[l]]=-666 ;take weighted mean
                    
                endfor
                
            endif else begin
                fluxout[i,j,*]=-666
                ferrout[i,j,*]=-666
                arcout[i,j,*]=-666            
                Npixout[i,j,*]=0
            endelse
            
        endfor
    endfor
endfor



;MAKS LOW EXPOSURE REGIONS
Nthresh=Nframes/Npointings/2.*3.
print, 'MASKING ANY SPECTRA WITH BLUE AND RED NPIX < ', NTHRESH
blue=where(lambdanew ge 3600 and lambdanew le 4300)
red=where(lambdanew ge 6000 and lambdanew le 6700)

for i=0, npixra-1 do begin
    for j=0, npixdec-1 do begin
        if (median(Npixout[i,j,blue]) lt Nthresh or median(Npixout[i,j,red]) lt Nthresh) then begin
                fluxout[i,j,*]=-666
                ferrout[i,j,*]=-666
                arcout[i,j,*]=-666            
                Npixout[i,j,*]=0
        endif        
    endfor
endfor




;CREATE MAP OF INSTRUMENTAL RESOLUTION
print, '---------- CREATING SPECTRAL RESOLUTION MAP ---------'

readcol, '/home/alicia/gblancm/idl/vengaidl/redarclines.txt', linelam, comment='#'
nlines=n_elements(linelam)
a0=fltarr(nlines, npixra, npixdec)
a1=fltarr(nlines, npixra, npixdec)
a2=fltarr(nlines, npixra, npixdec)
a3=fltarr(nlines, npixra, npixdec)

; fit each line with a gaussian
for j=0, npixra-1 do begin
    for k=0, npixdec-1 do begin
        for i=0, nlines-1 do begin
            xmin=min(where(lambdanew ge linelam[i]-10.))
            xmax=max(where(lambdanew le linelam[i]+10.))
            x=(lambdanew)[xmin:xmax]
            y=(arcout[j,k,*])[xmin:xmax]
            fit=jjgaussfit(x,y,param)
;        plot, x, y, title=strtrim(string(linelam[i]),2)
;        oplot, x, fit, color=250
;stop
            a0[i,j,k]=param[0]
            a1[i,j,k]=param[1]
            a2[i,j,k]=param[2]
            a3[i,j,k]=param[3]
        endfor
    endfor
endfor


FWHM0=2.355*a2
; fit FWHM as afunction of lambda with a 2nd order polynomial
coef=fltarr(3,npixra,npixdec)
for i=0, npixra-1 do begin &    for j=0, npixdec-1 do begin &        goodlines=where(a3[*,i,j] ne -666) &     if(goodlines ne [-1]) then   coef[*,i,j]=poly_fit(linelam[goodlines], FWHM0[goodlines,i,j], 2) &        fwhmout[i,j,*]=coef[0,i,j]+coef[1,i,j]*lambdanew+coef[2,i,j]*lambdanew^2 &    endfor & endfor


fwhmout[where(fluxout eq -666)]=0

;for i=0, npixra-1 do begin & for j=0, npixdec-1 do begin & plot, linelam, FWHM0[*,i,j], psym=2, yrange=[3.9,6], ystyle=1, title='PIXEL'+string(i)+' : '+string(j) & oplot, lambdanew, arcout[i,j,*]/max(arcout[i,j,*])+4 & oplot, lambdanew, coef[0,i,j]+coef[1,i,j]*lambdanew+coef[2,i,j]*lambdanew^2, linestyle=2 & pause & endfor & endfor

    


; transform flux per fiber to flux per new spatial pixel
sel=where(fluxout ne -666)
fluxout[sel]=fluxout[sel]/(!pi*(4.235/2.)^2)*delta_space^2
ferrout[sel]=ferrout[sel]/(!pi*(4.235/2.)^2)*delta_space^2

; ---------- WRITE OUTPUT FILES -------
print, '--------- WRITING OUTPUT FILES: DATACUBE FOR VISUALISATION---------'

if keyword_set(logbin) then filename=output+'_cube_log.fits' else filename=output+'_cube.fits' 


sxaddpar, hnew, 'HISTORY','------ Data cube made with IDL program mkvpcubereg.pro -------'
sxaddpar, hnew, 'OBJECT',strtrim(string(lambdanew[0], format='(f10.3)'),2)

sxaddpar, hnew, 'CTYPE1', 'RA---TAN', 'Coordinate Type'
sxaddpar, hnew, 'CTYPE2', 'DEC--TAN', 'Coordinate Type'
sxaddpar, hnew, 'CUNIT1', 'deg'
sxaddpar, hnew, 'CUNIT2', 'deg'
sxaddpar, hnew, 'CRPIX1', 1.0, 'Ref Pixel X'
sxaddpar, hnew, 'CRPIX2', 1.0, 'Ref Pixel X'
sxaddpar, hnew, 'CRVAL1', raarr[0,0], 'RA at Reference Pixel'
sxaddpar, hnew, 'CRVAL2', decarr[0,0], 'DEC at Reference Pixel'
sxaddpar, hnew, 'CD1_1', -1.0*delta_space/3600. , 'Linear projection matrix'
sxaddpar, hnew, 'CD1_2', 0.0 , 'Linear projection matrix'
sxaddpar, hnew, 'CD2_1', 0.0 , 'Linear projection matrix'
sxaddpar, hnew, 'CD2_2', delta_space/3600. , 'Linear projection matrix'


mwrfits ,fluxout[*,*,0], filename, hnew, /create
for i=1, npix-1 do begin &    sxdelpar,hnew,'OBJECT'   & sxaddpar,hnew,'OBJECT',strtrim(string(lambdanew[i], format='(f10.3)'),2) & mwrfits ,fluxout[*,*,i], filename, hnew & endfor    



print, '--------- WRITING OUTPUT FILES: RSS FILE FOR ANALYSIS ---------'

if keyword_set(logbin) then filename=output+'_log.fits' else filename=output+'.fits' 





objname=output ; get object name from header
;STORE WAVELENGTH INFORMATION IN THE HEADER:

mwrfits, transpose(reform(fluxout, npixra*npixdec, npix)), 'junk',/create
junk=mrdfits('junk',0,hdr)
sxaddpar, hdr, 'HISTORY','------ RSS file made with IDL program mkvpcubereg.pro -------'
sxaddpar,hdr,'CRVAL1', lambdanew[0], "the value at the reference pixel" ;starting wavelength (Ang)
sxaddpar,hdr,'CDELT1', delta_lambda , "the per pixel increment " ; increment (Ang) or bin size
sxaddpar,hdr,'CTYPE1', 'LINEAR',   "the type of transform" 
sxaddpar,hdr,'CRPIX1', 1.0, "use first pixel as reference pixel"
sxaddpar,hdr,'CD1_1', delta_lambda
sxaddpar,hdr,'CUNIT1','Angstrom'," Units of Wavelength"
sxdelpar,hdr,'OBJECT'      ;delete old object name
sxaddpar,hdr,'OBJECT',objname + " Flux","Object name" ; add object name + mask to header
sxaddpar,hdr,'BUNIT','ergs/s/cm^2/A'," Units of Flux"
mwrfits, transpose(reform(fluxout, npixra*npixdec, npix)), filename, hdr,/create


mwrfits, transpose(reform(ferrout, npixra*npixdec, npix)), 'junk'
junk=mrdfits('junk',1,hdr)
sxdelpar,hdr,'OBJECT'      ;
sxaddpar,hdr,'OBJECT',objname + " Error in Flux","Object name" 
sxaddpar,hdr,'CRVAL1', lambdanew[0], "the value at the reference pixel" ;starting wavelength (Ang)
sxaddpar,hdr,'CDELT1', delta_lambda , "the per pixel increment " ; increment (Ang) or bin size
sxaddpar,hdr,'CTYPE1', 'LINEAR',   "the type of transform" 
sxaddpar,hdr,'CRPIX1', 1.0, "use first pixel as reference pixel"
sxaddpar,hdr,'CD1_1', delta_lambda
sxaddpar,hdr,'CUNIT1','Angstrom'," Units of Wavelength"
sxdelpar,hdr,'BUNIT'
sxaddpar,hdr,'BUNIT','ergs/s/cm^2/A'," Units of Flux"
mwrfits, transpose(reform(ferrout, npixra*npixdec, npix)), filename, hdr


mwrfits, lambdanew, 'junk'
junk=mrdfits('junk',2,hdr)
sxdelpar,hdr,'OBJECT'      ;
sxaddpar,hdr,'OBJECT',objname + " Wavelength at pixel center","Object name" 
sxaddpar,hdr,'CRVAL1', lambdanew[0], "the value at the reference pixel" ;starting wavelength (Ang)
sxaddpar,hdr,'CDELT1', delta_lambda , "the per pixel increment " ; increment (Ang) or bin size
sxaddpar,hdr,'CTYPE1', 'LINEAR',   "the type of transform" 
sxaddpar,hdr,'CRPIX1', 1.0, "use first pixel as reference pixel"
sxaddpar,hdr,'CD1_1', delta_lambda
sxaddpar,hdr,'CUNIT1','Angstrom'," Units of Wavelength"
sxdelpar,hdr,'BUNIT'
sxaddpar,hdr,'BUNIT','Angstrom'," Units of Wavelength"
if keyword_set(logbin) then sxaddpar,hdr,'VELSCALE', velscale, "Pixel Velocity Scale in km/s"
mwrfits, lambdanew, filename, hdr


mwrfits, [[reform(raarr, npixra*npixdec)], [reform(decarr, npixra*npixdec)]], 'junk'
junk=mrdfits('junk',3,hdr)
sxdelpar,hdr,'OBJECT'      ;
sxaddpar,hdr,'OBJECT',objname + "Pixel Equatorial Coordinates","Object name" 
sxdelpar,hdr,'BUNIT'
sxaddpar,hdr,'BUNIT','Deg'," RA,DEC, decimal degrees"
mwrfits, [[reform(raarr, npixra*npixdec)], [reform(decarr, npixra*npixdec)]], filename, hdr

mwrfits, transpose(reform(fwhmout, npixra*npixdec, npix)), 'junk'
junk=mrdfits('junk',4,hdr)
sxdelpar,hdr,'OBJECT'      ;
sxaddpar,hdr,'OBJECT',objname + " Instrumental FWHM","Object name" 
sxaddpar,hdr,'CRVAL1', lambdanew[0], "the value at the reference pixel" ;starting wavelength (Ang)
sxaddpar,hdr,'CDELT1', delta_lambda , "the per pixel increment " ; increment (Ang) or bin size
sxaddpar,hdr,'CTYPE1', 'LINEAR',   "the type of transform" 
sxaddpar,hdr,'CRPIX1', 1.0, "use first pixel as reference pixel"
sxaddpar,hdr,'CD1_1', delta_lambda
sxaddpar,hdr,'CUNIT1','Angstrom'," Units of Wavelength"
sxdelpar,hdr,'BUNIT'
sxaddpar,hdr,'BUNIT','Angstrom'," Units of FWHM Instrumental Resolution"
mwrfits, transpose(reform(fwhmout, npixra*npixdec, npix)), filename, hdr


mwrfits, transpose(reform(Npixout, npixra*npixdec, npix)), 'junk'
junk=mrdfits('junk',5,hdr)
sxdelpar,hdr,'OBJECT'      ;
sxaddpar,hdr,'OBJECT',objname + " Number of combined pixels","Object name" 
sxaddpar,hdr,'CRVAL1', lambdanew[0], "the value at the reference pixel" ;starting wavelength (Ang)
sxaddpar,hdr,'CDELT1', delta_lambda , "the per pixel increment " ; increment (Ang) or bin size
sxaddpar,hdr,'CTYPE1', 'LINEAR',   "the type of transform" 
sxaddpar,hdr,'CRPIX1', 1.0, "use first pixel as reference pixel"
sxaddpar,hdr,'CD1_1', delta_lambda
sxaddpar,hdr,'CUNIT1','Angstrom'," Units of Wavelength"
sxdelpar,hdr,'BUNIT'
sxaddpar,hdr,'BUNIT','N'," Number of Fibers"
mwrfits, transpose(reform(Npixout, npixra*npixdec, npix)), filename, hdr



print, 'Data file named  ' + output + '  is made.'





stop
end

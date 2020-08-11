PRO calcRES, month_night

; RUN THIS CODE TO GET THE FITS FILES THAT RESMAKER_V3 NEEDS IN ORDER
; TO MAKE THE BINNAME.RES FILES USED FOR LOSVD EXTRACTION...

; MONTH_NIGHT: This is whatever hangs off the back of res.list.. So, if your
; list is res.list.march11, then you add just 'march11' YOU DO NOT
; NEED TO USE THE SAME THING YOU USED TO RUN VACCINE. IT'S
; GOING TO WORK FROM THAT LIST AND NOT CARE ABOUT WHETHER THESE NAMES MATCH.

; This code is essentially a wrapper for rescheckHgCd. It reads in the
; same things as rescheckHgCd (arc,ptow,mask,# of lines), runs
; rescheckHgCd, and then conducts some further analysis of the returned
; values. It's made to be flexible and changed as necessary to analyze
; the IR results

; COMPILE: rescheckhgcd   .r rescheckHgCd
; COMPILE: rescheckiron   .r rescheckIron - sometimes
; COMPILE: rescheckiron2  .r rescheckIron2 - sometimes
; COMPILE: calcRES        .r calcRES
; COMPILE: resmapf_v2        .r resmapF

; Any word in the Plot call will cause rescheckHgCd to plot.

; MODIFIED ON DEC 3rd, 2010: This code now looks for a list, called
; 'res.list.month' that are the arcs, ptow, mask files for in an
; observing run. YOU CAN EITHER RUN THIS NIGHT BY NIGHT, OR SEND IN
; THE VALUES FROM SEVERAL NIGHTS (I.E. OVER AN OBSERVING RUN) AND IT
; WILL TAKE THE MEDIAN OF EACH WAVELENGTH VALUE). THIS IS RECOMMENDED
; AS THE .RES VALUES CREATED IN LATER ROUTINES GET VERY CUMBERSOME
; WHEN YOU HAVE TOO MANY RES FILES.

; MODIFIED ON MAY 4, 2011: Made some general cleaning up and commenting.

; MODIFIED ON MAY 15, 2011 (OBSOLETE): There are now 3 different sets of values
; returned with this routine (both as fits and text files). They are
; described below:
; RAW: These are the individual gaussian fits. For an observing run,
; each line fit to each fiber is median-combined. So the 'raw' isn't
; as raw as a single night, but rather the median combination of every
; frame that you feed in via your res.list.month list.
; POLYFIT: This is a polynomial fit along a single fiber, AFTER the
; nights have been median combined. There's a switch now that allows
; you to switch to whether you're returning a 1st to 4th order
; polynomial fit.
; SMOOTH-[RF]: There are two 'smoothed' versions that are
; returned. The first is SMOOTH-RAW which is the smoothed version of
; the raw values as described above. The second is SMOOTH-FIT, which
; is the smoothed version of the poly fits to each fiber.

; MODIFIED ON SEPT 7TH, 2011: I have simplified things in this
; code. It used to put out several versions of the IR map. This has
; been trimmed to three fits files. The first, called
; VIRUS-P_IR_RAW_month_night.fits is the 10 (or so) raw Gaussian fits
; to individual arc lines. They are the median-combination if the
; res.list.month was any longer than 1 line. 
; These are then fit with a polynomial (2nd
; order seems the best) and the values at each of the 10 (or so)
; wavelength locations are returned as file
; VIRUS-P_IR_FIT_month_night.fits. Then, the polynomial coefficients
; are used to generate an IR value FOR ALL THE PIXELS IN THE
; RANGE. This is generally a 246x6144 fits file and is named
; VIRUS-P_IR_FITALL_month_night.fits. THERE IS NO LONGER ANY SMOOTHING
; STEP IN THE RESMAPF.PRO ROUTINE AS THE POLY FIT DOES ALL THE
; SMOOTHING THAT IS PHYSICALLY REASONABLE.
; THIS ALSO NOW USES VERSION 2 OF THE RESMAP ROUTINE.
;*********************************************************
;plot = 'off'
w1 = 3550.0 ;these three terms will be used to generate the poly-fit to the full VIRUS-P wavelength range.
;wd = 1.125
wd = 1.125
plotting  = 'on' ;plots the final polynomial fits
gaussplot = 'on' ;plots the gaussian fits for every 10th fiber
plot3d = 'on' ;plots the 3-D map
outnameR  = 'VIRUS-P_IR_RAW_'+month_night
outnameF  = 'VIRUS-P_IR_FIT_'+month_night 
outnameALL = 'VIRUS-P_IR_FITALL_'+month_night ;this is now the one to use for the subsequent reductions.
;outnameSR = 'VIRUS-P_IR_SMOOTH-RAW_'+month_night
;outnameSF  = 'VIRUS-P_IR_SMOOTH-FIT_'+month_night
;outnameSALL = 'VIRUS-P_IR_SMOOTH-ALL_'+month_night
smoothingF = 1 ;the size of the smoothing kernal used by filter_image in resmapF
boxsizeF = 1 ;the size of the median box (?) used by filter_image in resmapF.
smoothingR = 3 ;the size of the smoothing kernal used by filter_image in resmapF
boxsizeR = 3 ;the size of the median box (?) used by filter_image in resmapF.
;the 'F' is the smoothing parameters for the FIT version, and 'R' is
;for RAW.
pfit = 'second' ;this can be 'first', 'second', 'third', or 'fourth'. It's the order of the fitting polynomial.
factor = 1.0 ;THIS MUST BE SET TO 1.0 IF YOU WANT TO BUILD A PTOP FILE!
wdd = wd / factor
;*********************************************************

ans = ''
print,'Initial wavelength: '+strn(w1)
print,'Dispersion term: '+strn(wdd)
;print,'Is this correct? Hit "q" to quit or ENTER to continue...'
;read,ans
;if (ans eq 'q') then stop

readcol,'res.list.'+month_night,silent=1,f='a,a,a',arc,ptow,mask

n0 = n_elements(arc)            ; the number of nights

;wavelist = [4046.5469,4077.8403,4358.3262,4678.1558,4799.9038,$
;            4916.0962,5085.8110,5460.7397,5769.5972,5790.6348]
;wavelist = [4046.5469,4077.8403,4358.3262,4678.1558,4799.9038,$
;            4916.0962,5085.8110,5460.7397,5769.5972];for feb08
wavelist = [4046.5469,4077.8403,4358.3262,4678.1558,4799.9038,$
            5085.8110,5460.7397,5769.5972];for oct07
waveall = wavelist

for j=0,n0-1 do begin
    if (j eq 0) then begin
        readcol,ptow[0],silent=1,f='d,x,x,x,x',c0
        n1 = n_elements(c0)
        RESarr = dblarr(n_elements(wavelist),n1,5,n0)
    endif
;    print,'Fitting gaussians to '+arc[j]+'...'
    array = rescheckhgcd(arc[j],ptow[j],mask[j],n_elements(wavelist),Splot=gaussplot)
    RESarr[*,*,*,j] = array

endfor
i = where(RESarr eq 0.0)
RESarr[i] = !Values.F_NAN

if (plotting eq 'on') then begin
    window,2,retain=2
    device,decomposed=0
    for j=0,n0-1 do begin
        loadct,0,/silent
        plot,RESarr[*,0,0,j],RESarr[*,0,3,j],psym=sym(1),symsize=0.5,$
          xtitle='Wavelength',ytitle='FWHM (Ang)',yrange=[4,7],$
          title='FWHM values for '+arc[j]
        loadct,33,/silent
        for k=0,n1-1 do begin
            oplot,RESarr[*,k,0,j],RESarr[*,k,3,j],psym=sym(1),$
              symsize=0.5,color=k
;            wait,0.02
        endfor
        wait,2.0
    endfor
    if (n0 ge 2) then print,'The next ENTER deletes the plot(s)!'
    if (n_elements(looplist) eq 0) then wait,3
    if (n0 ge 2) then while !d.window ne -1 do wdelete, !d.window
endif

IR1 = dblarr(n_elements(wavelist),n1)

if (n0 ge 2) then for k=0,n1-1 do IR1[*,k] = median(RESarr[*,k,3,*],dim=4,/even) ;takes the median over all nights in the month
if (n0 eq 1) then for k=0,n1-1 do IR1[*,k] = RESarr[*,k,3,0] ;just the values are returned

;wave = median(RESarr[*,*,0,*],dim=2)

if (plotting eq 'on' and n0 ge 2) then begin
    window,2,retain=2
    device,decomposed=0
    loadct,0
    plot,wavelist,IR1[*,0],psym=sym(1),symsize=0.5,xtitle='Wavelength',$
      ytitle='FWHM (Ang)',yrange=[4,7],/nodata,xrange=[3500,6000],/xstyle,$
      title='Median values of all nights in the res.list list...'
    loadct,33
    
    for k=0,n1-1 do begin
        oplot,wavelist,IR1[*,k],psym=sym(1),symsize=0.5,color=k
        wait,0.03
    endfor
endif
    
;***********************************************************************
; Here the Iron lines taken in sept 2010 are used to increase the number
; of lines used in the fit. This is all hardwired. If you want to
; leave this out of the calculation of IR, then un-comment the
; goto,jump2 and jump2 lines. Then, a couple lines below this section,
; change IRall to equal just IR1, NOT IRall = [IR1,IR2,IR3]]
;***********************************************************************
goto,jump2
arc1 = 'iron_arcpe.fits'
ptow1 = 'ptow_sept10_n1.dat'
mask1 = 'mask_sept10_n1.dat'

out1 = rescheckiron(arc1,ptow1,mask1)
wavelist = [3592.32,3731.13,4556.23,4620.78,4669.06,4776.03]
waveall = [waveall,wavelist]
IR2 = out1[*,*,3]

wset,2
for j=0,n1-1 do begin
    oplot,wavelist,out1[*,j,3],psym=sym(1),symsize=0.5,color=j
;    wait,0.1
endfor

;This is a list of Iron arc frames
readcol,'res_IRONarc.list',silent=1,f='a',list

wavelist = [3530.37,3669.39,3994.69,4494.57,4559.12,4607.37,4714.33]
waveall = [waveall,wavelist]
RESarrFE = dblarr(7,246,5,6)
for j=0,5 do begin
    out = rescheckiron2(list[j],'ptow_IRONarc.dat','mask_IRONarc.dat')
    RESarrFE[*,*,*,j] = out
endfor

IR3 = dblarr(7,246)
for k=0,n1-1 do IR3[*,k] = median(RESarrFE[*,k,3,*],dim=4,/even)

wset,2
for j=0,n1-1 do begin
    oplot,wavelist,IR3[*,j],psym=sym(1),symsize=0.5,color=j
;    wait,0.05
endfor
jump2:
;***********************************************************************
;IRall = [IR1,IR2,IR3] ;comment this out (and un-comment the next line) to skip including the iron lines from sept 11
;***********************************************************************
IRall = IR1
iw = bsort(waveall)
wavestd = waveall[iw]
IRallS = IRall[iw,*]
nw = n_elements(wavestd)

IRout = [[IRall],[waveall]]
writefits,outnameR+'.fits',IRout ;the raw fits are written out.
;dirty trick installed here...
IRout = readfits(outnameR+'.fits',header)
sxaddpar,header,'COMMENT','This file is an instrumental resolution map.'
sxaddpar,header,'COMMENT','It is in units of Angstroms (FWHM)'
sxaddpar,header,'COMMENT','This was created by calcRES.pro'
writefits,outnameR+'.fits',IRout,header ;the raw fits are written out again, with a better header

IRout = fltarr(nw)              ;and now the variable is recycled

test = readfits(arc[0],/silent)
n0 = n_elements(test[*,0]) * factor
waveall2 = fltarr(n0)
n0 = round(n0)
IRoutALL = fltarr(n0)
for j=0,n0-1 do waveall2[j] = w1 + (j * wdd)

free_lun,55
openw,55,'ptop_'+month_night+'.dat'
for j=0,n1-1 do begin
    p1 = poly_fit(wavestd,IRallS[*,j],1)
    py1 = p1[0] + wavestd*p1[1]
    ppy1 = p1[0] + p1[1]*waveall2
    p2 = poly_fit(wavestd,IRallS[*,j],2)
    py2 = p2[0] + p2[1]*wavestd + p2[2]*wavestd^2
    ppy2 = p2[0] + p2[1]*waveall2 + p2[2]*waveall2^2

    kms = IRallS[*,j]
    kms = kms/wavestd*299792.458/2.35
    pkms = poly_fit(wavestd,kms,2)
    ppkms = pkms[0] + pkms[1]*waveall2 + pkms[2]*waveall2^2

    p3 = poly_fit(wavestd,IRallS[*,j],3)
    py3 = p3[0] + p3[1]*wavestd + p3[2]*wavestd^2 + p3[3]*wavestd^3
    ppy3 = p3[0] + p3[1]*waveall2 + p3[2]*waveall2^2 + p3[3]*waveall2^3
    p4 = poly_fit(wavestd,IRallS[*,j],4)
    py4 = p4[0] + p4[1]*wavestd + p4[2]*wavestd^2 + p4[3]*wavestd^3 + p4[4]*wavestd^4
    ppy4 = p4[0] + p4[1]*waveall2 + p4[2]*waveall2^2 + p4[3]*waveall2^3 + p4[4]*waveall2^4

; added to build the ptop files.
    ii = intarr(nw)
    for k=0,nw-1 do begin
        i = where(waveall2 ge wavestd[k])
        ii[k] = i[0]
    endfor
;    pindex = poly_fit(ii,IRallS[*,j],2)
;    ppindex = pindex[0] + pindex[1]*fltarr(n0*factor) + pindex[2]*fltarr(n0*factor)^2
;    printf,55,[pindex[0],pindex[1],pindex[2]]
    pindex = poly_fit(ii,kms,2)
    ppindex = pindex[0] + pindex[1]*(1.0+findgen(n0*factor)) + pindex[2]*(1.0+findgen(n0*factor))^2
    printf,55,[pindex[0],pindex[1],pindex[2]]

    if (plotting eq 'on') then begin
        if (j eq 0) then wset,2
        loadct,0,/silent
;        plot,wavestd,IRallS[*,j],psym=2,yrange=[4,7],title='Fiber # '+strn(j+1),$
;          ytitle='FWHM (A)',xtitle='Wavelength (A)',/ys,xrange=[3500,6000],/xs
;        plot,ii,IRallS[*,j],psym=2,yrange=[4,7],title='Fiber # '+strn(j+1),$
;          ytitle='FWHM (A)',xtitle='Pixel',/ys,xrange=[-50,2100*factor],/xs
        plot,ii,kms,psym=2,yrange=[80,220],title='Fiber # '+strn(j+1),$
          ytitle='km/sec',xtitle='Pixel',/ys,xrange=[-50,2100*factor],/xs
        loadct,4,/silent
        oplot,waveall2,ppy1,color=60
        oplot,waveall2,ppy2,color=110
        oplot,waveall2,ppy3,color=180
        oplot,waveall2,ppy4,color=150
        oplot,findgen(n0*factor)+1.0,ppindex
        wait,0.05
    endif

    if (pfit eq 'first')  then IRout = [[IRout],[py1]]
    if (pfit eq 'second') then IRout = [[IRout],[py2]]
    if (pfit eq 'third')  then IRout = [[IRout],[py3]]
    if (pfit eq 'fourth') then IRout = [[IRout],[py4]]
    if (pfit eq 'first')  then IRoutALL = [[IRoutALL],[ppy1]]
    if (pfit eq 'second') then IRoutALL = [[IRoutALL],[ppy2]]
    if (pfit eq 'third')  then IRoutALL = [[IRoutALL],[ppy3]]
    if (pfit eq 'fourth') then IRoutALL = [[IRoutALL],[ppy4]]
endfor
free_lun,55

if (plotting eq 'on') then begin
    print,'The next ENTER deletes the plot(s)...'
    while !d.window ne -1 do wdelete, !d.window
endif

IRout = IRout[*,1:*]            ;first row of 0's is dropped
IRout2 = [[IRout],[wavestd]]    ;the wavelength is added.
writefits,outnameF+'.fits',IRout2,header ;the poly fit values are written out.
; the resmapF routine
;if (plot3d eq 'on') then t = resmapF_v2(wavestd,IRout,smoothingF,boxsizeF,outnameF,outnameSF)
;if (plot3d eq 'on') then t = resmapF_v2(wavestd,IRallS,smoothingR,boxsizeR,outnameR,outnameSR)

IRoutALL = IRoutALL[*,1:*]      ;first row of 0's is dropped

writefits,outnameALL+'.fits',IRoutALL
IRoutALL = readfits(outnameALL+'.fits',header,/silent)
sxaddpar,header,'CRVAL1',w1,' The initial wavelength',after='DATE'
sxaddpar,header,'CDELT1',wdd,' The linear dispersion term',after='DATE'
sxaddpar,header,'COMMENT','This file is an instrumental resolution map.'
sxaddpar,header,'COMMENT','It is in units of A (FWHM)'
sxaddpar,header,'COMMENT','This was created by calcRES.pro'
writefits,outnameALL+'.fits',IRoutALL,header ;the poly fit for the entire wavelength regions are written out.

; the resmapF routine
if (plot3d eq 'on') then t = resmapf_v2(waveall2,IRoutALL,smoothingF,boxsizeF,outnameALL,outnameSALL,month_night)

stop
END

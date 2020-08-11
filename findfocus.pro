PRO findfocus, fitslist

; This is a stripped-down version of fiber_trans3, and should be used
; just to find the best focus. See the notes in fiber_trans3 and
; fiber_trans4 for more details on how to run this code. BUT, IN
; SHORT, THE INPUT IS THE SAME AS THE OTHER ROUTINES.

focustep = 0.030 ;this value sets the step size for the range of focus. If the horozontal lines exploring the spot diameter in the focus routine are too close, raise this number. If they are too far apart, lower this number.
fwhm = 52.0 ; (22.0 for the photometry bench) the FWHM in pixels...? If all your fibers are being found, this parameter has done it's job. It feeds into the find.pro routine.
hmin = 4000 ;The detection threshold, in pixel values. This is the counts above the background to consider a peak a peak. This is only used for the initial determination of peaks
sharplim = [0.001,5.0]
roundlim = [-5.0,5.0] ;To get the elongated fibers in the corners, the roundness criteria is loosened.
pitch = 109 ;MUST BE AN ODD INTEGER. An estimate of the pitch between fibers. This number is used for the fiber centriod finder and for sorting the fibers into order from 1 to nfibers (35 was the pitch for the Nikon lenses set up)
apertureR = 57.0 ;MUST BE AN ODD INTEGER. This aperture is the outer range when fitting circular apertures.
R1c = 12.0 ;used in the circular aperture photometry routine. The amount added to the major axis to include in the total counts. Also the inner radius for the background estimate
R2c = 18.0 ;used in the circular aperture photometry routine. The outer radius for the estimate of the background.
R1e = 5.0 ;used in the elliptical aperture photometry routine. The amount added to the major axis to include in the total counts. Also the inner radius for the background estimate
R2e = 12.0 ;used in the elliptical aperture photometry routine. The outer radius for the estimate of the background.

lop = 0.75 ;The percentage of the peak fiber value to set as the cutoff limit for flat-topping the fibers. As the fibers are not top-hat functions, this sets all the pixel values above this value to the same value, effectively flatten the fibers. This helps in determining the centers of the fibers as the ellipse-fitting routine used (fit_ellipse.pro) is a mass-weighting routine which would be skewed by lop-sided fiber intensities.
cutmagcent = 0.60 ; cutmagcent + (cutstepcent * 10) < 1.00 (This is the lower starting point when running the center-finding routine.)
cutstepcent = 0.02 ;the steps to take in the iteration of the center
cutmagell = 0.25 ;the starting point for the ellipse-fitting routine.
cutstepell = 0.005 ;the step size for the ellipse-fitting routine
wscale = 4.0       ;just scales plotting window sizes. Turn up for larger plots and down for smaller ones.
nf = 224 ;The number of fibers that should exist in the frame. IF THERE ARE DEAD FIBERS, DO NOT REDUCE THIS NUMBER! Just put their # into the 'deadfibers' array above.
flatfield = 'no' ;set to 'yes' to use the flat-fielding, set to 'no' to turn it off (I am dubious of using this...)
dxmin = 140.0 ;The minimum radial spacing between fibers. This # (and dymin) is used to determine which fibers are neighbors in order to estimate the plate scale of the imager.
dymin = 30.0
;*********************************************************************

APPTYPE = 'circle'
DARKSWITCH = 'on'
hp = floor(pitch/2.0)

readcol,fitslist,f='a,a,a',silent=1,fitsfile,darkfile,deadfile
nframes = n_elements(fitsfile)

counter = [0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1]

for all = 0,nframes-1 do begin    ;a loop through each individual fits frame
   print,''
   print,'Working on frame '+fitsfile[all]+'...'
   print,''
   
   readcol,deadfile[all],silent=1,f='i',deadfibers
   
   if (deadfibers[0] ne -1) then begin
       idead = deadfibers - 1.0
       ndead = n_elements(deadfibers)
       nfibers = nf - ndead
   endif else begin
       nfibers = nf
       ndead = 0
   endelse

   image = readfits(fitsfile[all],/silent,header)
   image = float(image)
   dark = readfits(darkfile[all],/silent,/noscale)
   dark = float(dark)
   image = image - dark
   if (max(image[100:200,100:200]) gt 30000.0) then image = image - 32768.0 ;a fudge to capture some funny-busniess which depends on the system you're running this on
   
   fiberarray = fltarr(6,nf)
   
   tempname = strsplit(fitsfile[all],'.',/extract)
   nameout4 = tempname[0]+'_focus.txt'

;*********************************************************************
;********* The beginning of the stripped find.pro code ***************
;*********************************************************************

   maxbox = 31                  ;Maximum size of convolution box in pixels 
   type = size(image)
   n_x  = type[1] & n_y = type[2]
   radius = 0.637 * FWHM > 2.001 ;Radius is 1.5 sigma
   radsq = radius^20
   nhalf = fix(radius) < (maxbox-1)/2 ;
   nbox = 2*nhalf + 1                 ;# of pixels inside of convolution box 
   middle = nhalf                     ;Index of central pixel
   lastro = n_x - nhalf
   lastcl = n_y - nhalf
   sigsq = ( fwhm/2.35482 )^2
   mask = bytarr( nbox, nbox )  ;Mask identifies valid pixels in convolution box 
   c = fltarr( nbox, nbox )     ;c will contain Gaussian convolution kernel
   
   dd = indgen(nbox-1) + 0.5 - middle ;Constants need to compute ROUND
   dd2 = dd^2
   w = 1. - 0.5*(abs(dd)-0.5) / (middle-.5)   
   ir = (nhalf-1) > 1
   
   row2 = (findgen(Nbox)-nhalf)^2
   
   for i = 0, nhalf do begin
      temp = row2 + i^2
      c[0,nhalf-i] = temp         
      c[0,nhalf+i] = temp                           
   endfor
   
   mask = fix(c LE radsq)       ;MASK is complementary to SKIP in Stetson's Fortran
   good = where( mask, pixels)  ;Value of c are now equal to distance to center
   
   c = c*mask               
   c[good] = exp(-0.5*c[good]/sigsq) ;Make c into a Gaussian kernel
   sumc = total(c)
   sumcsq = total(c^2) - sumc^2/pixels
   sumc = sumc/pixels
   c[good] = (c[good] - sumc)/sumcsq
   c1 = exp(-.5*row2/sigsq)
   sumc1 = total(c1)/nbox
   sumc1sq = total(c1^2) - sumc1
   c1 = (c1-sumc1)/sumc1sq
   sumc = total(w)              ;Needed for centroid computation
   
   h = convol(float(image),c)   ;Convolve image with kernel "c"
   
   h[0:nhalf-1,*] = 0 & h[n_x-nhalf:n_x-1,*] = 0
   h[*,0:nhalf-1] = 0 & h[*,n_y-nhalf:n_y-1] = 0
   
   mask[middle,middle] = 0      ;From now on we exclude the central pixel
   pixels = pixels -1           ;so the number of valid pixels is reduced by 1
   good = where(mask)           ;"good" identifies position of valid pixels
   xx= (good mod nbox) - middle	;x and y coordinate of valid pixels 
   yy = fix(good/nbox) - middle ;relative to the center
   offset = yy*n_x + xx
   
   SEARCH:                      ;Threshold dependent search begins here
   
   index = where(h GE hmin, nfound) ;Valid image pixels are greater than hmin
   if nfound lt nfibers then begin  ;Any maxima found?
      print,'WARNING! Not enough peaks found!'
      print,'Adjust your threshold value (hmin) and try again!'
      goto,FINISH    
   endif
   
   for i= 0L, pixels-1 do begin                             
      stars = where (h[index] GE h[index+offset[i]], nfound)
      if nfound lt nfibers then begin ;Do valid local maxima exist?
         print,'WARNING! Not enough peaks found!'
         print,'Adjust your threshold value (hmin) and try again!'
         goto,FINISH    
      endif
      index = index[stars]
   endfor 
   
   ix = index mod n_x           ;X index of local maxima
   iy = index/n_x               ;Y index of local maxima
   ngood = N_elements(index)       
   message,strtrim(ngood,2)+' local maxima located above threshold',/INF
   
   nstar = 0L                   ;NSTAR counts all stars meeting selection criteria
   badround = 0L & badsharp=0L  &  badcntrd=0L
   
   x = fltarr(ngood) & y = x
   
;  Loop over star positions; compute statistics
   for i = 0L,ngood-1 do begin   
      temp = float(image[ix[i]-nhalf:ix[i]+nhalf,iy[i]-nhalf:iy[i]+nhalf])
      d = h[ix[i],iy[i]]        ;"d" is actual pixel intensity        
      
;  Compute Sharpness statistic
      sharp1 = (temp[middle,middle] - (total(mask*temp))/pixels)/d
      if ( sharp1 LT sharplim[0] ) or ( sharp1 GT sharplim[1] ) then begin
         badsharp = badsharp + 1
         goto, REJECT           ;Does not meet sharpness criteria
      endif
      
;   Compute Roundness statistic
      dx = total( total(temp,2)*c1)   
      dy = total( total(temp,1)*c1)
      if (dx LE 0) or (dy LE 0) then begin
         badround = badround + 1
         goto, REJECT           ;Cannot compute roundness
      endif
      
      around = 2*(dx-dy) / ( dx + dy ) ;Roundness statistic
      if ( around LT roundlim[0] ) or ( around GT roundlim[1] ) then begin
         badround = badround + 1
         goto,REJECT            ;Does not meet roundness criteria
      endif

; Find X centroid
      derivat = shift(temp,-1,0) - temp
      derivat = total( derivat[0:nbox-2,middle-ir:middle+ir],2)
      sumd = total(w*derivat)
      sumxd = total(w*dd*derivat)
      sumxsq = total(w*dd2) 
      
      if ( sumxd GE 0. ) then begin
         badcntrd = badcntrd + 1
         goto,REJECT            ;Cannot compute X centroid
      endif
      
      dx =sumxsq*sumd/(sumc*sumxd)
      if abs(dx) GT nhalf then begin
         badcntrd = badcntrd + 1
         goto,REJECT            ;X centroid too far from local X maxima
      endif
      
      xcen = ix[i]-dx           ;Convert back to big image coordinates
      
; Find Y centroid                 
      derivat = shift(temp,0,-1) - temp 
      derivat = total( derivat[middle-ir:middle+ir,0:nbox-2], 1 )
      sumd = total( w*derivat )
      sumxd = total( w*dd*derivat )
      sumxsq = total( w*dd2 )
      if (sumxd GE 0) then begin
         badcntrd = badcntrd + 1
         goto, REJECT  
      endif
      
      dy = sumxsq*sumd/(sumc*sumxd)
      if ( abs(dy) GT nhalf ) then begin
         badcntrd = badcntrd + 1
         goto,REJECT 
      endif
      
      ycen = iy[i] - dy
      
;  This star has met all selection criteria.  Print out and save results
      x[nstar] = xcen & y[nstar] = ycen
      nstar = nstar+1
      
REJECT: 
      
   endfor
   
   nstar = nstar-1              ;NSTAR is now the index of last star found
   
   print,' No. of sources rejected by SHARPNESS criteria',badsharp
   print,' No. of sources rejected by ROUNDNESS criteria',badround
   print,' No. of sources rejected by CENTROID  criteria',badcntrd
   
   if nstar LT 0 then return    ;Any stars found?
   
   x = x[0:nstar]
   y = y[0:nstar]               ;these are still in index form, although they are not integers (###.### precision)
   
   nx = n_elements(x)
   
;*********************************************************************
;*************** The end of the stripped find.pro code ***************
;*********************************************************************

   window,0,retain=2,xsize=(n_x/wscale),ysize=(n_y/wscale)
   device,decomposed=0
   loadct,0,/silent
   plot,x,y,/nodata,xtitle='Pixels',ytitle='Pixels',title='Centers for '+fitsfile[all]
   for j=0,nx-1 do plots,x[j],y[j],psym=2
   loadct,4,/silent
   wait,1.0

   xout = fltarr(nfibers) & yout = xout
   cut = pitch * 0.4       ;used to pick off fibers found more than once by find.pro
   
   cntr = 0
   for j=0,nx-1 do begin ;the fibers are now sorted out to reject duplicate versions of the same fiber found earlier.
      xt = x[j]
      yt = y[j]
      
      ib = where(xt-cut le x and xt+cut gt x and yt-cut le y and yt+cut gt y,cb)
      if (cb eq 1) then begin
         xt = x[ib]
         yt = y[ib]
      endif else begin
         xt = median(x[ib],/even)
         yt = median(y[ib],/even)
      endelse
      test = where(xout gt xt-cut and xout lt xt+cut and yout gt yt-cut and yout lt yt+cut,count)
      if (count eq 0) then begin
         xout[cntr] = xt
         yout[cntr] = yt
         cntr = cntr + 1
      endif
   endfor
   
   for j=0,nfibers-1 do plots,xout[j],yout[j],psym=sym(1),color=150
   print,'Did all your fibers get found???'
   print,'(if not, change some parameters and try again!)'
   wait,2

;the fibers are sorted from 1 to nfibers (the dead fibers are not in this array yet...)
   iy = reverse(bsort(yout))
   ytemp = yout[iy]
   xtemp = xout[iy]
   youtS = 0.0
   xoutS = 0.0
   for i=0,nfibers-2 do begin
      iy1 = i
      if (i eq 0) then iy1l = iy1
      iy2 = i+1
      diff = abs(ytemp[iy1] - ytemp[iy2])
      if (diff gt pitch*0.5) or (i eq nfibers-2) then begin
         iy1 = iy1l
         iy1l = iy2
         if (i lt nfibers-2) then iy2 = iy2-1
         x1 = xtemp[iy1:iy2] ;the x values for a row of fibers is identified
         y1 = ytemp[iy1:iy2]
         ix = bsort(x1) ;they are sorted, from low to high, which places them in their proper fiber order
         youtS = [youtS,y1[ix]]
         xoutS = [xoutS,x1[ix]]
      endif
   endfor
   youtS = youtS[1:*] ;the leading zero's are dropped
   xoutS = xoutS[1:*] ;these ordered arrays now become the primary array so that all that follows is in order.
   
   if (ndead gt 0) then begin ;the dead fibers are re-inserted back into the xoutS and youtS arrays
      for j=0,ndead-1 do begin
         idead1 = idead[j]
         px1 = xoutS[0:idead1-1]
         px2 = xoutS[idead1:*]
         xoutS = [px1,-1.0,px2]
         py1 = youtS[0:idead1-1]
         py2 = youtS[idead1:*]
         youtS = [py1,-1.0,py2]
      endfor
   endif

;the radius array for this frame is created
   x2 = indgen(hp)+1
   x1 = reverse(indgen(hp+1))
   xline = [x1,x2]
   y2 = indgen(hp)+1
   y1 = reverse(indgen(hp+1))
   yline = transpose([y1,y2])
   
   Xarray = fltarr(pitch,pitch)
   Yarray = Xarray
   
   FOR i=0,pitch-1 DO Xarray(*,i) = xline
   FOR i=0,pitch-1 DO Yarray(i,*) = yline
   RADarr = sqrt(Xarray^2 + Yarray^2)
   iouter = where(RADarr gt apertureR-2.0)
   
   deadcntr = 0
   pcntr = 0.0
   for j=0,nf-1 do begin ;a loop through each of the fibers
      l = j + deadcntr
       print,'Working on fiber #'+strn(j+1)
      if ((ndead gt 0) and (ndead ne deadcntr)) then if (j eq idead[deadcntr]) then begin
          print,'Fiber #'+strn(idead[deadcntr]+1.0)+' is totally dead!!!'
         deadcntr = deadcntr + 1
         goto, deadjump
      endif
      x1 = uint(round(xoutS[j])-hp)
      x2 = uint(round(xoutS[j])+hp)
      y1 = uint(round(youtS[j])-hp)
      y2 = uint(round(youtS[j])+hp)
      chunk = image[x1:x2,y1:y2]
      chunkorig = chunk
      chunk[iouter] = 0.0 ;the values outside of the outer aperture are set to zero for the fiber location. This helps the fit_ellipse routine not get fooled by light from a neighboring fiber.
      if (j eq 0) then begin
         nn0 = n_elements(chunk[*,0])
         nn1 = n_elements(chunk[0,*])
         ssi = floor((nn0*nn1*lop)) ;the threshold for the top-hat of the fiber profile is set.
      endif
      schunk = chunk[bsort(chunk)]
      cut1 = schunk[ssi]
      ihigh = where(chunk ge cut1)
      chunk[ihigh] = cut1 ;the top is loped off the fiber to improve the center-finding routines.
      
      center = fltarr(2,11)
      pa = fltarr(11)
      axes = fltarr(2,11)
      for k=0,10 do begin       ;an iteration on the fit_ellipse routine.
         cutoff1 = (cut1 * (cutmagcent + (cutstepcent * float(k))))
         index1 = where(chunk gt cutoff1)
         cutoff2 = (cut1 * (cutmagell + (cutstepell * float(k))))
         index2 = where(chunk gt cutoff2)
         circle  = Fit_Ellipse(index1, XSize=nn0, YSize=nn0, CENTER=c1)
         ellipse = Fit_Ellipse(index2, XSize=nn0, YSize=nn0, orient=p1,axes=a1)
         center[0,k] = c1[0] & center[1,k] = c1[1]
         pa[k] = p1
         axes[0,k] = a1[0] & axes[1,k] = a1[1]
      endfor

      centerM = median(center,dim=2,/even)
      paM = median(pa,/even)
      axesM = median(axes,dim=2,/even)
; the [xy]out gives coordinates in the native frame. centerM is in the
; trimmed frame.
      fiberarray[0,j] = float(x1) + float(centerM[0]) + 1.0
      fiberarray[1,j] = float(y1) + float(centerM[1]) + 1.0
      fiberarray[2,j] = paM
      fiberarray[3,j] = axesM[0]
      fiberarray[4,j] = axesM[1]
      xx = fiberarray[0,j]
      yy = fiberarray[1,j]

      if (counter[all] eq 0) then begin
         Xfirst = fiberarray[0,0]
         Yfirst = fiberarray[1,0]
         fiberarray[5,j] = sqrt((xx-Xfirst)^2 + (yy-Yfirst)^2)
      endif
      if (counter[all] eq 1) then begin
         fiberarray[5,j] = sqrt((xx-Xfirst)^2 + (yy-Yfirst)^2)
      endif

; THE PRIMARY FOCUS ROUTINE

      schunk = chunkorig[bsort(chunkorig)]
      cut2 = median(schunk[n_elements(schunk)-1500,*],/even)
      diameter = fltarr(31)
      diaARR = fltarr(60)

      if (j eq 0) then begin
         window,5,retain=2,xsize=500,ysize=500,xpos=430,ypos=50
         window,15,retain=2,xsize=500,ysize=500,xpos=950,ypos=50         
      endif

      for k=1,31 do begin
         cutoff3 = cut2 - (cut2 * (focustep * k))
         index3 = where(chunkorig gt cutoff3)
         ellipse = Fit_Ellipse(index3, XSize=nn0, YSize=nn0,CENTER=c1)
         if (j eq pcntr) then begin ;plots the cross-sections and where the diameter line falls
            if (k eq 1) then begin
               xslice = median(chunkorig[*,c1[1]-3:c1[1]+3],dim=2)
               yslice = median(chunkorig[c1[0]-3:c1[0]+3,*],dim=1)
               wset,5
               plot,xslice,title='X Profile'
               oplot,[0,n_elements(xslice)],[cutoff3,cutoff3],thick=1
               wset,15
               plot,yslice,title='Y Profile'
               oplot,[0,n_elements(yslice)],[cutoff3,cutoff3],thick=1
            endif else begin
               wset,5
               oplot,[0,n_elements(xslice)],[cutoff3,cutoff3],thick=1
               wset,15
               oplot,[0,n_elements(yslice)],[cutoff3,cutoff3],thick=1
               wait,0.01
            endelse
         endif
         
         for ll=0,59 do begin ;makes an estimate of the diameter of the spot, based on the ellipse fitting.
            x1 = ellipse[0,ll]
            y1 = ellipse[1,ll]
            x2 = ellipse[0,ll+59]
            y2 = ellipse[1,ll+59]
            dx = double(abs(x1-x2))
            dy = double(abs(y1-y2))
            diaARR[ll] = sqrt(dx^2 + dy^2)
         endfor
         diameter[k-1] = median(diaARR,/even)
      endfor
      
      MEDdia = median(diameter)
      SDdia = stddev(diameter)
      if (j eq 0) then begin
         free_lun,56
         openw,56,nameout4
      endif 
      printf,56,MEDdia,' ',SDdia
      if (j eq pcntr) then pcntr = pcntr + 10.0
      
      deadjump:
   
   endfor
endfor                          ;the end of the loop over each fiber found
free_lun,56                     ;closing the focus estimate file

FINISH:
stop
END

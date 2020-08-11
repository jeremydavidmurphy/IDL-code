;THIS VERSION OF REALIGN IS FOR COLLAPSED DATA. AS THE NEW VERSION
;DOES ALL THE FIBER MASKING BEFORE ENTERING THIS FUNCTION, 

;-----------------------------------------------------------------
pro realignc, datain, ptowframe, wi, wd, forcefit, PLOT=p
;-----------------------------------------------------------------

;this code takes the place of runex. it's accepts collapsed data files
;and the corresponding ptow.dat files and makes a linear interpolation
;of the data so as to align the wavelengths of various files.

;DATAIN is an array the same size as the ptow.dat array and is the
;data that gets interpolated. The PTOWFRAME is a string and will search
;in the calling directory for the corresponding ptow.dat text file.
;WI is the initial wavelength and WD is the dispersion you want to use
;for your interpolation.

;FORCEFIT is a string that is either 'force' or 'free' 
;THIS IS SIGNIFICANT: When set to 'force' the length of the
;returned array matches that of the input array. Then set to 'free'
;the returned array will be whatever length is needed to complete the
;interpolation at that dispersion scale. The values outside of the
;interpolation are fixed at the same value as the last one.

;As the ptow.dat files skip all fibers that were labeled dead in the
;VACCINE reduction it's essential that the input data files have only
;-666 rows for the dead fibers.

;THE DATA HAS TO HAVE IT'S -666 FLAGS FOR THE DEAD FIBERS IN PLACE BUT
;NOT COMBINED WITH THE .R FILES (AS THIS CREATES MORE -666 ROWS AND
;WILL THROW OFF THE WAVELENGTH SOLUTION).

data = readfits(datain,/silent)
datain = data
n1 = n_elements(datain[*,0])
n2 = n_elements(datain[0,*])

tm = median(datain,dimension=1)
badfibers = where(tm eq -666)

;Here's a test to confirm the input files are aligned properly. As VP1
;is missing 6 fibers and VP2 none, if badfibers isn't 6 or 0 elements
;long then you receive a warning BUT THE ROUTINE CONTINUES! This can
;be changed with replacing "STOP" with "PAUSE" in the if statement
;below.

;if(n_elements(badfibers) ne 6) and (n_elements(badfibers) ne 1) then begin
;    print,'There is a problem with your input data!'
;    print,'There are too many fibers with -666 flags for frame '
;    print,'Fix this or the wavelength solution will not be correct...'
;    pause
;    stop
;endif

wr = dblarr(n1,n2)
readcol,ptowframe,F='D,D,D,D,D',c0,c1,c2,c3,c4
ptow = [[c0],[c1],[c2],[c3],[c4]];this is a 247x5 array, NOT a 5x247 array!

for k=0,n2-1 do begin
    for j=0,n1-1 do wr[j,k] = c0[k]+c1[k]*j+c2[k]*j*j+c3[k]*j*j*j+c4[k]*j*j*j*j
endfor

jump1:
if (forcefit eq 'free') then begin
    tup = max(wr)
    n3 = floor((tup-wi)/wd)
    wo = dblarr(n3)
    print,'The output frame will have '+strn(n3)+' elements.'
    dataout = dblarr(n3,n2)
    test = dataout
    for j=0,n3-1 do wo[j] = wi + (j*wd)
endif
if (forcefit eq 'force') then begin
    wo = dblarr(n1)
    dataout = dblarr(n1,n2)
    test = dataout
for j=0,n1-1 do wo[j] = wi + (j*wd)
endif

if (forcefit ne 'force') and (forcefit ne 'free') then begin
    ans=''
    print,'Force the interpolation to '+strn(n2)+' elements?'
    print,'Type FORCE or FREE:'
    read,ans
    goto, jump1
endif

n0 = n_elements(p)
if (n0 eq 0) then p = 'noplot'

cntr = 0
for k=0,n2-1 do begin
    if (k eq badfibers[cntr]) then begin
        dataout[*,k] = -666
        if (cntr+1 ne n_elements(badfibers)) then cntr = cntr+1
    endif else begin
;The only line that does anything meaningful...........
        dataout[*,k] = interpol(datain[*,k],wr[*,k],wo)
;*************************************************************
        if (p eq 'plot') then begin
            loadct,4
            window,0,retain=2,ysize=460,ypos=600
            device,decomposed=0
            plot,wr[*,k],datain[*,k],xrange=[3800,4000],xstyle=1,$
              title='Data for fiber # '+strn(k+1),$
              color=225,xtitle='Original data in red, interpolated data plotted in yellow.',charsize=1.5,/nodata
            oplot,wr[*,k],datain[*,k],color=150
            oplot,wo,dataout[*,k],color=255
            window,1,retain=2,ysize=460,ypos=50
            device,decomposed=0
            plot,wr[*,k],datain[*,k],xrange=[3800,5400],$
              xstyle=1,title='Data for fiber # '+strn(k+1),$
              color=225,xtitle='Original data in red, interpolated data plotted in yellow.',charsize=1.5,/nodata
            oplot,wr[*,k],datain[*,k],color=150
            oplot,wo,dataout[*,k],color=255
            pause
        endif
    endelse
endfor        
;********************************************************************
;The data has -666 flags from cosmic ray rejection. The output then
;will have negative values interpolated between -666 and some positive
;value. To correct for this, any negative value below some threshold
;is set equal to -666
cleaned = dataout
;cleaned[where(cleaned lt -30)] = -666
;********************************************************************

name = ''
print,'Enter the name of the output file:'
read,name
writefits,name,cleaned
stop
end

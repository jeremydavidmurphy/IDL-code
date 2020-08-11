;THIS ROUTINE IS USED TO COMBINE VIRUS-P DATA OF SEVERAL EXPOSURES AT
;THE SAME POINTING. IT READS IN A LIST OF FILE NAMES, THEN APPLIES THE
;WEIGHTED BIWEIGHT TO
;THE DATA FOR EACH FIBER. THE DATA IS NORMALIZED TO ACCOUNT FOR THE
;DIFFERENT LEVELS OF OVERALL INTENSITY BETWEEN EXPOSURES.

PRO combine_v3, WEIGHT=wgt

list=''
print,'Enter the list name:'
print,'(ENTER will look for "combine.list")' 
read,list
if (list eq '') then list='combine.list'
readcol,list,Format='A',files

;The size of the data set is determined
test = readfits(files[0])
n1 = n_elements(test[*,0]);wavelength
n2 = n_elements(test[0,*]);# of fibers
n3 = n_elements(files);# of files

data = dblarr(n1,n2,n3)
out = dblarr(n1,n2)
slice = dblarr(n1,n3)
wgtarray = dblarr(n3,n2)

if (n_elements(wgt) eq 0) then begin
    swch = 'on'
    wgt = dblarr(n3)
endif

if (n_elements(wgt) ne 0) then swch = 'off'

;the data is read in
for k=0,n_elements(files)-1 DO data[*,*,k] = readfits(files[k])

biwtout = dblarr(n1,n2)

ans=''
print,'Plot the data as it comes out? (y or n)'
print,'This shows plots of all the data being combined as well as the resulting output from the biweight.'
read,ans

for k=0,n2-1 do begin ;This loops on each fiber
    print,'Working on fiber '+strn(k+1)
    for l=0,n3-1 do begin
        slice[*,l] = data[*,k,l] ;all the exposures for one fiber are defined
        if (swch eq 'on') then begin
            wgt[l] = median(slice[n1-400:n1-100,l],/even)
            wgtarray[l,k] = wgt[l]
        endif
        if (swch eq 'off') then begin
            wgtarray[l,k] = wgt[l]
        endif
    endfor
    
    biwtout[*,k] = wbiweight(slice,wgt)

    if(ans eq 'y') then begin
        window,0,retain=2
        device,decomposed=0
        loadct,0
        plot,slice[30:n1-30,0],title='Fiber #: '+strn(k+1)+' The real data- the biweight is plotted in red.',$
          xtitle='Pixel',ytitle='ADU',/nodata
        loadct,4
        for j=0,n3-1 do oplot,slice[30:n1-30,j],color=(60+j*10)
        oplot,biwtout[30:n1-30,k],color=150,thick=2
        
        window,2,retain=2
        device,decomposed=0
        loadct,0
        plot,slice[30:n1-30,0],title='Fiber #: '+strn(k+1)+' The scaled data- the biweight is plotted in red.',$
          xtitle='Pixel',ytitle='ADU',/nodata
        loadct,4
        for j=0,n3-1 do oplot,normed[30:n1-30,j],color=60+j*10
        oplot,biwtout[30:n1-30,k],color=150,thick=2
        pause
        print,'The next ENTER deletes the plots.'
        pause
        wdelete,0,2
    endif

endfor

print,'Name the output frame (w/o .fits):'
read,ans
nam = ans+'.fits'
writefits,nam,biwtout
print,'The frame has been written to "combine.fits".'

STOP
END

;first, a slice of data is considered, composed of all the data that
;going to get combined for an individial fiber. this is labeled as
;'slice'.
;then 'boxed' is created. this is a smoothed version of each spectrum
;(it has the same dimensions as 'slice'). this is done with
;medbox.pro, which returns either a median or average boxcar.
;the data (slice) is then divided by 'boxed' and called 'normed'. this
;is done as there can be a wide variety of flux levels, with the
;biweight rejecting high and low values. in this version, the actual
;data is returned by the biweight. THIS IS THE PRIMARY DIFFERENCE
;BETWEEN THIS VERSION AND COMBINE_V2.
;If the biweight returns a value of 0.0, then this is put directly
;into the output file for that fiber. (this should occur only for dead
;fibers and the rare case where the data was masked at the same pixel
;for all data files.

;if the biweight is not equal to zero, it gets multiplied by the
;median (or average, depending upon which switch you flip) of the
;smoothed spectra. Perhaps this step is unnecessary. It's only job is
;to return the overall values to roughly their original ones, rather
;than being normalized to one.

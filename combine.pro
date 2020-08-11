PRO combine, med_ave

;This routine is used to combine data from several exposures. It uses
;a median boxcar smoothing (via medbox.pro) to weight, then combine
;each individual spectra. IT ASSUMES YOU HAVE ALIGNED THE WAVELENGTH
;OF THIS SPECTRA, THEN RESTACKED THEM.

;med_ave needs to be either 'median' or 'average' and dictates whether
;the 


boxsize=50;This determines the boxsize for the smoothing
list=''
num=intarr(1)
print,'Enter the number of frames you are combining:'
read,num
print,'Enter the list name:'
read,list
files = strarr(num)

openr,5,list
readf,5,files
free_lun,5

print,'The files to be combined are:'
print,files
pause

;The size of the data set is determined
test = readfits(files[0])
data = dblarr(n_elements(test[*,0]),n_elements(test[0,*]),n_elements(files))
out = dblarr(n_elements(test[*,0]),n_elements(test[0,*]))

;the data is read in
FOR k=0,n_elements(files)-1 DO data[*,*,k] = readfits(files[k])

slice = dblarr(n_elements(test[*,0]),n_elements(files));
boxed = slice
normed = boxed
meannorm = dblarr(n_elements(test[*,0]));wavelength
mn = meannorm

n1 = n_elements(data[0,*,0]);247
n2 = n_elements(normed[0,*]);# of files
n3 = n_elements(normed[*,0]);wavelength

unboxave = dblarr(n3)
unboxmed = unboxave

FOR k=0,n1-1 DO BEGIN ;This loops on each fiber
;FOR k=123,123 DO BEGIN ;This loops on each fiber
    print,'Working on fiber '+strn(k+1)
    FOR l=0,n2-1 DO slice[*,l] = data[*,k,l];all the exposures for one fiber are defined
    FOR l=0,n2-1 DO BEGIN ;the normalized spectra is calculated for fiber k

;********************************************************************************************************
        boxed[*,l] = medbox(slice[*,l],boxsize,2);1 is for median boxcar and 2 is for average boxcar
;********************************************************************************************************

        normed[*,l] = slice[*,l]/boxed[*,l] ;this is the normed spectra that needs biweight combination. It's all the spectra for a given fiber for one pointing.

    ENDFOR

    FOR j=0,n3-1 DO unboxave[j] = total(boxed[j,*])/n2
    FOR j=0,n3-1 DO unboxmed[j] = median(boxed[j,*],/EVEN)
    biwt = biweight(normed,'data') ;This is the biweight rejected, normalized spectra

;*********************************************************************************************
;SWITCH HERE TO EITHER TAKING THE AVERAGE (unboxave) OR THE MEDIAN (unboxmed) OF THE BOXCAR
    FOR j=0,n3-1 DO BEGIN
        IF (biwt[j] EQ 0.0) THEN BEGIN
            out[j,k] = 0.0 
;print,'-666 flag for pixel '+strn(j+1)+' at fiber '+strn(k+1)
;print,data[j,k,*]
        ENDIF ELSE BEGIN
            out[j,k]=biwt[j]*unboxmed[j]
        ENDELSE
    ENDFOR
;*********************************************************************************************
ENDFOR

writefits,'combine.fits',out
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

; THIS ROUTINE READS IN A LIST OF FITS FILES AND MEDIAN COMBINES THEM. THE
; HEADER FOR THE FIRST FITS FILE IS USED FOR ALL.

; The output is written out as median.fit unless the keyword is used

PRO medfits, list, OUTNAME=outname

readcol,list,silent=1,f='a',files
n0 = n_elements(files)
test = readfits(files[0],header,/silent)

nn1 = n_elements(test[*,0])
nn2 = n_elements(test[0,*])

array = dblarr(nn1,nn2,n0)

for j=0,n0-1 do begin
   frame = readfits(files[j],/silent)
   array[*,*,j] = frame
endfor

outarray = median(array,dim=3,/even)
sxaddpar,header,'COMMENT  ','= This frame made via medfits.pro'
if (n_elements(outname) eq 0) then writefits,'median.fits',outarray,header else $
writefits,outname,outarray,header

stop
END

; This routine is used to conduct a quick reduction and summation of
; the total counts over a given wavelength region for VIRUS-P data.

PRO sumover, datalist, wi, disp

readcol,datalist,format='a',files

test = readfits(files[0])
n0 = n_elements(test)
n1 = n_elements(files)

data = fltarr(n0,n1)
for j=0,n1-1 do data[*,j] = readfits(files[j],/silent)

wave = fltarr(n0)
for j=0,n0-1 do wave[j] = wi + j*disp

stop
END

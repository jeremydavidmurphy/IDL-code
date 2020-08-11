PRO rmsplot, rms, number='number'

;NUMBER: either 'one' or 'more' For ONE, rms is a single rms file. If
;MORE then rms is a list of rms files.

if (n_elements(number) eq 0) then begin
    number = ''
    print,'Enter ONE if a single file and MORE if a list of files:'
    read,number
endif

if (number eq 'one') then begin
    readcol,rms,format='a,f,x,x,f',name,rms,other
    goto,jump1
endif

if (number eq 'more') then readcol,rms,f='a',list
n0 = n_elements(list)
readcol,list[0],f='a,x,x,x,x',test
n1 = n_elements(test)
rmsarr = fltarr(n0,n1)

for j=0,n0-1 do begin
    readcol,rms,format='a,f,x,x,f',name,rms,other
    rmsarr[j,*] = rms
endfor


end

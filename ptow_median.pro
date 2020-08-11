PRO ptow_median, list

; This routine returns the median and mean of the c0 and c1
; coefficients for a list of ptow files. It returns the median and
; mean for each file, then a final median and mean of each of the values

readcol,list,f='a',silent=1,files
n0 = n_elements(files)

tr = fltarr(4,n0)

for j=0,n0-1 do begin
    file = files[j]
    readcol,file,f='d,d',c0,c1,silent=1
    izero = where(c0 eq 0.0,count)
    if (count gt 0) then begin
        igood = where(c0 ne 0.0)
        c0 = c0[igood]
        ci = c1[igood]
    endif
    tr[0,j] = median(c0)
    tr[1,j] = mean(c0)
    tr[2,j] = median(c1)
    tr[3,j] = median(c1)
    print,file
    print,'Median c0: '+strn(tr[0,j])
    print,'Mean c0: '+strn(tr[1,j])
    print,'Median c1: '+strn(tr[2,j])
    print,'Mean c1: '+strn(tr[3,j])
    print,''
endfor

print,'The median for all c0 values is: '+strn(median(tr[0,*],/even))
print,'The mean for all c0 values is: '+strn(mean(tr[1,*]))
print,'The median for all c1 values is: '+strn(median(tr[2,*],/even))
print,'The mean for all c1 values is: '+strn(median(tr[3,*]))

stop
END

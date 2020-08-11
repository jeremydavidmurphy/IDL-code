;This routine is used to calculate the mean and +/-S.D. from the
;name_region.fit files (the LOSVD profiles). This is used as an error
;estimate for the dynamical modeling.
pro vmsd,name,num

n1 = name+'.list'
n2 = name+'.out'
if (num eq 3) then begin
    readcol,n1,format='x,x,d,x,x,d,x,x,d',a,b,c
a = transpose(a)
b = tra
    openw,5,n2
    for j=0,28 do begin
        printf,5,a[j],b[j],c[j]
    endfor
    free_lun,5
endif
if (num eq 4) then begin
    readcol,n1,format='x,x,d,x,x,d,x,x,d,x,x,d',a,b,c,d
    openw,5,n2
    for j=0,28 do begin
        printf,5,a[j]
        printf,5,b[j]
        printf,5,c[j]
        printf,5,d[j]
    endfor
    free_lun,5
endif

end

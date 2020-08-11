PRO Pbinval

readcol,'losvd.list',format='a',files
n0 = n_elements(files)

for j=0,n0-1 do begin
    name = files[j]
    nameout = strsplit(name,'.',/extract)
    nameout = nameout[0]
    readcol,name,silent=1,vel,i1,i2,i3
    n1 = n_elements(vel)
    openw,5,nameout+'.fit'
    for k=0,n1-1 do printf,5,strn(k+1),vel[k],i1[k]
    free_lun,5
    openw,5,nameout+'.fitd'
    for k=0,n1-1 do printf,5,strn(k+1),vel[k],i2[k]
    free_lun,5
    openw,5,nameout+'.fitu'
    for k=0,n1-1 do printf,5,strn(k+1),vel[k],i3[k]
    free_lun,5
endfor

stop
END

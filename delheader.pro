PRO delheader, list

readcol,list,silent=1,f='a',files
n0 = n_elements(files)

for j=0,n0-1 do begin
    data = readfits(files[j],header,/silent)
    data = float(data)
    sxdelpar,header,'BZERO'
    sxaddpar,header,'BZERO',0.0,' a fudge to get Vaccine to work'
    sxdelpar,header,'BSCALE'
    sxaddpar,header,'BSCALE',1.0,' a fudge to get Vaccine to work'
    writefits,files[j],data,header
    print,'File '+files[j]+' complete!'
endfor

stop
END

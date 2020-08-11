pro alignall, list, ptow, mask, step

;This code essentially feeds the function REALIGN.PRO which does all
;the real work. It's purpose is to ALIGN INDIVIDUAL FRAMES, then get
;fed into combiwgt_v2.pro to get their final combination.

;The aligned data is returned by spliting the input name on either a
;'_' or a '.', then subscripted with a 'A.fits'

wi = 3530
wd = 1.12
readcol,list,format='A',files
readcol,ptow,format='A',ptows
readcol,mask,format='A',masks

n1 = n_elements(files)
if(n1 ne n_elements(ptows)) then begin
    print,'The two lists are not of the same length!'
    stop
endif

for j=0,n1-1 do begin
    file = files[j]
    data = readfits(file)
    test = n_elements(data)
    if (test eq 1) then begin
        print,'File '+file+' was not found! Skipping to next file.'
 ;       pause
        goto, jump1
    endif
    name = strsplit(file,'_.',/extract)
    print,''
    print,'File name: '+name[0]
    newname = name[0]+'_Asky.fits'
    out = realign(data,ptows[j],masks[j],wi,wd,step)
    writefits,newname,out
    jump1:
endfor

end

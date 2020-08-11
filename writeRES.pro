;under construction. this should take over from the cruder version of
;this, which is resmaker.pro
PRO writeRES, datalist

readcol,'bin.list',f='a',bins,silent=1
n0 = n_elements(bins)

readcol,datalist,f='x,a,a,a',pointings,arcs,dotR,silent=1
n1 = n_elements(pointings)

for j=0,n1-1 do begin
    arc = arcs[j]
    arc = strsplit(arc,'.',/extract)
    arcname = 'VIRUS-P_IR_FIT'+arc[0]+'.fits'
    if (j eq 0) then begin
        test = readfits(arcname,/silent)
        arcframes = fltarr(n_elements(test[*,0]),n_elements(test[0,*])-1,n1)
;        arcframes = fltarr[10,247,n1]
;        if (n_elements(test[0,*]) eq 247) 
        arcframes[*,*,0] = test[*,0:n_elements(test[0,*])-2]
        wave = test[*,n_elements(test[0,*])-1]
        readcol,dotR[j],silent=1,f='i',ifibR
        ibad = where(ifibR eq -1,ct)
        if (ct gt 0) then arcframes[*,ibad,j] = -1
    endif else begin
        test = readfits(arcname,/silent)
        if (n_elements(test[0,*]) eq 247) then test = [[test],[0,0,0,0,0,0,0,0,0,0]]
        arcframes[*,*,j] = test[*,0:n_elements(test[0,*])-2]
        readcol,dotR[j],silent=1,f='i',ifibR
        ibad = where(ifibR eq -1,ct)
        if (ct gt 0) then arcframes[*,ibad,j] = -1
    endelse
endfor ;the arc frames are now all contained in an array named "arcframes". The bad fibers have been IDed and masked.

for j=0,n0-1 do begin
    onebin = bins[j]
    readcol,onebin,silent=1,f='a',fibers
    igood = where(fibers ne -1,ct)
    if (ct gt 0) then fibers = fibers[igood]
    n2 = n_elements(fibers)
;    fibers = intarr(n2)
;    onepoint = strarr(n2)
    for k=0,n2-1 do begin
        onefiber = fibers[k]
        t = strsplit(onefiber,'_',/extract)
;        fibers[k] = uint(t[0])
;        onepoint[k] = t[1]
        fnum = t[0]
        ifiber =  uint(t[0]) - 1
        if (ifiber eq 246) then ifiber = 245
        onepoint = t[1]
        ione = where(pointings eq onepoint) ;index for which nights are needed
;        arctemp = arcframes[*,*,ione] ; just the relevant arc frames
        arcone = arcframes[*,ifiber,ione]
        if (k eq 0) then arcall = arcone else arcall = [[[arcall]],[[arcone]]]
    endfor ;the fibers and pointing for a given bin are in the array "fibers" and "onepoint"
    ibad = where(arcall eq -1,cntbad)
    if (cntbad gt 0) then arcall[ibad] = !VALUES.F_NAN
    ttt = median(arcall,dim=3,/even)
    cd, 'res_values'
    free_lun,55
    openw,55,onebin+'_median.res'
    for k=0,n_elements(wave)-1 do printf,55,wave[k],' ',ttt[k]
    free_lun,55
    for k=0,n_elements(arcall[0,0,*])-1 do begin
        openw,55,onebin+'_'+strn(k+1)+'.res'
        for l=0,n_elements(wave)-1 do printf,55,strn(wave[l]),'   ',strn(arcall[l,0,k]),'   ',strn(1)
        free_lun,55
    endfor
    cd,'../'
endfor

stop
END

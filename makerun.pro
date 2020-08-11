;this routine is used to generate the runall lists for the dynamical
;models. 

;MODIFIED ON Aug 13, 2010: To include a switch to NFW models.

PRO makerun, outputname, Limit=limit, Lists=lists, Override=oride,$
             Type=type

;the outputname is the name you want your output run files named. the
;limit is the size for each of your run files. They will subscript
;them each with letters.
;EX:
;outputname = 'run1'
;limit = 150
;each list will be written out as 'run1a','run1b',etc.

;If the limit callword isn't used, the entire list is written out in
;one piece.

; LISTS: If this callword is used and set to 'yes' then the code will
; look for a series of 4 lists in the calling directory that give the
; values for all the parameters that go into the dynamical models. If
; the keyword IS NOT USED then you will be prompted to enter values.

; If you enter the values by hand both the V_circ and R_scale will be
; integers unless you change your code...

if (n_elements(lists) eq 0) then lists = 'no'
if (n_elements(outputname)) eq 0 then begin
    outputname = ''
    print,'Enter a name for the run you are creating:'
    read,outputname
endif

if (n_elements(type) eq 0) then begin
    type = ''
    jump0:
    print,'NFW or Log DM halo?'
    print,'(NFW/log)'
    read,type
    print,type
    if (type ne 'NFW' and type ne 'log') then begin
        print,'Not a valid option! Try again!'
        goto, jump0
    endif
endif

if (lists eq 'yes') then begin
    readcol,'ml.list',format='a',ml
    readcol,'bhm.list',format='a',bhm
    readcol,'vcirc.list',format='a',vcirc
    readcol,'rscale.list',format='a',rscale
endif

restart:
if (type eq 'log') then begin
    bhu = fltarr(1)
    bhd = fltarr(1)
    bhi = fltarr(1)
    mlu = fltarr(1)
    mld = fltarr(1)
    mli = fltarr(1)
    vcu = fltarr(1)
    vcd = fltarr(1)
    vci = fltarr(1)
    rsu = fltarr(1)
    rsd = fltarr(1)
    rsi = fltarr(1)

    print,'Enter a lower black hole mass: (x10^9)'
    read,bhd
    print,'Enter an upper black hole mass:'
    read,bhu
    if (bhd eq bhu) then begin
        bhm = bhd
        bhi = 0.0
        goto,jump1
    endif
    print,'Enter a black hole increment:'
    read,bhi
    n0 = floor((bhu-bhd)/bhi)
    n0 = n0[0]
    bh = fltarr(n0)
    for j=0,n0-1 do begin
        bh[j] = bhd + bhi*j
    endfor
    bhm= [bh,bhu]
jump1:

    print,'Enter a lower M/L value:'
    read,mld
    print,'Enter an upper M/L value:'
    read,mlu
    if (mlu eq mld) then begin
        ml = mld
        mli = 0.0
        goto,jump2
    endif
    print,'Enter a M/L increment:'
    read,mli
    n0 = floor((mlu-mld)/mli)
    n0 = n0[0]
    ml = fltarr(n0)
    for j=0,n0-1 do begin
        ml[j] = mld + mli*j
    endfor
    ml = [ml,mlu]
jump2:

    print,'Enter a lower circular velocity value: (km/sec)'
    read,vcd
    print,'Enter an upper circular velocity value:'
    read,vcu
    if (vcd eq vcu) then begin
        vcirc = vcd
        vci = 0.0
        goto,jump3
    endif
    print,'Enter a circular velocity increment:'
    read,vci
    n0 = floor((vcu-vcd)/vci)
    n0 = n0[0]
    vc = fltarr(n0)
    for j=0,n0-1 do begin
        vc[j] = vcd + vci*j
    endfor
    vcirc = [vc,vcu]
jump3:

    print,'Enter a lower DM scale radius: (kpc)'
    read,rsd
    print,'Enter an upper DM scale radius:'
    read,rsu
    if (rsd eq rsu) then begin
        rscale = rsd
        rsi = 0.0
        goto,jump4
    endif
    print,'Enter a DM scale radius increment:'
    read,rsi
    n0 = floor((rsu-rsd)/rsi)
    n0 = n0[0]
    rs = fltarr(n0)
    for j=0,n0-1 do begin
        rs[j] = rsd + rsi*j
    endfor
    rscale = [rs,rsu]
jump4:

endif

if (type eq 'NFW') then begin
    bhu = fltarr(1)
    bhd = fltarr(1)
    bhi = fltarr(1)
    mlu = fltarr(1)
    mld = fltarr(1)
    mli = fltarr(1)
    vcu = fltarr(1)
    vcd = fltarr(1)
    vci = fltarr(1)
    rsu = fltarr(1)
    rsd = fltarr(1)
    rsi = fltarr(1)

    print,'Enter a lower black hole mass: (x10^9)'
    read,bhd
    print,'Enter an upper black hole mass:'
    read,bhu
    if (bhd eq bhu) then begin
        bhm = bhd
        bhi = 0.0
        goto,jump5
    endif
    print,'Enter a black hole increment:'
    read,bhi
    n0 = floor((bhu-bhd)/bhi)
    n0 = n0[0]
    bh = fltarr(n0)
    for j=0,n0-1 do begin
        bh[j] = bhd + bhi*j
    endfor
    bhm= [bh,bhu]
jump5:

    print,'Enter a lower M/L value:'
    read,mld
    print,'Enter an upper M/L value:'
    read,mlu
    if (mlu eq mld) then begin
        ml = mld
        mli = 0.0
        goto,jump6
    endif
    print,'Enter a M/L increment:'
    read,mli
    n0 = floor((mlu-mld)/mli)
    n0 = n0[0]
    ml = fltarr(n0)
    for j=0,n0-1 do begin
        ml[j] = mld + mli*j
    endfor
    ml = [ml,mlu]
jump6:

    print,'Enter a lower concentration:'
    read,vcd
    print,'Enter an upper concentration:'
    read,vcu
    if (vcd eq vcu) then begin
        vcirc = vcd
        vci = 0.0
        goto,jump7
    endif
    print,'Enter a concentration increment:'
    read,vci
    n0 = floor((vcu-vcd)/vci)
    n0 = n0[0]
    vc = fltarr(n0)
    for j=0,n0-1 do begin
        vc[j] = vcd + vci*j
    endfor
    vcirc = [vc,vcu]
jump7:

    print,'Enter a lower NFW DM scale radius: (kpc)'
    read,rsd
    print,'Enter an upper NFW DM scale radius:'
    read,rsu
    if (rsd eq rsu) then begin
        rscale = rsd
        rsi = 0.0
        goto,jump8
    endif
    print,'Enter a NFW DM scale radius increment:'
    read,rsi
    n0 = floor((rsu-rsd)/rsi)
    n0 = n0[0]
    rs = fltarr(n0)
    for j=0,n0-1 do begin
        rs[j] = rsd + rsi*j
    endfor
    rscale = [rs,rsu]
jump8:


endif

form1 = '(a13,1x,f10.5,2x,f10.5,2x,f9.5)'
free_lun,5
openw,5,outputname+'_params.txt'
printf,5,'black hole: ',bhd,bhu,bhi,format=form1
printf,5,'M/L: ',mld,mlu,mli,format=form1
if (type eq 'log') then printf,5,'V_circ: ',vcd,vcu,vci,format=form1
if (type eq 'NFW') then printf,5,'C: ',vcd,vcu,vci,format=form1
printf,5,'R_scale: ',rsd,rsu,rsi,format=form1
free_lun,5

n1 = n_elements(ml)
n2 = n_elements(bhm)
n3 = n_elements(vcirc)
n4 = n_elements(rscale)

n5 = n1*n2*n3*n4
backup:
print,'You have '+strn(n5)+' elements in your modeling run.'
print,'(C)ontinue or (R)edo?'
ans = ''
read,ans
if (ans eq 'R') then goto,restart
if (ans eq 'C') then print,'Very well...' else goto,backup

outfile=strarr(5,n5)
cntr = 0
for j=0,n1-1 do begin
    for k=0,n2-1 do begin
        for l=0,n3-1 do begin
            for m=0,n4-1 do begin
                bhpiece = strn(bhm[k],length=5)+'E+9'
                if (type eq 'log') then line = ['runallLOG',strn(ml[j],length=5),bhpiece,$
                        strn(vcirc[l],length=5),strn(rscale[m],length=5)]
                if (type eq 'NFW') then line = ['runallNFW',strn(ml[j],length=5),bhpiece,$
                        strn(vcirc[l],length=5),strn(rscale[m],length=5)]
                outfile[*,cntr] = line
                cntr = cntr + 1
            endfor
        endfor
    endfor
endfor

n6 = n_elements(outfile[0,*])
boost = 0

if (n_elements(oride) ne 0) then begin
    readcol,oride,silent=1,format='f,f,i,i,x,x,x',mlO,bhO,vcO,rsO
;    imlO = intarr
    for j=0,n1-1 do begin
        temp = where(mlO eq ml[j])
    endfor       
endif

if (n_elements(limit) eq 0) then begin
    openw,5,outputname
    for j=0,n6-1 do printf,5,outfile[*,j]
    free_lun,5
endif else begin
    addon = ['a','b','c','d','e','f','g','h','i','j','k','l','m','n','o',$
            'p','q','r','s','t','u','v','x','y','z','aa','bb','cc','dd','ee',$
            'ff','gg','hh','ii','jj','kk','ll','mm','nn','oo','pp','qq','rr',$
            'ss','tt','uu','vv','xx','yy','zz']
    n7 = floor(n6/limit)
    for j=0,n7 do begin
        name = outputname+addon[j]
        openw,5,name
        if (boost+limit lt n6) then begin
            for k=0,limit-1 do printf,5,outfile[*,k+boost]
            boost = boost + k
            free_lun,5
        endif else begin
            diff = n6 - boost
            for k=0,diff-1 do printf,5,outfile[*,k+boost]
            free_lun,5
        endelse
    endfor
endelse

STOP
END

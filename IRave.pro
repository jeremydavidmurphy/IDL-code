PRO IRave, irvalues, binlist

readcol,irvalues,silent=1,skipline=1,fib,w1,w2,w3,w4,w5,w6,w7,w8,w9,w10

IRarr = [[w1],[w2],[w3],[w4],[w5],[w6],[w7],[w8],[w9],[w10]]
IRarr = transpose(IRarr)

readcol,binlist,silent=1,f='a',bins
n0= n_elements(bins)

readcol,'wave.list',silent=1,wave

for j=0,n0-1 do begin
    bin = bins[j]
    print,'Working on bin '+bin
    readcol,bin,silent=1,f='a',fibers
    fibers = fibers[1:*]
    n1 = n_elements(fibers)

    irtemp = fltarr(10)
    for k=0,n1-1 do begin
        fiber = fibers[k]
        s = strsplit(fiber,'_',/extract)
        fn = s[0]
        fn = float(fn)
        ifn = uint(fn - 1.0)
        irtemp = [[irtemp],[IRarr[*,ifn]]]
    endfor
    irtemp = irtemp[*,1:*]
    if (n1 eq 1) then begin
        irout0 = irtemp
        irout1 = irtemp
        irout2 = irtemp
    endif else begin
        irout1 = median(irtemp,dim=2,/even)
        irout0 = fltarr(10) & irout2 = irout0
        for k=0,9 do irout0[k] = min(irtemp[k,*])
        for k=0,9 do irout2[k] = max(irtemp[k,*])
    endelse
    free_lun,5
    openw,5,bin+'.res'
    for k=0,9 do printf,5,wave[k],'  ',irout0[k],'  ',irout1[k],'  ',irout2[k]
    free_lun,5
endfor

stop
END        

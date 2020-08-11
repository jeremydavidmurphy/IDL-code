; This routine is just a wrapper. It feeds realigncF.
PRO shift, list, wi, wd

plotting = 'on'

readcol,list,f='a,a',files,ptows
n0 = n_elements(files)

for j=0,n0-1 do begin ; a loop through each file (the length of the datalist.
    print,''
    print,'Reading in file '+files[j]
    frame = readfits(files[j],/silent,h)

    filename = strsplit(files[j],'.',/extract)
    filename = filename[0]+'A.fits'
    
    cntr = 0
    out = 'lost'
    repeat begin
        hh = strsplit(h[cntr],' ',/extract)
        if (hh[0] eq 'WAVEZP') then out = 'found' else cntr = cntr + 1
    endrep until (out eq 'found')
    wavezp = float(hh[2])
    if (wavezp gt 1.0) or (wavezp lt -1.0) then begin
        print,'You have a crazy big wavelength ZP!!!'
        print,wavezp[j]
        stop
    endif
    print,'The wavelength offset for frame '+files[j]+' is '+hh[2]
    cntr = 0
    out = 'lost'
    repeat begin
        hh = strsplit(h[cntr],' ',/extract)
        if (hh[0] eq 'VHELIO') then out = 'found' else cntr = cntr + 1
    endrep until (out eq 'found')
    helioc = float(hh[2])
    print, 'The heliocentric correction is '+hh[2]

    ptowframe = ptows[j]

    out = realigncF(frame, ptowframe, wi, wd, WZP=wavezp, HELIOC=helioc)
    H1 = 'W-INIT  = '+strn(wi)
    H2 = 'W-DISP  = '+strn(wd)
    n1 = n_elements(h)
    hh = [h[0:n1-2],'HISTORY = Aligned frames',H1,H2,'END     ']
    writefits,filename,out,hh

endfor

stop
end

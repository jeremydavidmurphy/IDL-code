; This is the function version of Rmajor. It returns the Rmajor values
; from either fiber(bin) RA and Dec positions (e.g. from the
; finder_code coordinate files) or dRA and dDec offsets.

FUNCTION FRmajor, RA, Dec, PA, axisratio, galRA=galRA, galDec=galDec

; RA and DEC: The RA and Dec (either offsets or coordinate positions)

; PA: position angle of the galaxy, Eastward of North
; AXISRATIO: R_major/R_minor

; GALRA and GALDEC: the galaxy RA and Dec ('09:19:43.6')

; The galoffsets array is what the fiber dRA and dDec is IN THE
; REFERENCE FRAME OF THE GALAXY. 

theta = 90.0 - pa
theta = theta / (360.0/(2*!pi))

fRA = RA
fDec = Dec
n0 = n_elements(fRA)

fcoords = dblarr(2,n0) ;fiberRA, fiberDec (in decimal coords)
foffsets = dblarr(3,n0) ;deltaRA, deltaDec, radius
galoffsets = dblarr(3,n0) ;deltaRA, deltaDec, Rmajor (rotated coord. system)

test = fRA[0]
t = strsplit(test,':')
if (n_elements(t) eq 1) then test = t[0]

if (test eq 0) then begin ;the coords are assumed to be offsets!
    for j=0,n0-1 do begin
        foffsets[0,j] = double(fRA[j])
        foffsets[1,j] = double(fDec[j])
        foffsets[2,j] = sqrt(foffsets[0,j]^2 + foffsets[1,j]^2)
    endfor
    goto, jump1
endif else begin ;the coordinates are positions, NOT offsets!
    if (n_elements(galRA) eq 0) then begin
        galRA = ''
        print,'Enter the galaxy RA (00:00:00.0):'
        read,galRA
    endif
    if (n_elements(galDec) eq 0) then begin
        galDec = ''
        print,'Enter the galaxy Dec (00:00:00.0):'
        read,galDec
    endif
    RA = strn(galRA) & Dec = strn(galDec)
    RA = strsplit(RA,':',/EXTRACT)
    Dec = strsplit(Dec,':',/EXTRACT)
    galRA = double((RA[0]+RA[1]/60.+RA[2]/3600.)*15)
    if (Dec[0] lt 0) then galDec = double(Dec[0]-Dec[1]/60.-Dec[2]/3600.)
    if (Dec[0] ge 0) then galDec = double(Dec[0]+Dec[1]/60.+Dec[2]/3600.)

    for j=0,n0-1 do begin
        fc = fra[j]
        fd = fdec[j]
        fc = strsplit(fc,':',/EXTRACT)
        fd = strsplit(fd,':',/EXTRACT)
        fcoords[0,j] = double((fc[0]+fc[1]/60.+fc[2]/3600.)*15.)
        if (offsetcorrect eq 'on') then fcoords[0,j] = fcoords[0,j] - raoffset
        if (fd[0] lt 0) then fcoords[1,j] = double(fd[0]-fd[1]/60.-fd[2]/3600.)
        if (fd[0] ge 0) then fcoords[1,j] = double(fd[0]+fd[1]/60.+fd[2]/3600.)
        if (offsetcorrect eq 'on') then fcoords[1,j] = fcoords[1,j] - decoffset
    endfor
endelse

for j=0,n0-1 do begin
    dx = double((fcoords[0,j]-galRA))
    dy = double((fcoords[1,j]-galDec))
    foffsets[0,j] = dx*3600.0
    foffsets[1,j] = dy*3600.0
    foffsets[2,j] = sqrt(dx^2 + dy^2)*3600.0
endfor

jump1:

for j=0,n0-1 do begin
    galoffsets[0,j] = foffsets[0,j]*cos(theta) - foffsets[1,j]*sin(theta)
    galoffsets[1,j] = foffsets[0,j]*sin(theta) + foffsets[1,j]*cos(theta)
    galoffsets[2,j] = sqrt(galoffsets[0,j]^2 + (galoffsets[1,j]^2/axisratio^2))
endfor

return,galoffsets[2,*]
END

; This routine reads in the output of pipe1.pro (the galaxy_D#.coord)
; and generates bins based on the radial positions of each fiber.

; MODIFIED ON DEC 8, 2010: The code now accepts a list of binning
; regions, in arcsecs, and uses this to determine the bins.

pro radial_bin, coordlist

;readcol,coordlist,f='a',files
;n0 = n_elements(files)

readcol,coordlist,f='a,f',fiber,radius
n1 = n_elements(fiber)

readcol,'radial.bnn',f='f',Rbins
;readcol,'radial.bin',f='f',Rbins

n2 = n_elements(Rbins)

free_lun,5
for j=0,n2-2 do begin
    m1 = Rbins[j]
    m2 = Rbins[j+1]
    i = where(radius gt m1 and radius lt m2,count)
    if (count gt 0) then begin
        name = 'bnn'+strn(uint(m1))+'_'+strn(uint(m2))
;        name = 'bin'+strn(uint(m1))+'_'+strn(uint(m2))
        print,fiber[i],' ',radius[i]
        openw,5,name
        for k=0,count-1 do printf,5,fiber[i[k]]
        free_lun,5
    endif
    print,count
endfor

stop 
end

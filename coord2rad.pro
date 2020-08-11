; This routine reads in a list of name.coord files (output from
; pipe2.pro) and returns a luminosity-weighted radius, meaning, the
; radial positions of each fiber are weighted by the flux in each
; fiber.

; This code assumes the bin.coord files are in the new format:
; EX:  F #             dRA          dDec          Rad          NFlux
;      103       -78.55200     103.69700    130.09100       28.63715
;       64       -45.59300     122.19900    130.42700       31.26864
; etc.

PRO coord2rad

readcol,'coord.list',f='a',list,silent=1

n1 = n_elements(list)

radius = fltarr(n1)
names = strarr(n1)

for j=0,n1-1 do begin
    file = list[j]
    readcol,file,silent=1,format='x,x,x,f,f',skipline=1,rad,flux
    i = where(flux ne -1.0)
    rad = rad[i]
    flux = flux[i]
    weights = flux / max(flux)
    top = total(weights * rad)
    bot = total(weights)
    ave = top / bot
    radius[j] = ave
    t = strsplit(file,'.',/extract)
    names[j] = t[0]
endfor

free_lun,5
openw,5,'bin.radius'
for j=0,n1-1 do printf,5,names[j],radius[j],format='(A16,2x,F10.5)'
free_lun,5

stop
END

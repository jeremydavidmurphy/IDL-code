;THIS CODE HAS, TO SOME DEGREE, BEEN MADE OBSOLETE BY
;BINOFFSET.PRO. Binoffset.pro creates the binning regions, but does so
;with weights based on the intensity of the data in the various fibers
;being combined.

PRO radbin

;This code accepts the bin#### regions (which fibers are included) and
;the name.std.rad file and returns the weighted radial position of the
;bin in arcseconds from the center of the galaxy.
name=''
print,'Enter the .std.rad file name (w/o the .std.rad):'
read,name
name = name+'.std.rad'

readcol,name,FORMAT='I,F',fiber,radius
ind = bsort(fiber)
fiberstd = fiber(ind)
radiusstd = radius(ind)

print,'Enter the list of bin#### files:'
read,name

readcol,name,FORMAT='A',bnames
n = n_elements(bnames)

for k=0,n-1 do begin
    nm = bnames[k]
    readcol,nm,FORMAT='I',fibnum
    print,'Fibers in bin '+nm+' is : '+strn(fibnum)
    bin = radiusstd(fibnum-1)
    print,bin
endfor

STOP
END

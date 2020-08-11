PRO split

;This code reads in a fits file and breaks it up into individual
;spectra

name =''
ans1=''
ix=''

print,'Enter the name of the file (w/o the fits):'
read,name
name = name+'.fits'

data = readfits(name)

n1 = n_elements(data[*,0])
n2 = n_elements(data[0,*])

for k=0,n2-1 do begin
    spec = data[*,k]
    writefits,'F'+strn(k+1)+'.fits',spec
endfor

end

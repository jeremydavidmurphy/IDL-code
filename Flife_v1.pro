PRO Flife_v1

readcol,'trim.list',f='a',list
n0 = n_elements(list)

test = readfits(list[0],h)
n1 = n_elements(test[*,0])
n2 = n_elements(test[0,*])

array = fltarr(n1,n2,25)
for j=0,24 do array[*,*,j] = readfits(list[j],/silent)

window,0,retain=2,xsize=n1/2.0,ysize=n2/2.0
loadct,27

for j=12,n0-13 do begin
    print,list[j]
    frame = readfits(list[j],/silent)
    med = median(array,dim=3)
    div = frame / med
    print,stddev(div)
    tvscale,div
    add = readfits(list[j+13],/silent)
    array = [[[array[*,*,1:24]]],[[add]]]
endfor

stop
END

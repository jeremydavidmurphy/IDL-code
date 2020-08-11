PRO stability

readcol,'all.list',f='a',list
n0 = n_elements(list)

totals = fltarr(n0)

for j=0,n0-1 do begin
    frame = readfits(list[j],/silent)
    totals[j] = total(frame)
    print,j,totals[j]
endfor

window,0,retain=2
plot,totals,/ynozero,title='Time (sec)',ytitle='Total Counts'

window,1,retain=2
plot,totals/max(totals),/ynozero,title='Time (sec)',ytitle='Total Counts'
stop
END

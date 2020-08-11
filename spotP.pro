PRO spotP, list

winsize = 8.0

readcol,list,format='a',slist
n0 = n_elements(slist)
test = readfits(slist[0],/silent)
n1 = n_elements(test[*,0])
n2 = n_elements(test[0,*])

window,0,retain=2,xsize=n1/winsize,ysize=n2/winsize
device,decomposed=0
for j=0,n0-1 do begin
    name = strsplit(slist[j],'.',/extract)
    name = name[0]
    frame = readfits(slist[j],/silent)
    loadct,37
    TVImage, BytScl(frame, Top=!D.Table_Size-3)
    loadct,0
    xyouts,0.1,0.1,name,charsize=1.5,/normal
    pause
endfor

stop
END

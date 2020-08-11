pro Pmcfit2loop

; Looks for the list mcfit2.list of the format
;bin1001_HK.mcfit2
;bin1001_GB.mcfit2
;bin1001_HB.mcfit2
;bin1001_MGw.mcfit2
;bin1001_MGwo.mcfit2
;bin1001_FE.mcfit2
;*
;bin1105_HK.mcfit2
;bin1105_GB.mcfit2
;bin1105_HB.mcfit2
;bin1105_MGw.mcfit2
;bin1105_MGwo.mcfit2
;bin1105_FE.mcfit2
;*
;etc...

; This routine lets you plot simultaneously all the losvds (the
; mcfit2's) for a given bin so as to compare spectral regions.

readcol,'mcfit2.list',silent=1,format='a',list
n0 = n_elements(list)
count = 0.0
for j=0,n0-1 do begin
    temp = list[j]
    if (temp ne '*') then count = count + 1.0 else goto, jump1
endfor
jump1:
n1 = uint(count)

; the list is now read in again, this time stepping over the
; delimiters...
readcol,'mcfit2.list',silent=1,format='a',list,comment='*'
nall = n_elements(list)

nfiles = nall/n1
print,'The number of bins to be plotted is '+strn(nfiles)

window,2,retain=2,xsize=1200,ysize=900
device,decomposed=0
loadct,0
    
n2 = 0
ans = ''
repeat begin
    piece = list[n2:n2+n1-1]
    print,'The files to be plotted...'
    print,transpose(piece)
    
    for j=0,n1-1 do begin
        !p.multi = [j,3,2,0,1]
        readcol,silent=1,piece[j],format='f,f,f,f',v,i,d,u
        temp = strsplit(piece[j],'_.',/extract)
        plot,v,i,xtitle='Velocity (km/sec)',title=temp[0]+'  '+temp[1],charsize=1.5
        oplot,v,u,linestyle=2
        oplot,v,d,linestyle=2
    endfor
    n2 = n2 + n1
    if (n2 ge nfiles*n1) then begin
        ans = 'done'
        goto,jump2
    endif
    print,'Enter "q" to quit or RETURN to continue...'
    read,ans
    jump2:
endrep until (ans eq 'q' or ans eq 'done')
pause
wdelete
stop
end

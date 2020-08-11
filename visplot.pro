;---------------------------------------------------------------
;
; PURPOSE:
; This routine is used to plot the velocity and velocity dispersion
; profiles for VIRUS-P data. It can be used when multiple pointings of
; the galaxy exist.
;
; REQUIREMENTS:
; "bin.list" must exist in the calling directory. "bin.list" is a list
; of all the bins, each with the format ###_a, ###_b, etc. with each
; letter corresponding to a specific pointing on the galaxy.
;
; "dotR.list" must exist in the calling directory. This is a list of
; the master .R files for each pointing ('master' meaning that no
; fibers are rejected in the .R file). THESE MUST BE IN ORDER,
; CORRESPONDING TO THE LETTERS OF THE BINS. So, the 'A' pointing must
; go first, then the 'B' pointing, and so on.
;
; "pfitlov.out" must exist in the calling directory. This should
; contain the final values to be plotted. 

;---------------------------------------------------------------
PRO visplot
;---------------------------------------------------------------

color=[60,110,180,150]
readcol,'bin.list',format='a',binlist
readcol,'dotR.list',format='a',dotR
;readcol,'pfitlov.out',format='a,f,f,x,f,f,x,x,x',wholename,vel,sigma,h3,h4

n1 = n_elements(binlist)
binarr = strarr(n1,250,5)

for j=0,n1-1 do begin
    readcol,binlist[j],format='a',bin
    n2 = n_elements(bin)
    for k=0,n2-1 do begin
        binarr[j,k,0:1] = strsplit(bin[k],'_',/extract)
    endfor
endfor

n3 = n_elements(dotR)
dotRarr = fltarr(5,250,n3)
binloc = strarr(n1,4)
for j=0,n3-1 do begin
    readcol,dotR[j],f,r,dra,ddec,i
    n4 = n_elements(f)
    extra = fltarr(250-n4)
    print,n_elements(extra)
    dotRarr[0,*,j] = [f,extra]
    dotRarr[1,*,j] = [r,extra]
    dotRarr[2,*,j] = [dra,extra]
    dotRarr[3,*,j] = [ddec,extra]
    dotRarr[4,*,j] = [i,extra]
endfor

for j=0,n1-1 do begin
    ia = where(binarr[j,*,1] eq 'a')
    ib = where(binarr[j,*,1] eq 'b')
    ic = where(binarr[j,*,1] eq 'c')
    id = where(binarr[j,*,1] eq 'd')
    ie = where(binarr[j,*,1] eq 'e')
    if (ia[0] ne -1) then begin
        fa = uint(binarr[j,ia,0]); the fibers are found
        dRA_ia = intarr(n_elements(fa))
        for k=0,n_elements(fa)-1 do dRA_ia[k] = where(uint(dotRarr[0,*,0]) eq fa[k])
        binarr[j,ia,2] = dotRarr[2,dRA_ia,0]
        binarr[j,ia,3] = dotRarr[3,dRA_ia,0]
        binarr[j,ia,4] = dotRarr[1,dRA_ia,0]
    endif
    if (ib[0] ne -1) then begin
        fb = uint(binarr[j,ib,0])
        dRA_ib = intarr(n_elements(fb))
        for k=0,n_elements(fb)-1 do dRA_ib[k] = where(uint(dotRarr[0,*,1]) eq fb[k])
        binarr[j,ib,2] = dotRarr[2,dRA_ib,1]
        binarr[j,ib,3] = dotRarr[3,dRA_ib,1]
        binarr[j,ib,4] = dotRarr[1,dRA_ib,1]
    endif
    if (ic[0] ne -1) then begin
        fc = uint(binarr[j,ic,0])
        dRA_ic = intarr(n_elements(fc))
        for k=0,n_elements(fc)-1 do dRA_ic[k] = where(uint(dotRarr[0,*,2]) eq fc[k])
        binarr[j,ic,2] = dotRarr[2,dRA_ic,2]
        binarr[j,ic,3] = dotRarr[3,dRA_ic,2]
        binarr[j,ic,4] = dotRarr[1,dRA_ic,2]
    endif
    if (id[0] ne -1) then begin
        fd = uint(binarr[j,id,0])
        dRA_id = intarr(n_elements(fd))
        for k=0,n_elements(fd)-1 do dRA_id[k] = where(uint(dotRarr[0,*,3]) eq fd[k])
        binarr[j,id,2] = dotRarr[2,dRA_id,3]
        binarr[j,id,3] = dotRarr[3,dRA_id,3]
        binarr[j,id,4] = dotRarr[1,dRA_id,3]
    endif
    if (ie[0] ne -1) then begin
        fe = uint(binarr[j,ie,0])
        dRA_ie = intarr(n_elements(fe))
        for k=0,n_elements(fe)-1 do dRA_ie[k] = where(uint(dotRarr[0,*,4]) eq fe[k])
        binarr[j,ie,2] = dotRarr[2,dRA_ie,4]
        binarr[j,ie,3] = dotRarr[3,dRA_ie,4]
        binarr[j,ie,4] = dotRarr[1,dRA_ie,4]
    endif
    ind = where(binarr[j,*,0] ne 0)
    binloc[j,0] = binlist[j]
    binloc[j,1] = mean(float(binarr[j,ind,2]))
    binloc[j,2] = mean(float(binarr[j,ind,3]))
    binloc[j,3] = mean(float(binarr[j,ind,4]))
endfor    

yup = max(dotRarr[3,*,*])+10
ydown = min(dotRarr[3,*,*])-10
xup = max(dotRarr[2,*,*])+10
xdown = min(dotRarr[2,*,*])-10
window,1,retain=2,xsize=800,ysize=800,ypos=50
device,decomposed=0
loadct,0
plot,dotRarr[2,where(dotRarr[2,*,0] ne 0),0],dotRarr[3,where(dotRarr[3,*,0] ne 0),0],$
  psym=3,yrange=[ydown,yup],xrange=[xdown,xup],/nodata,xstyle=1,ystyle=1

loadct,4
for j=0,n3-1 do begin
    oplot,dotRarr[2,where(dotRarr[2,*,0] ne 0),j],dotRarr[3,where(dotRarr[3,*,0] ne 0),j],psym=3,color=color[j]
endfor
loadct,27
for j=0,n1-1 do begin
    plots,binloc[j,1],binloc[j,2],psym=1
    xyouts,binloc[j,1]+1,binloc[j,2]+5,binloc[j,0],orientation=30,charsize=1.0,color=j*5
    xyouts,binloc[j,1]+1,binloc[j,2]+5,binloc[j,3],orientation=60,charsize=1.0,color=j*5
endfor

f1 = '(A10,2X,A10,2X,A10,2X,A10)'
f2 = '(A10,2X,F10.4,2X,F10.4,2X,F10.4)'
title = ['bin','radius','deltaRA','deltaDec']
openw,5,'binloc.txt'
printf,5,title,format=f1
for j=0,n1-1 do begin
    printf,5,binloc[j,0],binloc[j,3],binloc[j,1],binloc[j,2],format=f2
endfor
free_lun,5

stop
END

;----------------------------------------------------------------------
;The code was modifed on 11-11-08 to accept a list of recipe files and
;is now flexible enough to accept any list of ranges.
;It expects to find a list called recipe.list
;----------------------------------------------------------------------

PRO recipe

;This code plots the stellar recipe used to make the velocity fits
;from pfitlov.  It looks for 'Tlist' and all the '.recipe' files in the
;'results' directory.  This routine is hard-coded to expect 4 regions-
;H+K,Gband,Mg and All


ans=''
ans1=''

readcol, 'recipe.list',Format='a',rlist
readcol,'Tlist',Format='a,a',junk,tlist
n1 = n_elements(rlist)
nn = n_elements(tlist)

names = strarr(4,n1)

starmix = fltarr(n,nn)
xaxis = findgen(nn)+1

FOR k=0,n-1 DO BEGIN
    readcol,files[k],F='F,A',percent,junk
    starmix[k,*] = percent
    names[*,k] = strsplit(files[k],'._',/EXTRACT)
ENDFOR

radius = strarr(n_elements(names[1,*]))
FOR j=0,n_elements(radius)-1 DO radius[j] = strsplit(names[1,j],'f',/EXTRACT)

window,2,retain=2
device,decomposed=0
loadct,0
plot,starmix[0,*],xrange=[0,nn+1],xstyle=1,yrange=[-0.2,1.0],ystyle=1,title='!3Stellar Recipe for '+names[0,0],$
  ytitle='Relative Contribution',xtitle='Stellar Type',/nodata,charsize=1.5
FOR j=0,nn-1 DO xyouts,xaxis[j]+0.4,-0.18,tlist[j],orientation=90,charsize=1.5
loadct,4

ctr = 0
FOR k=0,n/4-1 DO BEGIN
    xyouts,xaxis,starmix[k+ctr,*],radius[k+ctr],color=255,charsize=1.5
    ctr = ctr + 3
ENDFOR
ctr = 0
FOR k=0,n/4-1 DO BEGIN
    xyouts,xaxis,starmix[k+ctr+1,*],radius[k+ctr+1],color=110,charsize=1.5
    ctr = ctr + 3
ENDFOR
ctr = 0
FOR k=0,n/4-1 DO BEGIN
    xyouts,xaxis,starmix[k+ctr+2,*],radius[k+ctr+2],color=60,charsize=1.5
    ctr = ctr + 3
ENDFOR
ctr = 0
FOR k=0,n/4-1 DO BEGIN
    xyouts,xaxis,starmix[k+ctr+3,*],radius[k+ctr+3],color=150,charsize=1.5
    ctr = ctr + 3
ENDFOR

xyouts,0.3,0.9,names[2,0]+' region',/normal,color=255,charsize=1.5
xyouts,0.3,0.85,names[2,1]+' region',/normal,color=110,charsize=1.5
xyouts,0.6,0.9,names[2,2]+' region',/normal,color=60,charsize=1.5
xyouts,0.6,0.85,names[2,3]+' region',/normal,color=150,charsize=1.5

print,'Save the plot?'
read,ans

IF (ans EQ 'y') THEN BEGIN
    print,'Add a note to the plot?'
    read,ans1
    IF (ans1 EQ 'y') THEN BEGIN
        note = strarr(1)
        print,'Type away...'
        read,note
        note = 'NOTE: '+note
    ENDIF

    set_plot,'ps'
    device,file='recipe.ps',/color
loadct,0
plot,starmix[0,*],xrange=[0,nn+1],xstyle=1,yrange=[-0.2,1.0],ystyle=1,title='!6Stellar Recipe for '+names[0,0],$
  ytitle='Relative Contribution',xtitle='Stellar Type',/nodata,charsize=1.3
FOR j=0,nn-1 DO xyouts,xaxis[j]+0.4,-0.18,tlist[j],orientation=90,charsize=1.2
loadct,4

ctr = 0
FOR k=0,n/4-1 DO BEGIN
    xyouts,xaxis,starmix[k+ctr,*],radius[k+ctr],color=180,charsize=1
    ctr = ctr + 3
ENDFOR
ctr = 0
FOR k=0,n/4-1 DO BEGIN
    xyouts,xaxis,starmix[k+ctr+1,*],radius[k+ctr+1],color=110,charsize=1
    ctr = ctr + 3
ENDFOR
ctr = 0
FOR k=0,n/4-1 DO BEGIN
    xyouts,xaxis,starmix[k+ctr+2,*],radius[k+ctr+2],color=60,charsize=1
    ctr = ctr + 3
ENDFOR
ctr = 0
FOR k=0,n/4-1 DO BEGIN
    xyouts,xaxis,starmix[k+ctr+3,*],radius[k+ctr+3],color=150,charsize=1
    ctr = ctr + 3
ENDFOR

xyouts,0.3,0.88,names[2,0]+' region',/normal,color=180,charsize=1.2
xyouts,0.3,0.83,names[2,1]+' region',/normal,color=60,charsize=1.2
xyouts,0.6,0.88,names[2,2]+' region',/normal,color=110,charsize=1.2
xyouts,0.6,0.83,names[2,3]+' region',/normal,color=150,charsize=1.2

IF (ans1 EQ 'y') THEN BEGIN
    loadct,0
    xyouts,0.19,0.18,note,charsize=0.7,ORIENTATION=90,/normal
ENDIF

device,/close_file
set_plot,'x'
ENDIF ELSE print,'Word. No plot for you!'

print,'Next ENTER deletes the plot.'
pause
wdelete,2
cd,'../'
STOP
END

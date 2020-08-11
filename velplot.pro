PRO velplot, arp

IF (n_elements(arp) EQ 0) THEN BEGIN
    arp=''
    print,'Enter an ARP file:'
    read,arp
ENDIF
;THIS ROUTINE IS CALLED FROM REDUC.  IT EXPECTS A LIST, CALLED
;'PFITLOV.OUT', TO EXIST IN THE 'RESULTS' DIRECTORY.

;The 'arp' file comes from the makebinall.pro routine.  It's the
;average radial position of each bin

ptnnum = 4 ;The number of regions you have in your pfitlov.out file (H&K, All, etc)
ans=''
ans1=''
band=''
color = [60,110,150,255]

cd,'results'

readcol,'pfitlov.out',Format='A,F,F,F,F,F',temp,vel,disp,junk,h3,h4

n = n_elements(temp)
names = strarr(4,n)
junk = strarr(2,n)
FOR j=0,n-1 DO names[*,j] = strsplit(temp[j],'._',/EXTRACT)
FOR j=0,n-1 DO junk[*,j] = strsplit(names[1,j],'f',/EXTRACT) ;This just dumps the 'f' from the name
names = [names[0,*],junk[0,*],names[2,*]]

cd,'../lists'
radfile = strn(arp)+'.arp'
radius = fltarr(n/ptnnum) 
openr,5,radfile
readf,5,radius
free_lun,5
cd,'../results'
xmin = min(radius)-15
xmax = max(radius)+15

ptn = names[0,0]
indy = where(names[0,*] EQ ptn)
n = n_elements(indy)

plotnames = strarr(3,n)
plotnames = names[0:2,indy]
data = fltarr(4,n)
data[0,*] = vel[indy]
data[1,*] = disp[indy]
data[2,*] = h3[indy]
data[3,*] = h4[indy]
ymin = min(data[1,*])-20
ymax = max(data[1,*])+20

window,2,retain=2
device,decomposed=0

REPEAT BEGIN
loadct,0
plot,radius[0:1],data[1,0:1],title='LOSVD measurements for '+plotnames[0,0],ytitle='!7r!6(!ukm!n/!dsec!n)',/nodata,/ynozero,$
  charsize=2,xrange=[xmin,xmax],xstyle=1,yrange=[ymin,ymax],ystyle=1,xtitle='Radius (arcsec)'
loadct,4
cntr = 0
    
FOR k=0,n/ptnnum-1 DO BEGIN
    plots,radius[k],data[1,k+cntr],psym=2,symsize=2,thick=2,color=color[3]
    xyouts,radius[k]+2,data[1,k+cntr]+2,plotnames[2,k+cntr],charsize=1.5,color=color[3],orientation=45
    plots,radius[k],data[1,k+cntr+1],psym=2,symsize=2,thick=2,color=color[0] 
    xyouts,radius[k]+2,data[1,k+cntr+1]+2,plotnames[2,k+cntr+1],charsize=1.5,color=color[0],orientation=45
    plots,radius[k],data[1,k+cntr+2],psym=2,symsize=2,thick=2,color=color[1]
    xyouts,radius[k]+2,data[1,k+cntr+2]+2,plotnames[2,k+cntr+2],charsize=1.5,color=color[1],orientation=45
    plots,radius[k],data[1,k+cntr+3],psym=2,symsize=2,thick=2,color=color[2]
    xyouts,radius[k]+2,data[1,k+cntr+3]+2,plotnames[2,k+cntr+3],charsize=1.5,color=color[2],orientation=45
    cntr = cntr + 3
ENDFOR

print,'Limit the Y-axis?'
read,ans
IF (ans EQ 'y') THEN BEGIN
    print,'Enter a new Y upper limit'
    read,ymax
ENDIF

ENDREP UNTIL (ans EQ 'n')

print,'Save the plot?'
read, ans

IF (ans EQ 'y') THEN BEGIN
    print,'Add a note to the plot?'
    read,ans1
    IF (ans1 EQ 'y') THEN BEGIN
        note = strarr(1)
        print,'Type away...'
        read,note
        note = 'NOTE: '+note
    ENDIF
    color = [60,110,150,180]
    set_plot,'ps'
    device,file=plotnames[0,0]+'_losvd.ps',/color
    loadct,0
    plot,radius[0:1],data[1,0:1],title='LOSVD measurements for '+plotnames[0,0],ytitle='!7r!6(!ukm!n/!dsec!n)',/nodata,/ynozero,$
      charsize=1.5,xrange=[xmin,xmax],xstyle=1,yrange=[ymin,ymax],ystyle=1,xtitle='Radius (arcsec)'
    loadct,4
    cntr = 0
    cntr = 0
    FOR k=0,n/ptnnum-1 DO BEGIN
        plots,radius[k],data[1,k+cntr],psym=2,symsize=1.2,thick=2,color=color[3]
        xyouts,radius[k]+2,data[1,k+cntr]+2,plotnames[2,k+cntr],charsize=1,color=color[3],orientation=45
        plots,radius[k],data[1,k+cntr+1],psym=2,symsize=1.2,thick=2,color=color[0] 
        xyouts,radius[k]+2,data[1,k+cntr+1]+2,plotnames[2,k+cntr+1],charsize=1,color=color[0],orientation=45
        plots,radius[k],data[1,k+cntr+2],psym=2,symsize=1.2,thick=2,color=color[1]
        xyouts,radius[k]+2,data[1,k+cntr+2]+2,plotnames[2,k+cntr+2],charsize=1,color=color[1],orientation=45
        plots,radius[k],data[1,k+cntr+3],psym=2,symsize=1.2,thick=2,color=color[2]
        xyouts,radius[k]+2,data[1,k+cntr+3]+2,plotnames[2,k+cntr+3],charsize=1,color=color[2],orientation=45
        cntr = cntr + 3
    ENDFOR
    IF (ans1 EQ 'y') THEN BEGIN
        loadct,0
        xyouts,0.22,0.19,note,charsize=0.7,/normal
    ENDIF

    device,/close
    set_plot,'x'
ENDIF

print,'Next ENTER deletes the plot.'
pause

wdelete,2
cd,'../'
STOP
END


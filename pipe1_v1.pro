;This is the first piece of the pipeline, which ties together
;radius.pro and radplot.pro, and works from lists. Each pointing will
;be sent through this piece, followed by PIPE2.PRO

;CHANGE MADE ON 01-17-09
;AS THE NEW PIPELINE WORKS FROM THE UNCOLLAPSED DATA THIS NEW VERSION
;OF PIPE1 USES THE PEFSM.FITS FILES AND CORRESPONDING PEFW.FITS FILES
;TO DO A QUICK COLLAPSE OF THE DATA BEFORE BEING FED INTO
;RADPLOT_FUNC.PRO. 

;The quick collapse is of the form
;         collapsed = sum_over_i(pefsm_i * pefw_i)/sum_over_i(pefw_i)

;Another alteration is eliminating the .rad files. Now, all the
;values are kept in the .R file with the following format.

; fiber #   Radius   dRA    dDec   Flux (normalized)   Flux (fiber median)

; The dRA and dDec values (deltaRA, deltaDec) is the offset from the
; center of the galaxy as entered by you in the procedure call.
;****************************************************************************
PRO pipe1, pointing, galRA, galDec, Data=datalist 
;****************************************************************************
;EX: 'n1600cent','hh:mm:ss','+dd:mm:ss','n1600cent.list'

;POINTING: The name used in finder_code_gamma, WITHOUT the "_D1_coords.txt"
;GALRA AND GALDEC: The galaxy RA and Dec (feed in as strings)
;DATALIST: A list that contains all the frames FOR A GIVEN POINTING. In the
;case of multiple sky-subtractions, only one needs to be sent through
;as the name.R file will hold for all of them.
;(NOTE: If the datalist keyword isn't used, the routine looks for a
;list named pointing.list.)

;NOTE THAT THEY SHOULD ALL BE THE SAME LEVEL OF SUBTRACTION AS THE
;UN-NORMALIZED FLUX IS USED AS A WEIGHT TO DETERMINE THE BINNING
;LOCATIONS IN THE VISPLOT AND OTHER ROUTINES...

;The format of the list is the same as used in pipe2.pro, so some of
;the columns are unnecessary.
;WARNING: YOU WILL NEED TO PAY ATTENTION TO THE LEVEL OF SKY-SUBTRACTION IN
;THE LIST. IF THESE LEVELS VARY MUCH OTHER ROUTINES THAT USE THE
;ESTIMATED FLUX WILL BE THROWN OFF.

;The format of the list should be as follows:
;   jm0302cc_   a   _jan08_n2.dat    jm0302.R


;****************************************************************************
;-----------------------------------------------------------------
nprofile = 5 ;this is the number of elements in your fiber profile
;-----------------------------------------------------------------

ans=''
print,'Write out the collapsed files as name_pefsmC.fits files? (y or n)'
read,ans

if (n_elements(datalist) eq 0) then datalist = pointing+'.list'
readcol,datalist,format='a,x,a,a',filenames,dat,dotR
print,'The files in the list are:'
print,transpose(filenames)
pause 

file = pointing+'_D1_coords.txt'
RA = strn(galRA) & Dec = strn(galDec)
RA = strsplit(RA,':',/EXTRACT)
Dec = strsplit(Dec,':',/EXTRACT)
RA = double((RA[0]+RA[1]/60.+RA[2]/3600.)*15)
;The + or - sign of the declination is accounted for
IF (Dec[0] LT 0) THEN Dec = double(Dec[0]-Dec[1]/60.-Dec[2]/3600.)
IF (Dec[0] GE 0) THEN Dec = double(Dec[0]+Dec[1]/60.+Dec[2]/3600.)

;The RA and Dec for each fiber is read in.
readcol,file,format='I,A,A',fib,fra,fdec
n1 = n_elements(fib)
rad = dblarr(4,n1) ;The radius array (fiber#, rad, deltaRA, deltaDec)
cor = dblarr(3,n1) ;The coordinate array (fiber#, RA, Dec)

FOR j=0,n1-1 DO BEGIN ;A loop through each fiber
    cor[0,j] = fib[j]
    fc = fra[j]
    fd = fdec[j]
    fc = strsplit(fc,':',/EXTRACT)
    fd = strsplit(fd,':',/EXTRACT)
;The coordinates are converted to decimal degrees
    cor[1,j] = double((fc[0]+fc[1]/60.+fc[2]/3600.)*15.)
    IF (fd[0] LT 0) THEN cor[2,j] = double(fd[0]-fd[1]/60.-fd[2]/3600.)
    IF (fd[0] GT 0) THEN cor[2,j] = double(fd[0]+fd[1]/60.+fd[2]/3600.)
ENDFOR

;The radial distance from the galaxy center, deltaRA and deltaDec are
;calculated, then converted back to arcseconds.
FOR j=0,n1-1 DO BEGIN
    dx = double((cor[1,j]-RA))
    dy = double((cor[2,j]-Dec))
    rad[0,j] = j+1
    rad[1,j] = sqrt(dx^2+dy^2)*3600.
    rad[2,j] = dx*3600.
    rad[3,j] = dy*3600.
;    Print, 'The radius for fiber '+strn(j+1)+ ' is: ',strn(rad[1,j])
ENDFOR

xmin = min(cor[1,*])-0.001
xmax = max(cor[1,*])+0.003
ymin = min(cor[2,*])-0.002
ymax = max(cor[2,*])+0.002

;A visual confirmation of the positions of each fiber. This should
;always look just like a plot of the IFU head.

window,0,retain=2,xsize=550,ysize=450,xpos=0,ypos=700
device,decomposed=0
loadct,0
plot,cor[1,*],cor[2,*],psym=3,xtitle='RA',ytitle='Dec',/nodata,$
  title='Radial positions of each fiber relative to galaxy center',$
  charsize=1.2,xrange=[xmin,xmax],yrange=[ymin,ymax],xstyle=1,ystyle=1
loadct,4
for j=0,n1-1 do xyouts,cor[1,j],cor[2,j],strn(round(rad[1,j])),color=150

openw,6,pointing+'.rad'
for j=0,n1-1 do begin
    printf,6,uint(rad[0,j]),rad[1,j],rad[2,j],rad[3,j]
endfor
free_Lun,6
print,'The file '+pointing+'.rad has been written.'

;Now that the "name.rad" file is created, it is fed into a modified
;version of radplot. Here, individual fibers are rejected and the
;"name.R" files are generated.

for j=0,n_elements(filenames)-1 do begin ;a loop through each file
    data = readfits(filenames[j]+'pefsm.fits',/silent)
    weight = readfits(filenames[j]+'pefw.fits',/silent)
    maskname = 'mask'+dat[j]

;the data is collapsed
    datac = wgtcollapse(data,weight,maskname,nprofile)
    if (ans eq 'y') then writefits,filenames[j]+'pefsmC.fits',datac

;the following is all just for plotting purposes.
    values = dblarr(n1)
    for k=0,n1-1 do values[k] = median(datac[*,k])
    n666 = where(values eq -666)
    if (n666[0] ne -1) then values[n666] = 0.0
    values = double(values/max(values))
    values = round(values*250)

    ans1=''
    print,'Plot contour lines? (y/n)'
    read,ans1

    lvl = findgen(20)*15+min(values[where(values gt 0)])

    window,2,retain=2,xsize=480,ysize=450,xpos=800,ypos=700
    device,decomposed=0
    loadct,0
    plot,cor[1,*],cor[2,*],psym=3,xtitle='RA',ytitle='Dec',/nodata,$
      title='Fiber Number (colored by flux level)',$
      charsize=1.2,xrange=[xmin,xmax],yrange=[ymin,ymax],xstyle=1,ystyle=1
    loadct,27
    if (ans1 eq 'y') then contour,values,cor[1,*],cor[2,*],$
      /irregular,/overplot,levels=lvl,min_value=1
    for k=0,n1-1 do xyouts,cor[1,k],cor[2,k],strn(k+1),color=values[k]
    print,''
    print,'Working on file '+dotR[j]
    print,''

    window,1,retain=2,xsize=480,ysize=450,xpos=340,ypos=700
    device,decomposed=0
    loadct,0
    plot,cor[1,*],cor[2,*],psym=3,xtitle='RA',ytitle='Dec',/nodata,$
      title='Relative Flux in Each Fiber (normalized to 250)',$
      charsize=1.2,xrange=[xmin,xmax],yrange=[ymin,ymax],xstyle=1,ystyle=1
    loadct,27
    if (ans1 eq 'y') then contour,values,cor[1,*],cor[2,*],$
      /irregular,/overplot,levels=lvl,min_value=1
    for k=0,n1-1 do xyouts,cor[1,k],cor[2,k],strn(values[k]),color=values[k]
    print,''
    print,'Working on file '+dotR[j]
    print,''

;Radplot is now a function. The string can be either 'm' for median
;or 't' for total. This is simply what radplot.pro plots to aid in the
;visual rejection. i.e. the median or total of the fiber vs. radius 
    nameout = radplot_func(pointing+'.rad',datac,dotR[j],'m')

endfor

wdelete,0,1,2

END


;----------------------------------------------------------------------
;The code was modifed on 11-11-08 to accept a list of recipe files and
;is now flexible enough to accept any list of ranges.
;It expects to find a list called recipe.list, or use the name keyword
;call to enter your own name.
;This code will use the first parts of the recipe files as the names
;for the bins

;  |
;  |                     %
;T |                         %
;Y |
;P |
;E |%
;  | %              %                   %
;  |   %   %    %                %
;  -----------------------------------------------------
;              BINNING REGION
;
;This is the general format of the plots.
;----------------------------------------------------------------------

PRO recipe, NAME= recipelist

ans=''
ans1=''

if (n_elements(recipelist) eq 0) then recipelist = 'recipe.list'
readcol,recipelist,Format='a',rlist


n1 = n_elements(rlist)
readcol,rlist[0],format='d,x',temp
n3 = n_elements(temp) 

readcol,'Tlist',Format='x,a,a,x',tlist,metalicity
n2 = n_elements(tlist)
tlist = reverse(tlist)

if (n2 ne n3) then begin
    print,'Your Tlist does not match your recipes in length!'
    print,'Try again.'
    stop
endif

starmix = fltarr(n1,n2,4)
; The star mix has 4 layers, the actual recipe values, the star index
; (as referenced from Tlist) the rounded % of the recipe (this is what
; gets plotted). The fourth layer is the plotting color
names = strarr(n1)

;---------------------------------------------------------------
sss = ''
print,''
print,rlist[0]
print,'What element do you want to split the name on?'
read,sss
sss = sss+'.'
plotname = strsplit(rlist[0],sss,/extract)
print,plotname
pause

;---------------------------------------------------------------

for j=0,n1-1 do begin
    readcol,rlist[j],format='d,i',temp1,temp2
    idx = bsort(temp1)
    starmix[j,*,0] = reverse(temp1)
    starmix[j,*,1] = reverse(temp2)
    starmix[j,*,2] = round(starmix[j,*,0]*100)
    starmix[j,*,3] = starmix[j,*,0]*255
    temp = strsplit(rlist[j],sss,/extract)
    names[j] = temp[0]
endfor

xaxis = findgen(n1)+1

window,0,retain=2,xsize=1100,ysize=600
device,decomposed=0
shiftx = 0.0
shifty = 0.0

jump1:
loadct,0
plot,xaxis,starmix[*,0],/nodata,xrange=[-1,n1+4],yrange=[-10,n2+1],$
  xstyle=1,ystyle=1,title='Stellar Recipe',ytitle='Stellar Type'
     
loadct,27

for j=0,n1-1 do begin
    xyouts,j+3+shiftx,-9.5,names[j],orientation=90,charsize=1.0,color=0
   for k=0,n2-1 do begin
      xyouts,j+3,k-1.5,string(uint(starmix[j,k,2])),color=starmix[j,k,3],orientation=90
      xyouts,-0.5,k+0.35+shifty,tlist[k],charsize=1.0,color=0
   endfor
endfor

print,'Make a X or Y shift?'
print,'(Enter "yx", "yy" or "n")'
read,ans
if (ans eq 'yx') then begin
    print,'Shift the bin names by how much?'
    read,shiftx
    goto, jump1
endif
if (ans eq 'yy') then begin
    print,'Shift the stellar names by how much?'
    read,shifty
    goto, jump1
endif

print,'Save the plot?'
read,ans
;ans='n'

IF (ans EQ 'y') THEN BEGIN
;   print,'Add a note to the plot?'
;   read,ans1
;   IF (ans1 EQ 'y') THEN BEGIN
;      note = strarr(1)
;      print,'Type away...'
;      read,note
;      note = 'NOTE: '+note
;   ENDIF
   
   set_plot,'ps'
   device,file='recipe.ps',/color,/landscape
   loadct,0
   
   plot,xaxis,starmix[*,0],/nodata,xrange=[-1,n1+6],yrange=[-10,n2+1],$
     xstyle=1,ystyle=1,title='Stellar Recipe',ytitle='Stellar Type',$
     xthick=2,ythick=2,charthick=2
   loadct,27
   
   for j=0,n1-1 do begin
      xyouts,j+5+shiftx,-9.5,names[j],orientation=90,charsize=0.8,color=150,charthick=1.5
      for k=0,n2-1 do begin
         xyouts,j+5,k-1.3,string(uint(starmix[j,k,2])),color=starmix[j,k,3],charsize=0.8,orientation=90
         xyouts,-0.5,k+0.35+shifty,tlist[k],charsize=1.0,color=150,charthick=1.5
      endfor
   endfor
;   IF (ans1 EQ 'y') THEN BEGIN
;      loadct,0
;      xyouts,0.15,0.95,note,charsize=1.0,/normal,charthick=1.5
;   ENDIF
   
   device,/close_file
   set_plot,'x'
   
ENDIF ELSE print,'Word. No plot for you!'

print,'Next ENTER deletes the plot.'
pause
wdelete,0

STOP
END

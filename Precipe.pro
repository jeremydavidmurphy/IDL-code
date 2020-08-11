PRO Precipe, list

; This routine reads in a list of name.recipe files. It is generally
; the case that these recipes will be from one spectral region.

; It will look for the Tlist file...

readcol,'Tlist',format='a,a,f,i',star,type,metal,x
readcol,list,format='a',files
n0 = n_elements(files)

readcol,files[0],silent=1,format='f',test
n1 = n_elements(test)

if (n_elements(star) ne n1) then stop

recipe = fltarr(n0,n1)
fnames = strarr(n0)

for j=0,n0-1 do begin
    file = files[j]
    readcol,silent=1,file,format='f',temp
    recipe[j,*] = temp
    file = strsplit(file,'_',/extract)
    fnames[j] = file[0]
endfor

rcolors = intarr(n0)+60
incre = 195./n0
for j=0,n0-1 do rcolors[j] = rcolors[j]+(incre*j)
scolors = intarr(n1)+60
incre = 195./n1
for j=0,n1-1 do scolors[j] = scolors[j]+(incre*j)

ans = ''
print,'Plot to the screen? (y)'
print,'(A postscript file, named recipe.ps, will be generated automatically...)'

read,ans

if (ans eq 's') then begin
    window,0,retain=2
    device,decomposed=0
    loadct,0
    plot,recipe[0,*],psym=2,title='Recipe for '+list,yrange=[0,1],$
      xrange=[-1,n1+1],xstyle=1,/nodata,ytitle='Percentage',xtitle='Star'
    loadct,4
    for j=0,n0-1 do begin
        oplot,recipe[j,*],psym=2,color=rcolors[j]
    endfor
    
    for j=0,n1-1 do begin
        xyouts,j,0.87,type[j],color=scolors[j],orientation=90,charsize=1.5
    endfor
endif

set_plot,'ps'
device,filename=list+'.ps',/color
loadct,0
plot,recipe[0,*],psym=2,title='Recipe for '+list,yrange=[0,1],$
  xrange=[-1,n1+1],xstyle=1,/nodata,ytitle='Percentage',xtitle='Star',$
  xthick=3,ythick=3,charthick=3

loadct,4
for j=0,n0-1 do begin
    oplot,recipe[j,*],psym=2,color=rcolors[j],thick=2
endfor

for j=0,n1-1 do begin
    xyouts,j,0.87,type[j],color=scolors[j],orientation=90,charsize=1.2,$
      charthick=2
endfor
device,/close_file
set_plot,'x'


STOP 
END

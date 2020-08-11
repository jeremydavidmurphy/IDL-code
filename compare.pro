;This routine is used to visually compare the quality of spectra from
;different exposures of the same fiber. It simply overplots the
;spectra to allow for visual inspection. It works from a list of fits
;files that are assumed to be 1D and of the same dimension and
;pointing on a galaxy

PRO compare, list

wi = 3530
wd = 1.12

color = [60,110,150,180,255]
yup = intarr(1) & ydown = intarr(1)

readcol,list,format='A',files

test = readfits(files[0])
n1 = n_elements(test[*,0])
n2 = n_elements(test[0,*]) 
n3 = n_elements(files)

data = dblarr(n1,n2,n3)
data[*,*,0] = test
for l=1,n3-1 do data[*,*,l] = readfits(files[l])
data[where(data eq -666)] = 0.0

wave = dblarr(n1)
for j=0,n1-1 do wave[j] = wi + (j*wd)

w1 = wave[1]
w2 = wave[n1-1]
decision=''

repeat begin
   fiber=intarr(1)
   print,'Choose a fiber to plot:'
   read,fiber
   fiber = fiber-1

jump1:
   firsttitle = 'Spectra for fiber '+strn(fiber+1)
   
   window,0,retain=2
   device,decomposed=0
   loadct,0
   plot,wave,data[*,fiber,0],xtitle='Wavelength (A)',ytitle='Flux',$
     title=firsttitle,charsize=1.2,xrange=[w1,w2],/nodata,/ynozero
   oplot,wave,data[*,fiber,0]
   xyouts,0.7,0.35,files[0],/normal,charsize=1.5
   loadct,4
   if (n_elements(color) le n3) then begin
       for l=1,n3-1 do begin
           oplot,wave,data[*,fiber,l],color=color[l-1]
           xyouts,0.7,0.35-(l*0.05),files[l],/normal,charsize=1.5,color=color[l-1]
       endfor
       loadct,0
   endif else begin
       for l=1,n3-1 do begin
           oplot,wave,data[*,fiber,l],color=60+(20*l)
           xyouts,0.7,0.35-(l*0.05),files[l],/normal,charsize=1.5,color=60+(20*l)
       endfor
       loadct,0
   endelse

jump2:
   print,'What next?'
   print,'("n" for new fiber, "w" for wavelength range, "y" for y-axis change, "d" for done)'
   read,decision
   if (decision eq 'w') then begin
       print,'Enter a new starting wavelength:'
       read,w1
       print,'Enter a new ending wavelength:'
       read,w2
       goto, jump1
   endif
   if (decision eq 'y') then begin
       print,'Enter a new lower limit for the Y-axis:'
       read,ydown
       print,'Enter a new upper limit for the Y-axis:'
       read,yup
       plot,wave,data[*,fiber,0],xtitle='Wavelength (A)',ytitle='Flux',$
         title=firsttitle,charsize=1.2,xrange=[w1,w2],yrange=[ydown,yup],$
         /nodata,ystyle=1
       oplot,wave,data[*,fiber,0]
       xyouts,0.7,0.35,files[0],/normal,charsize=1.5
       loadct,4
       if (n_elements(color) le n3) then begin
           for l=1,n3-1 do begin
               oplot,wave,data[*,fiber,l],color=color[l-1]
               xyouts,0.7,0.35-(l*0.05),files[l],/normal,charsize=1.5,color=color[l-1]
           endfor
           loadct,0
       endif else begin
           for l=1,n3-1 do begin
               oplot,wave,data[*,fiber,l],color=60+(20*l)
               xyouts,0.7,0.35-(l*0.05),files[l],/normal,charsize=1.5,color=60+(20*l)
           endfor
           loadct,0
       endelse
       goto, jump2
   endif

endrep until (decision eq 'd')
   
ans=''
print,'Save a hardcopy of this plot? (y or n)'
read,ans

if (ans eq 'y') then begin
    name = ''
    print,'Name the output file (w/o the .ps):'
    read,name
    name = name+'.ps'
   set_plot,'ps'
   device,file=name,/color
   loadct,0
   if (yup ne 0) then begin
       plot,wave,data[*,fiber,0],xtitle='Wavelength (A)',ytitle='Flux',$
         title=firsttitle,xrange=[w1,w2],xstyle=1,yrange=[ydown,yup],$
         ystyle=1,charsize=1.2,/nodata,/ynozero,xthick=2,ythick=2,charthick=2
   endif else begin
       plot,wave,data[*,fiber,0],xtitle='Wavelength (A)',ytitle='Flux',$
         title=firsttitle,xrange=[w1,w2],xstyle=1,charsize=1.2,/nodata,/ynozero,$
         xthick=2,ythick=2,charthick=2
   endelse
   oplot,wave,data[*,fiber,0],color=60
   xyouts,0.7,0.35,files[0],/normal,charsize=1.2
   loadct,4
   if (n_elements(color) le n3) then begin
       for l=1,n3-1 do begin
           oplot,wave,data[*,fiber,l],color=color[l-1]
           xyouts,0.7,0.35-(l*0.05),files[l],/normal,charsize=1.2,$
             color=color[l-1],charthick=2
       endfor
       loadct,0
   endif else begin
       for l=1,n3-1 do begin
           oplot,wave,data[*,fiber,l],color=60+(20*l)
           xyouts,0.7,0.35-(l*0.05),files[l],/normal,charsize=1.2,$
             color=60+(20*l),charthick=2
       endfor
       loadct,0
   endelse
   device,/close
   set_plot,'x'
   print,'The file has been saved as compare.ps.'
endif

wdelete,0

STOP
END

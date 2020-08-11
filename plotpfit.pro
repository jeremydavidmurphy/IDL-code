PRO plotpfit, pfitlovlist, radiusfile, scale=scale, galname=galname

; This routine plots the output of pfitlov. It reads in a list of
; pfitlov files (pfitlovlist). Each of these must be the same length
; as the radiusfile. If the median switches are turned on below
; (either pmedian or jmedian) then a median of all the pfitlov files
; in the list is plotted. The colors indicate the different angular
; bins (black,blue,green,orange,red), while symbols indicate the
; different pfitlov.out files being used.

; These set the plotting ranges...
yv = [-300,300]
yd = [100,350]
yh3 = [-0.2,0.2]
yh4 = [-0.2,0.2]

; PFITLOVLIST: a list of pfitlov files you want plotted. THIS IS A
; LIST, EVEN IF YOU WANT TO PLOT A SINGLE FILE. So, make a list that
; has a single file in it in that case.

; RADIUSFILE: A list of 2 columns of the following format:
; 01     -7.002
; 05     -6.417
; 02    -11.912
; 05     -8.919
; 01    -14.254
; 04    -12.814

; Where the first column is the angular bin and the second is the
; radial distance, in arcseconds.
 
; SCALE: This is the scale of the galaxy (parsecs/arcseconds). If you
; include this, the 2nd x-axis is plotted along the top, given in
; units of kpc.

; The output is a postscript figure showing the four moments. The
; postscript will get named, based on the name of your
; pfitlovlist. EX: if your list is pfitlov.list, your ps file will be
; named pfitlov_moments.ps. If you give a galaxy name, you will get
; galaxyname_moments.ps

; MODIFIED ON FEB 18, 2012: This code now accepts in a LIST of
; pfitlov.out files. If the list is just a single pfitlov file, then
; it follows the same routine as before, where the bin.radius file is
; used to color the data points depending on their angular
; position. If there are more than one (i.e. for different spectra
; regions, then it uses a change in symbols to denote the different
; pfitlov.out files, keeping the color for the angular bins.

velshift = 'on' ;set to 'on' and the velocity will be redshifted to force median(velocity) = 0.0
pmedian = 'off' ;set to 'on' to OVERplot the median ON TOP OF all the different pfitlov.out files.
jmedian = 'on' ;set to 'on' to plot JUST THE MEDIAN of all the different spectral regions.
;**********************************************************************

readcol,pfitlovlist,silent=1,f='a',list
n0 = n_elements(list)

if (n0 eq 1) then begin ;the switches are overridden
   pmedian = 'off'
   jmedian = 'off'
endif

cstep = floor(255.0/float(n0))
colors = floor((indgen(n0)+1)*cstep)
   
for j=0,n0-1 do begin ;a loop through each pfitlov.out file in the list
   
   if (j eq 0) then begin
      readcol,radiusfile,silent=1,f='i,f',bin2,radius
      n1 = n_elements(radius)
      moments = fltarr(4,n1,n0)
      pmoments = moments
      if (min(radius) lt 0) then x1 = min(radius)+(0.15*min(radius))
      if (min(radius) gt 0) then x1 = min(radius)-(0.15*min(radius))
      x2 = max(radius)+(0.15*max(radius))
      i1 = where(bin2 eq 01)
      i2 = where(bin2 eq 02)
      i3 = where(bin2 eq 03)
      i4 = where(bin2 eq 04)
      i5 = where(bin2 eq 05)
      r1 = 0 & r2 = 0 & r3 = 0 & r4 = 0
   endif
   
   pfitfile = list[j]
   readcol,pfitfile,silent=1,f='a,f,f,x,f,f,x,x,x',bin1,vel,disp,h3,h4
      
   if (velshift eq 'on') then begin
      shift = median(vel)
      vel = vel - shift
      if (j eq 0) then offset = shift else offset = [offset,shift]
   endif
      
   moments[0,*,j] = vel
   pmoments[0,*,j] = vel
   ibad = where(vel gt yv[1] or vel lt yv[0],c1)
   if (c1 gt 0) then pmoments[0,ibad,j] = 1e6

   moments[1,*,j] = disp
   pmoments[1,*,j] = disp
   ibad = where(disp gt yd[1] or disp lt yd[0],c2)
   if (c2 gt 0) then pmoments[1,ibad,j] = 1e6

   moments[2,*,j] = h3
   pmoments[2,*,j] = h3
   ibad = where(h3 gt yh3[1] or h3 lt yh3[0],c3)
   if (c3 gt 0) then pmoments[2,ibad,j] = 1e6

   moments[3,*,j] = h4
   pmoments[3,*,j] = h4
   ibad = where(h4 gt yh4[1] or h4 lt yh4[0],c4)
   if (c4 gt 0) then pmoments[3,ibad,j] = 1e6

   r1 = c1+r1 & r2 = c2+r2 & r3 = c3+r3 & r4 = c4+r4
      
endfor
    
if (velshift eq 'on' and n0 gt 1) then moffset = median(offset,/even)
if (velshift eq 'on' and n0 eq 1) then moffset = offset

if (pmedian eq 'on' or jmedian eq 'on') then begin
   mvel = median(moments[0,*,*],dim=3)
   mdisp = median(moments[1,*,*],dim=3)
   mh3 = median(moments[2,*,*],dim=3)
   mh4 = median(moments[3,*,*],dim=3)
endif

if (n_elements(galname) gt 0) then nameout = galname+'_moments.ps' else begin
   s = strsplit(pfitlovlist,'.',/extract)
   nameout = s[0]+'_moments.ps'
endelse

set_plot,'ps'
device,file=nameout,/color
!p.region = [0.1,0.1,0.9,0.9]
!p.multi = [0,1,4,0,1]
!y.omargin = [2,4]
      
;VELOCITY
if (n_elements(scale) ne 0) then begin
   loadct,0,/silent
   plot,radius,pmoments[0,*,0],xthick=2,ythick=2,charthick=2,psym=sym(6),$
        xrange=[x1,x2],xstyle=9,yrange=[yv[0],yv[1]],/ys,pos=[0.15,0.675,0.85,0.85],$
        xtickformat='(A1)',xticklen=0.07,yminor=2,/nodata,symsize=1.5
   xyouts,0.105,0.7,'Velocity (km/s)',/normal,orientation=90,charsize=0.6,$
          charthick=1.5
   xyouts,0.88,0.755,'('+strn(r1)+')',/normal,charsize=0.8,charthick=1.5
   if (velshift eq 'on') then xyouts,0.86,0.80,'Offset = '+strn(moffset),/normal,charsize=0.6,charthick=1.5
   oplot,[x1,x2],[0.0,0.0],linestyle=1
endif else begin
   loadct,0,/silent
   plot,radius,pmoments[0,*,0],xthick=2,ythick=2,charthick=2,psym=sym(6),$
        xrange=[x1,x2],xstyle=1,yrange=[yv[0],yv[1]],/ys,pos=[0.15,0.675,0.85,0.85],$
        xtickformat='(A1)',xticklen=0.07,yminor=2,/nodata,symsize=1.5
   xyouts,0.105,0.7,'Velocity (km/s)',/normal,orientation=90,charsize=0.6,$
          charthick=1.5
   xyouts,0.88,0.755,'('+strn(r1)+')',/normal,charsize=0.8,charthick=1.5
   if (velshift eq 'on') then xyouts,0.86,0.80,'Offset = '+strn(moffset),/normal,charsize=0.6,charthick=1.5
   oplot,[x1,x2],[0.0,0.0],linestyle=1
endelse
loadct,33,/silent
      
if (jmedian eq 'on') then goto,jumpv
      
for j=0,n0-1 do begin
   plots,radius[i2],pmoments[0,i2,j],psym=sym(7),color=colors[j],thick=3,symsize=0.9
   plots,radius[i3],pmoments[0,i3,j],psym=sym(8),color=colors[j],thick=3,symsize=0.9
   plots,radius[i4],pmoments[0,i4,j],psym=sym(16),color=colors[j],thick=3,symsize=0.9
   plots,radius[i5],pmoments[0,i5,j],psym=sym(17),color=colors[j],thick=3,symsize=0.9
   plots,radius[i1],pmoments[0,i1,j],psym=sym(1),color=colors[j],thick=3,symsize=0.9
endfor
      
if (pmedian eq 'on') then begin
   jumpv:
   loadct,0,/silent
   plots,radius[i2],mvel[i2],psym=sym(10),thick=3,symsize=0.9
   plots,radius[i3],mvel[i3],psym=sym(10),thick=3,symsize=0.9
   plots,radius[i4],mvel[i4],psym=sym(10),thick=3,symsize=0.9
   plots,radius[i5],mvel[i5],psym=sym(13),thick=3,symsize=0.9
   plots,radius[i1],mvel[i1],psym=sym(1),thick=3,symsize=0.9
endif
      
if (n_elements(scale) ne 0) then begin ;(the second xaxis is put in kpc)
   irad = bsort(radius)
   radkpc = (radius * scale) / 1000.0
   radkpc = radkpc[irad]
   x12 = min(radkpc)-((25.0*scale)/1000.0)
   x22 = max(radkpc)+((25.0*scale)/1000.0)
   loadct,0,/silent
   axis,xaxis=1,xrange=[x12,x22],/save,xtitle='Radius (kpc)',$
        charthick=2,xthick=2,charsize=1.3,xstyle=1,xticklen=0.07
endif
      
;VELOCITY DISPERSION
loadct,0,/silent
plot,radius,pmoments[1,*,0],xthick=2,ythick=2,charthick=2,psym=sym(6),$
     xrange=[x1,x2],/xs,yrange=[yd[0],yd[1]],/ys,pos=[0.15,0.50,0.85,0.675],$
     xtickformat='(A1)',xticklen=0.07,yminor=2,/nodata,symsize=1.5
xyouts,0.105,0.51,'Dispersion (km/s)',/normal,orientation=90,charsize=0.6,$
       charthick=1.5
xyouts,0.88,0.58,'('+strn(r2)+')',/normal,charsize=0.8,charthick=1.5
loadct,33,/silent

if (jmedian eq 'on') then goto,jumpd

for j=0,n0-1 do begin
   plots,radius[i2],pmoments[1,i2,j],psym=sym(7),color=colors[j],thick=3,symsize=0.9
   plots,radius[i3],pmoments[1,i3,j],psym=sym(8),color=colors[j],thick=3,symsize=0.9
   plots,radius[i4],pmoments[1,i4,j],psym=sym(16),color=colors[j],thick=3,symsize=0.9
   plots,radius[i5],pmoments[1,i5,j],psym=sym(17),color=colors[j],thick=3,symsize=0.9
   plots,radius[i1],pmoments[1,i1,j],psym=sym(1),color=colors[j],thick=3,symsize=0.9
endfor

if (pmedian eq 'on') then begin
   jumpd:
   loadct,0,/silent
   plots,radius[i2],mdisp[i2],psym=sym(10),thick=3,symsize=0.9
   plots,radius[i3],mdisp[i3],psym=sym(10),thick=3,symsize=0.9
   plots,radius[i4],mdisp[i4],psym=sym(10),thick=3,symsize=0.9
   plots,radius[i5],mdisp[i5],psym=sym(13),thick=3,symsize=0.9
   plots,radius[i1],mdisp[i1],psym=sym(1),thick=3,symsize=0.9
endif

;H3
loadct,0,/silent
plot,radius,pmoments[2,*,0],xthick=2,ythick=2,charthick=2,psym=sym(6),$
     xrange=[x1,x2],/xs,yrange=[yh3[0],yh3[1]],/ys,pos=[0.15,0.325,0.85,0.50],$
     xtickformat='(A1)',xticklen=0.07,yminor=2,/nodata,symsize=1.5
xyouts,0.105,0.412,'H3',/normal,orientation=90,charsize=0.6,$
       charthick=1.5
oplot,[x1,x2],[0.0,0.0],linestyle=1
xyouts,0.88,0.41,'('+strn(r3)+')',/normal,charsize=0.8,charthick=1.5
loadct,33,/silent

if (jmedian eq 'on') then goto,jumph3

for j=0,n0-1 do begin
   plots,radius[i2],pmoments[2,i2,j],psym=sym(7),color=colors[j],thick=3,symsize=0.9
   plots,radius[i3],pmoments[2,i3,j],psym=sym(8),color=colors[j],thick=3,symsize=0.9
   plots,radius[i4],pmoments[2,i4,j],psym=sym(16),color=colors[j],thick=3,symsize=0.9
   plots,radius[i5],pmoments[2,i5,j],psym=sym(17),color=colors[j],thick=3,symsize=0.9
   plots,radius[i1],pmoments[2,i1,j],psym=sym(1),color=colors[j],thick=3,symsize=0.9
endfor

if (pmedian eq 'on') then begin
   jumph3:
   loadct,0,/silent
   plots,radius[i2],mh3[i2],psym=sym(10),thick=3,symsize=0.9
   plots,radius[i3],mh3[i3],psym=sym(10),thick=3,symsize=0.9
   plots,radius[i4],mh3[i4],psym=sym(10),thick=3,symsize=0.9
   plots,radius[i5],mh3[i5],psym=sym(13),thick=3,symsize=0.9
   plots,radius[i1],mh3[i1],psym=sym(1),thick=3,symsize=0.9
endif

;H4
loadct,0,/silent
plot,radius,pmoments[3,*,0],xthick=2,ythick=2,charthick=2,psym=sym(6),$
     xrange=[x1,x2],/xs,yrange=[yh4[0],yh4[1]],/ys,pos=[0.15,0.15,0.85,0.325],$
     xticklen=0.07,yminor=2,/nodata,symsize=1.5,xtitle='Radius (arcsec)',$
     charsize=1.0,xcharsize=1.3
xyouts,0.105,0.237,'H4',/normal,orientation=90,charsize=0.6
oplot,[x1,x2],[0.0,0.0],linestyle=1
charthick=1.5
xyouts,0.88,0.24,'('+strn(r4)+')',/normal,charsize=0.8,charthick=1.5
loadct,33,/silent

if (jmedian eq 'on') then goto,jumph4

for j=0,n0-1 do begin
   plots,radius[i2],pmoments[3,i2,j],psym=sym(7),color=colors[j],thick=3,symsize=0.9
   plots,radius[i3],pmoments[3,i3,j],psym=sym(8),color=colors[j],thick=3,symsize=0.9
   plots,radius[i4],pmoments[3,i4,j],psym=sym(16),color=colors[j],thick=3,symsize=0.9
   plots,radius[i5],pmoments[3,i5,j],psym=sym(17),color=colors[j],thick=3,symsize=0.9
   plots,radius[i1],pmoments[3,i1,j],psym=sym(1),color=colors[j],thick=3,symsize=0.9
endfor

if (pmedian eq 'on') then begin
   jumph4:
   loadct,0,/silent
   plots,radius[i2],mh4[i2],psym=sym(10),thick=3,symsize=0.9
   plots,radius[i3],mh4[i3],psym=sym(10),thick=3,symsize=0.9
   plots,radius[i4],mh4[i4],psym=sym(10),thick=3,symsize=0.9
   plots,radius[i5],mh4[i5],psym=sym(13),thick=3,symsize=0.9
   plots,radius[i1],mh4[i1],psym=sym(1),thick=3,symsize=0.9
endif

device,/close_file
set_plot,'x'

stop
END

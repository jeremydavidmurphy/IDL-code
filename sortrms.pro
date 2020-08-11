; This is a different routine than plotrms_v1.pro. This routine works
; from the output of filter.pro. It also requires a list of the bin
; names. 

;----------------------------------------------------------------
pro sortrms
;----------------------------------------------------------------

;This routine requires three files to exist in the calling directory:
;PFITLOV.ALL, RUN.FINAL, and BIN.LIST

;The run.final file isn't used for anything. It just comes
;along for the ride so as to get sorted in the same fashion that
;pfitlov.all is. THE OUTPUT OF THIS ROUTINE IS PFITLOV.ALL.STD AND
;RUN.FINAL.STD.

;NOTE THAT THE RMS VALUES ARE NOW CORRECTED FOR THEIR DIFFERENCE IN
;EXPOSURE TIMES.

readcol,silent=1,'pfitlov.all',format='a,f,f,f,f,f,f,f,i,f',$
        pfile,vel,sig1,sig2,h3,h4,r1,r2,r3,rms
n1 = n_elements(pfile)

readcol,silent=1,'run.final',format='a,a,i,i,f,i,i,a',$
        rfit,rfile,w1,w2,z,v,sm,reg

readcol,silent=1,'bin.list',format='a',bins

if (n_elements(pfile) ne n_elements(rfile)) then stop
n0 = n_elements(bins)

;the pfitlov.all names are split up
pfilearr = strarr(3,n1)
for j=0,n1-1 do begin
   temp = strsplit(pfile[j],'_.',/extract)
   pfilearr[*,j] = temp[0:2]
endfor

; The different levels of sky subtraction are picked off.
temp = pfilearr[1,*]
mag = strarr(n1)
pmag = strarr(n1)
for j=0,n1-1 do begin
    t = strsplit(temp[j],'f',/extract)
    mag[j] = t[0]
endfor

i1 = where(mag eq 'aa') ;4/3.5
i2 = where(mag eq 'ab' or mag eq 'ba') ;4/3.625
i3 = where(mag eq 'ac' or mag eq 'ca' or mag eq 'bb') ;4/3.75
i4 = where(mag eq 'ad' or mag eq 'da' or mag eq 'bc' or mag eq 'cb') ;4/3.875
i5 = where(mag eq 'ae' or mag eq 'ea' or mag eq 'bd' or mag eq 'db' or mag eq 'cc') ;4/4
i6 = where(mag eq 'cd' or mag eq 'dc' or mag eq 'be' or mag eq 'eb');4/4.125
i7 = where(mag eq 'ce' or mag eq 'ec' or mag eq 'dd') ;4/4.25
i8 = where(mag eq 'de' or mag eq 'ed') ;4/4.375
i9 = where(mag eq 'ee') ;4/4.5

rmsc = fltarr(n1)
if (i1[0] ne -1) then rmsc[i1] = rms[i1]*4.0/3.5
if (i2[0] ne -1) then rmsc[i2] = rms[i2]*4.0/3.625
if (i3[0] ne -1) then rmsc[i3] = rms[i3]*4.0/3.75
if (i4[0] ne -1) then rmsc[i4] = rms[i4]*4.0/3.875
if (i5[0] ne -1) then rmsc[i5] = rms[i5]*4.0/4.0
if (i6[0] ne -1) then rmsc[i6] = rms[i6]*4.0/4.125
if (i7[0] ne -1) then rmsc[i7] = rms[i7]*4.0/4.25
if (i8[0] ne -1) then rmsc[i8] = rms[i8]*4.0/4.375
if (i9[0] ne -1) then rmsc[i9] = rms[i9]*4.0/4.5

; a loop through each bin. this creates 'indy' which is the index for
; the order of the output frames. They are sorted by bin first, then,
; within each bin, by the corrected rms values.
indy = intarr(1)
steps = intarr(n0)
for j=0,n0-1 do begin
   bin = bins[j]
   piece = where(pfilearr[0,*] eq bin)
   if (piece[0] ne -1) then begin
       steps[j] = n_elements(piece)
       rmscpiece = rmsc[piece]
       piece = piece[bsort(rmscpiece)]
       indy = [indy,piece]
   endif
endfor
indy = indy[1:*]

fp = '(A27,2x,F10.3,2x,F8.3,2x,F8.3,2x,F6.3,2x,F6.3,2x,F8.3,2x,F8.3,2x,I5,2x,D12.10)'
fr = '(A7,1x,A16,2x,I4,2x,I4,2x,F8.6,2x,I3,1x,I3,1x,A6)'

openw,5,'pfitlov.all.std'
for j=0,n1-1 do begin
   i = indy[j]
   printf,5,pfile[i],vel[i],sig1[i],sig2[i],$
          h3[i],h4[i],r1[i],r2[i],r3[i],rmsc[i],format=fp
endfor
free_lun,5

openw,5,'run.final.std'
for j=0,n1-1 do begin
   i = indy[j]
   printf,5,rfit[i],rfile[i],w1[i],w2[i],z[i],v[i],sm[i],reg[i]
endfor
free_lun,5

;they've been written, but we all like plots...

ans=''
print,'Plot to screen? (y or n)'
read,ans
if (ans eq 'y') then begin
    set_plot,'x'
   window,0, retain=2,xsize=900,ysize=900
   device,decomposed=0
   
   for j=0,n0-1 do begin
      bin = bins[j]
      piece = where(pfilearr[0,*] eq bin)
      if (piece[0] ne -1) then begin
         prms = rmsc[piece]
         i = piece[bsort(prms)]
         n = n_elements(i)
         pm = mag[i]+' '+pfilearr[2,i]
         kolor = intarr(n)
         ii = where(pfilearr[2,i] eq 'HK15')
         if (ii[0] ne -1) then kolor[ii] = 60
         ii = where(pfilearr[2,i] eq 'GB15')
         if (ii[0] ne -1) then kolor[ii] = 110
         ii = where(pfilearr[2,i] eq 'HB15')
         if (ii[0] ne -1) then kolor[ii] = 180
         ii = where(pfilearr[2,i] eq 'MG15')
         if (ii[0] ne -1) then kolor[ii] = 150
         ii = where(pfilearr[2,i] eq 'MG15b')
         if (ii[0] ne -1) then kolor[ii] = 150
         pv = vel[i]
         ps = sig1[i]
         ph3 = h3[i]
         ph4 = h4[i]
         prms = rmsc[i]
         yaxis = indgen(n)+1
         y2 = max(yaxis)+2
         loadct,0
         !p.multi = [0,2,2,0,0]
         x1 = min(pv)-15 & x2 = max(pv)+5
         if (x1 lt -150) then x1 = -150
         if (x2 gt 150) then x2 = 150
         plot,pv,yaxis,title=bin,xtitle='Velocity',ytitle='RMS Range:  '+strn(min(prms))+'        '+strn(max(prms)),$
              yrange=[0,y2],xrange=[x1,x2],xstyle=1,ystyle=1,charsize=1.2,/nodata
         loadct,4
         for k=0,n-1 do plots,pv[k],yaxis[k],psym=sym(1),color=kolor[k]
         for k=0,n-1 do xyouts,x1+2,yaxis[k],pm[k],charsize=1.0,color=kolor[k]
         
         loadct,0
         !p.multi = [2,2,2,0,0]
         x1 = min(ps)-25 & x2 = max(ps)+5
;         if (x1 lt 0) then x1 = 0
;         if (x2 gt 400) then x2 = 400
         plot,ps,yaxis,title=bin,xtitle='Velocity Dispersion',ytitle='RMS Range:  '+strn(min(prms))+'        '+strn(max(prms)),$
              yrange=[0,y2],xrange=[x1,x2],xstyle=1,ystyle=1,charsize=1.2,/nodata
         loadct,4
         for k=0,n-1 do plots,ps[k],yaxis[k],psym=sym(1),color=kolor[k]
         for k=0,n-1 do xyouts,x1+2,yaxis[k],pm[k],charsize=1.0,color=kolor[k]
         
         loadct,0
         !p.multi = [3,2,2,0,0]
         x1 = min(ph3)-0.04 & x2 = max(ph3)+0.01
         if (x1 lt -0.5) then x1 = -0.5
         if (x2 gt 0.7) then x2 = 0.7
         plot,ph3,yaxis,title=bin,xtitle='H!d3!n',ytitle='RMS Range:  '+strn(min(prms))+'        '+strn(max(prms)),$
              yrange=[0,y2],xrange=[x1,x2],xstyle=1,ystyle=1,charsize=1.2,/nodata
         loadct,4
         for k=0,n-1 do plots,ph3[k],yaxis[k],psym=sym(1),color=kolor[k]
         for k=0,n-1 do xyouts,x1+0.005,yaxis[k],pm[k],charsize=1.0,color=kolor[k]
         
         loadct,0
         !p.multi = [1,2,2,0,0]
         x1 = min(ph4)-0.04 & x2 = max(ph4)+0.01
         if (x1 lt -0.5) then x1 = -0.5
         if (x2 gt 0.7) then x2 = 0.7
         plot,ph4,yaxis,title=bin,xtitle='H!d4!n',ytitle='RMS Range:  '+strn(min(prms))+'        '+strn(max(prms)),$
              yrange=[0,y2],xrange=[x1,x2],xstyle=1,ystyle=1,charsize=1.2,/nodata
         loadct,4
         for k=0,n-1 do plots,ph4[k],yaxis[k],psym=sym(1),color=kolor[k]
         for k=0,n-1 do xyouts,x1+0.005,yaxis[k],pm[k],charsize=1.0,color=kolor[k]
         pause
         
      endif
   endfor
endif

print,'Generate postscript files? (y or n):'
read,ans

if (ans eq 'y') then begin
;Now .ps files are created.
   set_plot,'ps'
;   !p.region = [0.1, 0.1, 0.9, 0.9]
   for j=0,n0-1 do begin
      bin = bins[j]
      piece = where(pfilearr[0,*] eq bin)
      if (piece[0] ne -1) then begin
         prms = rmsc[piece]
         i = piece[bsort(prms)]
         n = n_elements(i)
         pm = mag[i]+' '+pfilearr[2,i]
         kolor = intarr(n)
         ii = where(pfilearr[2,i] eq 'HK15')
         if (ii[0] ne -1) then kolor[ii] = 60
         ii = where(pfilearr[2,i] eq 'GB15')
         if (ii[0] ne -1) then kolor[ii] = 110
         ii = where(pfilearr[2,i] eq 'HB15')
         if (ii[0] ne -1) then kolor[ii] = 180
         ii = where(pfilearr[2,i] eq 'MG15')
         if (ii[0] ne -1) then kolor[ii] = 150
         ii = where(pfilearr[2,i] eq 'MG15b')
         if (ii[0] ne -1) then kolor[ii] = 150
         pv = vel[i]
         ps = sig1[i]
         ph3 = h3[i]
         ph4 = h4[i]
         prms = rmsc[i]
         yaxis = indgen(n)+1
         y2 = max(yaxis)+2

         device,file=bin+'_vel.ps',/color,/landscape
         loadct,0
;         !p.multi = [0,2,2,0,0]
         !x.margin = [4,2]
         !y.margin = [2,2]
;         !p.region = [0.1,0.5,0.5,0.9]
         x1 = min(pv)-15 & x2 = max(pv)+5
         if (x1 lt -150) then x1 = -150
         if (x2 gt 150) then x2 = 150
         plot,pv,yaxis,title=bin,xtitle='Velocity',ytitle='RMS Range:  '+strn(min(prms))+' to '+strn(max(prms)),$
           yrange=[0,y2],xrange=[x1,x2],xstyle=1,ystyle=1,charsize=2.0,/nodata,$
           xthick=3,ythick=3,charthick=3
         loadct,4
         for k=0,n-1 do plots,pv[k],yaxis[k],psym=sym(1),color=kolor[k],thick=2,symsize=1.5
         for k=0,n-1 do xyouts,x1+2,yaxis[k],pm[k],charsize=1.0,color=kolor[k],charthick=3
         device,/close_file

         device,file=bin+'_disp.ps',/color,/landscape
         loadct,0
         !x.margin = [4,2]
         !y.margin = [2,2]
;         !p.multi = [2,2,2,0,0]
;         !p.region = [0.1,0.1,0.5,0.5]
         x1 = min(ps)-25 & x2 = max(ps)+5
;         if (x1 lt 0) then x1 = 0
;         if (x2 gt 400) then x2 = 400
         plot,ps,yaxis,title=bin,xtitle='Velocity Dispersion',ytitle='RMS Range:  '+strn(min(prms))+' to '+strn(max(prms)),$
           yrange=[0,y2],xrange=[x1,x2],xstyle=1,ystyle=1,charsize=2.0,/nodata,$
           xthick=3,ythick=3,charthick=3
         loadct,4
         for k=0,n-1 do plots,ps[k],yaxis[k],psym=sym(1),color=kolor[k],thick=2,symsize=1.5
         for k=0,n-1 do xyouts,x1+2,yaxis[k],pm[k],charsize=1.0,color=kolor[k],charthick=3
         device,/close_file
;         
         device,file=bin+'_h3.ps',/color,/landscape
         loadct,0
         !x.margin = [4,2]
         !y.margin = [2,2]
;         !p.multi = [3,2,2,0,0]
;         !p.region = [0.5,0.5,0.9,0.9]
         x1 = min(ph3)-0.04 & x2 = max(ph3)+0.01
         if (x1 lt -0.5) then x1 = -0.5
         if (x2 gt 0.7) then x2 = 0.7
         plot,ph3,yaxis,title=bin,xtitle='H!d3!n',ytitle='RMS Range:  '+strn(min(prms))+' to '+strn(max(prms)),$
           yrange=[0,y2],xrange=[x1,x2],xstyle=1,ystyle=1,charsize=2.0,/nodata,$
           xthick=3,ythick=3,charthick=3
         loadct,4
         for k=0,n-1 do plots,ph3[k],yaxis[k],psym=sym(1),color=kolor[k],thick=2,symsize=1.5
         for k=0,n-1 do xyouts,x1+0.005,yaxis[k],pm[k],charsize=1.0,color=kolor[k],charthick=3
         device,/close_file
;         
         device,file=bin+'_h4.ps',/color,/landscape
         loadct,0
         !x.margin = [4,2]
         !y.margin = [2,2]
;         !p.multi = [1,2,2,0,0]
;         !p.region = [0.5,0.1,0.9,0.5]
         x1 = min(ph4)-0.04 & x2 = max(ph4)+0.01
         if (x1 lt -0.5) then x1 = -0.5
         if (x2 gt 0.7) then x2 = 0.7
         plot,ph4,yaxis,title=bin,xtitle='H!d4!n',ytitle='RMS Range:  '+strn(min(prms))+' to '+strn(max(prms)),$
           yrange=[0,y2],xrange=[x1,x2],xstyle=1,ystyle=1,charsize=2.0,/nodata,$
           xthick=3,ythick=3,charthick=3
         loadct,4
         for k=0,n-1 do plots,ph4[k],yaxis[k],psym=sym(1),color=kolor[k],thick=2,symsize=1.5
         for k=0,n-1 do xyouts,x1+0.005,yaxis[k],pm[k],charsize=1.0,color=kolor[k],charthick=3
         device,/close_file
      endif
   endfor
   set_plot,'x'
endif

stop
end

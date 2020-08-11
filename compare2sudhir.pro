PRO compare2sudhir, me, sudhir

;this routine is used to compare your kinematics with
;Sudhir's. The "H" difference is just that now the h3 and h4
;terms are added

;ME: a list of pfitlov.out files, to be median combined
;SUDHIR: his single file of the form:
; radius,vel,velerr,sigma,sigmaerr,h3,h3err,h4,h4err

Pvel = [-400,400]
Pdisp = [0,400]

readcol,sudhir,f='f,f,f,f,f',rad,v1,ve1,d1,de1
n1 = n_elements(rad)

readcol,me,f='a',myfiles
n0 = n_elements(myfiles)

if n0 gt 1 then begin
   array = fltarr(2,n1,n0)
   error = fltarr(2,n1)

   for j=0,n0-1 do begin
      file = myfiles[j]
      readcol,file,f='x,f,f,x,x,x,x,x,x',a,b
;   readcol,file,f='x,f,x,f,x,x,x,x,x',a,b ;gaussian
      array[0,*,j] = a
      array[1,*,j] = b
   endfor
   
   for j=0,n1-1 do begin
      error[0,j] = stddev(array[0,j,*])
      error[1,j] = stddev(array[1,j,*])
   endfor

   myvel = median(array[0,*,*],dim=3)
   mydisp = median(array[1,*,*],dim=3)
   ms = median(v1)
   mm = median(myvel)
   v1s = v1 - ms

;velocity
   set_plot,'ps'
   device,file='velocity.ps',/color
   loadct,0
   plot,rad,v1s,psym=2,symsize=2,yrange=[Pvel[0],Pvel[1]],$
        xtitle='Radius (arcsec)',ytitle='Velocity (km/s)',$
        charthick=2,ythick=2,xthick=2,thick=3
;   oplot,rad,array[0,*,5]-mm,psym=sym(1)
   oplot,rad,array[0,*,0]-mm,psym=4,color=180
   oplot,rad,array[0,*,1]-mm,psym=4,color=90
   oplot,rad,array[0,*,2]-mm,psym=4,color=120
   oplot,rad,array[0,*,3]-mm,psym=4,color=60
   oplot,rad,array[0,*,4]-mm,psym=4,color=150
   loadct,4
   oplot,rad,myvel-mm,psym=4,symsize=2,color=150
   oplot,rad,myvel-mm,psym=4,symsize=1.9,color=150
   oplot,rad,myvel-mm,psym=4,symsize=1.8,color=150
   device,/close
   
;dispersion
   device,file='dispersion.ps',/color
   loadct,0
   plot,rad,d1,psym=2,symsize=2,yrange=[Pdisp[0],Pdisp[1]],$
        xtitle='Radius (arcsec)',ytitle='Velocity Dispersion(km/s)',$
        charthick=2,ythick=2,xthick=2,thick=3
;   oplot,rad,array[1,*,5],psym=sym(1)
   oplot,rad,array[1,*,0],psym=4,color=180
   oplot,rad,array[1,*,1],psym=4,color=90
   oplot,rad,array[1,*,2],psym=4,color=120
   oplot,rad,array[1,*,3],psym=4,color=60
   oplot,rad,array[1,*,4],psym=4,color=150
   loadct,4
   oplot,rad,mydisp,psym=4,symsize=2,color=150
   oplot,rad,mydisp,psym=4,symsize=1.9,color=150
   oplot,rad,mydisp,psym=4,symsize=1.8,color=150
   device,/close

endif else begin
   readcol,myfiles,f='x,f,f,x,x,x,x,x,x',a,b
   ms = median(v1)
   mm = median(a)
   v1s = v1 - ms
   myvel = a
   mydisp = b

   set_plot,'ps'
   device,file='velocity.ps',/color
   loadct,0
   plot,rad,v1s,psym=2,symsize=2,yrange=[Pvel[0],Pvel[1]],$
        xtitle='Radius (arcsec)',ytitle='Velocity (km/s)',$
        charthick=2,ythick=2,xthick=2,thick=3
   loadct,4
   oplot,rad,myvel-mm,psym=4,symsize=2,color=150
   oplot,rad,myvel-mm,psym=4,symsize=1.9,color=150
   oplot,rad,myvel-mm,psym=4,symsize=1.8,color=150
   device,/close
   
   device,file='dispersion.ps',/color
   loadct,0
   plot,rad,d1,psym=2,symsize=2,yrange=[Pdisp[0],Pdisp[1]],$
        xtitle='Radius (arcsec)',ytitle='Velocity Dispersion(km/s)',$
        charthick=2,ythick=2,xthick=2,thick=3
   loadct,4
   oplot,rad,mydisp,psym=4,symsize=2,color=150
   oplot,rad,mydisp,psym=4,symsize=1.9,color=150
   oplot,rad,mydisp,psym=4,symsize=1.8,color=150
   device,/close
endelse

stop
END

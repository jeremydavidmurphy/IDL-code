PRO skymodelcheck, month, skylist, FILE=finalspectra

; This routine is a modification to skycheck.pro. The difference here
; is that it's looking specifically at the models for
; individual fibers rather than the global average. This is just a
; simplified plotting routine that will overplot the individual model
; sky fibers on top of the median (or biwt) combination of all of them.

; MONTH: The month of the observing (same as what you put into vaccine.param
; SKYLIST: The list of sky nods and model skies. 
;example:

;* *
;venga0149_
;venga0150_
;venga0152_
;* *
;venga0825_
;venga0826_

; (you can just trim down the science_n#.list.2 frames to get these)
; The nights are separated by * *. This will lead the routine to
; iterate on the night. NOTE YOU HAVE TO START WITH THE * *.

wi = 3550.0
biwtormed = 'med' ;set to either 'med' or 'biwt' to make the collapse with a median or biweight. the median is much quicker.
xd = [4050,4250] ;the plotted wavelength range
yd = [15,35] ;the plotted y-range (for the data)
yr = [-3,3] ;the plotted y-range (for the residuals)
;************************************************************************

wave = findgen(2048)*1.125 + wi
readcol,skylist,f='a,a',modelframes,silent=1
n0 = n_elements(modelframes)
ncntr = 1
cntr = 0
framenames = ''

for j=0,n0-1 do begin
    if (modelframes[j] eq '*') then begin
        mask = 'mask_'+month+'_n'+strn(ncntr)+'.dat'
        ptow = 'ptow_'+month+'_n'+strn(ncntr)+'.dat'
        ncntr = ncntr + 1
        goto, jumpend
    endif

    framenames = [framenames,modelframes[j]]
    model = readfits(modelframes[j]+'pefy.fits',/silent)
    wgt = readfits(modelframes[j]+'pefw.fits',/silent)
    modelc = wgtcollapsef(model,wgt,mask,5)
    modela = realigncf(modelc,ptow,wi,1.125,1)
    n1 = n_elements(modela[0,*])

    if (biwtormed eq 'biwt') then begin
        bmodel = biweightf(modela)
    endif
    if (biwtormed eq 'med') then begin
        bmodel = median(modela,dim=2)
    endif

    set_plot,'ps'
    device,file=modelframes[j]+'model_comparison.ps',/color
    loadct,0,/silent
    !p.multi = [0,1,2,0,0]
    plot,wave,bmodel,yrange=[yd[0],yd[1]],xthick=2,ythick=2,charthick=2,$
         position=[0.1,0.3,0.9,0.9],xtickformat='(A1)',ytitle='Pixel Counts',$
         xrange=[xd[0],xd[1]], ystyle=1,xstyle=1
    loadct,4,/silent
    for k=0,n1-1 do oplot,wave,modela[*,k],color=k
    xyouts,0.15,0.81,'Model: '+modelframes[j],charsize=1.3,/normal,$
           color=150,charthick=2.0
    loadct,0,/silent
    oplot,wave,bmodel,thick=3
    oplot,[4100,4100],[0,100]
    oplot,[4210,4210],[0,100]
    oplot,[4090,4090],[0,100],linestyle=5
    oplot,[4200,4200],[0,100],linestyle=5
    oplot,[4110,4110],[0,100],linestyle=5
    oplot,[4220,4220],[0,100],linestyle=5
    plot,wave,bmodel,yrange=[yr[0],yr[1]],xthick=2,ythick=2,charthick=2,$
         position=[0.1,0.1,0.9,0.3],ytitle='Model-Sky',$
         xtitle='Wavelength (A)',ystyle=1,xrange=[xd[0],xd[1]],/nodata
    loadct,4,/silent
    for k=0,n1-1 do begin
       diff = bmodel - modela[*,k]
       oplot,wave,diff,color=k
    endfor
    loadct,0,/silent
    oplot,[4100,4100],[yr[0],yr[1]]
    oplot,[4210,4210],[yr[0],yr[1]]
    oplot,[4090,4090],[yr[0],yr[1]],linestyle=5
    oplot,[4200,4200],[yr[0],yr[1]],linestyle=5
    oplot,[4110,4110],[yr[0],yr[1]],linestyle=5
    oplot,[4220,4220],[yr[0],yr[1]],linestyle=5
    device,/close


    device,file=modelframes[j]+'model_comparisonALL.ps',/color
    loadct,0,/silent
    !p.multi = [0,1,2,0,0]
    plot,wave,bmodel,yrange=[10,50],xthick=2,ythick=2,charthick=2,$
         position=[0.1,0.3,0.9,0.9],xtickformat='(A1)',ytitle='Pixel Counts',$
         xrange=[3700,5500], ystyle=1,xstyle=1
    loadct,4,/silent
    for k=0,n1-1 do oplot,wave,modela[*,k],color=k
    xyouts,0.15,0.81,'Model: '+modelframes[j],charsize=1.3,/normal,$
           color=150,charthick=2.0
    loadct,0,/silent
    oplot,wave,bmodel,thick=2
    plot,wave,bmodel,yrange=[-3,3],xthick=2,ythick=2,charthick=2,$
         position=[0.1,0.1,0.9,0.3],ytitle='% (Model-Sky)/Model',$
         xtitle='Wavelength (A)',ystyle=1,xrange=[3700,5500],/nodata
    loadct,4,/silent
    for k=0,n1-1 do begin
       diff = bmodel - modela[*,k]
       oplot,wave,diff,color=k
    endfor
    loadct,0,/silent
    oplot,[4100,4100],[yr[0],yr[1]]
    oplot,[4210,4210],[yr[0],yr[1]]
    oplot,[4090,4090],[yr[0],yr[1]],linestyle=5
    oplot,[4200,4200],[yr[0],yr[1]],linestyle=5
    oplot,[4110,4110],[yr[0],yr[1]],linestyle=5
    oplot,[4220,4220],[yr[0],yr[1]],linestyle=5
    device,/close
    set_plot,'x'
    !p.multi = [0,1,1,0,0]
    jumpend:
endfor
free_lun,5

stop
END

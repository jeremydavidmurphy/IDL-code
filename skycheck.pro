PRO skycheck, month, skylist, FILE=finalspectra

; This routine generates .ps files to explore how well vaccine is
; handling the sky modeling. It reads in a list of all of the skies
; you want looked at. YOU CAN FEED IN AS MANY NIGHTS AS YOU WANT.

; MONTH: The month of the observing (same as what you put into
; vaccine.param). The routine will then look for the flat, ptow and
; mask files for that night!

; SKYLIST: The list of sky nods and model skies. 
;example:

;* *
;venga0149_ sky0148c_0151c_
;venga0150_ sky0148c_0151c_
;venga0152_ sky0151c_0154c_
;* *
;venga0825_ sky0824c_0827c_
;venga0826_ sky0824c_0827c_

; (you can just trim down the science_n#.list.2 frames to get these)
; The nights are separated by * *. This will lead the routine to
; iterate on the night. NOTE YOU HAVE TO START WITH THE * *.

; MODIFIED ON DEC 19TH, 2011: If a "FILE" is added then that frame
; will get plotted along with versions of all the skies subtracted
; from it (sorta like Figure #15 in your M87 paper). This works well
; for understanding where the rapidly changing emission lines are at.

wi = 3550.0
skyfiber = 100 ;the fiber in the model sky frame to plot when plotting sky residuals. This is just for the residual plot. For the frame-to-frame comparison the median (or biweight) of the entire sky is completed.
biwtormed = 'med' ;set to either 'med' or 'biwt' to make the collapse with a median or biweight. the median is much quicker.
xd = [3500,6000] ;the plotted wavelength range
yd = [0,100] ;the plotted y-range (for the data)
yr = [-5,5] ;the plotted y-range (for the residuals)
;************************************************************************

wave = findgen(2048)*1.125 + wi
readcol,skylist,f='a,a',modelframes,skyframes,silent=1
n0 = n_elements(skyframes)
ncntr = 1
cntr = 0
framenames = ''
skyfiber = skyfiber - 1

if (n_elements(file) eq 1) then begin
    skyarray = fltarr(2048,n0)
    skycompare = 'on'
endif else skycompare = 'off'

for j=0,n0-1 do begin ;a loop over each frame
    if (skyframes[j] eq '*') then begin
        name1 = 'flat_'+month+'_n'+strn(ncntr)+'en.fits'
        flat = readfits(name1,h)
        mask = 'mask_'+month+'_n'+strn(ncntr)+'.dat'
        ptow = 'ptow_'+month+'_n'+strn(ncntr)+'.dat'
        ncntr = ncntr + 1
        goto, jumpend
    endif
    framenames = [framenames,modelframes[j]]
    sky = readfits(skyframes[j]+'pe.fits',header,/silent)
    if n_elements(sky) eq 0 or sky[0] eq -1 then goto, jumpend
    skyf = sky / flat
    model = readfits(modelframes[j]+'pefy.fits',/silent)
    if n_elements(model) eq 0 or model[0] eq -1 then goto, jumpend
    wgt = readfits(modelframes[j]+'pefw.fits',/silent)
    if n_elements(wgt) eq 0 or wgt[0] eq -1 then goto, jumpend
    skyc = wgtcollapsef(skyf,wgt,mask,5)
    skya = realigncf(skyc,ptow,wi,1.125,1)
    modelc = wgtcollapsef(model,wgt,mask,5)
    modela = realigncf(modelc,ptow,wi,1.125,1)
    if (biwtormed eq 'biwt') then begin
        bsky = biweightf(skya)
        bmodel = biweightf(modela)
    endif
    if (biwtormed eq 'med') then begin
;        bsky = median(skya,dim=2)
        bskya = median(skya[*,0:50],dim=2)
        bsky = median(skya,dim=2)
        bmodel = median(modela,dim=2)
        bmodela = median(modela[*,0:50],dim=2)
    endif
    diff = bsky - bmodel
    div = ((bmodel-bsky) / bmodel) * 100.0
    diffa = bskya - bmodela
    diva = ((bmodela-bskya) / bmodela) * 100.0
;    div = ((bsky-bmodel) / bmodel) * 100.0
    if (skycompare eq 'on') then skyarray[*,j] = modela[*,skyfiber]

;    v1 = median(diff)
;    v2 = stddev(diff)
;    v1 = float(v1)
;    v2 = float(v2)
;    v1b = median(diff[71:271])
;    v1r = median(diff[1351:1551])
;    v2b = stddev(diff[71:271])
;    v2r = stddev(diff[1351:1551])

    v1 = median(div)
    v2 = stddev(div)
    v1 = float(v1)
    v2 = float(v2)
    v1b = median(div[30:330]) ;~3584A to 3921A
    v1r = median(div[1189:1489]);~4887A to 5225A
    v2b = stddev(div[80:380])
    v2r = stddev(div[1189:1489])

    if (cntr eq 0) then begin
;        openw,5,month+'_difference_values.txt'
        openw,5,month+'_percentdifference_values.txt'
        printf,5,['med_diff','std_diff','med_blue','std_blue','med_red','std_red']
        cntr = 1
    endif
    printf,5,[v1,v2,v1b,v2b,v1r,v2r],$
      format='(f7.4,2x,f7.4,2x,f7.4,2x,f7.4,2x,f7.4,2x,f7.4)'

    set_plot,'ps'
    device,file=modelframes[j]+'model_v_datasky.ps',/color
    loadct,0,/silent
    !p.multi = [0,1,2,0,0]
    plot,wave,bsky,yrange=[yd[0],yd[1]],xthick=2,ythick=2,charthick=2,$
         position=[0.1,0.3,0.9,0.9],xtickformat='(A1)',ytitle='Pixel Counts',$
         xrange=[xd[0],xd[1]], ystyle=1,xstyle=1
    xyouts,0.15,0.85,'Nods: '+skyframes[j],charsize=1.3,/normal,charthick=2.0
;    oplot,[4100,4100],[0,100]
;    oplot,[4210,4210],[0,100]
;    oplot,[4090,4090],[0,100],linestyle=5
;    oplot,[4200,4200],[0,100],linestyle=5
;    oplot,[4110,4110],[0,100],linestyle=5
;    oplot,[4220,4220],[0,100],linestyle=5
    loadct,4,/silent
;    oplot,wave,bskya,color=60
    xyouts,0.15,0.81,'Model: '+modelframes[j],charsize=1.3,/normal,$
      color=150,charthick=2.0
    oplot,wave,bmodel,color=150
;    oplot,wave,bmodela,color=110
    loadct,0,/silent
    xyouts,0.4,0.44,'Difference Median: '+strn(v1),/normal,$
      charthick=2,charsize=1.3
    xyouts,0.4,0.40,'Difference Stddev: '+strn(v2),/normal,$
      charthick=2,charsize=1.3
    loadct,4,/silent
    xyouts,0.4,0.36,'Difference Blue: '+strn(v1b),/normal,$
      charthick=2,charsize=1.3,color=60
    xyouts,0.4,0.32,'Difference Red: '+strn(v1r),/normal,$
      charthick=2,charsize=1.3,color=150
    loadct,0,/silent
;    plot,wave,diff,yrange=[-2,2],xthick=2,ythick=2,charthick=2,$
;      position=[0.1,0.1,0.9,0.3],ytitle='Model - Nods',xtitle='Wavelength (A)'
    plot,wave,div,yrange=[yr[0],yr[1]],xthick=2,ythick=2,charthick=2,$
         position=[0.1,0.1,0.9,0.3],ytitle='% (Model-Sky)/Model',$
         xtitle='Wavelength (A)',ystyle=1,xrange=[xd[0],xd[1]],xstyle=1
;    oplot,[4100,4100],[yr[0],yr[1]]
;    oplot,[4210,4210],[yr[0],yr[1]]
    loadct,4,/silent
;    oplot,wave,diva,color=60
    oplot,[xd[0],xd[1]],[0,0],color=150
    loadct,0,/silent
    device,/close
    set_plot,'x'
    !p.multi = [0,1,1,0,0]
    jumpend:
endfor
free_lun,5

goto, jumpend2
;now, all the values are plotted.
framenames = framenames[1:*]

readcol,month+'_difference_values.txt',f='f,f,f,f,f,f',am,as,bm,bs,rm,rs
files = framenames
n1 = n_elements(am)

x = findgen(n1)

y1 = max([bm,rm])+max([bs,rs])
y2 = min([bm,rm])-max([bs,rs])

set_plot,'ps'
device,file='blue_2_red_comparison_'+month+'.ps',/color

plot,x,bm,xrange=[-2,n1+2],/xs,yrange=[y2,y1],/ys,$
  xthick=2,ythick=2,charthick=2,xtitle='Files',/nodata,$
  ytitle='Median of the difference (model - realsky)'

for j=0,n1-1 do begin
    xyouts,x[j]+0.3,-0.55,files[j],orientation=90,charsize=0.8
endfor
loadct,4,/silent
oplot,x,bm,psym=sym(6),color=60
oplot,x+0.1,rm,psym=sym(6),color=150

for j=0,n1-1 do begin
    plots,[x[j],x[j]],[bm[j]-bs[j]*0.5,bm[j]+bs[j]*0.5],color=60
    plots,[x[j]+0.1,x[j]+0.1],[rm[j]-rs[j]*0.5,rm[j]+rs[j]*0.5],color=150
endfor

device,/close_file

if (skycompare eq 'on') then begin ;the sky variance is plotted
    wave1 = findgen(6144)*0.375+3550.0
    realframe = readfits(file,header)
    t = cont_normf(realframe)
    rfC = t[*,0]
    skyarrayC = skyarray
    var = fltarr(2048)
    for j=0,n0-1 do begin
        t = cont_normf(skyarray[*,j])
        skyarrayC[*,j] = t[*,0]
     endfor
    for j=0,2047 do begin
        t = moment(skyarrayC[j,*],/Nan)
        var[j] = t[1]
    endfor

    device,file='skyvariance.ps',/color
    !p.region = [0.05,0.05,0.95,0.95]
    !p.multi = [0,1,2,0,1]
    loadct,0,/silent
    plot,wave1,rfC,xrange=[3800,5800],/xs,/nodata,charthick=2,xthick=2,ythick=2,$
      ytitle='Continuum Normalized Counts',position=[0.1,0.4,0.9,0.9],xtickformat='(A1)'
    loadct,4,/silent
    for j=0,n0-1 do begin
        oplot,wave,skyarrayC[*,j],color=110
    endfor
    loadct,0,/silent
    oplot,wave1,rfC,thick=2

    plot,wave,var,xrange=[3800,5800],/xs,charthick=2,xthick=2,ythick=2,$
      xtitle='Observed Wavelength (A)',ytitle='Sky Variance',thick=3,$
      position=[0.1,0.1,0.9,0.4]
    device,/close_file
endif

jumpend2:

set_plot,'x'

stop
END

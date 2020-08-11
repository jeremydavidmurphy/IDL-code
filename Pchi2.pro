pro Pchi2, list, SORTON=sorton

; This routine accepts a list of SORTED CHI-SQUARED VALUES (the output
; of Pchi.pro) and then overplots the lowest chi-squared values
; against the parameter of interest.

; SORTON: Use this callword to sort the chi-squared values on either
; 'stars', 'gc' or 'both'

; WARNING: If you are sending several cres.#### files through then the
; number of different parameter values explored has to be either the
; same for all, or the largest in the first one. Otherwise, the array
; created will not be large enough and the routine will crash.

if (n_elements(sorton) eq 0) then begin
    sorton = ''
    print,'Enter which chi-squared values to plot (stars, gs, both):'
    read,sorton
endif

;bhr = ((findgen(1000)/20.)+1.)*1e9
;mlr = ((findgen(100)/10.)+1.)
;vcr = ((findgen(300)/0.2)+100.)
;rsr = ((findgen(50))+1.)
;vcr = long(vcr)
;rsr = long(rsr)

readcol,list,format='a',silent=1,files
n0 = n_elements(files)

for j=0,n0-1 do begin ;a loop through the number of cres files you've got

    readcol,files[j],silent=1,format='f,f,i,i,f,f,f',ml,bh,vc,rs,stars,gc,both
    if (sorton eq 'stars') then isort = bsort(stars)
    if (sorton eq 'gc') then isort = bsort(gc)
    if (sorton eq 'both') then isort = bsort(both)

    Sstars = stars[isort]
    Sgc = gc[isort]
    Sboth = both[isort]

    Sstarsf = 0
    Sgcf = 0
    Sbothf = 0
    Sbhf = 0
    Sbh = bh[isort]
    for k=0,n_elements(Sbh)-1 do begin
        t = where(Sbhf eq Sbh[k],c)
        if (c eq 0) then begin
            Sbhf = [Sbhf,Sbh[k]]
            Sstarsf = [Sstarsf,Sstars[k]]
            Sgcf = [Sgcf,Sgc[k]]
            Sbothf = [Sbothf,Sboth[k]]
        endif
    endfor
    bh1f = Sstarsf[1:*]
    bh2f = Sgcf[1:*]
    bh3f = Sbothf[1:*]
    Sbhf = Sbhf[1:*]

    Sstarsf = 0
    Sgcf = 0
    Sbothf = 0
    Smlf = 0
    Sml = ml[isort]
    for k=0,n_elements(Sml)-1 do begin
        t = where(Smlf eq Sml[k],c)
        if (c eq 0) then begin
            Smlf = [Smlf,Sml[k]]
            Sstarsf = [Sstarsf,Sstars[k]]
            Sgcf = [Sgcf,Sgc[k]]
            Sbothf = [Sbothf,Sboth[k]]
        endif
    endfor
    ml1f = Sstarsf[1:*]
    ml2f = Sgcf[1:*]
    ml3f = Sbothf[1:*]
    Smlf = Smlf[1:*]

    Sstarsf = 0
    Sgcf = 0
    Sbothf = 0
    Svcf = 0
    Svc = vc[isort]
    for k=0,n_elements(Svc)-1 do begin
        t = where(Svcf eq Svc[k],c)
        if (c eq 0) then begin
            Svcf = [Svcf,Svc[k]]
            Sstarsf = [Sstarsf,Sstars[k]]
            Sgcf = [Sgcf,Sgc[k]]
            Sbothf = [Sbothf,Sboth[k]]
        endif
    endfor
    vc1f = Sstarsf[1:*]
    vc2f = Sgcf[1:*]
    vc3f = Sbothf[1:*]
    Svcf = Svcf[1:*]

    Sstarsf = 0
    Sgcf = 0
    Sbothf = 0
    Srsf = 0
    Srs = rs[isort]
    for k=0,n_elements(Srs)-1 do begin
        t = where(Srsf eq Srs[k],c)
        if (c eq 0) then begin
            Srsf = [Srsf,Srs[k]]
            Sstarsf = [Sstarsf,Sstars[k]]
            Sgcf = [Sgcf,Sgc[k]]
            Sbothf = [Sbothf,Sboth[k]]
        endif
    endfor

    rs1f = Sstarsf[1:*]
    rs2f = Sgcf[1:*]
    rs3f = Sbothf[1:*]
    Srsf = Srsf[1:*]

    if (j eq 0) then begin
        tml = n_elements(Smlf)
        if (tml gt 1) then mlarr = fltarr(tml,3,n0)
        tbh = n_elements(Sbhf)
        if (tbh gt 1) then bharr = fltarr(tbh,3,n0)
        tvc = n_elements(Svcf)
        if (tvc gt 1) then vcarr = fltarr(tvc,3,n0)
        trs = n_elements(Srsf)
        if (trs gt 1) then rsarr = fltarr(trs,3,n0)
    endif

    if (tml gt 1) then begin
        mlarr[*,0,j] = Smlf
        if (sorton eq 'stars') then mlarr[*,1,j] = ml1f
        if (sorton eq 'gc') then mlarr[*,1,j] = ml2f
        if (sorton eq 'both') then mlarr[*,1,j] = ml3f
    endif

    if (tbh gt 1) then begin
        bharr[*,0,j] = Sbhf
        if (sorton eq 'stars') then bharr[*,1,j] = bh1f
        if (sorton eq 'gc') then bharr[*,1,j] = bh2f
        if (sorton eq 'both') then bharr[*,1,j] = bh3f
    endif

    if (tvc gt 1) then begin
        vcarr[*,0,j] = Svcf
        if (sorton eq 'stars') then vcarr[*,1,j] = vc1f
        if (sorton eq 'gc') then vcarr[*,1,j] = vc2f
        if (sorton eq 'both') then vcarr[*,1,j] = vc3f
    endif

    if (trs gt 1) then begin
        rsarr[*,0,j] = Srsf
        if (sorton eq 'stars') then rsarr[*,1,j] = rs1f
        if (sorton eq 'gc') then rsarr[*,1,j] = rs2f
        if (sorton eq 'both') then rsarr[*,1,j] = rs3f
    endif

endfor

; The data is now in a usable array...
if (tml gt 1) then begin
    for j=0,n0-1 do begin
        mlarr[*,2,j] = mlarr[*,1,j] + (max(mlarr[0,1,*]) - mlarr[0,1,j])
        i = bsort(mlarr[*,0,j])
        mlarr[*,*,j] = mlarr[i,*,j]
    endfor
endif

if (tbh gt 1) then begin
    for j=0,n0-1 do begin
        bharr[*,2,j] = bharr[*,1,j] + (max(bharr[0,1,*]) - bharr[0,1,j])
        i = bsort(bharr[*,0,j])
        bharr[*,*,j] = bharr[i,*,j]
    endfor
endif

if (tvc gt 1) then begin
    for j=0,n0-1 do begin
        vcarr[*,2,j] = vcarr[*,1,j] + (max(vcarr[0,1,*]) - vcarr[0,1,j])
        i = bsort(vcarr[*,0,j])
        vcarr[*,*,j] = vcarr[i,*,j]
    endfor
endif

if (trs gt 1) then begin
    for j=0,n0-1 do begin
        rsarr[*,2,j] = rsarr[*,1,j] + (max(rsarr[0,1,*]) - rsarr[0,1,j])
        i = bsort(rsarr[*,0,j])
        rsarr[*,*,j] = rsarr[i,*,j]
    endfor
endif

; THE PLOTTING STAGES...

if (n0 eq 1) then colors = 150 else colors = intarr(n0-1)
step = floor(200/(n0-1.))
for j=0,n0-2 do colors[j] = (j+1)*step 

if (tml gt 1) then begin
    set_plot,'x'
    window,0,retain=2
    device,decomposed=0
    xup = max(mlarr[*,0,*])+0.25 & xdn = min(mlarr[*,0,*])-0.25
    yup = max(mlarr[*,2,*])+25 & ydn = min(mlarr[*,2,*])-25
    loadct,0
    
    plot,mlarr[*,0,0],mlarr[*,2,0],psym=-1,/ynozero,xrange=[xdn,xup],/xstyle,$
      title='DM Halo M/L Chi-squared Values ('+sorton+')',$
      xtitle='M/L',ytitle='scaled chi-squared',charsize=1.2,$
      yrange=[ydn,yup],/ystyle
    xyouts,0.3,0.87,files[0],/normal,charsize=1.3
    xyouts,0.45,0.87,strn(mlarr[0,1,0]),/normal,charsize=1.3
    xyouts,0.45,0.91,'Min chi2 Value',/normal,charsize=1.3
    xyouts,0.30,0.91,'File',/normal,charsize=1.3
    
    loadct,27
    for j=1,n0-1 do begin
        oplot,mlarr[*,0,j],mlarr[*,2,j],psym=-1,color=colors[j-1]
        xyouts,0.3,0.87-(j*0.03),files[j],color=colors[j-1],/normal,charsize=1.3
        xyouts,0.45,0.87-(j*0.03),strn(mlarr[0,1,j]),/normal,charsize=1.3,color=colors[j-1]
    endfor
    
    ans = ''
    print,'Change the y-axis plotting range? (y/n)'
    read,ans
    if (ans eq 'y') then begin
        repeat begin
            yup = intarr(1)
            ydn = intarr(1)
            print,'Enter an upper Y range:'
            read,yup
            print,'Enter a lower Y range:'
            read,ydn
            loadct,0
            plot,mlarr[*,0,0],mlarr[*,2,0],psym=-1,/ynozero,xrange=[xdn,xup],/xstyle,$
              title='DM Halo M/L Chi-squared Values ('+sorton+')',$
              xtitle='M/L',ytitle='scaled chi-squared',yrange=[ydn,yup],/ystyle,charsize=1.2
            xyouts,0.3,0.87,files[0],/normal,charsize=1.0
            xyouts,0.45,0.87,strn(mlarr[0,1,0]),/normal,charsize=1.3
            xyouts,0.45,0.91,'Min chi2 Value',/normal,charsize=1.3
            xyouts,0.30,0.91,'File',/normal,charsize=1.3
            loadct,27
            for j=1,n0-1 do begin
                oplot,mlarr[*,0,j],mlarr[*,2,j],psym=-1,color=colors[j-1]
                xyouts,0.3,0.87-(j*0.03),files[j],color=colors[j-1],/normal,charsize=1.3
                xyouts,0.45,0.87-(j*0.03),strn(mlarr[0,1,j]),/normal,charsize=1.3,color=colors[j-1]
            endfor
            print,'Another change of the Y-axis? (y)'
            read,ans
        endrep until (ans ne 'y')
    endif
    
    set_plot,'ps'
    device,file='X2_'+sorton+'_ml.ps',/color
    
    loadct,0
    plot,mlarr[*,0,0],mlarr[*,2,0],psym=-1,/ynozero,xrange=[xdn,xup],/xstyle,$
      title='DM Halo M/L Chi-squared Values ('+sorton+')',$
      xtitle='M/L',ytitle='scaled chi-squared',charsize=1.0,$
      xthick=3,ythick=3,charthick=2,thick=3,yrange=[ydn,yup],/ystyle
    xyouts,0.3,0.84,files[0],/normal,charsize=1.1,charthick=3
    xyouts,0.45,0.84,strn(mlarr[0,1,0]),/normal,charsize=1.1,charthick=3
    xyouts,0.45,0.88,'Min chi2 Value',/normal,charsize=1.1,charthick=3
    xyouts,0.30,0.88,'File',/normal,charsize=1.1,charthick=3
    
    loadct,4
    for j=1,n0-1 do begin
        oplot,mlarr[*,0,j],mlarr[*,2,j],psym=-1,color=colors[j-1],thick=3
        xyouts,0.3,0.84-(j*0.03),files[j],color=colors[j-1],/normal,charsize=1.1,charthick=3
        xyouts,0.45,0.84-(j*0.03),strn(mlarr[0,1,j]),/normal,charsize=1.1,color=colors[j-1],charthick=3
    endfor
    device,/close_file
    pause
endif


if (tbh gt 1) then begin
    set_plot,'x'
    window,1,retain=2
    device,decomposed=0
    xup = max(bharr[*,0,*])+0.25 & xdn = min(bharr[*,0,*])-0.25
    yup = max(bharr[*,2,*])+25 & ydn = min(bharr[*,2,*])-25
    loadct,0
    
    plot,bharr[*,0,0],bharr[*,2,0],psym=-1,/ynozero,xrange=[xdn,xup],/xstyle,$
      title='DM Halo Black Hole Mass Chi-squared Values ('+sorton+')',$
      xtitle='BH Mass (e9 M_solar)',ytitle='scaled chi-squared',charsize=1.2,$
      yrange=[ydn,yup],/ystyle
    xyouts,0.3,0.87,files[0],/normal,charsize=1.3
    xyouts,0.45,0.87,strn(bharr[0,1,0]),/normal,charsize=1.3
    xyouts,0.45,0.91,'Min chi2 Value',/normal,charsize=1.3
    xyouts,0.30,0.91,'File',/normal,charsize=1.3
    
    loadct,27
    for j=1,n0-1 do begin
        oplot,bharr[*,0,j],bharr[*,2,j],psym=-1,color=colors[j-1]
        xyouts,0.3,0.87-(j*0.03),files[j],color=colors[j-1],/normal,charsize=1.3
        xyouts,0.45,0.87-(j*0.03),strn(bharr[0,1,j]),/normal,charsize=1.3,color=colors[j-1]
    endfor
    
    ans = ''
    print,'Change the y-axis plotting range? (y/n)'
    read,ans
    if (ans eq 'y') then begin
        repeat begin
            yup = intarr(1)
            ydn = intarr(1)
            print,'Enter an upper Y range:'
            read,yup
            print,'Enter a lower Y range:'
            read,ydn
            loadct,0
            plot,bharr[*,0,0],bharr[*,2,0],psym=-1,/ynozero,xrange=[xdn,xup],/xstyle,$
              title='DM Halo Black Hole Mass Chi-squared Values ('+sorton+')',$
              xtitle='BH Mass (e9 M_solar)',ytitle='scaled chi-squared',yrange=[ydn,yup],/ystyle,charsize=1.2
            xyouts,0.3,0.87,files[0],/normal,charsize=1.0
            xyouts,0.45,0.87,strn(bharr[0,1,0]),/normal,charsize=1.3
            xyouts,0.45,0.91,'Min chi2 Value',/normal,charsize=1.3
            xyouts,0.30,0.91,'File',/normal,charsize=1.3
            loadct,27
            for j=1,n0-1 do begin
                oplot,bharr[*,0,j],bharr[*,2,j],psym=-1,color=colors[j-1]
                xyouts,0.3,0.87-(j*0.03),files[j],color=colors[j-1],/normal,charsize=1.3
                xyouts,0.45,0.87-(j*0.03),strn(bharr[0,1,j]),/normal,charsize=1.3,color=colors[j-1]
            endfor
            print,'Enter "y" for another y-axis iteration'
            read,ans
        endrep until (ans ne 'y')
    endif
    
    set_plot,'ps'
    device,file='X2_'+sorton+'_bh.ps',/color
    
    loadct,0
    plot,bharr[*,0,0],bharr[*,2,0],psym=-1,/ynozero,xrange=[xdn,xup],/xstyle,$
      title='DM Halo Black Hole Mass Chi-squared Values ('+sorton+')',$
      xtitle='BH Mass (e9 M_solar)',ytitle='scaled chi-squared',charsize=1.0,$
      xthick=3,ythick=3,charthick=2,thick=3,yrange=[ydn,yup],/ystyle
    xyouts,0.3,0.84,files[0],/normal,charsize=1.1,charthick=3
    xyouts,0.45,0.84,strn(bharr[0,1,0]),/normal,charsize=1.1,charthick=3
    xyouts,0.45,0.88,'Min chi2 Value',/normal,charsize=1.1,charthick=3
    xyouts,0.30,0.88,'File',/normal,charsize=1.1,charthick=3
    
    loadct,4
    for j=1,n0-1 do begin
        oplot,bharr[*,0,j],bharr[*,2,j],psym=-1,color=colors[j-1],thick=3
        xyouts,0.3,0.84-(j*0.03),files[j],color=colors[j-1],/normal,charsize=1.1,charthick=3
        xyouts,0.45,0.84-(j*0.03),strn(bharr[0,1,j]),/normal,charsize=1.1,color=colors[j-1],charthick=3
    endfor
    device,/close_file
   pause
endif


if (tvc gt 1) then begin
    set_plot,'x'
    window,2,retain=2
    device,decomposed=0
    xup = max(vcarr[*,0,*])+25 & xdn = min(vcarr[*,0,*])-25
    yup = max(vcarr[*,2,*])+25 & ydn = min(vcarr[*,2,*])-25
    loadct,0
    
    plot,vcarr[*,0,0],vcarr[*,2,0],psym=-1,/ynozero,xrange=[xdn,xup],/xstyle,$
      title='DM Halo Circular Velocity Chi-squared Values ('+sorton+')',$
      xtitle='Velocity (km/sec)',ytitle='scaled chi-squared',charsize=1.2,$
      yrange=[ydn,yup],/ystyle
    xyouts,0.3,0.87,files[0],/normal,charsize=1.3
    xyouts,0.45,0.87,strn(vcarr[0,1,0]),/normal,charsize=1.3
    xyouts,0.45,0.91,'Min chi2 Value',/normal,charsize=1.3
    xyouts,0.30,0.91,'File',/normal,charsize=1.3
    
    loadct,27
    for j=1,n0-1 do begin
        oplot,vcarr[*,0,j],vcarr[*,2,j],psym=-1,color=colors[j-1]
        xyouts,0.3,0.87-(j*0.03),files[j],color=colors[j-1],/normal,charsize=1.3
        xyouts,0.45,0.87-(j*0.03),strn(vcarr[0,1,j]),/normal,charsize=1.3,color=colors[j-1]
    endfor
    
    ans = ''
    print,'Change the y-axis plotting range? (y/n)'
    read,ans
    if (ans eq 'y') then begin
        repeat begin
            yup = intarr(1)
            ydn = intarr(1)
            print,'Enter an upper Y range:'
            read,yup
            print,'Enter a lower Y range:'
            read,ydn
            loadct,0
            plot,vcarr[*,0,0],vcarr[*,2,0],psym=-1,/ynozero,xrange=[xdn,xup],/xstyle,$
              title='DM Halo Circular Velocity Chi-squared Values ('+sorton+')',$
              xtitle='Velocity (km/sec)',ytitle='scaled chi-squared',yrange=[ydn,yup],/ystyle,charsize=1.2
            xyouts,0.3,0.87,files[0],/normal,charsize=1.0
            xyouts,0.45,0.87,strn(vcarr[0,1,0]),/normal,charsize=1.3
            xyouts,0.45,0.91,'Min chi2 Value',/normal,charsize=1.3
            xyouts,0.30,0.91,'File',/normal,charsize=1.3
            loadct,27
            for j=1,n0-1 do begin
                oplot,vcarr[*,0,j],vcarr[*,2,j],psym=-1,color=colors[j-1]
                xyouts,0.3,0.87-(j*0.03),files[j],color=colors[j-1],/normal,charsize=1.3
                xyouts,0.45,0.87-(j*0.03),strn(vcarr[0,1,j]),/normal,charsize=1.3,color=colors[j-1]
            endfor
            print,'Enter "y" for another y-axis iteration'
            read,ans
        endrep until (ans ne 'y')
    endif
    
    set_plot,'ps'
    device,file='X2_'+sorton+'_vc.ps',/color
    
    loadct,0
    plot,vcarr[*,0,0],vcarr[*,2,0],psym=-1,/ynozero,xrange=[xdn,xup],/xstyle,$
      title='DM Halo Circular Velocity Chi-squared Values ('+sorton+')',$
      xtitle='Velocity (km/sec)',ytitle='scaled chi-squared',charsize=1.0,$
      xthick=3,ythick=3,charthick=2,thick=3,yrange=[ydn,yup],/ystyle
    xyouts,0.3,0.84,files[0],/normal,charsize=1.1,charthick=3
    xyouts,0.45,0.84,strn(vcarr[0,1,0]),/normal,charsize=1.1,charthick=3
    xyouts,0.45,0.88,'Min chi2 Value',/normal,charsize=1.1,charthick=3
    xyouts,0.30,0.88,'File',/normal,charsize=1.1,charthick=3
    
    loadct,4
    for j=1,n0-1 do begin
        oplot,vcarr[*,0,j],vcarr[*,2,j],psym=-1,color=colors[j-1],thick=3
        xyouts,0.3,0.84-(j*0.03),files[j],color=colors[j-1],/normal,charsize=1.1,charthick=3
        xyouts,0.45,0.84-(j*0.03),strn(vcarr[0,1,j]),/normal,charsize=1.1,color=colors[j-1],charthick=3
    endfor
    device,/close_file
    pause
endif


if (trs gt 1) then begin
    set_plot,'x'
    window,3,retain=2
    device,decomposed=0
    xup = max(rsarr[*,0,*])+1 & xdn = min(rsarr[*,0,*])-1
    yup = max(rsarr[*,2,*])+25 & ydn = min(rsarr[*,2,*])-25
    loadct,0
    
    plot,rsarr[*,0,0],rsarr[*,2,0],psym=-1,/ynozero,xrange=[xdn,xup],/xstyle,$
      title='DM Halo Scale Radius Chi-squared Values ('+sorton+')',$
      xtitle='Radius',ytitle='scaled chi-squared',charsize=1.2,$
      yrange=[ydn,yup],/ystyle
    xyouts,0.3,0.87,files[0],/normal,charsize=1.3
    xyouts,0.45,0.87,strn(rsarr[0,1,0]),/normal,charsize=1.3
    xyouts,0.45,0.91,'Min chi2 Value',/normal,charsize=1.3
    xyouts,0.30,0.91,'File',/normal,charsize=1.3
    
    loadct,27
    for j=1,n0-1 do begin
        oplot,rsarr[*,0,j],rsarr[*,2,j],psym=-1,color=colors[j-1]
        xyouts,0.3,0.87-(j*0.03),files[j],color=colors[j-1],/normal,charsize=1.3
        xyouts,0.45,0.87-(j*0.03),strn(rsarr[0,1,j]),/normal,charsize=1.3,color=colors[j-1]
    endfor
    
    ans = ''
    print,'Change the y-axis plotting range? (y/n)'
    read,ans
    if (ans eq 'y') then begin
        repeat begin
            yup = intarr(1)
            ydn = intarr(1)
            print,'Enter an upper Y range:'
            read,yup
            print,'Enter a lower Y range:'
            read,ydn
            loadct,0
            plot,rsarr[*,0,0],rsarr[*,2,0],psym=-1,/ynozero,xrange=[xdn,xup],/xstyle,$
              title='DM Halo Scale Radius Chi-squared Values ('+sorton+')',$
              xtitle='Radius',ytitle='scaled chi-squared',yrange=[ydn,yup],/ystyle,charsize=1.2
            xyouts,0.3,0.87,files[0],/normal,charsize=1.0
            xyouts,0.45,0.87,strn(rsarr[0,1,0]),/normal,charsize=1.3
            xyouts,0.45,0.91,'Min chi2 Value',/normal,charsize=1.3
            xyouts,0.30,0.91,'File',/normal,charsize=1.3
            loadct,27
            for j=1,n0-1 do begin
                oplot,rsarr[*,0,j],rsarr[*,2,j],psym=-1,color=colors[j-1]
                xyouts,0.3,0.87-(j*0.03),files[j],color=colors[j-1],/normal,charsize=1.3
                xyouts,0.45,0.87-(j*0.03),strn(rsarr[0,1,j]),/normal,charsize=1.3,color=colors[j-1]
            endfor
            print,'Enter "y" for another y-axis iteration'
            read,ans
        endrep until (ans ne 'y')
        pause
    endif
    
    set_plot,'ps'
    device,file='X2_'+sorton+'_rs.ps',/color
    
    loadct,0
    plot,rsarr[*,0,0],rsarr[*,2,0],psym=-1,/ynozero,xrange=[xdn,xup],/xstyle,$
      title='DM Halo Scale Radius Chi-squared Values ('+sorton+')',$
      xtitle='Radius',ytitle='scaled chi-squared',charsize=1.0,$
      xthick=3,ythick=3,charthick=2,thick=3,yrange=[ydn,yup],/ystyle
    xyouts,0.3,0.84,files[0],/normal,charsize=1.1,charthick=3
    xyouts,0.45,0.84,strn(rsarr[0,1,0]),/normal,charsize=1.1,charthick=3
    xyouts,0.45,0.88,'Min chi2 Value',/normal,charsize=1.1,charthick=3
    xyouts,0.30,0.88,'File',/normal,charsize=1.1,charthick=3
    
    loadct,4
    for j=1,n0-1 do begin
        oplot,rsarr[*,0,j],rsarr[*,2,j],psym=-1,color=colors[j-1],thick=3
        xyouts,0.3,0.84-(j*0.03),files[j],color=colors[j-1],/normal,charsize=1.1,charthick=3
        xyouts,0.45,0.84-(j*0.03),strn(rsarr[0,1,j]),/normal,charsize=1.1,color=colors[j-1],charthick=3
    endfor
    device,/close_file
endif
set_plot,'x'
print,'An ENTER deletes the windows...'
pause
wdelete & wdelete & wdelete & wdelete
stop
END

FUNCTION biweightEF, array, plotF=plotF

; The difference here is that this routine putputs two values for each
; position. The second is now the uncertainty as defined by the Beers
; paper.

compile_opt idl2
on_error,2

;This code has been simplifed from the original biweight which
;conducted the normalization WITHIN the routine. Now, the data sent in
;is not manipulated within the program. Thus, if any normalization is
;required, it must happen before the data is sent in.

;PLOTBWF: Set this to anything other than 'on' to supress the plotting
;routine. If not used, the default is to not plot.

;The format of the biweight was taken from the Beers 1990 paper.

;THE CONCENTRATION PARAMETER****************************************
;The HIGHER the number the WEAKER the rejection
c = 6.0 ;The concentration (tension) parameter.
cE = 9.0 ;The concentration for the uncertainty parameter
;plotF = 'on'
speed = 10e-6 ;10e-4 is very fast, 10e-6 is slow, and slightly more 
              ;robust for arrays with small numbers
;*******************************************************************

print,'Running biweightEF...'

if (n_elements(plotF) eq 0) then plotF = 'off'
;The data array is a MxN array, where the biweight is run over EACH
;COLUMN.

nn1 = n_elements(array[*,0])
nn2 = n_elements(array[0,*])
biwt = dblarr(nn1,2) ;the new array, including error

;MODIFIED ON DEC 27, 2010: The old method (which sucked) has been
;superceded as the median function happily steps around the NaN flag
ibadF = where(array le -666,countF)
;print,countF
if (countF gt 0) then array[ibadF] = !values.F_NAN

for jj=0,nn1-1 do begin ;a loop through wavelength
;    print,'Biweight on X-pixel '+strn(jj+1)
    temp = array[jj,*]
    pp = finite(temp)
    ii = where(pp eq 1,countF);the non-Nan values are indexed
;    print,'The number of good values is '+strn(countF)

    if (countF gt 1) then begin ;if there's more than 1 # to work on
        temp = temp[ii]
;        plot,temp,title=strn(jj)
        nn3 = n_elements(temp)
        nn3E = float(nn3)
        M = median(temp,/even,/double)
        goto,Fjump1
    endif
    if (countF eq 1) then begin ;if there's only 1 #, the routine returns that #
        biwt[jj,0] = temp[ii]
        biwt[jj,1] = 0.0
        goto,Fjumpend
    endif
    if (countF eq 0) then begin ;if no real #'s exist, you get a -666 returned
        biwt[jj,0] = -666
        biwt[jj,1] = -666
        goto,Fjumpend
    endif

    Fjump1:
    biarr1 = dblarr(nn3)
    biarr2 = dblarr(nn3)
    biarr1E = dblarr(nn3)
    biarr2E = dblarr(nn3)
    madarr = dblarr(nn3)
    ui = dblarr(nn3)
    uiE = dblarr(nn3)
    cntrF = 0
    delta = 1.0

;------------------------------------------------------------------
;THE BIWEIGHT
;------------------------------------------------------------------
    while (abs(delta) gt speed and cntrF lt 20) do begin
        madarr = abs(temp - M) ;the MAD array, BEFORE the median is taken
        MAD = median(madarr,/even,/double)
        if (MAD eq 0.0) then MAD = 0.0000000001

        kk = 0.0
        repeat begin
            ui[kk] = (temp[kk]-M)/(c*MAD)
            uiE[kk] = (temp[kk]-M)/(cE*MAD)
            t = abs(ui[kk])
            te = abs((uie[kk]))
            if (t le 1.0) then begin ;if the point is included
                biarr2[kk] = (1 - (ui[kk]^2))^2
                biarr1[kk] = (temp[kk]-M) * biarr2[kk]
                biarr1E[kk] = ((temp[kk]-M)^2) * ((1-(uiE[kk])^2)^4)
                biarr2E[kk] = (1-(uiE[kk]^2)) * (1-(5*(uiE[kk]^2)))
            endif else begin    ;if the point is rejected
;                print,'Pixel #'+strn(kk+1)+' at column #'+strn(jj+1)+' was rejected!' 
                biarr1[kk]  = 0.0
                biarr2[kk]  = 0.0
                biarr1E[kk] = 0.0
                biarr2E[kk] = 0.0
            endelse
            kk = kk + 1.0
        endrep until (kk eq nn3)
        
        t1 = total(biarr1)
        t2 = total(biarr2)
        t1E = sqrt(total(biarr1E)) * sqrt(nn3E)
        t2E = abs(total(biarr2E))

        if (t2 eq 0.0) then begin
            biwt[jj,0] = -666
            biwt[jj,1] = -666
            goto,Fjumpend
        endif else delta = t1/t2
        M = M + delta
        cntrF = cntrF + 1
    endwhile
    
    biwt[jj,0] = M
    biwt[jj,1] = double(t1E/t2E)

;------------------------------------------------------------------
;THE BIWEIGHT
;------------------------------------------------------------------
Fjumpend:
endfor

if (plotF eq 'on') then begin

    xup=nn1-floor(nn1*0.1)
    xdown=50
    yup = max(biwt[xdown:xup]) + (0.1*max(biwt[xdown:xup]))
    ydown = min(biwt[xdown:xup]) - (0.1*min(biwt[xdown:xup]))

    step = floor(255.0/nn2)
    colors = intarr(nn2)
    if (step eq 0) then step = 1.0
    for jj=0,nn2-1 do colors[jj] = step*jj

    set_plot,'x'
    window,/free,retain=2
    device,decomposed=0
    loadct,0
    plot,array[*,0],/nodata,yrange=[ydown,yup],xrange=[xdown,xup],$
      xstyle=1,title='The biweight value is plotted in white',ystyle=1
    loadct,27
    for jj=0,nn2-1 do begin
        oplot,array[*,jj],color=colors[jj]
    endfor
    loadct,0
    oplot,biwt

endif

RETURN,biwt

end

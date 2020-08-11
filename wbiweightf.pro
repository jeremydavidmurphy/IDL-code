FUNCTION wbiweightf, array, weights
ON_ERROR,2

;This code has been simplifed from the original biweight which
;conducted the normalization WITHIN the routine. Now, the data sent in
;is not manipulated within the program. Thus, if any normalization is
;required, it must happen before the data is sent in.

;The format of the biweight was taken from the Beers 1990 paper. The
;one significant modification is that now an external weight is
;included that combines with the inherent weight from the
;biweight. This allows for inherently higher quality data to receive
;more weight in the final spectrum.

;THE CONCENTRATION PARAMETER****************************************
;The HIGHER the number the WEAKER the rejection
c=6.0 ;The concentration (tension) parameter.
;*******************************************************************

;The data array is a MxN array, where the biweight is run over EACH
;COLUMN. The weight array must be of the same size and be composed of
;values between 0 and 1 in order not to overwhelm the natural
;weighting of the biweight.

if (n_elements(array) ne n_elements(weights)) then begin
    print,'TRY AGAIN!!!'
    print,'The size of the weights array does not match the data array!'
    stop
endif

n1 = n_elements(array[*,0])
n2 = n_elements(array[0,*])
biwt = dblarr(n1)

;---------------------------------------------------------------
;Plotting limits, in pixels
xup=n1-floor(n1*0.1)
xdown=50
yup = max(array[xdown:xup,*])
ydown = min(array[xdown:xup,*])
;---------------------------------------------------------------
cntr = 0
for j=0,n1-1 do begin ;a loop through wavelength
    ind = where(array[j,*] eq -666)
    n3 = n_elements(ind)

;-----------------------------------------------------------------
;THE -666 FLAGS ARE DEALT WITH
;-----------------------------------------------------------------
;Here, TEMP is the array of values the biweight will ultimately run
;on. This is all the values, for a given wavelength, that are being
;combined.

    if (n3 eq n2) then begin    ;if all data carries the -666 flag
        biwt[j] = -666
        goto,jumpend
    endif
    if (n3 eq 1) then begin
        if (ind eq -1) then begin ;if no -666 flags were found
            temp = array[j,*]
            wgts = weights[j,*]
            M = median(temp,/even)
            goto, jump1
        endif else begin        ;if one -666 flag is found
            temp = dblarr(n2-n3)
            goodi = where(array[j,*] ne -666)
            temp = array[j,goodi]
            wgts = weights[j,goodi]
            M = median(temp,/even)
            goto,jump1
        endelse
    endif else begin            ;if more than one -666 flag is found
        temp = dblarr(n2-n3)
        goodi = where(array[j,*] ne -666)
        temp = array[j,goodi]
        wgts = weights[j,goodi]
        M = median(temp,/even)
        goto,jump1
    endelse

;********************************************************************* 
;At this point, all the median values (M) are set.
;The 'temp' array are the remaining data points, and the size of this
;array varies
    
jump1:
;------------------------------------------------------------------
;for the oddball chance that there's only 1 surviving data point
;(i.e. all the rest were masked in VACCINE) then a -666 flag is returned
    n4 = n_elements(temp)
    if (n4 eq 1) then begin
        biwt[j] = -666
        goto, jumpend
    endif
;------------------------------------------------------------------
    biarr1 = dblarr(n4)
    biarr2 = dblarr(n4)
    madarr = dblarr(n4)
    ui = dblarr(n4)
    cntr2 = 0
    delta = 1.0
;--------------------------------------------------------------------------------------
;THE WEIGHTED BIWEIGHT
;--------------------------------------------------------------------------------------
    while (abs(delta) gt 10e-6 and cntr2 lt 20) do begin
        for k=0,n4-1 do begin
            madarr[k] = abs(temp[k]-M) ;the MAD array, BEFORE the median is taken
        endfor
        MAD = median(madarr,/even)
        
        for k=0,n4-1 do begin ;a loop through each element in the TEMP array
            ui[k] = (temp[k]-M)/(c*MAD)
            t = abs(ui[k])
            if (t le 1.0) then begin ;if the point receives a weight
                biarr2[k] = (1-ui[k]^2)^2+(wgts[k]^2)
                biarr1[k] = (temp[k]-M)*biarr2[k]
            endif else begin    ;if the point is rejected
                cntr=cntr+1
                biarr2[k] = 0.0
                biarr1[k] = 0.0
            endelse
        endfor
        
        t0 = total(biarr1)
        t1 = total(biarr2)
        delta = double(t0/t1)
;        if (t1 eq 0) then delta = 0.0 else delta = t0/t1
        M = M + delta
        cntr2 = cntr2 + 1

    endwhile
    
    biwt[j] = M
;    print,'The number of biweight iterations was '+strn(cntr2)
;    print,'Delta equals '+strn(delta)
;    if (t1 eq 0) then biwt[j] = -666 else biwt[j] = M

jumpend:
endfor
;print,'The number of rejected data points was '+strn(cntr)


;ans='n'
;print,'Plot the results to screen? (y or n):'
;read,ans

;if (ans eq 'y') then begin
;    step = floor(255/n2)
;    color = intarr(n2)
;    for j=0,-1 do color[j] = step*j
    
;    window,0,retain=2
;    device,decomposed=0
;    loadct,0
;    plot,array[*,0],/nodata,yrange=[ydown,yup],xrange=[xdown,xup],xstyle=1
;    loadct,27
;    for j=0,n2-1 do begin
;        oplot,array[*,j],color=color[j]
;    endfor
;    loadct,0
;    pause
;    oplot,biwt,thick=2

;    print,'The next ENTER deletes the plot:'
;    pause
;    wdelete,0
;endif

RETURN,biwt

end

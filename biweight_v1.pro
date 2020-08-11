; THIS VERSION OF THE BIWEIGHT ALLOWS FOR THE NORMALIZATION OF THE
; INCOMING DATA TO BE MADE WITHIN THE ROUTINE. THE NEWER ROUTINE, USED
; IN PIPE2.PRO, 


FUNCTION biweight_v1, array, NORM=normvalues
ON_ERROR,2



;NOTE: The normvalues get DIVIDED by the data array, NOT
;MULTIPLIED. So, they are better thought of as weights. They are used
;to normalize the data via the equation,

;          normed = array / normvalues

;Normvalues must be an array with the same number of values as there
;are n-elements in the data array. This modification from the original
;biweight code allows the raw data (rather than normalized data) to be
;fed into the biweight, avoiding the confusion with the -666 flags.

;If the normvalues call isn't activated, then the straight biweight is
;run.

;THE CONCENTRATION PARAMETER****************************************
;The HIGHER the number the WEAKER the rejection
c=6.0 ;The concentration (tension) parameter.
;*******************************************************************

if (n_elements(normvalues) ne 0) then begin
    if (n_elements(array[0,*]) ne n_elements(normvalues)) then begin
        print,'TRY AGAIN!!!'
        print,'The size of the normvalues array does not match the data array!'
        stop
    endif
    test = 'n'
endif else test = 'a'

n1 = n_elements(array[*,0])
n2 = n_elements(array[0,*])
biwt = dblarr(n1)


if (test eq 'n') then begin
    normed = dblarr(n1,n2)
    for k=0,n2-1 do begin
        normed[*,k] = array[*,k] / normvalues[k]
    endfor
endif

cntr=0

for j=0,n1-1 do begin ;a loop through wavelength
    ind = where(array[j,*] eq -666)
    n3 = n_elements(ind)

;-----------------------------------------------------------------
;THE -666 FLAGS ARE DEALT WITH
;-----------------------------------------------------------------
;Here, TEMP is the array of values the biweight will ultimately run
;on. This is all the values, for a given wavelength, that are being
;combined.

;THE NORMALIZED ARRAY
    if(test eq 'n') then begin
        if (n3 eq n2) then begin ;if all data carries the -666 flag
            biwt[j] = -666
            goto,jumpend
        endif
        if (n3 eq 1) then begin
            if (ind eq -1) then begin ;if no -666 flags were found
                temp = normed[j,*]
                M = median(temp,/even)
                goto, jump1
            endif else begin    ;if one -666 flag is found
                temp = dblarr(n2-n3)
                cntr1=0
                for l=0,n2-1 do begin ;the truncated data array (w/o -666 flags) is created
                    if (array[j,l] ne -666) then begin
                        temp[cntr1] = normed[j,l]
                        cntr1=cntr1+1
                    endif
                endfor
                M = median(temp,/even)
                goto,jump1
            endelse
        endif else begin ;if more than one -666 flag is found
            temp = dblarr(n2-n3)
            cntr1=0
            for l=0,n2-1 do begin
                if (array[j,l] ne -666) then begin
                    temp[cntr1] = normed[j,l]
                    cntr1=cntr1+1
                endif
            endfor
            M = median(temp,/even)
            goto,jump1
        endelse
    endif

;THE REGULAR ARRAY. THIS RUNS IF NO NORMALIZATION IS INCLUDED IN THE
;FUNCTION CALL    
    if(test eq 'a') then begin
        if (n3 eq n2) then begin ;if all data carries the -666 flag
            biwt[j] = -666
            goto,jumpend
        endif
        if (n3 eq 1) then begin
            if (ind eq -1) then begin ;if no -666 flags were found
                temp = array[j,*]
                M = median(temp,/even)
                goto, jump1
            endif else begin    ;if one -666 flag is found
                temp = dblarr(n2-n3)
                cntr1=0
                for l=0,n2-1 do begin ;the truncated data array (w/o -666 flags) is created
                    if (array[j,l] ne -666) then begin
                        temp[cntr1] = array[j,l]
                        cntr1=cntr1+1
                    endif
                endfor
                M = median(temp,/even)
                goto,jump1
            endelse
        endif else begin        ;if more than one -666 flag is found
            temp = dblarr(n2-n3)
            cntr1=0
            for l=0,n2-1 do begin
                if (array[j,l] ne -666) then begin
                    temp[cntr1] = array[j,l]
                    cntr1=cntr1+1
                endif
            endfor
            M = median(temp,/even)
            goto,jump1
        endelse
    endif
    
;********************************************************************* 
;At this point, all the median values (M) are set, coming from either
;the data array or normed array, depending upon whether normvalues
;were included in the function call.
;The 'temp' array are the remaining data points, and the size of this
;array varies
    
jump1:
;------------------------------------------------------------------
;for the oddball chance that there's only 1 surviving data point
;(i.e. all the rest were masked in VACCINE) then a -666 flag is returned
    n4 = n_elements(temp)
    if (n4 eq 1) then begin
        ind = where(array[j,*] ne -666)
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
;THE BIWEIGHT
;--------------------------------------------------------------------------------------
    while (abs(delta) gt 10e-6 and cntr2 lt 20) do begin
        for k=0,n4-1 do begin
            madarr[k] = abs(M-temp[k]) ;the MAD array, BEFORE the median is taken
        endfor
        MAD = median(madarr,/even)
        
        for k=0,n4-1 do begin   ; a loop through the number of spectra
            ui[k] = (temp[k]-M)/(c*MAD)
            t = abs(ui[k])
            if (t le 1.0) then begin ;If the point receives a weight
                biarr2[k] = (1-ui[k]^2)^2
                biarr1[k] = (temp[k]-M)*biarr2[k]
            endif else begin    ;If the point is rejected
                cntr=cntr+1
                biarr2[k] = 0.0
                biarr1[k] = 0.0
            endelse
        endfor
        
        t0 = total(biarr1)
        t1 = total(biarr2)
        if (t1 eq 0) then delta = 0.0 else delta = t0/t1
        M = M + delta
        cntr2 = cntr2 + 1

    endwhile
;    print,'The number of biweight iterations was '+strn(cntr2)
;    print,'Delta equals '+strn(delta)
    if (t1 eq 0) then biwt[j] = -666 else biwt[j] = M

jumpend:
endfor
print,'The number of rejected data points was '+strn(cntr)


;print,'Plot the results to screen? (y or n):'
;read,ans
ans='n'
if (ans eq 'y') then begin
    if (n2 le 10) then begin
        color = [60,80,100,120,140,160,180,200,220,240]
    endif else begin
        step = floor(195/n2)
        color = intarr(n2)
        for j=0,-1 do color[j] = 60+(step*j)
    endelse

    window,0,retain=2
    device,decomposed=0
    if(test eq 'n') then begin
        loadct,0
        plot,normed[*,0],/nodata,yrange=[0,40],xrange=[200,1500]
        loadct,4
        for j=0,n2-1 do begin
            oplot,normed[*,j],color=60+(20*j)
        endfor
        loadct,0
        oplot,biwt
    endif
    if(test eq 'a') then begin
        high = max(array[1600:2000,*])
        loadct,0
;        plot,array[*,0],yrange=[0,high],/nodata
        plot,array[*,0],yrange=[0,40],/nodata,xrange=[200,1500]
        loadct,4
        for j=0,n2-1 do begin
            oplot,array[*,j],color=60+(20*j)
        endfor
        loadct,0
        oplot,biwt
    endif
;    print,'The next ENTER deletes the plot:'
;    pause
;    wdelete,0
endif

RETURN,biwt

end


;The new biweight has done away with the dori junk

FUNCTION biweight_dori, array, NORM=normvalues, DATAorINDEX=dori
ON_ERROR,2

;NOTE: The normvalues get DIVIDED by the data array, NOT
;MULTIPLIED. So, they are better thought of as weights. They are used
;to normalize the data via the equation,

;          normed = array / normvalues

;This code accepts an array (mxn) and returns EITHER THE BIWEIGHT OR
;THE INDEX OF REJECTED VALUES length m. (As this is rarely used, it's
;been relegated to a keyword call.)

;The 'dori' needs to be either "data" or "index". When set for index,
;the array returned indicates which values were rejected: -1 means the
;data point was kept while another number is the index of the value
;rejected.

;Normvalues must be an array with the same number of values as there
;are n-elements in the data array. This modification from the original
;biweight code allows the raw data (rather than normalized data) to be
;fed into the biweight, avoiding the confusion with the -666 flags.

;If the normvalues call isn't activated, then the straight biweight is
;run.

;THE CONCENTRATION PARAMETER****************************************
;The HIGHER the number the WEAKER the rejection
c=6.0 ;The concentration(tension) parameter.
;*******************************************************************

if (n_elements(dori) ne 0) then begin
    if (dori ne 'index') and (dori ne 'data') then begin
        print,'TRY AGAIN!!!'
        print,'The calling syntax is: out = biweight(array, "data" or "index")'
        stop
    endif
    test2 = 'i'
endif else test2 = 'd'

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
indexarray = intarr(n1,n2)
biwt = dblarr(n1)
cntr=0


if (test eq 'n') then begin
    normed = dblarr(n1,n2)
    for k=0,n2-1 do begin
        normed[*,k] = array[*,k] / normvalues[k]
    endfor
endif

for k=0,n1-1 do begin ;a loop through wavelength
    ind = where(array[k,*] eq -666)
    n3 = n_elements(ind)

;THE -666 FLAGS ARE DEALT WITH
;*********************************************************************    
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
                temp = dblarr(n2-1)
                cntr1=0
                for l=0,n2-1 do begin ;the truncated data array (w/o -666 flags) is created
                    if (array[j,l] ne -666) then begin
                        temp[cntr1] = normed[j,l]
                        cntr1=cntr1+1
                    endif else indexarray[j,l] = -666
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
                endif else indexarray[j,l] = -666
            endfor
            M = median(temp,/even)
            goto,jump1
        endelse
    endif
    
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
                temp = dblarr(n2-1)
                cntr1=0
                for l=0,n2-1 do begin ;the truncated data array (w/o -666 flags) is created
                    if (array[j,l] ne -666) then begin
                        temp[cntr1] = array[j,l]
                        cntr1=cntr1+1
                    endif else indexarray[j,l] = -666
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
                endif else indexarray[j,l] = -666
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

;for the oddball chance that there's only 1 surviving data point
;(i.e. all the rest were masked in VACCINE) then a -666 flag is returned
    n4 = n_elements(temp)
    if (n4 eq 1) then begin
        ind = where(array[j,*] ne -666)
        biwt[j] = -666
        goto, jumpend
    endif

    biarr1 = dblarr(n4)
    biarr2 = dblarr(n4)
    madarr = dblarr(n4)
    ui = dblarr(n4)
    cntr2 = 0
    delta = 1.0

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
                indexarray[j,k] = -1
            endif else begin    ;If the point is rejected
                cntr=cntr+1
                biarr2[k] = 0.0
                biarr1[k] = 0.0
                indexarray[j,k] = j
            endelse
        endfor
        
        t0 = total(biarr1)
        t1 = total(biarr2)
        if (t1 eq 0) then delta = 0.0 else begin
            delta = t0/t1
            M = M + delta
            cntr2 = cntr2 + 1
        endelse

    endwhile
    print,'The number of biweight iterations was '+strn(cntr2)
;    print,'Delta equals '+strn(delta)
    if (t1 eq 0) then biwt[j] = -666 else biwt[j] = M

jumpend:
endfor
print,'The number of data points rejected by the biweight was :'+strn(cntr)

;biwt[where(biwt eq -1)] = -666
;biwt[where(biwt lt -50)] = -666

ans='n'
;print,'Plot the results to screen? (y or n):'
;read,ans
if (ans eq 'y') then begin
    window,0,retain=2
    device,decomposed=0
    if(test eq 'n') then begin
        loadct,0
        plot,normed[*,0],/nodata
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
        plot,array[*,0],yrange=[0,high],/nodata
        loadct,4
        for j=0,n2-1 do begin
            oplot,array[*,j],color=60+(20*j)
        endfor
        loadct,0
        oplot,biwt
    endif
    print,'The next ENTER deletes the plot:'
    pause
    wdelete,0
endif

if (test2 eq 'd') then begin
;    print,'The actual data array will be returned.'
    RETURN,biwt
endif

if (test2 eq 'i') then begin
;    print,'The index of rejected points will be returned.'
;    print,'(-1 entries indicate good data)'
    RETURN,indexarray
endif

end


;The only change from version 1 is that now this routine accepts a
;list of lists of files so it doesn't need to be babysat when reducing
;several frames...also changed is the introduction of a second list of
;the .R files. This allows you to mask the stars and other
;interstellar junk. 

;THIS ROUTINE IS USED TO COMBINE VIRUS-P DATA OF SEVERAL EXPOSURES AT
;THE SAME POINTING. IT READS IN A LIST OF FILE NAMES, THEN APPLIES THE
;BIWEIGHT (A SWITCH TO THE WEIGHTED BIWEIGHT IS BEING DEVELOPED) TO
;THE DATA FOR EACH FIBER. THE DATA IS NORMALIZED TO ACCOUNT FOR THE
;DIFFERENT LEVELS OF OVERALL INTENSITY BETWEEN EXPOSURES.

PRO combiwgt_v2, listolists, WEIGHT=wgt, Mask=dotR

;The 'list0lists' is a list of names of lists of files.
;ex: m87b1_aa.list
;    m87b2_aa.list
;    etc.
;with each list containing the lists of frames to have combined.
;The output fits files are then named the same name as the list, only
;with '.list' replaced with '.fits'

;This is a change from combine_v2.pro. The difference is that now
;instead of the normalized values being fed into the biweight, the raw
;data is fed into the weighted biweight, along with weights either
;determined by this routine, or fed into combine_v3 as an array, 'wgt'

;IF THE WEIGHT ISN'T FED DIRECTLY IN AS AN ARRAY, THIS ROUTINE
;DETERMINES THE WEIGHT BY TAKING THE MEDIAN OF A REGION OF 300 PIXELS
;IN THE RED REGION OF THE SPECTRUM IN QUESTION.

;This routine will look for a list named 'combine.list' of the data
;files to be combined.

readcol,listolists,format='A',lists

n0 = n_elements(lists)
for m=0,n0-1 do begin ;the main loop THROUGH EACH LIST is set up
    readcol,lists[m],format='A',files
    test = readfits(files[0])
    n1 = n_elements(test[*,0])  ;wavelength
    n2 = n_elements(test[0,*])  ;# of fibers
    n3 = n_elements(files)      ;# of files

    data = dblarr(n1,n2,n3)
    out = dblarr(n1,n2)
    slice = dblarr(n1,n3)
    wgtarray = dblarr(n3,n2)

;THIS IS GETTING SHUT OFF FOR NOW = A WEIGHT FILE CAN NOT BE FED IN...
;if (n_elements(wgt) ne 0) then swch = 'off'
;if (n_elements(wgt) eq 0) then begin
    swch = 'on'
    wgt = dblarr(n3)
    nwgt = dblarr(n3)
;endif

;the data is read in
    for k=0,n3-1 DO data[*,*,k] = readfits(files[k])
    tube = where(data[500,150,*] eq -1)
    n4 = n_elements(tube)
    if(n4 gt 1) then begin ;not all data frames in the list were found
        print,''
        print,'Data frames in the list '+lists[m]+' were not found!'
        print,'The missing frames are:'
        print,files[tube]
        stop
    endif else begin
        if (tube ne -1) then begin
            print,''
            print,'Data frames in the list '+lists[m]+' were not found!'
            print,'The missing frames are:'
            print,files[tube]
            stop
        endif
    endelse
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;this section can be used to just continue the reduction with only the
;frames found. IT'S NOT COMPLETE!
;        newfiles = strarr(n3-n_elements(tube))
;        cntr=0
;        print,cntr
;        pause
;        for p=0,n3-1 do begin
;            if (tube[cntr] ne p) then newfiles[p-cntr] = files[p] else cntr=cntr+1
;        endfor
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

    biwtout = dblarr(n1,n2)

    ans='n'
;print,'Plot the data as it comes out? (y or n)'
;print,'This shows plots of all the data being combined as well as the resulting output from the biweight.'
;read,ans

    for k=0,n2-1 do begin       ;This loops on each fiber
        print,'Working on fiber '+strn(k+1)
        for l=0,n3-1 do begin
            slice[*,l] = data[*,k,l] ;all the exposures for one fiber are defined
            if (swch eq 'on') then begin
                wgt[l] = median(slice[n1-400:n1-100,l],/even)
            endif
        endfor
        
;        if (swch eq 'off') then wgtarray[*,k] = wgt
        if (swch eq 'on') then begin
            for l=0,n3-1 do nwgt[l] = wgt[l]/median(wgt,/even)
            wgtarray[*,k] = nwgt
        endif else nwgt = wgt
        
        

        biwtout[*,k] = biweight(slice,NORM=nwgt)
        
        if(ans eq 'y') then begin
            window,0,retain=2
            device,decomposed=0
            loadct,0
            plot,slice[30:n1-30,0],title='Fiber #: '+strn(k+1)+' The real data- the biweight is plotted in red.',$
              xtitle='Pixel',ytitle='ADU',/nodata
            loadct,4
            for j=0,n3-1 do oplot,slice[30:n1-30,j],color=(60+j*10)
            oplot,biwtout[30:n1-30,k],color=150,thick=2
            
            window,2,retain=2
            device,decomposed=0
            loadct,0
            plot,slice[30:n1-30,0],title='Fiber #: '+strn(k+1)+' The scaled data- the biweight is plotted in red.',$
              xtitle='Pixel',ytitle='ADU',/nodata
            loadct,4
            for j=0,n3-1 do oplot,normed[30:n1-30,j],color=60+j*10
            oplot,biwtout[30:n1-30,k],color=150,thick=2
            pause
            print,'The next ENTER deletes the plots.'
            pause
            wdelete,0,2
        endif
        
    endfor

;    print,''
;    print,'The first file in the list is:'
;    print,files[0]
;    print,''
;    print,'Name the output frame (w/o .fits):'
;    read,ans
    n = strsplit(lists[m],'.',/extract)
    name1 = n[0]+'.fits'
    name2 = n[0]+'_w.txt'
    writefits,name1,biwtout
    openw,5,name2
    printf,5,wgtarray
    free_lun, 5

endfor

STOP
END


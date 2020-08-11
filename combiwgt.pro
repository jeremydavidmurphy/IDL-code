;THIS ROUTINE IS USED TO COMBINE VIRUS-P DATA OF SEVERAL EXPOSURES AT
;THE SAME POINTING. IT READS IN A LIST OF FILE NAMES, THEN APPLIES THE
;BIWEIGHT (A SWITCH TO THE WEIGHTED BIWEIGHT IS BEING DEVELOPED) TO
;THE DATA FOR EACH FIBER. THE DATA IS NORMALIZED TO ACCOUNT FOR THE
;DIFFERENT LEVELS OF OVERALL INTENSITY BETWEEN EXPOSURES.

PRO combiwgt, List=list, WEIGHT=wgt, Mask=dotR

;LIST is the list of files. The code will prompt you if this keyword
;isn't used.

;WEIGHT is a weight file, equal in length to the number of files being
;combined. If none is used, the code calculates its own weight based
;on the median value of 300 pixels in the red region of the spectra.

;MASK is the list of the .R files created in radplot.pro. This list
;needs to be of the same length as the list of files and is used to
;mask the data (with a -666 flag) for all fibers that shouldn't be
;used in the combination.

;This is a change from combine_v2.pro. The difference is that now
;instead of the normalized values being fed into the biweight, the raw
;data is fed into the weighted biweight, along with weights either
;determined by this routine, or fed into combine_v3 as an array, 'wgt'

if (n_elements(list) eq 0) then begin
    list=''
    print,'Enter the list name:'
    print,'(ENTER will look for "combine.list")' 
    read,list
    if (list eq '') then list='combine.list'
endif

readcol,list,Format='A',files

if (n_elements(dotR) ne 0) then begin
    print,'The .R files are being used to mask bad fibers.'
    readcol,dotR,Format='A',fibers
    swchdotR = 'on'
endif

;The size of the data set is determined
test = readfits(files[0])
n1 = n_elements(test[*,0]);wavelength
n2 = n_elements(test[0,*]);# of fibers
n3 = n_elements(files);# of files

data = dblarr(n1,n2,n3)
out = dblarr(n1,n2)

;the data is read in
cntr = 0
for k=0,n_elements(files)-1 do begin
    file = readfits(files[k])
    test = n_elements(file)
    if (test eq 1) then begin
        print,'File '+files[k]+' not found. It was skipped in the combination.'
        pause
        cntr = cntr + 1
        goto,jump1
    endif else data[*,*,k] = file
jump1:
endfor

;the data array is corrected for missing files.
if (cntr ne 0) then begin
    n3 = n3-cntr
    data = data[*,*,0:n3-1]
    print,'n3: '+strn(n3)
    pause
endif

slice = dblarr(n1,n3)
wgtarray = dblarr(n3,n2)

if (n_elements(wgt) ne 0) then swch = 'off'
if (n_elements(wgt) eq 0) then begin
    swch = 'on'
    wgt = dblarr(n3)
    nwgt = dblarr(n3)
endif
    
biwtout = dblarr(n1,n2)

;MASKING...
if (swchdotR eq 'on') then begin
    for l=0,n3-1 do begin
        readcol,fibers[l],format='I,X,X',gdfibs
        gdfibs = gdfibs[bsort(gdfibs)]-1
        n4 = n_elements(gdfibs)
        cntr = 0
        for k=0,n2-1 do begin
            if (gdfibs[cntr] eq k) then begin
                cntr = cntr + 1
                if (cntr eq n4) then goto, jump2
            endif else begin
                print,'Fiber # '+strn(k+1)+' was rejected from frame '+strn(l+1)
                data[*,k,l] = -666
            endelse
        endfor
jump2:
    endfor
endif

ans='n'
;print,'Plot the data as it comes out? (y or n)'
;print,'This shows plots of all the data being combined as well as the resulting output from the biweight.'
;read,ans

for k=0,n2-1 do begin ;This loops on each fiber
    print,'Working on fiber '+strn(k+1)
    for l=0,n3-1 do begin
        slice[*,l] = data[*,k,l] ;all the exposures for one fiber are defined
        if (swch eq 'on') then begin
            wgt[l] = median(slice[n1-400:n1-100,l],/even)
        endif
    endfor
    
    if (swch eq 'off') then wgtarray[*,k] = wgt
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

biwtout[where(biwtout lt -10)] = -666

print,''
print,'The first file in the list is:'
print,files[0]
jumpback:
print,''
print,'Name the output frame (w/o .fits):'
read,ans
if(ans eq '') then goto,jumpback
name1 = ans+'.fits'
name2 = ans+'_w.txt'
writefits,name1,biwtout
openw,5,name2
printf,5,wgtarray
free_lun, 5


END


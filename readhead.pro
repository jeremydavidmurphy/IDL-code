PRO readhead, list, ELE=headELE, MULTI=many

;This routine is used to locate when the Fracker wasn't running.

;headELE is the header element you are looking for

;If the keyword MULTI is used, then the routine will look for a list
;of keywords. So, feed the routine with a list of header words, and
;give the name of the list to the MULTI keyword

readcol,list,silent=1,f='a',files

n0 = n_elements(files)

free_lun,5
openw,5,'readhead.txt'

cntr = 0
repeat begin
    j = cntr
    file = files[j]
    h = headfits(file,/silent)
    if (n_elements(headELE) gt 0) then out1 = sxpar(h,headELE)
    if (n_elements(many) gt 0) then begin
        readcol,many,f='a',elements,/silent
        n1 = n_elements(elements)
        out1 = strarr(n1)
        for j=0,n1-1 do begin
            out1[j] = sxpar(h,elements[j])
        endfor
    endif
    print,file,' ',out1
    printf,5,file,' ',out1
cntr = cntr + 1
endrep until (cntr eq n0)

free_lun,5

stop
END

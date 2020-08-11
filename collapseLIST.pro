pro collapseLIST, list

; compile: collapseF

; This routine feeds collapseF with a series of data and weight
; files. They are collapsed, using the weight files, and written out
; as fits files. I give them a capital C (as opposed to a lowercase)
; so as to distinguish them from a Vaccine collapsed file.

;**********************************************************************
aperture = 5.0
;**********************************************************************

readcol,list,format='a,a,a,i',Dlist,Wlist,Mlist,Alist ;data,wgt,mask,aperture

n0 = n_elements(Dlist)

for j=0,n0-1 do begin ;a loop through each frame
    print,'Collapsing '+Dlist[j]
    data = readfits(Dlist[j],header,/silent)
    weight = readfits(Wlist[j],/silent)
    mask = Mlist[j]
    aperture = Alist[j]
    out = wgtcollapseF(data,weight,mask,aperture)
    t = strsplit(Dlist[j],'.',/extract)
    tt = t[0]
    outname = tt+'C.fits'
    writefits,outname,out,header
endfor

stop
END

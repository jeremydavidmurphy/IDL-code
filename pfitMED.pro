PRO pfitMED, list

; This routine is a modification of pfitave.pro and now works from a
; list. This list contains pfitlov.out files. The routine will read
; these in, in turn, and return median values for the 4 moments, and
; SD estimates for each.

;This routine calculates the meand and SD for a pfitlov.out
;frame. These output values are written to a text file.

;The pfitlov frame needs a line of the form
;0  0  0  0  0  0  0  0  0 (9 of them) in between each bin region

readcol,list,f='a',silent=1,files
n0 = n_elements(files)

readcol,files[0],format='a,f,f,x,f,f,x,x,x',binname,v,d,h3,h4,silent=1
n1 = n_elements(v)

data = fltarr(4,n1,n0)

for j=0,n0-1 do begin
    readcol,files[j],format='a,f,f,x,f,f,x,x,x',binname,v,d,h3,h4
    data[0,*,j] = v
    data[1,*,j] = d
    data[2,*,j] = h3
    data[3,*,j] = h4
endfor

;final = fltarr(4,n1)
;for j=0,n1-1 do begin
;    final[*,j] = 

for j=0,n1-1 do begin
    trim = strsplit(binname[j],'_',/extract)
    binname[j] = trim[0]
endfor

;outname = ['NAME','NOF','MEAN SIG','SD SIG','MEAN V','SD V','MEAN H3','SD H3','MEAN H4','SD H4']
;f0 = '(a10,2x,a3,2x,a10,2x,a10,2x,a10,2x,a10,2x,a10,2x,a10,2x,a10,2x,a10,2x)'
;f1 = '(a10,2x,i3,2x,f10.5,2x,f10.5,2x,f10.5,2x,f10.5,2x,f10.5,2x,f10.5,2x,f10.5,2x,f10.5)'

;free_lun,5
;openw,5,'pfitave.txt'
;printf,5,outname,format=f0

stop
END


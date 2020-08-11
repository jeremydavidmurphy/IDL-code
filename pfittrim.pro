;This routine is used to cull down a pfitlov.out file. It will search
;for dispersion values above the cutoff value and toss them out. It
;returns a text file name TRIMMED.OUT.

PRO pfittrim,pfitfile

;-------------------------------------------------------------------------
cutoff = 400.0
;-------------------------------------------------------------------------

readcol,pfitfile,format='a,f,f,f,f,f,f,f,f',name,vel,sig1,sig2,h3,h4,r1,r2,r3

ind = where(sig1 lt cutoff)
n1 = n_elements(ind)

f1 = '(A27,2x,F8.3,2x,F8.3,2x,F8.3,2x,F6.3,2x,F6.3,2x,F8.3,2x,F7.3,2x,I4)'

openw,5,'trimmed.out'
for j=0,n1-1 do begin
    i = ind[j]
    printf,5,name[i],vel[i],sig1[i],sig2[i],h3[i],h4[i],r1[i],r2[i],r3[i],format=f1
endfor
free_lun,5

stop
END

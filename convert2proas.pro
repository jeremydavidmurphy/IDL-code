pro convert2proas, list
;this code is used to take Jenny's list and convert to H:M:S in
;the form proas likes

;readcol,list,f='a,f,f,x,x,x,f,x,x,i,i,x,a',name,ra,dec,sig,r1,r2,p
readcol,list,f='a,f,f',name,ra,dec
;readcol,list,f='f,f,a',ra,dec,name

n0 = n_elements(ra)

s = strsplit(list,'.',/extract)
sout = s[0]+'PROAS.txt'

openw,5,'proas.txt'

for j=0,n0-1 do begin
   radec,ra[j],dec[j],a,b,c,d,e,f
   a = strtrim(a,2)
   if b lt 10 then b = '0'+strtrim(b,2) else b = strtrim(b,2)
   c = round(c)
   if c lt 10 then c = '0'+strtrim(c,2) else c = strtrim(c,2)
   d = strtrim(d,2)
   if e lt 10 then e = '0'+strtrim(e,2) else e = strtrim(e,2)
   f = round(f)
   if f lt 10 then f = '0'+strtrim(f,2) else f = strtrim(f,2)
   out = name[j]+',  '+a+'.'+b+c+',  '+d+'.'+e+f
   printf,5,out
endfor

free_lun,5

stop
end

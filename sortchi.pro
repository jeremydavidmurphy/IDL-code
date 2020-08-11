;OBSOLETE: This routine has been superceded by Pchi.pro

PRO sortchi, file, OUTNAME=outname
;this routine is used to sort the results of a modeling run. After the
;cres files are catted together, these are sorted in order from best
;to worst chi-squared.

readcol,file,format='f,f,i,i,f,f,f',ml,bhm,vc,rs,x2stars,x2gc,x2tot

n1 = n_elements(ml)

sindex = bsort(x2stars)
sml = ml(sindex)
sbhm = bhm(sindex)
svc = vc(sindex)
srs = rs(sindex)
sx2stars = x2stars(sindex)
sx2gc = x2gc(sindex)
sx2tot = x2tot(sindex)

f1 = '(f5.2,2x,e12.5,2x,i4,2x,i3,2x,f9.3,2x,f9.3,2x,f9.3)'

if(n_elements(outname) eq 0) then outname = 'chi2.std.txt'
openw,5,outname
for j=0,n1-1 do printf,5,sml[j],sbhm[j],svc[j],srs[j],$
  sx2stars[j],sx2gc[j],sx2tot[j],format=f1
free_lun,5

STOP
END

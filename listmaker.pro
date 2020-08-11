;This routine is used to generate lists of all the possible
;combinations of sky subtraction for a given pointing. The list is a
;list of all the frames for a given pointing (ex: jm0324). THE NAME OF
;THE LIST IS USED TO GENERATE THE OUTPUT LIST NAMES. For example, if
;the list is named M87a.list then each output file will be named
;M87a_1.list, M87a_2.list, etc.
 
PRO listmaker, list

readcol,list,format='A',files

tag = ['aa','ab','ac','ba','bb','bc','ca','cb','cc']

n1 = n_elements(files)
n2 = n_elements(tag)

namearr = strarr(n2,n1)
outarray = strarr(n2,n1)
for j=0,n1-1 do namearr[*,j] = files[j] + tag

for j=0,n1-1 do begin
    outarray

stop
end

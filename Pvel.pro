; This routine is used to plot all the final model values for, in this
; case, the M87 data.

PRO Pvel

readcol,'pfitlov.out',/silent,format='a,f,f,x,f,f,x,x,x',bins1,vel,disp,h3,h4
n0 = n_elements(bins1)
readcol,'velbin.sav.stars.all',/silent,format='a,f,f,x,x,x,x,x',bins2,ra,dec
n1 = n_elements(bins2)

bin2trim = strarr(n1,2)
for j=0,n1-1 do begin
    temp = strsplit(bins2[j],'.',/extract)
    bin2trim[j,*] = temp
endfor

bin1trim = strarr(n0,2)
for j=0,n0-1 do begin
    temp = strsplit(bins1[j],'.',/extract)
    bin1trim[j,*] = temp
endfor

isauron = where(bin2trim[*,1] eq 'datP')
ivp = where(bin2trim[*,1] eq 'MlosvdP')

; now the bins i used in the modeling are IDed
readcol,'velbin.sav.stars',/silent,format='a,x,x,x,x,x,x,x',mybins2
n2 = n_elements(mybins2)

mybin2trim = strarr(n2,2)
for j=0,n2-1 do begin
    temp = strsplit(mybins2[j],'.',/extract)
    mybin2trim[j,*] = temp
endfor

imodel = [where(mybin2trim[*,1] eq 'dat') and where(mybin2trim[*,1] eq 'Mlosvd')]

stop

END

PRO pfitave, NAME=name

;This routine calculates the meand and SD for a pfitlov.out
;frame. These output values are written to a text file.

;The pfitlov frame needs a line of the form
;0  0  0  0  0  0  0  0  0 (9 of them) in between each bin region

if (n_elements(name) eq 0) then begin
    readcol,'pfitlov.out',format='a,f,f,x,f,f,x,x,x',binname,v,d,h3,h4
endif else readcol,name,format='a,f,f,x,f,f,x,x,x',binname,v,d,h3,h4

n1 = n_elements(v)
for j=0,n1-1 do begin
    trim = strsplit(binname[j],'_',/extract)
    binname[j] = trim[0]
endfor

outname = ['NAME','NOF','MEAN SIG','SD SIG','MEAN V','SD V','MEAN H3','SD H3','MEAN H4','SD H4']
f0 = '(a10,2x,a3,2x,a10,2x,a10,2x,a10,2x,a10,2x,a10,2x,a10,2x,a10,2x,a10,2x)'
f1 = '(a10,2x,i3,2x,f10.5,2x,f10.5,2x,f10.5,2x,f10.5,2x,f10.5,2x,f10.5,2x,f10.5,2x,f10.5)'

free_lun,5
openw,5,'pfitave.txt'
printf,5,outname,format=f0

cntr = 0
n2 = 0
repeat begin
cntr = cntr + n2
if (cntr ge n1) then goto, jumpend
first = binname[cntr]
i = where(binname eq first)
n2 = n_elements(i)
varr = v[i]
darr = d[i]
h3arr = h3[i]
h4arr = h4[i]
mv = mean(varr)
sv = stddev(varr)
md = mean(darr)
sd = stddev(darr)
mh3 = mean(h3arr)
sh3 = stddev(h3arr)
mh4 = mean(h4arr)
sh4 = stddev(h4arr)

printf,5,binname[cntr],n2,md,sd,mv,sv,mh3,sh3,mh4,sh4,format=f1

endrep until (n2 ge n1)
jumpend:
free_lun,5

END


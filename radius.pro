; This routine generates the luminosity-weighted radius file for a
; galaxy. (bin.radius.galaxyname). It uses the binname.coord files
; generated during the run of pipe1.pro

PRO radius, list, galname=galname

readcol,list,/silent,f='a',files
n0 = n_elements(files)

radiusout = fltarr(n0)
names = strarr(n0)

for j=0,n0-1 do begin
    name = files[j]
    readcol,name,f='i,f,f,f,f',/silent,skipline=1,$
      fiber,dra,ddec,rad,flux
    wgts = flux / max(flux)
;    for k=0,n_elements(wgts)-1 do print,strn(wgts[k])+' '+strn(flux[k])
    ibad = where(flux eq -1,count)
    if (count gt 0) then wgts[ibad] = 0.0
    top = total(rad * wgts)
    bot = total(wgts)
    radiusout[j] = top / bot
    n = strsplit(name,'.',/extract)
    names[j] = n[0]
endfor

if (n_elements(galname) eq 0) then begin
    galname = ''
    print,'Enter the name of your output file (bin.radius will be added)'
    read,galname
endif

nameout = 'bin.radius.'+galname
free_lun,5
openw,5,nameout
for j=0,n0-1 do printf,5,names[j],radiusout[j],format='(a12,2x,f9.3)'
free_lun,5

stop
END


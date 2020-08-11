; This routine generates the luminosity-weighted radius file for a
; galaxy. (bin.radius.galaxyname). It uses the binname.coord files
; generated during the run of pipe1.pro

PRO radius2D, list, galname=galname

; The change to this code is that now it returns the dRA and dDec
; offsets rather than the radial offsets.

readcol,list,/silent,f='a',files
n0 = n_elements(files)

offsetout = fltarr(2,n0)
names = strarr(n0)

for j=0,n0-1 do begin
    name = files[j]
    readcol,name,f='i,f,f,f,f',/silent,skipline=1,$
      fiber,dra,ddec,rad,flux
    wgts = flux / max(flux)
;    for k=0,n_elements(wgts)-1 do print,strn(wgts[k])+' '+strn(flux[k])
    ibad = where(flux eq -1,count)
    if (count gt 0) then wgts[ibad] = 0.0
    top1 = total(dra * wgts)
    top2 = total(ddec * wgts)
    bot = total(wgts)
    offsetout[0,j] = top1 / bot
    offsetout[1,j] = top2 / bot
    n = strsplit(name,'.',/extract)
    names[j] = n[0]
endfor

if (n_elements(galname) eq 0) then begin
    galname = ''
    print,'Enter the name of your output file (bin.radius will be added)'
    read,galname
endif

nameout = 'bin.offset.'+galname
free_lun,5
openw,5,nameout
for j=0,n0-1 do printf,5,names[j],offsetout[0,j],$
  offsetout[1,j],format='(a12,2x,f9.3,f9.3)'
free_lun,5

stop
END


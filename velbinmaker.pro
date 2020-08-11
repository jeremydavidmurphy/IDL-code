; This routine generates the velbin.dat file that goes into the
; dynamical modeling.

pro velbinmaker, losvdlist

; LOSVDLIST is the name of the losvds you want to go into
; velbin.dat. It has the following format:
; bin1001.Mlosvd  10  01
; bin1002.Mlosvd  10  02

readcol,'bin_r.out',skipline=2,format='i,x,f,x,x',radM,rmid
readcol,'bin_v.out',skipline=2,format='i,x,f,x,x',angM,amid

; degrees are converted to radians
amid = amid * 0.0174532925

readcol,losvdlist,format='a,i,i',losvdnames,radi,angi
n0 = n_elements(radi)
xout = fltarr(n0)
yout = fltarr(n0)

for j=0,n0-1 do begin
    r = rmid[where(radM eq radi[j])]
    a = amid[where(angM eq angi[j])]
    xout[j] = r * cos(a)
    yout[j] = r * sin(a)
;   print,'The file is '+losvdnames[j]
;   print,''
;   print,r,a,xout[j],yout[j]
;   pause
endfor

f1 = '(a15,2x,f10.4,f10.4,2x,i1,2x,i1,2x,i1,2x,i2,2x,i1)'
openw,5,'velbin.dat.new'
for j=0,n0-1 do printf,5,losvdnames[j],xout[j],yout[j],2,1,0,-1,0,format=f1
free_lun,5

stop
end

PRO estimateastrometryradius

readcol,'radiusmatch.txt',f1,f2
n0 = n_elements(f1)

if1 = f1 - 1
if2 = f2 - 1

readcol,'VP2.7.txt',f='i,f,f,x,x,x,x,x',skipline=1,fiber,x,y

radius = fltarr(n0)

for j=0,n0-1 do begin
   x1 = x[if1[j]]
   x2 = x[if2[j]]
   y1 = y[if1[j]]
   y2 = y[if2[j]]
   dx = x2-x1
   dy = y2-y1
   r = sqrt(dx^2+dy^2)
   radius[j] = r
endfor

stop
END

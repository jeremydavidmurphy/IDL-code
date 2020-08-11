;This program adds back in the dead fibers with values of -666 or whatever...

PRO addline, name, Newname=out, Mask=number
; The default replacement value is -666
;the fibers assumed dead are 
cd, 'data'
file = name+'.fits'
data = readfits(file)

IF (N_Elements(number) EQ 0) THEN number=-666
line = fltarr(1024)+number

newdata = [[data[*,0:14]], [line], [data[*,15:138]], [line], [data[*,139:196]], [line],$
 [data[*,197:210]], [line], [data[*,211:227]], [line], [line], [data[*,228:240]]]

IF (N_Elements(out) EQ 0) THEN out=file ELSE out=out+'.fits'

cd, '../combo'
writefits,out,newdata
cd, '../'

END


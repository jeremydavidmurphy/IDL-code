PRO flat_field_maker

; This routine is used to generate the flat-field for the re-imager
; bench. The routine now needs only an estimate of the center as it
; uses the ellipse fitting routine to improve the center finding part.

; As there are three fibers in each frame that are fully illuminated I
; use all three to determine the relative flux. To alter this, change
; the dx and dy values below (specifying, in pixels, the offset
; between the first fiber and the others) and the 'l' forloop. 
;***************************************************************
wscale = 3.0
slow = 'off'
major = 3.0
r1 = 48.0 ;the inner radii for the circ app phot
r2 = 52.0 ;the outer radii for the circ app phot
r3 = 50.0
dx = [0.0,104.0,207.0] ;the typical offset in X
dy = [0.0,1.0,1.0] ;the typical offset in Y

cut0 = 3500.0
plotting = 'off' ;just the surface plotting...

plotting2 = 'off' ;this plotting may not reflect the real ellipse fitting. the plotting can make the ellipse offset from the spot, but this is just a scaling issue and not real, meaning, just leave this plotting routine off.

; THERE IS ANOTHER PLOTTING ROUTINE IN CIRC_APPHOT. YOU CAN TURN THAT
; ON OR OFF ON LINE ~#77. This routine plots where the circular
; appertures for background subtraction are falling and is worth
; having on so that you know where the background is being taken from.

;***************************************************************
readcol,'centers.txt',f='a,a,i,i',silent=1,files,bgfiles,x,y
n0 = n_elements(files)
x = x - 1.0
y = y - 1.0
flux= fltarr(3,n0)
bg = fltarr(3,n0)

center = fltarr(2,11)

xreal = intarr(n0,3)
yreal = intarr(n0,3)

if (plotting2 eq 'on') then $
window,0,retain=2,xsize=(r3*2+1)*8,ysize=(r3*2+1)*8
free_lun,5
openw,5,'values.txt'
for j=0,n0-1 do begin ;a loop over each frame...
    filename = files[j]
    print,'Working on frame: '+filename
    file = readfits(files[j],/silent)
    bgfile = readfits(bgfiles[j],/silent)
    file = float(file) - float(bgfile)
    nn0 = n_elements(file[*,0])
    nn1 = n_elements(file[0,*])

    for l=0,2 do begin ;a loop through the number of fibers fully illuminated
        x1 = x[j]+dx[l]
        y1 = y[j]+dy[l]

        chunk = file[x1-r3:x1+r3,y1-r3:y1+r3]
        for k=0,10 do begin     ;an iteration to refine the centers
            cutoff1 = cut0 + 200.0*k
            index1 = where(chunk gt cutoff1)
            circle  = Fit_Ellipse(index1, XSize=n_elements(chunk[*,0]), YSize=n_elements(chunk[0,*]), CENTER=c1)
            center[0,k] = c1[0] & center[1,k] = c1[1]
            if (plotting2 eq 'on') then begin
                if (k eq 0) then begin
                    wset,0
                    loadct,37,/silent
                    TVImage, BytScl(chunk, Top=!D.Table_Size-3)
                endif
                plots,circle*8,/Device,Color=FSC_Color('white'),thick=2
            endif
        endfor
        if (slow eq 'on') then wait,1.0
        xc = median(center[0,*],/even)
        yc = median(center[1,*],/even)
        
        ixreal = x1 + (xc - r3)
        iyreal = y1 + (yc - r3)
        xreal[j,l] = round(ixreal)
        yreal[j,l] = round(iyreal)

        xc = xc + 1.0
        yc = yc + 1.0
        out = circ_appphotf(file,ixreal+1,iyreal+1,major,r1,r2,win=wscale,plot='off')
        flux[l,j] = out[0]
        bg[l,j] = out[1]
        printf,5,xreal[j,l], yreal[j,l], out[0], out[1]
        print,xreal[j,l], yreal[j,l], out[0], out[1]
    endfor
endfor
free_lun,5

nflux1 = flux[0,*]/max(flux[0,*])
nflux2 = flux[1,*]/max(flux[1,*])
nflux3 = flux[2,*]/max(flux[2,*])

; notes on triangulate: TRIANGULATE conducts a triangulation on the
; data. X and Y and the x and y points to be triangulated. triangles
; here is a 3xN array, where the 0,1 and 2 values of array element 'i'
; is the indices for the 3 triangle coordinates that define the
; corners. Here, 'b' is optional and gives a list of the indicies of
; the boundry points. This is needed if you want to extrapolate beyond
; the surface you're working on.

; notes on trigrid: TRIGRID takes the output of triangulate. here, x,
; y, and nflux are the 3 dimensions you are trying to interpolate
; between. X and Y being self-explanatory, and nflux# being the Z
; component. 

myxout = findgen(4007)
myyout = findgen(2671)

triangulate,x,y,triangles,b
gridspace = [1.0,1.0]
result1 = trigrid(xreal[*,0],yreal[*,0],nflux1,triangles,gridspace,/quintic,xout=myxout,yout=myyout)
result2 = trigrid(xreal[*,1],yreal[*,1],nflux2,triangles,gridspace,/quintic,xout=myxout,yout=myyout)
result3 = trigrid(xreal[*,2],yreal[*,2],nflux3,triangles,gridspace,/quintic,xout=myxout,yout=myyout)

i = where(result1 lt 0.1)
result1[i] = !values.f_nan
s = smooth(result1,550,/edge_truncate,/nan)
s[i] = 10000.0
writefits,'flat1.fits',s
s1 = s

i = where(result2 lt 0.1)
result2[i] = !values.f_nan
s = smooth(result2,550,/edge_truncate,/nan)
s[i] = 10000.0
writefits,'flat2.fits',s
s2 = s

i = where(result3 lt 0.1)
result3[i] = !values.f_nan
s = smooth(result3,550,/edge_truncate,/nan)
s[i] = 10000.0
writefits,'flat3.fits',s
s3 = s

ss = [[[s1]],[[s2]],[[s3]]]
mss = median(ss,dim=3)
writefits,'flatMED.fits',mss

if (plotting eq 'on') then begin
   window,0,retain=2
   device,decomposed=0
   surface,result1,xreal[*,0],yreal[*,0],az=120,zrange=[0.92,1.0]
   pause
   surface,result2,xreal[*,1],yreal[*,1],az=120,zrange=[0.92,1.0]
   pause
   surface,result3,xreal[*,2],yreal[*,2],az=120,zrange=[0.92,1.0]
   pause
   while !d.window ne -1 do wdelete, !d.window
endif

stop
END

;----------------------------------------------------------------------
pro display_bins, xbin, ybin, vel, x, y, RANGE=range, $
		 PA=pa, PIXELSIZE=pixelsize, _EXTRA=extra

compile_opt idl2
on_error, 2

; HISTORY:
;   V1.0: Written by Michele Cappellari, Leiden, 18 February 2003
;   V2.0: Written by Jesus Falcon, Univ. of Nottingham, 01 April 2003
;   	Addition of more arguments in the procedure.
;

; Perform a Voronoi tessellation starting from the coordinates
; of the generators and the coordinates of the original pixels
;
out = x*0
FOR j=0L,n_elements(x)-1 DO BEGIN
    tmp = min((x[j]-xbin)^2+(y[j]-ybin)^2,index)
    out[j] = vel[index]
ENDFOR
display_pixels, x, y, out, RANGE=range, $
	    PIXELSIZE=pixelsize, PA=pa, _EXTRA=extra

end
;----------------------------------------------------------------------

; Used to display a fits image to screen

PRO display_fits, fitsfile
COMPILE_OPT IDL2

image = readfits(fitsfile,/silent)

swx = !D.X_VSIZE    ; Window size in X in device units
swy = !D.Y_VSIZE    ; Size in Y
s = SIZE(image)
six = FLOAT(s[1])   ; Image size in X in pixels
siy = FLOAT(s[2])   ; Size in Y

ratio = siy/swy > six/swx
TVSCL, POLY_2D(image, [[0,0],[ratio,0]],[[0,ratio],[0,0]], 0, six/ratio, siy/ratio)

END

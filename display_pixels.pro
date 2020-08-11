;######################################################################
;+
; NAME:
;     DISPLAY_PIXELS
;
; AUTHOR:
;       Michele Cappellari, Leiden Observatory, The Netherlands
;       cappellari@strw.leidenuniv.nl
;
; PURPOSE:
;       Plot a set of square pixels, given a set of (X,Y) coordinates
;       and of corresponding values to display. This is uesful when the pixels
;       are not arranged in an full array and cannot be displayed as an image.
;       This is often the case for integral-field spectroscopic data.
;
; CALLING SEQUENCE:
;       DISPLAY_PIXELS, x, y, vel, PA=pa, PIXELSIZE=pixelSize, RANGE=[minVel,maxVel]
;
; INPUTS:
;       X: Vector containing the X coordinate of the center of each pixel.
;       Y: same as XBIN for the Y coordinate of each bin.
;       VEL: Vector containing the quantity to plot (e.g. velocity)
;           associated to each pixel of coordinates (X,Y).
;
; KEYWORDS:
;       PA: position agle of the side of each pixel.
;       PIXELSIZE: vector with the size of each pixel in x and y,
;            in the same units as the coordinates.
;           If this is not given then its value is computed automatically.
;       RANGE: two elements vector [minVel,maxVel] defining the range of
;           VEL values to plot. Any value outside this range will be
;           visualized with the maximum or minimum color in the colormap.
;           By default RANGE=[min(VEL),max(VEL)].
;       NOTISO: Does not plot with /ISO 
;       _EXTRA: any additional keyword can be passed to PLOT via the _EXTRA
;           mechanism (e.g. TITLE, XTITLE, XRANGE,...).
;
; MODIFICATION HISTORY:
;       V1.0: Michele Cappellari, Leiden, January 2000
;       V1.1: added RANGE keyword and cuts for numbers out of range.
;           MC, Leiden, 3 December 2002
;       V1.2: Added more input checking. MC, Leiden, 12 December 2004
;       V1.21: Shortened loop. MC, Leiden, 10 September 2005
;       V1.22: Perform robust but slow pixel size calculation by default.
;           Prevent tick marks to be covered by pixels. Updated documentation.
;           MC, Leiden, 29 September 2005
;       V1.23: Allowed for unqual pixel sizes 
;              Added NOTISO keyword
;       V1.24: fixed YRANGE scaling for non-squre pixels, RvdB, Leiden
;              Okt 2006
;-
;----------------------------------------------------------------------------
pro display_pixels, x, y, vel, PA=pa, PIXELSIZE=pixelSize, RANGE=range, $
                    NOTISO=notiso, _EXTRA=ex
compile_opt idl2
on_error, 2
 
n = N_ELEMENTS(vel)
if n ne n_elements(x) or n ne n_elements(y) then $
    message, 'The vectors (X, Y, VEL) must have the same size'

; For each point, find the distance to all other points and select the minimum.
; This is a robust but slow way of determining the pixel size of unbinned data.
;
if n_elements(pixelSize) eq 0 then begin
    dx = 1e30
    for j=0,n-2 do dx = min((x[j]-x[j+1:*])^2 + (y[j]-y[j+1:*])^2) < dx
    pixelSize = sqrt(dx)
endif
if n_elements(pixelSize) eq 1 then begin
    pixelSize = [pixelSize,pixelSize]
endif

if n_elements(pa) eq 0 then begin
    cpa = 1.0
    spa = 0.0
endif else begin
    cpa = COS(pa/!RADEG)
    spa = SIN(pa/!RADEG)
endelse
if n_elements(range) eq 0 then begin
    velMax = MAX(vel,MIN=velMin)
    velRange = velMax - velMin
endif else begin
    velMin = range[0]
    velRange = range[1] - range[0]
endelse

maxx = max(x, MIN=minx)
maxy = max(y, MIN=miny)
xaper  = [-0.5, -0.5, +0.5, +0.5, -0.5] * pixelSize[0]
yaper  = [+0.5, -0.5, -0.5, +0.5, +0.5] * pixelSize[1]
x1 = xaper*cpa - yaper*spa
y1 = xaper*spa + yaper*cpa
color = 255.0/velRange*(vel-velMin) < 255 > 0

; Plot an empty frame first and then overplot the axes on top of it
; to prevent the axes tick marks from being covered by the pixels.
;
dx = pixelSize[0]/sqrt(2) ; Make sure pixels always lie inside the plot
dy = pixelSize[1]/sqrt(2)

if (keyword_set(NOTISO)) then $
    PLOT, [minx-dx, maxx+dx], [miny-dy, maxy+dy], /NODATA, XSTYLE=1, YSTYLE=1 , _EXTRA=ex

if (not keyword_set(NOTISO)) then $
    PLOT, [minx-dx, maxx+dx], [miny-dy, maxy+dy], /NODATA, XSTYLE=5, YSTYLE=5 , /ISO,_EXTRA=ex

FOR j=0L, n-1L DO POLYFILL, x[j]+x1, y[j]+y1, COLOR=color[j]
AXIS, XAXIS=0, XRANGE=[minx-dx, maxx+dx], /XSTYLE, _EXTRA=ex
AXIS, XAXIS=1, XRANGE=[minx-dx, maxx+dx], /XSTYLE, XTICKFORMAT='(A1)'
AXIS, YAXIS=0, XRANGE=[miny-dx, maxy+dx], /YSTYLE, _EXTRA=ex
AXIS, YAXIS=1, XRANGE=[miny-dx, maxy+dx], /YSTYLE, YTICKFORMAT='(A1)'

END
;----------------------------------------------------------------------

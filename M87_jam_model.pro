;#############################################################################
;
; In this example/test we fit a one dimensional MGE to the density 
; of a Hernquist (1990, ApJ, 356, 359; hereafter H90) model.
; - IMPORTANT: One needs the MGE_FIT_1D routine from the MGE_FIT_SECTORS
; package available from http://purl.org/cappellari/idl. 
; We then compute the following quantities for a spherical model:
; (i) The circular velocity with MGE_CIRCULAR_VELOCITY;
; (ii) The sigma of an isotropic model, with both the JAM_SPHERICAL_RMS 
;     and JAM_AXISYMMETRIC_RMS routines;
; (iii) The sigma of a fully tangential anisotropic model with both the 
;     JAM_SPHERICAL_RMS and JAM_AXISYMMETRIC_RMS routines.
; In all cases the solutions are compared with the analytic results by H90.   
;
; Michele Cappellari, Oxford, 28 November 2008
; 
;#############################################################################
PRO M87_jam_model

M = 1d11          ; Total mass in Solar Masses
a = 1d3           ; Break radius in pc
distance = 16.5   ; Assume Virgo distance in Mpc (Mei et al. 2007)
pc = distance*!dpi/0.648d ; Constant factor to convert arcsec --> pc
G = 0.00430237d    ; (km/s)^2 pc/Msun [6.674e-11 SI units]
mbh = 0d

n = 100 ; Number of values to sample the H90 profile for the fit
r = range(a/100,a*100,n,/LOG) ; logarithically spaced radii in pc
rho = M*a/(2d*!dpi*r*(r+a)^3) ; Density in Msun/pc^3 (H90 equation 2)

mge_fit_1d, r, rho, NGAUSS=16, SOL=sol

surf = sol[0,*] ; Surface density in Msun/pc^2
sigma = sol[1,*]/pc ; Gaussian dispersion in arcsec
qObs = surf*0+1 ; Assume spherical model
inc = 90d      ; Edge-on view
rad = range(0.01,50,100) ; desired output radii in arcsec (avoid R=0)
meshgrid, rad, rad, x, y ; Create regular grid of coordinates 


;################# Circular Velocity #################

mge_circular_velocity, surf, sigma, qObs, inc, mbh, distance, rad, vcirc
loadct, 12
plot, rad, vcirc, XTITLE='R (arcsec)', YTITLE='Vc, Vrms (km/s)', THICK=8

stop
; Compare with analytic result
;
vc = sqrt(G*M*r)/(r+a) ; H90 equation (16)
oplot, r/pc, vc, THICK=3, COLOR=200
xyouts, 30, 310, 'Circular Velocity'

;################### Isotropic Vrms ###################

; Spherical isotropic H90 model
;
jam_spherical_rms, surf, sigma, surf, sigma, mbh, distance, rad, sigp
oplot, rad, sigp, THICK=8

; Axisymmetric isotropic model in the spherical limit. 
; This is plotted on the major axis, but the Vrms has circular symmetry
;
jam_axisymmetric_rms, surf, sigma, qObs, surf, sigma, qObs, $
    inc, mbh, distance, x, y, vrms
oplot, rad, vrms[*,0], THICK=5, COLOR=100
    
; Analytic surface brightness from H90
;
s = r/a
xs = r*0
w1 = where(s lt 1, COMPLEMENT=w2)
acosh = alog( 1d/s[w1] + sqrt(1d/s[w1]^2 - 1d) )
xs[w1] = acosh / sqrt(1d - s[w1]^2) ; H90 equation (33)
xs[w2] =  acos(1d/s[w2]) / sqrt(s[w2]^2 - 1d) ; H90 equation (34)
IR = M*((2 + s^2)*xs - 3) / (2d*!dpi*a^2*(1d - s^2)^2) ; H90 equation (32)

; Projected second moments of isotropic model from H90
;
sigp = sqrt( G*M^2/(12*!dpi*a^3*IR) $ H90 equation (41)
     * ( 0.5/(1-s^2)^3 $
     * ( -3*s^2*xs* (8*s^6 - 28*s^4 + 35*s^2 - 20) $
     - 24*s^6 + 68*s^4 - 65*s^2 + 6) - 6*!dpi*s) )
oplot, r/pc, sigp, THICK=2, COLOR=200        
xyouts, 30, 110, 'Isotropic Vrms'

;################## Anisotropic Vrms ##################

; Projected second moments for a H90 model with sigma_R=0.
; This implies beta=-Infinity but I adopt as an approximation 
; below a large negative beta.
;
beta = -20
jam_spherical_rms, surf, sigma, surf, sigma, mbh, distance, rad, sigp, BETA=beta
oplot, rad, sigp, THICK=8

; Axisymmetric anisotropic model in the spherical limit.
; The spherical Vrms is the quadratic average of major 
; and minor axes of the axisymmetric model.
;
jam_axisymmetric_rms, surf, sigma, qObs, surf, sigma, qObs, $
    inc, mbh, distance, x, y, vrms, BETA=beta ; Major axis
oplot, rad, sqrt((vrms[0,*]^2+vrms[*,0]^2)/2d), THICK=5, COLOR=100

; Projected second moments of fully tangential model from H90
;
sigp = sqrt( G*M^2*r^2/(2*!dpi*a^5*IR) $ H90 equation (42)
     * ( 1d/(24*(1-s^2)^4) $
     * ( -xs*(24*s^8 - 108*s^6 + 189*s^4 - 120*s^2 + 120) $
     - 24*s^6 + 92*s^4 - 117*s^2 + 154) + 0.5*!dpi/s) )
oplot, r/pc, sigp, THICK=2, COLOR=200        
xyouts, 30, 180, 'Fully Tangential Vrms'

END
;----------------------------------------------------------------------------

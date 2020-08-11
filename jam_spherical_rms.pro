;#############################################################################
;
; Copyright (C) 2004-2010, Michele Cappellari
; E-mail: cappellari_at_astro.ox.ac.uk
;
; Updated versions of the software are available from my web page
; http://purl.org/cappellari/idl
;
; If you have found this software useful for your research,
; I would appreciate an acknowledgment to the use of the
; "JAM modelling method of Cappellari (2008)"
;
; This software is provided as is without any warranty whatsoever.
; Permission to use, for non-commercial purposes is granted.
; Permission to modify for personal or internal use is granted,
; provided this copyright and disclaimer are included unchanged
; at the beginning of the file. All other rights are reserved.
;
;#############################################################################
;+
; NAME:
;   JAM_SPHERICAL_RMS
;
; PURPOSE:
;    This procedure calculates a prediction for the projected second
;    velocity moments V_RMS = sqrt(V^2 + sigma^2), or for a non-rotating
;    galaxy V_RMS = sigma, for an anisotropic spherical galaxy model.
;    It implements the solution of the anisotropic Jeans equations
;    presented in equation (50) of Cappellari (2008, MNRAS, 390, 71).
;    PSF convolution is done as described in the Appendix of that paper.
;    http://adsabs.harvard.edu/abs/2008MNRAS.390...71C
;
; CALLING SEQUENCE:
;    JAM_SPHERICAL_RMS, $
;        surf_lum, sigma_lum, surf_pot, sigma_pot, mbh, distance, rad, rmsModel, $
;        BETA=beta, CHI2=chi2, ERMS=erms, NORMPSF=normPsf, ML=ml, NRAD=nrad, $
;        PIXSIZE=pixSize, /PLOT, /QUIET, RMS=rms, SIGMAPSF=sigmaPsf, STEP=step
;
; INPUT PARAMETERS:
;   SURF_LUM: vector of length N containing the peak surface brightness of the
;       MGE Gaussians describing the galaxy surface brightness in units of
;       Lsun/pc^2 (solar luminosities per parsec^2).
;   SIGMA_LUM: vector of length N containing the dispersion in arcseconds of
;       the MGE Gaussians describing the galaxy surface brightness.
;   SURF_POT: vector of length M containing the peak value of the MGE Gaussians
;       describing the galaxy surface density in units of Msun/pc^2 (solar
;       masses per parsec^2). This is the MGE model from which the model
;       potential is computed.
;     - In a common usage scenario, with a self-consistent model, one has
;       the same Gaussians for both the surface brightness and the potential.
;       This implies SURF_POT = SURF_LUM, SIGMA_POT = SIGMA_LUM.
;       The M/L, by which SURF_POT has to be multiplied to best match the
;       data, is fitted by the routine when passing the RMS and ERMS
;       keywords with the observed kinematics.
;   SIGMA_POT: vector of length M containing the dispersion in arcseconds of
;       the MGE Gaussians describing the galaxy surface density.
;   MBH: Mass of a nuclear supermassive black hole in solar masses.
;     - VERY IMPORTANT: The model predictions are computed assuming SURF_POT
;       gives the total mass. In the common self-consistent case one has
;       SURF_POT = SURF_LUM and if requested (keyword ML) the program can scale
;       the output RMSMODEL to best fit the data. The scaling is equivalent to
;       multiplying *both* SURF_POT and MBH by a factor M/L. To avoid mistakes,
;       the actual MBH used by the output model is printed on the screen.
;   DISTANCE: distance of the galaxy in Mpc.
;   RAD: Vector of length P with the (positive) radius from the galaxy center
;       in arcseconds of the bins (or pixels) at which one wants to compute
;       the model predictions.
;     - When no PSF/pixel convolution is performed (SIGMAPSF=0 or PIXSIZE=0)
;       there is a singularity at RAD=0 which should be avoided.
;
; KEYWORDS:
;   BETA: Vector of length N with the anisotropy
;       beta = 1 - (sigma_theta/sigma_R)^2 of the individual MGE Gaussians. A
;       scalar can be used if the model has constant anisotropy.
;   CHI2: Reduced chi^2 describing the quality of the fit
;       chi^2 = total( ((rms-rmsModel)/erms)^2 ) / n_elements(rms)
;   ERMS: Vector of length P with the 1sigma errors associated to the RMS
;       measurements. From the error propagation
;       ERMS = sqrt((dVel*velBin)^2 + (dSig*sigBin)^2)/RMS,
;       where velBin and sigBin are the velocity and dispersion in each bin
;       and dVel and dSig are the corresponding errors
;       (Default: constant errors ERMS=0.05*MEDIAN(RMS)).
;   ML: Mass-to-light ratio to multiply the values given by SURF_POT.
;       Setting this keyword is completely equivalent to multiplying the
;       output RMSMODEL by SQRT(M/L) after the fit. This implies that the
;       BH mass becomes MBH*(M/L).
;     - If this keyword is set to a negative number in input, the M/L is
;       fitted from the data and the keyword returns the best-fitting M/L
;       in output. The BH mass of the best-fitting model is MBH*(M/L).
;   NORMPSF: Vector of length Q with the fraction of the total PSF flux
;       contained in the various circular Gaussians describing the PSF of the
;       observations. It has to be total(NORMPSF) = 1. The PSF will be used for
;       seeing convolution of the model kinematics.
;   NRAD: Number of logarithmically spaced radial positions for which the
;       models is evaluated before interpolation and PSF convolution. One may
;       want to increase this value if the model has to be evaluated over many
;       orders of magnitutes in radius (default: NRAD=50).
;   PIXSIZE: Size in arcseconds of the (square) spatial elements at which the
;       kinematics is obtained. This may correspond to the size of the spaxel
;       or lenslets of an integral-field spectrograph. This size is used to
;       compute the kernel for the seeing and aperture convolution.
;     - If this is not set, or PIXSIZE = 0, then convolution is not performed.
;   /PLOT: Set this keyword to produce a plot at the end of the calculation.
;   /QUIET: Set this keyword not to print values on the screen.
;   RMS: Vector of length P with the input observed stellar
;       V_RMS=sqrt(velBin^2 + sigBin^2) at the coordinates positions
;       given by the vector RAD.
;     - If RMS is set and ML is negative or not set, then the model is fitted to
;       the data, otherwise the adopted ML is used and just the chi^2 is returned.
;   SIGMAPSF: Vector of length Q with the dispersion in arcseconds of the
;       circular Gaussians describing the PSF of the observations.
;     - If this is not set, or SIGMAPSF = 0, then convolution is not performed.
;     - IMPORTANT: PSF convolution is done by creating a 2D image, with pixels
;       size given by STEP=MAX(SIGMAPSF,PIXSIZE/2)/4, and convolving it with
;       the PSF + aperture. If the input radii RAD are very large with respect
;       to STEP, the 2D image may require a too large amount of memory. If this
;       is the case one may compute the model predictions at small radii
;       separately from those at large radii, where PSF convolution is not
;       needed.
;   STEP: Spatial step for the model calculation and PSF convolution in arcsec.
;       This value is automatically computed by default as
;       STEP=MAX(SIGMAPSF,PIXSIZE/2)/4. It is assumed that when PIXSIZE or
;       SIGMAPSF are big, high resolution calculations are not needed. In some
;       cases however, e.g. to accurately estimate the central Vrms in a very
;       cuspy galaxy inside a large aperture, one may want to override the
;       default value to force smaller spatial pixels using this keyword.
;     - Use this keyword to set the desired scale of the model when no PSF or
;       pixel convolution is performed (SIGMAPSF=0 or PIXSIZE=0).
;
; OUTPUT PARAMETER:
;   RMSMODEL: Vector of length P with the model predictions for the velocity
;       second moments (sigma in the spherical non-rotating case) of each bin.
;
; USAGE EXAMPLE:
;    A simple usage example is given in the procedure TEST_JAM_SPHERICAL_RMS at
;    the end of this file.
;
; REQUIRED ROUTINES:
;       By M. Cappellari (included in the JAM distribution):
;       - ANY
;       - DIFF
;       - IBETAM
;       - MESHGRID
;       - QUADVA
;       - RANGE
;
; MODIFICATION HISTORY:
; V1.0: Written and tested isotropic case.
;    Michele Cappellari, Vicenza, 10 August 2004
; V2.0: Included anisotropic case with 1D integral. MC, Oxford, 4 April 2008
; V3.1: First released version. MC, Oxford, 12 August 2008
; V3.2: Updated documentation. MC, Oxford, 14 August 2008
; V4.0: Implemented PSF convolution using interpolation on polar grid.
;     Dramatic speed-up of calculation. Further documentation.
;     MC, Oxford, 11 September 2008
; V4.01: Included keyword STEP. MC, Windhoek, 29 September 2008
; V4.02: Added keywords NRAD. Thanks to Michael Williams for reporting possible
;     problems with too coarse interpolation. MC, Oxford, 21 November 2008
; V4.1: Added keywords CHI2, ERMS, ML, /PRINT, /QUIET, RMS as in the
;     JAM_AXISYMMETRIC_RMS routine. Updated the usage example routine
;     TEST_JAM_SPHERICAL_RMS. MC, Oxford, 04 February 2010
; -
;#############################################################################
FUNCTION sphani_integrand_spherical_jeans, r, $
    SIG_L=sig_l, SIG_M=sig_m, LUM=lum, MASS=mass, MBH=Mbh, RMIN=rmin, BETA=beta
compile_opt idl2, hidden
;
; This function implements the integrand of equation (50) of Cappellari (2008).
; The routine tries to speed up the calculation by treating differently the three
; cases: (i) isotropic, (ii) constant-anisotropy and (iii) variable-anisotropy.

mass_r = Mbh
for j=0L,n_elements(mass)-1 do begin
    h = r/(SQRT(2d)*sig_m[j])
    mass_r += mass[j]*(ERF(h) - 2d/SQRT(!dpi)*h*EXP(-h^2)) ; equation (49)
endfor

if array_equal(beta,beta[0]) then begin   ; Faster constant-anisotropy model
    if beta[0] eq 0 then $                ; Isotropic case
        er = sqrt(r^2-rmin^2) $           ; equation (44)
    else begin                            ; Anisotropic case
        rat = (rmin/r)^2
        er = 0.5d*rmin/rat^beta[0]* $     ; equation (43)
            ( beta[0]*ibetam(0.5d + beta[0],0.5d,rat) - ibetam(beta[0] - 0.5d,0.5d,rat) $
            + sqrt(!dpi)*(1.5d - beta[0])*gamma(beta[0] - 0.5d)/gamma(beta[0]) )
    endelse
    lum_dens = 0d
    FOR j=0L,n_elements(lum)-1 DO $      ; equation (47)
        lum_dens += lum[j]*EXP(-0.5d*(r/sig_l[j])^2)/(SQRT(2d*!dpi)*sig_l[j])^3
    fun = er * lum_dens
endif else begin                          ; Slower variable-anisotropy model
    fun = 0d
    rat = (rmin/r)^2
    for j=0,n_elements(lum)-1 DO begin
        if beta[j] eq 0 then $            ; Isotropic case
            er = sqrt(r^2-rmin^2) $       ; equation (44)
        else $                            ; Anisotropic case
            er = 0.5d*rmin/rat^beta[j]* $ ; equation (43)
                ( beta[j]*ibetam(0.5d + beta[j],0.5d,rat) - ibetam(beta[j] - 0.5d,0.5d,rat) $
                + sqrt(!dpi)*(1.5d - beta[j])*gamma(beta[j] - 0.5d)/gamma(beta[j]) )
        fun += er * lum[j]*EXP(-0.5d*(r/sig_l[j])^2)/(SQRT(2d*!dpi)*sig_l[j])^3
    endfor
endelse

; This routine returns a vector of values computed at different values of r
;
G = 0.00430237d    ; (km/s)^2 pc/Msun [6.674e-11 SI units]

return, 2d*G*fun*mass_r/r^2
END
;------------------------------------------------------------------
FUNCTION sphani_weighted_sigma2, R, s
compile_opt idl2, hidden
;
; Integration of equation (50)

s.rmin = R
rmax = 3d*max(s.sig_l)
if R ge rmax then message, 'R > rmax'
quadva, 'sphani_integrand_spherical_jeans', [R, rmax], int, $
    RELTOL=1d-5, ABSTOL=0d, FUNCTARGS=s

return, int
END
;------------------------------------------------------------------
function shpani_second_moment, R, s
compile_opt idl2, hidden
;
; This routine gives the second V moment after convolution with a PSF.
; The convolution is done using interpolation of the model on a
; polar grid, as described in Appendix A of Cappellari (2008).

if max(s.sigmaPsf) gt 0 and s.pixSize gt 0 then begin ; PSF convolution

    ; Kernel step is 1/4 of largest value between sigma(min) and 1/2 pixel side.
    ; Kernel half size is the sum of 3*sigma(max) and 1/2 pixel diagonal.
    ;
    if s.step gt 0 then step = s.step $
        else step = (s.pixSize/2d > min(s.sigmaPsf))/4d
    mx = 3d*max(s.sigmaPsf) + s.pixSize/sqrt(2d)

    ; Make grid linear in log of radius RR
    ;
    rmax = max(R) + mx ; Radius of circle containing all data + convolution
    nrad = s.NRAD ; radial elements
    logRad = range(alog(step),alog(rmax),nrad) ; Linear grid in log(RR)
    rr = exp(logRad)

    ; The model Vrms computation is only performed on the radial grid
    ; which is then used to interpolate the values at any other location
    ;
    wm2Pol = rr*0
    mgePol = wm2Pol
    for j=0,n_elements(rr)-1 do begin
        wm2Pol[j] = sphani_weighted_sigma2(rr[j],s)
        mgePol[j] = total( s.surf_l * exp(-0.5d*(rr[j]/s.sig_l)^2) )
    endfor

    nx = ceil(rmax/step/64d)*64  ; Make 2*nx dimensions divisible by 2^7
    x1 = range(-nx,nx,2*nx)*step ; for much faster FFT computation.
    meshgrid, x1, x1, xCar, yCar ; Cartesian grid for convolution

    ; Interpolate MGE model and Vrms over cartesian grid
    ;
    r1 = 0.5d*alog(xCar^2 + yCar^2) ; Log radius of cartesian grid
    indx = (r1 - logRad[0])/(logRad[1] - logRad[0]) ; Convert to indices for INTERPOLATE
    wm2Car = interpolate(wm2Pol,indx,CUBIC=-0.5)
    mgeCar = interpolate(mgePol,indx,CUBIC=-0.5)

    nk = ceil(mx/step)
    kgrid = range(-nk,nk,2*nk)*step
    meshgrid, kgrid, kgrid, xgrid, ygrid ; Kernel is square

    ; Compute kernel with equation (A6) of Cappellari (2008).
    ; Normaliztion is irrelevant here as it cancels out.
    ;
    kernel = 0d
    dx = s.pixSize/2d
    sp = sqrt(2d)*s.sigmaPsf
    for j=0,n_elements(s.sigmapsf)-1 do $
        kernel += s.normPsf[j] $
            * (erf((dx-xgrid)/sp[j]) + erf((dx+xgrid)/sp[j])) $
            * (erf((dx-ygrid)/sp[j]) + erf((dx+ygrid)/sp[j]))

    ; Seeing and aperture convolution with equation (A3)
    ;
    muCar = sqrt(convolve(wm2Car,kernel)/convolve(mgeCar,kernel))

    ; Interpolate convolved image at observed apertures.
    ; Aperture integration was already included in the kernel.
    ;
    indx = (R/sqrt(2) - x1[0])/(x1[1] - x1[0]) ; Convert to indices for INTERPOLATE
    mu = interpolate(muCar,indx,indx)

endif else begin ; No PSF convolution: just compute values

    mu = R*0
    for j=0,n_elements(R)-1 do begin
        wm2Pol = sphani_weighted_sigma2(R[j],s)
        mgePol = total( s.surf_l * exp(-0.5d*(R[j]/s.sig_l)^2) )
        mu[j] = sqrt(wm2Pol/mgePol)
    endfor

endelse

return, mu
end
;----------------------------------------------------------------------
pro jam_spherical_rms, $
    surf_lum, sigma_lum, surf_pot, sigma_pot, mbh, distance, rad, rmsModel, $
    BETA=beta, NORMPSF=normPsf, PIXSIZE=pixSize1, SIGMAPSF=sigmaPsf1, $
    STEP=step, NRAD=nrad, RMS=rms, ERMS=erms, CHI2=chi2, ML=ml, $
    QUIET=quiet, PLOT=plot
compile_opt idl2
on_error, 2

if n_elements(beta) eq 0 then beta = 0d
if n_elements(sigmaPsf1) eq 0 then sigmaPsf1 = 0d
if n_elements(normPsf) eq 0 then normPsf = 1d
if n_elements(pixSize1) eq 0 then pixSize1 = 0d
if n_elements(erms) eq 0 && n_elements(rms) gt 0 then erms = rms*0+median(rms)*0.05 ; Constant 5% errors
if n_elements(mbh) eq 0 then mbh = 0d
if n_elements(step) eq 0 then step = 0
if n_elements(nrad) eq 0 then nrad = 50

pc = distance*!dpi/0.648d ; Constant factor to convert arcsec --> pc

sigmaPsf = sigmaPsf1*pc
pixSize = pixSize1*pc
step_pc = step*pc

sigma_lum_pc = sigma_lum*pc     ; Convert from arcsec to pc
lum = 2d*!dpi*surf_lum*sigma_lum_pc^2

sigma_pot_pc = sigma_pot*pc     ; Convert from arcsec to pc
mass = 2d*!dpi*surf_pot*sigma_pot_pc^2

s = {SIG_L:sigma_lum_pc, SIG_M:sigma_pot_pc, LUM:lum, MASS:mass, MBH:mbh, $
    RMIN:0d, SURF_L:surf_lum, SIGMAPSF:sigmapsf, NORMPSF:normpsf, $
    PIXSIZE:pixSize, BETA:beta, STEP:step_pc, NRAD:nrad}

t = systime(1)

rmsModel = shpani_second_moment(rad*pc, s)

if ~keyword_set(quiet) then print, 'Elapsed time sec:', systime(1) - t

;###### Output and optional M/L fit
; If RMS keyword is not given all this section is skipped

if n_elements(rms) gt 0 then begin

    if n_elements(ml) eq 0 || ml le 0 then begin

        ; y1 = rms; dy1 = erms (y1 are the data, y2 the model)
        ; scale = total(y1*y2/dy1^2)/total(y2^2/dy1^2)  (equation 51)
        ;
        scale = total(rms*rmsModel/erms^2) / total((rmsModel/erms)^2)
        ml = scale^2

    endif else scale = sqrt(ml)

    rmsModel *= scale
    chi2 = total( ((rms-rmsModel)/erms)^2 ) / n_elements(rms)

    if ~keyword_set(quiet) then begin
        print, chi2, ml, mbh*ml, beta[0], $
            FORMAT='("chi^2/DOF: ", g0.3, "; M/L: ", g0.3, "; M_BH: ", e0.2, "; beta: ", f0.2)'
        print, 'Total mass MGE:', total(mass*ml)
    endif

    if keyword_set(plot) then begin
        str = '!7b!X=' + STRING(beta[0],FORMAT='(f0.2)') + $
              ' M/L=' + STRING(ml,FORMAT='(g0.3)') + $ ; Note multiplication of BH by M/L
              ' BH=' + STRING(mbh*ml,FORMAT='(e0.1)')  ; as calculation is done for M/L=1

        LOADCT, 12, /SILENT
        mx = max([rms+erms,rms-erms],MIN=mn)
        plot, rad>(0.38*pixSize1), rms, TITLE=str, /YNOZERO, /NODATA, $
            YRANGE=[mn,mx], XTITLE='R (arcsec)', YTITLE='!4r!X (km/s)', /XLOG
        oploterr, rad, rms, erms, 4
        oplot, rad, rmsModel, COLOR=200, THICK=3
        wait, 0.1
    endif

endif

END
;--------------------------------------------------------------------
pro test_jam_spherical_rms
;
; This example takes 7s on a 2GHz computer

; Realistic MGE galaxy surface brightness.
; The surface brightness is in L_sun/pc^2 and the sigma in arcsec
;
surf_pc = [6229., 3089., 5406., 8443., 4283., 1927., 708.8, 268.1, 96.83]
sigma_arcsec = [0.0374, 0.286, 0.969, 2.30, 4.95, 8.96, 17.3, 36.9, 128.]

; Realistic observed stellar kinematics. It comes from AO observations
; at R<2" and seeing-limited long slit observations at larger radii.
; The galaxy has negligible rotation and we can use sigma as V_RMS
;
sig = [395.,390.,387.,385.,380.,365.,350.,315.,310.,290.,260.] ; km/s
rad = [0.15, 0.2, 0.3, 0.4, 0.5, 1, 1.5, 3, 5, 9, 15] ; arcsec

; Realistic anisotropy profile from a Schwarzschild model.
; The anisotropy varies smoothly between the following three regimes:
; 1. beta = -1 for R < 1"
; 2. beta = 0.3 for 1" < R < 30"
; 3. beta = -0.2 for R > 30"
;
beta = sigma_arcsec*0
beta[where(sigma_arcsec le 1)] = -1.0
beta[where(sigma_arcsec gt 1 and sigma_arcsec le 30)] = 0.3
beta[where(sigma_arcsec gt 30)] = -0.2

; Compute V_RMS profiles and optimize M/L to best fit the data.
; Assume self-consistency: same MGE for luminosity and potential.
;
pixSize = 0.1   ; Spaxel size in arcsec
psf = 0.2/2.355 ; sigma of the PSF in arcsec from AO observations
mbh = 1.5e8 ; Black hole mass in solar masses before multiplication by M/L
distance = 20d    ; Mpc

jam_spherical_rms, surf_pc, sigma_arcsec, surf_pc, sigma_arcsec, mbh, distance, rad, sigp, $
    BETA=beta, SIGMAPSF=psf, PIXSIZE=pixSize, RMS=sig, /PLOT

END
;--------------------------------------------------------------------

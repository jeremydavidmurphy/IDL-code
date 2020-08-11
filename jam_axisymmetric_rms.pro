;##############################################################################
;
; Copyright (C) 2003-2010, Michele Cappellari
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
;##############################################################################
;+
; NAME:
;   JAM_AXISYMMETRIC_RMS
;
; PURPOSE:
;    This procedure calculates a prediction for the projected second velocity
;    moments V_RMS = sqrt(V^2 + sigma^2) for an anisotropic axisymmetric galaxy
;    model. It implements the solution of the anisotropic Jeans equations
;    presented in equation (28) of Cappellari (2008, MNRAS, 390, 71). PSF
;    convolution in done as described in the Appendix of that paper.
;    http://adsabs.harvard.edu/abs/2008MNRAS.390...71C
;
; CALLING SEQUENCE:
;    JAM_AXISYMMETRIC_RMS, $
;        surf_lum, sigma_lum, qObs_lum, surf_pot, sigma_pot, qObs_pot, $
;        inc_deg, mbh, distance, xbin, ybin, rmsModel, $
;        BETA=beta, CHI2=chi2, ERMS=erms, FLUX=flux, GOODBINS=goodBins, $
;        ML=ml, NORMPSF=normPsf, NANG=nang, NRAD=nrad, PIXANG=pixAng, $
;        PIXSIZE=pixSize, /PLOT, /QUIET, RBH=rbh, RMS=rms, $
;        SIGMAPSF=sigmaPsf, STEP=step
;
; INPUT PARAMETERS:
;   SURF_LUM: vector of length N containing the peak surface brightness of the
;       MGE Gaussians describing the galaxy surface brightness in units of
;       Lsun/pc^2 (solar luminosities per parsec^2).
;   SIGMA_LUM: vector of length N containing the dispersion in arcseconds of
;       the MGE Gaussians describing the galaxy surface brightness.
;   QOBS_LUM: vector of length N containing the observed axial ratio of the MGE
;       Gaussians describing the galaxy surface brightness.
;   SURF_POT: vector of length M containing the peak value of the MGE Gaussians
;       describing the galaxy surface density in units of Msun/pc^2 (solar
;       masses per parsec^2). This is the MGE model from which the model
;       potential is computed.
;     - In a common usage scenario, with a self-consistent model, one has
;       the same Gaussians for both the surface brightness and the potential.
;       This implies SURF_POT = SURF_LUM, SIGMA_POT = SIGMA_LUM and
;       QOBS_POT = QOBS_LUM. The global M/L of the model is fitted by the
;       routine when passing the RMS and ERMS keywords with the observed kinematics.
;   SIGMA_POT: vector of length M containing the dispersion in arcseconds of
;       the MGE Gaussians describing the galaxy surface density.
;   QOBS_POT: vector of length M containing the observed axial ratio of the MGE
;       Gaussians describing the galaxy surface density.
;   INC_DEG: inclination in degrees (90 being edge-on).
;   MBH: Mass of a nuclear supermassive black hole in solar masses.
;     - VERY IMPORTANT: The model predictions are computed assuming SURF_POT
;       gives the total mass. In the common self-consistent case one has
;       SURF_POT = SURF_LUM and if requested (keyword ML) the program can scale
;       the output RMSMODEL to best fit the data. The scaling is equivalent to
;       multiplying *both* SURF_POT and MBH by a factor M/L. To avoid mistakes,
;       the actual MBH used by the output model is printed on the screen.
;   DISTANCE: distance of the galaxy in Mpc.
;   XBIN: Vector of length P with the X coordinates in arcseconds of the bins
;       (or pixels) at which one wants to compute the model predictions. The
;       X-axis is assumed to coincide with the galaxy projected major axis. The
;       galaxy center is at (0,0).
;     - When no PSF/pixel convolution is performed (SIGMAPSF=0 or PIXSIZE=0)
;       there is a singularity at (0,0) which should be avoided by the input
;       coordinates.
;   YBIN: Vector of length P with the Y coordinates in arcseconds of the bins
;       (or pixels) at which one wants to compute the model predictions. The
;       Y-axis is assumed to concide with the projected galaxy symmetry axis.
;
; KEYWORDS:
;   BETA: Vector of length N with the anisotropy
;       beta_z = 1 - (sigma_z/sigma_R)^2 of the individual MGE Gaussians.
;       A scalar can be used if the model has constant anisotropy.
;   CHI2: Reduced chi^2 describing the quality of the fit
;        chi^2 = total( ((rms[goodBins]-rmsModel[goodBins])/erms[goodBins])^2 )
;              / n_elements(goodBins)
;   ERMS: Vector of length P with the 1sigma errors associated to the RMS
;        measurements. From the error propagation
;        ERMS = sqrt((dVel*velBin)^2 + (dSig*sigBin)^2)/RMS,
;        where velBin and sigBin are the velocity and dispersion in each bin
;        and dVel and dSig are the corresponding errors.
;        (Default: constant errors ERMS=0.05*MEDIAN(RMS))
;   FLUX: In output this contains a vector of length P with the unconvolved MGE
;       surface brightness of each bin, used to plot the isophotes on the model
;       results.
;   GOODBINS: Vector of length <=P with the indices of the bins which have to
;        be included in the fit (if requested) and chi^2 calculation.
;        (Default: fit all bins).
;   ML: Mass-to-light ratio to multiply the values given by SURF_POT.
;       Setting this keyword is completely equivalent to multiplying the
;       output RMSMODEL by SQRT(M/L) after the fit. This implies that the
;       BH mass becomes MBH*(M/L).
;     - If this keyword is set to a negative number in input, the M/L is
;       fitted from the data and the keyword returns the best-fitting M/L
;       in output. The BH mass of the best-fitting model is MBH*(M/L).
;   NORMPSF: Vector of length Q with the fraction of the total PSF flux
;       contained in the circular Gaussians describing the PSF of the
;       observations. It has to be total(NORMPSF) = 1. The PSF will be used
;       for seeing convolution of the model kinematics.
;   NRAD: Number of logarithmically spaced radial positions for which the
;       models is evaluated before interpolation and PSF convolution. One may
;       want to increase this value if the model has to be evaluated over many
;       orders of magnitutes in radius (default: NRAD=50). The computation time
;       scales as NRAD*NANG.
;   NANG: Same as for NRAD, but for the number of angular intervals
;       (default: NANG=10).
;   PIXANG: angle between the observed spaxels and the galaxy major axis X.
;   PIXSIZE: Size in arcseconds of the (square) spatial elements at which the
;       kinematics is obtained. This may correspond to the side of the spaxel
;       or lenslets of an integral-field spectrograph. This size is used to
;       compute the kernel for the seeing and aperture convolution.
;     - If this is not set, or PIXSIZE = 0, then convolution is not performed.
;   /PLOT: Set this keyword to produce a plot at the end of the calculation.
;   /QUIET: Set this keyword not to print values on the screen.
;   RBH: This scalar gives the sigma in arcsec of the Gaussian representing the
;       central black hole of mass MBH (See Section 3.1.2 of Cappellari 2008).
;       The gravitational potential is indistinguishable from a point source
;       for radii > 2*RBH, so the default RBH=0.01 arcsec is appropriate in
;       most current situations.
;     - RBH should not be decreased unless actually needed!
;   RMS: Vector of length P with the input observed stellar
;       V_RMS=sqrt(velBin^2 + sigBin^2) at the coordinates positions given by
;       the vectors XBIN and YBIN.
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
;
; OUTPUT PARAMETER:
;   RMSMODEL: Vector of length P with the model predictions for the velocity
;       second moments V_RMS ~ sqrt(vel^2 + sig^2) for each bin.
;
; USAGE EXAMPLE:
;    A simple usage example is given in the procedure TEST_JAM_AXISYMMETRIC_RMS
;    at the end of this file.
;
; REQUIRED ROUTINES:
;       By M. Cappellari (included in the JAM distribution):
;       - ANY
;       - DIFF
;       - MESHGRID
;       - QUADVA
;       - RANGE
;       - ROTATE_POINTS
;       - SYMMETRIZE_VELFIELD
;
;       By M. Cappellari available at
;       http://www-astro.physics.ox.ac.uk/~mxc/idl/#jam
;       - PLOT_VELFIELD
;       - SAURON_COLORMAP
;
; MODIFICATION HISTORY:
; V1.0: Written and tested by Michele Cappellari, Vicenza, 19 November 2003
; V2.0: Introduced new solution of the MGE Jeans equations with constant
;   anisotropy sig_R = b*sig_z. MC, Oxford, 20 September 2007
; V3.13: First released version. MC, Oxford, 12 August 2008
; V3.2: Updated documentation. MC, Oxford, 14 August 2008
; V4.0: Implemented PSF convolution using interpolation on polar grid.
;     Dramatic speed-up of calculation. Further documentation.
;     MC, Oxford, 11 September 2008
; V4.01: Bug fix: when ERMS was not given, the default was not properly set.
;     Included keyword STEP. The keyword FLUX is now only used for output:
;     the surface brightness for plotting is computed from the MGE model.
;     MC, Windhoek, 29 September 2008
; V4.02: Added keywords NRAD and NANG. Thanks to Michael Williams for
;     reporting possible problems with too coarse interpolation.
;     MC, Oxford, 21 November 2008
; V4.03: Added keyword RBH. MC, Oxford 4 April 2009
; V4.04: Compute FLUX even when not plotting. MC, Oxford, 29 May 2009
; V4.05: Skip unnecessary interpolation when computing few points
;     without PSF convolution. After feedback from Eric Emsellem.
;     MC, Oxford, 6 July 2009
; V4.06: Updated documentation. The routine TEST_JAM_AXISYMMETRIC_RMS with
;     the usage example now adopts a more realistic input kinematics.
;     MC, Oxford, 08 February 2010
; V4.07: Forces q_lum && q_pot < 1. MC, Oxford, 01 March 2010
;-
;#############################################################################
function janis2_jeans_mge_integrand, u, $
    DENS_lum=dens_lum, SIGMA_lum=sigma_lum, Q_lum=q_lum, $
    DENS_pot=dens_pot, SIGMA_pot=sigma_pot, Q_pot=q_pot, $
    X1=x1, Y1=y1, INC=inc, BETA=beta
compile_opt idl2, hidden
;
; This routine computes the integrand of Eq.(28) of Cappellari (2008)
; for a model with constant anisotropy sigma_R^2 = b*sigma_z^2 and <V_R*V_z> = 0.

kani = 1d/(1d - beta) ; Anisotropy ratio b = (sig_R/sig_z)^2
si2 = sin(inc)^2
ci2 = cos(inc)^2
x2 = x1^2
y2 = y1^2
u2 = u^2

s2_lum = sigma_lum^2
q2_lum = q_lum^2
e2_lum = 1d - q2_lum
s2q2_lum = s2_lum*q2_lum

s2_pot = sigma_pot^2
e2_pot = 1d - q_pot^2

; Double summation over (j,k) of eq.(28) vectorized over integration variable u.
; The j-index refers to the Gaussians describing the total mass,
; from which the potential is derived, while the k-index is used
; for the MGE components describing the galaxy stellar luminosity.
;
sum = 0d
for j=0,n_elements(dens_pot)-1 do begin ; loop over mass Gaussians
    e2u2_pot = e2_pot[j]*u2
    for k=0,n_elements(dens_lum)-1 do begin ; loop over luminous Gaussians
        a = 0.5d*(u2/s2_pot[j] + 1d/s2_lum[k]) ; equation (29)
        b = 0.5d*(e2u2_pot*u2/(s2_pot[j]*(1d - e2u2_pot)) + e2_lum[k]/s2q2_lum[k]) ; equation (30)
        c = e2_pot[j] - s2q2_lum[k]/s2_pot[j] ; equation (22)
        d = 1d - kani[k]*q2_lum[k] - ((1d - kani[k])*c + e2_pot[j]*kani[k])*u2 ; equation (23)
        e = a + b*ci2
        sum += dens_lum[k]*q_pot[j]*dens_pot[j]*u2 $  ; sum vector has the size of u
            * (s2q2_lum[k]*(ci2 + kani[k]*si2) + x2*si2*d) $
            * exp(-a*(x2 + y2*(a + b)/e)) / ((1d - c*u2)*sqrt((1d - e2u2_pot)*e))
    endfor
endfor

G = 0.00430237d    ; (km/s)^2 pc/Msun [6.674e-11 SI units]

return, 4d*!dpi^1.5d*G*sum
end
;----------------------------------------------------------------------
function janis2_weighted_second_moment_squared, x, y, inc_deg, s
compile_opt idl2, hidden
;
; This routine gives the projected non-centered
; weighted second moment squared \Sigma*<V_los^2>

; Axisymmetric deprojection of both luminous and total mass.
; See equation (12)-(14) of Cappellari (2008)
;
inc = inc_deg/!radeg

qintr_lum = s.qobs_lum^2 - cos(inc)^2
if any(qintr_lum le 0d) then message, 'Inclination too low q < 0'
qintr_lum = sqrt(qintr_lum)/sin(inc)
if any(qintr_lum lt 0.05d) then message, 'q < 0.05 components'
dens_lum = s.surf_lum*s.qobs_lum / (s.sigma_lum*qintr_lum*sqrt(2d*!dpi))

qintr_pot = s.qobs_pot^2 - cos(inc)^2
if any(qintr_pot le 0d) then message, 'Inclination too low q < 0'
qintr_pot = sqrt(qintr_pot)/sin(inc)
if any(qintr_pot lt 0.05d) then message, 'q < 0.05 components'
dens_pot = s.surf_pot*s.qobs_pot / (s.sigma_pot*qintr_pot*sqrt(2d*!dpi))

functargs = {DENS_lum:dens_lum, SIGMA_lum:s.sigma_lum, Q_lum:qintr_lum, $
             DENS_pot:dens_pot, SIGMA_pot:s.sigma_pot, Q_pot:qintr_pot, $
             X1:x, Y1:y, INC:inc, BETA:s.beta}
quadva, 'janis2_jeans_mge_integrand', [0d,1d], sb_mu2, RELTOL=1d-5, ABSTOL=0d, FUNCTARGS=functargs

return, sb_mu2
end
;----------------------------------------------------------------------
function janis2_second_moment, x, y, inc_deg, s
compile_opt idl2, hidden
;
; This routine gives the second V moment after convolution with a PSF.
; The convolution is done using interpolation of the model on a
; polar grid, as described in Appendix A of Cappellari (2008).

; Define parameters of polar grid for interpolation
;
w = where(s.sigma_lum lt max(x),m) ; Characteristic MGE axial ratio in observed range
if m lt 3 then qmed = median(s.qobs_lum) else qmed = median(s.qobs_lum[w])
rell = sqrt(x^2 + (y/qmed)^2) ; Elliptical radius of input (x,y)

; Kernel step is 1/4 of largest value between sigma(min) and 1/2 pixel side.
; Kernel half size is the sum of 3*sigma(max) and 1/2 pixel diagonal.
;
if max(s.sigmaPsf) gt 0 && s.pixSize gt 0 then begin ; PSF convolution
    if s.step gt 0 then step = s.step $
        else step = (s.pixSize/2d > min(s.sigmaPsf))/4d
    mx = 3d*max(s.sigmaPsf) + s.pixSize/sqrt(2d)
endif else begin ; No convolution
    nbins = n_elements(x)
    if s.NRAD*s.NANG ge nbins then begin
        mu = x*0
        for j=0,nbins-1 do begin
            wm2 = janis2_weighted_second_moment_squared(x[j], y[j], inc_deg, s)
            mge = total(s.surf_lum * exp(-0.5d/s.sigma_lum^2 * (x[j]^2 + (y[j]/s.qobs_lum)^2)))
            mu[j] = sqrt(wm2/mge)
        endfor
        return, mu ; Skip the interpolation when computing just a few points
    endif
    step = min(rell>1d) ; Minimum radius of 1pc
    mx = 0d
endelse

; Make linear grid in log of elliptical radius RAD and eccentric anomaly ANG
; See Appendix A
;
rmax = max(rell) + mx ; Major axis of ellipse containing all data + convolution
nrad = s.NRAD ; radial elements
nAng = s.NANG ; angular sectors
logRad = range(alog(step),alog(rmax),nrad) ; Linear grid in log(rell)
ang = range(0d,!dpi/2d,nAng) ; Linear grid in eccentric anomaly
meshgrid, exp(logRad), ang, radGrid, angGrid
xPol = radGrid*cos(angGrid)
yPol = radGrid*sin(angGrid) * qmed

; The model Vrms computation is only performed on the polar grid
; which is then used to interpolate the values at any other location
;
wm2Pol = xPol*0
mgePol = wm2Pol
for j=0,n_elements(xPol)-1 do begin
    wm2Pol[j] = janis2_weighted_second_moment_squared(xPol[j], yPol[j], inc_deg, s)
    mgePol[j] = total(s.surf_lum * exp(-0.5d/s.sigma_lum^2 * (xPol[j]^2 + (yPol[j]/s.qobs_lum)^2)))
endfor

if max(s.sigmaPsf) gt 0 && s.pixSize gt 0 then begin ; PSF convolution

    nx = ceil(rmax/step/64d)*64       ; Make 2*nx dimensions divisible by 2^7
    ny = ceil(rmax*qmed/step/64d)*64  ; for much faster FFT computation
    x1 = range(-nx,nx,2*nx)*step
    y1 = range(-ny,ny,2*ny)*step
    meshgrid, x1, y1, xCar, yCar ; Cartesian grid for convolution

    ; Interpolate MGE model and Vrms over cartesian grid
    ;
    r1 = 0.5d*alog(xCar^2 + (yCar/qmed)^2) ; Log elliptical radius of cartesian grid
    e1 = atan(abs(yCar/qmed),abs(xCar))    ; Eccentric anomaly of cartesian grid
    indx = (r1 - logRad[0])/(logRad[1] - logRad[0]) ; Convert to indices for INTERPOLATE
    indy = (e1 - ang[0])/(ang[1] - ang[0])
    wm2Car = interpolate(wm2Pol,indx,indy,CUBIC=-0.5)
    mgeCar = interpolate(mgePol,indx,indy,CUBIC=-0.5)

    nk = ceil(mx/step)
    kgrid = range(-nk,nk,2*nk)*step
    meshgrid, kgrid, kgrid, xx, yy ; Kernel is square
    rotate_points, xx, yy, s.pixAng, xgrid, ygrid

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
    kernel /= total(kernel)

    ; Seeing and aperture convolution with equation (A3)
    ;
    tmp = convolve(wm2Car,kernel)/convolve(mgeCar,kernel)
    tmp >= min(tmp[where(tmp gt 0)]) ; Fix possible rounding errors
    muCar = sqrt(tmp)

    ; Interpolate convolved image at observed apertures.
    ; Aperture integration was already included in the kernel.
    ;
    indx = (x - x1[0])/(x1[1] - x1[0]) ; Convert to indices for INTERPOLATE
    indy = (y - y1[0])/(y1[1] - y1[0])
    mu = interpolate(muCar,indx,indy)

endif else begin ; No PSF convolution: just interpolate values

    tmp = wm2Pol/mgePol
    tmp >= min(tmp[where(tmp gt 0)]) ; Fix possible rounding errors
    muPol = sqrt(tmp)
    r1 = 0.5d*alog(x^2 + (y/qmed)^2) ; Log elliptical radius of input (x,y)
    e1 = atan(abs(y/qmed),abs(x))    ; Eccentric anomaly of input (x,y)
    indx = (r1 - logRad[0])/(logRad[1] - logRad[0]) ; Convert to indices for INTERPOLATE
    indy = (e1 - ang[0])/(ang[1] - ang[0])
    mu = interpolate(muPol,indx,indy,CUBIC=-0.5)

endelse

return, mu
end
;----------------------------------------------------------------------
pro jam_axisymmetric_rms, $
    surf_lum, sigma_lum, qobs_lum, surf_pot, sigma_pot, qobs_pot1, $
    inc, mbh, distance, xbin, ybin, rmsModel, $
    ML=ml, NORMPSF=normPsf, PIXANG=pixAng, PIXSIZE=pixSize, PLOT=plot, $
    RMS=rms, ERMS=erms, CHI2=chi2, SIGMAPSF=sigmaPsf, GOODBINS=goodBins, $
    QUIET=quiet, FLUX=flux, BETA=beta1, STEP=step, NRAD=nrad, NANG=nang, $
    RBH=rbh
compile_opt idl2
on_error, 2

if n_elements(beta1) eq 0 then beta = sigma_lum*0 ; Anisotropy parameter beta = 1 - (sig_z/sig_R)^2
if n_elements(beta1) eq 1 then beta = sigma_lum*0+beta1 ; All components have the same anisotropy
if n_elements(beta1) eq n_elements(surf_lum) then beta = beta1

if n_elements(sigmaPsf) eq 0 then sigmaPsf = 0d
if n_elements(normPsf) eq 0 then normPsf = 1d
if n_elements(pixAng) eq 0 then pixAng = 0d
if n_elements(pixSize) eq 0 then pixSize = 0d
if n_elements(erms) eq 0 && n_elements(rms) gt 0 then erms = rms*0+median(rms)*0.05 ; Constant 5% errors
if n_elements(step) eq 0 then step = 0
if n_elements(nrad) eq 0 then nrad = 50
if n_elements(nang) eq 0 then nang = 10
if n_elements(rbh) eq 0 then rbh = 0.01d ; arcsec

pc = distance*!dpi/0.648d ; Constant factor to convert arcsec --> pc

surf_lum_pc = surf_lum
surf_pot_pc = surf_pot
sigma_lum_pc = sigma_lum*pc         ; Convert from arcsec to pc
sigma_pot_pc = sigma_pot*pc         ; Convert from arcsec to pc
xbin_pc = xbin*pc                   ; Convert all distances to pc
ybin_pc = ybin*pc
pixSize_pc = pixSize*pc
sigmaPsf_pc = sigmaPsf*pc
step_pc = step*pc

; Add a Gaussian with small sigma and the same total mass as the BH.
; The Gaussian provides an excellent representation of the second moments
; of a point-like mass, to 1% accuracy out to a radius 2*sigmaBH.
; The error increses to 14% at 1*sigmaBH; Independently of the BH mass.
;
if mbh gt 0 then begin
    sigmaBH_pc = rbh*pc ; Adopt for the BH just a very small size
    surfBH_pc = mbh/(2d*!dpi*sigmaBH_pc^2)
    surf_pot_pc = [surfBH_pc,surf_pot_pc[*]] ; Add Gaussian to potential only!
    sigma_pot_pc = [sigmaBH_pc,sigma_pot_pc[*]]
    qobs_pot = [1d,qobs_pot1[*]]  ; Make sure vectors do not have extra dimensions
endif else qobs_pot = qobs_pot1  ; Important: do not change input vector qobs_pot!

s = {SURF_LUM:surf_lum_pc, SIGMA_LUM:sigma_lum_pc, QOBS_LUM:qobs_lum<0.99d, $
     SURF_POT:surf_pot_pc, SIGMA_POT:sigma_pot_pc, QOBS_POT:qobs_pot<0.99d, $
     SIGMAPSF:sigmaPsf_pc, NORMPSF:normPsf, STEP:step_pc, BETA:beta, $
     PIXSIZE:pixSize_pc, PIXANG:pixAng, NRAD:nrad, NANG:nang}

t = systime(1)

rmsModel = janis2_second_moment(xbin_pc, ybin_pc, inc, s)

if ~keyword_set(quiet) then print, 'Elapsed time sec:', systime(1) - t

flux = 0d ; Total MGE surface brightness for plotting
for j=0,n_elements(surf_lum_pc)-1 do $
    flux += surf_lum_pc[j]*exp(-0.5d/sigma_lum[j]^2*(xbin^2 + (ybin/qObs_lum[j])^2))

;###### Output and optional M/L fit
; If RMS keyword is not given all this section is skipped

if n_elements(rms) gt 0 then begin

    ; Only consider the good bins for the chi^2 estimation
    ;
    if n_elements(goodBins) eq 0 then goodBins = indgen(n_elements(xbin))

    if n_elements(ml) eq 0 || ml le 0 then begin

        ; y1 = rms; dy1 = erms (y1 are the data, y2 the model)
        ; scale = total(y1*y2/dy1^2)/total(y2^2/dy1^2)  (equation 51)
        ;
        scale = total(rms[goodBins]*rmsModel[goodBins]/erms[goodBins]^2) $
              / total((rmsModel[goodBins]/erms[goodBins])^2)
        ml = scale^2

    endif else scale = sqrt(ml)

    rmsModel *= scale
    chi2 = total( ((rms[goodBins]-rmsModel[goodBins])/erms[goodBins])^2 ) / n_elements(goodBins)

    if ~keyword_set(quiet) then begin
        print, inc, chi2, ml, mbh*ml, beta[0], $
            FORMAT='("inc: ", f0.1, "; chi^2/DOF: ", g0.3, "; M/L: ", g0.3, "; M_BH: ", e0.2, "; beta: ", f0.2)'
        mass = 2d*!dpi*surf_pot*qobs_pot1*sigma_pot_pc^2
        print, 'Total mass MGE:', total(mass*ml)
    endif

    if keyword_set(plot) then begin
        str = 'i=' + STRING(inc,FORMAT='(f0.1)') + $
              ' !7b!X!Dz!N=' + STRING(beta[0],FORMAT='(f0.2)') + $
              ' M/L=' + STRING(ml,FORMAT='(g0.3)') + $ ; Note multiplication of BH by M/L
              ' BH=' + STRING(mbh*ml,FORMAT='(e0.1)') ; as calculation is done for M/L=1

        ; The X-axis was aligned with the major axis as input
        ;
        dx = randomn(seed,n_elements(xbin))*0.001 ; Avoids co-linear points IDL bug in TRIGRID
        symmetrize_velfield, xbin[goodBins]+dx, ybin[goodBins]+dx, rms[goodBins], tmp, SYM=2, PA=90d
        rms1 = rms
        rms1[goodBins] = tmp ; Only symmetrize good bins

        !p.multi = [0,1,2]
        sauron_colormap
        mx = max(sigrange(rms1[goodBins]),MIN=mn)
        range = [mn,mx]
        plot_velfield, xbin, ybin, rms1, RANGE=range, FLUX=flux

        ; Overplot bad bins on the data
        ;
        tmp = xbin*0+1
        tmp[goodBins] = 0
        badBins = where(tmp,m)
        if m gt 0 then begin
            oplot, xbin[badBins], ybin[badBins], PSYM=4, SYMSIZE=0.5, COLOR=0, THICK=1
            oplot, xbin[badBins], ybin[badBins], PSYM=4, SYMSIZE=0.25, COLOR=255, THICK=1
        endif
        plot_velfield, xbin, ybin, rmsModel, RANGE=range, TITLE=str, FLUX=flux
        !p.multi = 0
        wait, 0.1
    endif

endif

end
;----------------------------------------------------------------------
pro test_jam_axisymmetric_rms
;
; This example takes 5.5s on a 2GHz computer

; Assume as input a simple but realistic toy model for the
; integral-field stellar velocity and velocity dispersion.
; In a real case the user has to provide four vectors
; with the observed values for XBIN, YBIN, VEL, and SIG.
;
inc = 60d ; assumed galaxy inclination
xbin = randomu(s,1000)*100-50 ; Random X sky coordinates in [-50,50] arcsec
ybin = randomu(s,1000)*100-50 ; Random Y sky coordinates in [-50,50] arcsec
r = sqrt(xbin^2 + (ybin/cos(inc*!dtor))^2) ; Radius in the plane of the disk
a = 40d                          ; Scale length in arcsec
vr = 2000*sqrt(r)/(r+a)          ; Assumed velocity profile
vel = vr * sin(inc*!dtor)*xbin/r ; Projected velocity field
sig = 8700/(r+a)                 ; Assumed velocity dispersion field
rms = sqrt(vel^2 + sig^2)        ; Vrms field in km/s

; Realistic MGE galaxy surface brightness
; surf in LSun/pc^2; sigma in arcsec
;
surf = [39483, 37158, 30646, 17759, 5955.1, 1203.5, 174.36, 21.105, 2.3599, 0.25493]
sigma = [0.153, 0.515, 1.58, 4.22, 10, 22.4, 48.8, 105, 227, 525]
qObs = [0.57, 0.57, 0.57, 0.57, 0.57, 0.57, 0.57, 0.57, 0.57, 0.57]

; Assume self-consistency: same MGE for both light and mass
;
surf_lum = surf
sigma_lum = sigma
qobs_lum = qObs
surf_pot = surf
sigma_pot = sigma
qobs_pot = qObs

distance = 16.5   ; Assume Virgo distance in Mpc (Mei et al. 2007)
mbh = 1e8 ; Black hole mass in solar masses

; The model is similar but not the same as the adopted kinematics!

jam_axisymmetric_rms, $
    surf_lum, sigma_lum, qobs_lum, surf_pot, sigma_pot, qobs_pot, $
    inc, mbh, distance, xbin, ybin, rmsModel, $
    BETA=0.3, RMS=rms, /PLOT, SIGMAPSF=0.6, PIXSIZE=0.8

end
;----------------------------------------------------------------------

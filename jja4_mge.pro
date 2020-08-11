;A wrapper to call jam with DM halo parameters
;This code produces an alpha vs. chisq curve with all other 
;variables marginalized over.
;----------------------------------------------------------------------------
FUNCTION wrap_jam, p_in, XVAL=x_in, YVAL=y_in, ERRVAL=err_in, $
SUR_LUM=surf_lum,SIG_LUM=sigma_lum,QOB_LUM=qobs_lum,RAVAL=ra_in, $
DECVAL=dec_in, PC_IN=pc_in, DIST_IN=distance
compile_opt idl2, hidden
;create MGE model of DM halo
n_pl = 1000
r_pl = range(1d-1,5000d,n_pl,/LOG) ; in pc
rho_pl = p_in[3]*(r_pl)^(-p_in[4]) ; in Msun/pc^3
mge_fit_1d, r_pl, rho_pl, NGAUSS=16, SOL=sol_pl
tmp_var1 = sol_pl[0,*]
tmp_var2 = sol_pl[1,*]/pc_in
tmp_var3 = tmp_var1*0+1
size_tmp = Size(tmp_var1, /DIMENSIONS)
surf_pot=[surf_lum*p_in[0],reform(tmp_var1,size_tmp[1])]
sigma_pot=[sigma_lum,reform(tmp_var2,size_tmp[1])]
qobs_pot=[qobs_lum,reform(tmp_var3,size_tmp[1])]
print, "Formed"
print, p_in
jam_axisymmetric_rms, surf_lum, sigma_lum, qobs_lum, surf_pot, $
   sigma_pot, qobs_pot, p_in[1], 0d0, $
   distance, ra_in, dec_in, jam_model, ML=1., BETA=p_in[2], $
   SIGMAPSF=0.8, PIXSIZE=16d, PIXANG=-53d, /QUIET
return, (y_in-jam_model)/err_in
end
;----------------------------------------------------------------------------
PRO jja4_mge
compile_opt idl2
;read x_in, y_in, err_in, ra_in, dec_in
nlines=file_lines('mindata.txt')
ra_in=dblarr(nlines)
dec_in=dblarr(nlines)
y_in=dblarr(nlines)
err_in=dblarr(nlines)
x_in=dindgen(nlines)
covar_out=dblarr(5,5)
readcol, 'mindata.txt', ra_in, dec_in, y_in, err_in

lum_R = 2.95E8 ; My numerical integration for the MGE fit in R band L_sol within 1.5 kpc
distance = 3.2d ; Mpc
pc = distance*!dpi/0.648d ; Constant factor to convert arcsec --> pc
surf_lum = [2298.3470,68.851126,27.769026,117.76896,19.550785]
sigma_lum = [0.280440,0.986237,6.35659,49.8893,108.260]
qobs_lum = [0.990000,0.990000,0.990000,0.480000,0.542514]
;ML_0 = 0.662122*3.
nbins=13
xbin=DBLARR(nbins)
ybin=DBLARR(nbins)
stelM_bin=DBLARR(nbins)
DMM_bin=DBLARR(nbins)
for ix=0, nbins-1, 1 do begin
   xbin[ix]=double(ix)*1.3d/double(nbins-1)+0.1d
ML_0 = 1.0
inc_0 = 66.5900
beta_0 = 0.400000
rhoDM_0 = 0.1
alphaDM_0 = 0.5
;putting limits on fit
parinfo = replicate({limited:[0,0],limits:[0d0,0d0],step:0d,fixed:0}, 5)
parinfo[0].limited[0] = 1
parinfo[0].limited[1] = 1
parinfo[0].limits[0] = 5d-4
parinfo[0].limits[1] = 1d2
parinfo[0].step = 1d-3
parinfo[0].fixed = 0
parinfo[1].limited[0] = 1
parinfo[1].limited[1] = 1
parinfo[1].limits[0] = 61.45d0
parinfo[1].limits[1] = 80d0
parinfo[1].step = 1d-3
parinfo[1].fixed = 0
parinfo[2].limited[0] = 1
parinfo[2].limited[1] = 1
parinfo[2].limits[0] = -0.5d0
parinfo[2].limits[1] = 1.0d0
parinfo[2].step = 1d-2
parinfo[2].fixed = 0
parinfo[3].limited[0] = 1
parinfo[3].limited[1] = 1
parinfo[3].limits[0] = 1d-8
parinfo[3].limits[1] = 5d2
parinfo[3].step = 1d-2
parinfo[3].fixed = 0
parinfo[4].limited[0] = 1
parinfo[4].limited[1] = 1
parinfo[4].limits[0] = 0d0
parinfo[4].limits[1] = 2d0
parinfo[4].step = 1d-3
parinfo[4].fixed = 1
;p0 = [ML_0,inc_0,beta_0,rhoDM_0,alphaDM_0]
p0 = [ML_0,inc_0,beta_0,rhoDM_0,xbin[ix]]
fa = {XVAL:x_in, YVAL:y_in, ERRVAL:err_in,SUR_LUM:surf_lum,SIG_LUM:sigma_lum,$
QOB_LUM:qobs_lum,RAVAL:ra_in,DECVAL:dec_in,PC_IN:pc,DIST_IN:distance}
fin_p = mpfit('wrap_jam', p0, functargs=fa, STATUS=status_out, COVAR=covar_out,$
   BESTNORM=chisq_out, DOF=DOF_out,PARINFO=parinfo)
   ybin[ix]=chisq_out
   stelM_bin[ix]=fin_p[0]*lum_R ; in solar masses
   DMM_bin[ix]=4d0*!dpi*fin_p[3]*(1.5d3)^(3d0-fin_p[4])/(3d0-fin_p[4]) ; in solar masses
;print, "Status"
;print, status_out
;print, chisq_out, DOF_out
;print, fin_p
;print, covar_out
endfor
;Also now printing ratio of DM mass to stellar mass within 1.5kpc
openw, 28, "marg_chisq.txt"
for iy=0, nbins-1, 1 do begin
   printf, 28, xbin[iy], ybin[iy], stelM_bin[iy], DMM_bin[iy], format='(f12.5,1x,f12.5,1x,f12.5,1x,f12.5)'
endfor
END

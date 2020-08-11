specfile = 'speclist'
errfile = 'errlist'
zspecfile = 'zspeclist'
velfile = 'moments_w_errorsMOD.txt'

z = 0.00426d
wave0 = 3545d
dlambda = 1.125d

readcol, specfile, specname, format='a'
readcol, errfile, errname, format='a'
nspecs = N_elements(specname)
readcol, zspecfile, zspecname, format='a'
readcol, velfile, bin, vel, vel_err, vdisp, vdisp_err, $
         format='a,f,f,f,f'

for i=0, nspecs - 1 do begin

   spec = readfits(specname[i], header)
   err = readfits(errname[i])

   newheader = header

   Sxaddpar, newheader, 'CRVAL1', wave0 / (1. + (z + vel[i]/3e5))
   Sxaddpar, newheader, 'CD1_1', dlambda / (1. + (z + vel[i]/3e5))

   mwrfits, spec, zspecname[i], newheader
   mwrfits, err, zspecname[i], newheader

endfor


end

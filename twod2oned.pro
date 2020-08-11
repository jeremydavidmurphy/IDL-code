PRO twod2oned, list

;this routine accepts the data cubes from Sudhir and converts them to
;1D spectra. The initial wavelength solution is log. This routine
;converts it to a linear solution. It uses the log2lin.pro routine.

;NOTE: Sudhir's data comes in 3D fits arrays. They are
;wavelength, the different radial bins; the 3rd dimesion is data,
;noise and wavelength.

wi = 3550.0 ;the linear wavelength 0th order term
wd = 0.375 ;the linear wavelength dispersion term
trim = [50,5700] ;what the spectra gets trimmed to
add = ['a01','a02','a03','a04','a05','b01','b02','b03','b04','b05']
;*******************************************

witrim = wi+trim[0]*wd

readcol,list,f='a',files

n0 = n_elements(files)

for j=0,n0-1 do begin
   file = files[j]
   data = readfits(file,h)
;   logwave1 = sxpar(h,'CRVAL1')
;   logwavedisp = sxpar(h,'CDELT1')
;   logwave1 = 3541.11525008
   logwave1 = double(data[0,0,2])
   logwavedisp = 3.45879379893e-5
   s = size(data)
   if j eq 0 then logwave = logwave1 * 10^(logwavedisp * dindgen(s[1]))
   ss = strsplit(file,'.',/extract)
   ss = ss[0]
   sxaddpar,h,'CRVAL1',witrim,' The linear, zeroth-order term'
   sxaddpar,h,'CDELT1',wd,' The linear, dispersion term'
   sxaddpar,h,'CTYPE1','LINEAR',' the type of transform'
   for k=0,s[2]-1 do begin
      out = log2lin(data[*,k,0],logwave1,logwavedisp,wi,wd)
      if (s[0] eq 3) then outN = log2lin(data[*,k,1],logwave1,logwavedisp,wi,wd)
      out = out[trim[0]:trim[1]]
      if (s[0] eq 3) then outN = outN[trim[0]:trim[1]]
;      ssout = ss+add[j]+'.fits'
      ssout = ss+'_'+strn(k+1)+'.fits'
      if (s[0] eq 3) then ssoutN = ss+'_'+strn(k+1)+'N.fits'
      writefits,ssout,out,h
      if (s[0] eq 3) then writefits,ssoutN,outN,h
      wout = findgen(n_elements(out))*wd + witrim
      loadct,0,/silent
      plot,wout,out,xrange=[3500,5400]
      wait,0.1
      loadct,4,/silent
      oplot,data[*,k,2],data[*,k,0],color=150
;      oplot,logwave,data[*,k]
      wait,0.2
   endfor
endfor

stop

END

; This is a different routine than plotrms_v1.pro. This routine works
; from the output of filter.pro. It also requires a list of the bin
; names. 

;----------------------------------------------------------------
pro Prms
;----------------------------------------------------------------

;This routine requires three files to exist in the calling directory:
;PFITLOV.ALL, RUN.FINAL, and BIN.LIST

;The run.final file isn't used for anything. It just comes
;along for the ride so as to get sorted in the same fashion that
;pfitlov.all is. THE OUTPUT OF THIS ROUTINE IS PFITLOV.ALL.STD AND
;RUN.FINAL.STD.

readcol,silent=1,'pfitlov.all',format='a,f,f,f,f,f,f,f,i,f',$
        pfile,vel,sig1,sig2,h3,h4,r1,r2,r3,rms
n1 = n_elements(pfile)

readcol,silent=1,'run.final',format='a,a,i,i,f,i,i,a',$
        rfit,rfile,w1,w2,z,sigg,smooth,region

readcol,silent=1,'bin.list',format='a',bins

if (n_elements(pfile) ne n_elements(rfile)) then stop
n0 = n_elements(bins)

;the pfitlov.all names are split up
pfilearr = strarr(3,n1)
for j=0,n1-1 do begin
   temp = strsplit(pfile,'_.',/extract)
   pfilearr[*,j] = temp[0:2]
endfor

; a loop through each bin
indy = intarr(1)
for j=0,n0-1 do begin
   bin = bins[j]
   temp = where(pfilearr[0,*] eq bin)
   indy = [indy,temp]
endfor
indy = indy[1:*]

stop
end

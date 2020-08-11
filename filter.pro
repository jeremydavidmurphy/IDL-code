; The purpose of this code is to hone the kinematic extractions
; through alternations of the wavelength range of a specific
; wavelength region.

;This routine is used in an iterative manner to hone the
;kinematics. The routine reads in PFITLOV.ALL, PFITLOV.OUT RUN0.FINAL and RUN0.
;Pfitlov.out is then filtered for values above the specified
;cutoff. Any values above this cutoff are written out into a new run
;file (called run1, unless the keyword is used)

;The rms values are also fed in and attached to the pfitlov.all file
;in the last column.

;************************************************************************

; If a pfitlov.all and run.final frame exists, they get appended
; to. If not, they are created. If the keywords are not used then the
; following files are looked for.
; PFITOUT = pfitlov.out
; RZERO = run0
; RMS = rmsall.out
; ROUT = run1


;************************************************************************

PRO filter, PFA=pfitall, RA=runall,$
            PFO=pfitout, RZERO=run0, RMS=rmsall, ROUT=rout

;PFA=pfitall: USE THIS WHEN YOU ARE RUNNING LATER ITERATIONS. This is
;critical; if you run this a second time without using this keyword
;then the pfitlov.all file is overwritten. (As a safety feature, a
;file named pfitlov.all.first is also created.)

;RA=runall: Same as above. This must be used in order to avoid
;overwriting the run.final frame that's getting created. Again, a
;"run.final.first" frame is created on the first run.

;PFO=pfitout: This only gets used if your new pfitlov.out file is
;named something other than "pfitlov.out". If this keyword isn't used
;then the code will search for and use "pfitlov.out"

;RZERO=run0: Only used if your run file is named something other than
;"run0"

;RMSALL=rmsall: Same as above. Without this keyword "rmsall.out" will
;be searched for.

;ROUT=rout; Use this keyword to name the output run file something
;other than the default, "run1"
;-------------------------------------------------------------------------
cutoffu = 560.0
cutoffl = 200.0
;-------------------------------------------------------------------------

;First, the various files are searched for, then read in.

if (n_elements(pfitout) eq 0) then pfitout = 'pfitlov.out'
if (n_elements(run0) eq 0) then run0 = 'run0'
if (n_elements(rmsall) eq 0) then rmsall = 'rmsall.out'
readcol,pfitout,silent=1,format='a,f,f,f,f,f,f,f,f',nameO,velO,sig1O,sig2O,h3O,h4O,r1O,r2O,r3O
readcol,run0,silent=1,format='a,a,i,i,f,i,i,a',ro1,nameo1,wo1,wo2,zo,so,smo,subo
readcol,rmsall,silent=1,format='d,x,x,x',rmso
n1 = n_elements(sig1o)

;to check for rfitlovc's in the run list and get rid of them...
i = where(ro1 eq 'rfitlovc')
if (i[0] ne -1) then begin
    i = where(ro1 eq 'rfitlov')
    ro1 = ro1[i]
    nameo1 = nameo1[i]
    wo1 = wo1[i]
    wo2 = wo2[i]
    zo = zo[i]
    so = so[i]
    smo = smo[i]
    subo = subo[i]
endif
    
;a check to see that the files are the correct size
if (n_elements(nameo) ne n_elements(ro1) or n_elements(nameo) ne n_elements(rmso)) then begin
    print,''
    print,'Your '+pfitout+' and '+runo+' files are not the same length!'
    stop
endif

fp = '(A27,2x,F10.3,2x,F8.3,2x,F8.3,2x,F6.3,2x,F6.3,2x,F8.3,2x,F8.3,2x,I5,2x,D12.10)'
fr = '(A7,1x,A16,2x,I4,2x,I4,2x,F8.6,2x,I3,1x,I3,1x,A6)'

;the rejected values are picked out
ikeep = where(sig1o le cutoffu and sig1o ge cutoffl)
i1 = where(sig1o gt cutoffu)
i2 = where(sig1o lt cutoffl)
if (i1[0] ne -1 and i2[0] ne -1) then itoss = [i1,i2]
if (i1[0] ne -1 and i2[0] eq -1) then itoss = i1
if (i1[0] eq -1 and i2[0] ne -1) then itoss = i2
if (i1[0] eq -1 and i2[0] eq -1) then begin
    print,''
    print,'NO VALUES WERE THROWN OUT!'
    print,'No improvement made. Try again or call it good...'
    stop
endif

;by sorting the array the i2 rejects (low values) are put in order
itoss = itoss[bsort(itoss)]

nk = n_elements(ikeep)
nt = n_elements(itoss)

print,''
print,'The number of filtered bins: '+strn(n1)
print,'The number of bins that were kept: '+strn(nk)
pause

; if the pfitall and runall keywords are NOT used then it writes NEW
; pfitlov.all and run.final files
if (n_elements(pfitall) eq 0) then begin
    openw,5,'pfitlov.all'
    for j=0,nk-1 do begin
        i = ikeep[j]
        printf,5,nameo[i],velo[i],sig1o[i],sig2o[i],h3o[i],h4o[i],$
          r1o[i],r2o[i],r3o[i],rmso[i],format=fp
    endfor
    free_lun,5
; a backup is created
    openw,5,'pfitlov.all.first'
    for j=0,nk-1 do begin
        i = ikeep[j]
        printf,5,nameo[i],velo[i],sig1o[i],sig2o[i],h3o[i],h4o[i],$
          r1o[i],r2o[i],r3o[i],rmso[i],format=fp
    endfor
    free_lun,5
endif

if (n_elements(runall) eq 0) then begin
    openw,5,'run.final'
    for j=0,nk-1 do begin
        i = ikeep[j]
        printf,5,ro1[i],nameo1[i],wo1[i],wo2[i],zo[i],so[i],smo[i],subo[i],format=fr
    endfor
    free_lun,5
; a backup is created
    openw,5,'run.final.first'
    for j=0,nk-1 do begin
        i = ikeep[j]
        printf,5,ro1[i],nameo1[i],wo1[i],wo2[i],zo[i],so[i],smo[i],subo[i],format=fr
    endfor
    free_lun,5
endif

;if the pfitall keyword is used then the old pfitlov.all and run.final
;files are read in, then appended to.
if (n_elements(pfitall) ne 0) then begin

    readcol,pfitall,silent=1,format='a,f,f,f,f,f,f,f,f,d',nameA,velA,sig1A,sig2A,h3A,h4A,r1A,r2A,r3A,rmsA

    nameAout = [nameA,nameo[ikeep]]
    velAout = [velA,velo[ikeep]]
    sig1Aout = [sig1A,sig1o[ikeep]]
    sig2Aout = [sig2A,sig2o[ikeep]]
    h3Aout = [h3A,h3o[ikeep]]
    h4Aout = [h4A,h4o[ikeep]]
    r1Aout = [r1A,r1o[ikeep]]
    r2Aout = [r2A,r2o[ikeep]]
    r3Aout = [r3A,r3o[ikeep]]
    rmsAout = [rmsA,rmso[ikeep]]
    n2 = n_elements(nameAout)

    openw,5,pfitall
    for j=0,n2-1 do begin
        printf,5,nameAout[j],velAout[j],sig1Aout[j],sig2Aout[j],h3Aout[j],h4Aout[j],$
          r1Aout[j],r2Aout[j],r3Aout[j],rmsAout[j],format=fp
    endfor
    free_lun,5
endif

;same thing, only for the run.final file
if (n_elements(runall) ne 0) then begin

    readcol,runall,silent=1,format='a,a,i,i,f,i,i,a',rA1,nameA1,wA1,wA2,zA,sA,smA,subA

    rA1out = [rA1,ro1[ikeep]]
    nameA1out = [nameA1,nameo1[ikeep]]
    wA1out = [wA1,wo1[ikeep]]
    wA2out = [wA2,wo2[ikeep]]
    zAout = [zA,zo[ikeep]]
    sAout = [sA,so[ikeep]]
    smAout = [smA,smo[ikeep]]
    subAout = [subA,subo[ikeep]]
    n3 = n_elements(rA1out)

    openw,5,runall
    for j=0,n3-1 do begin
        printf,5,ra1out[j],namea1out[j],wa1out[j],wa2out[j],zaout[j],saout[j],$
          smaout[j],subaout[j],format=fr
    endfor
    free_lun,5

endif

;the good frames have now been written out. now a new run script is
;created for the bad ones.

if (n_elements(rout) eq 0) then rout = 'run1'
openw,5,rout
for j=0,nt-1 do begin
    i = itoss[j]
    printf,5,ro1[i],nameo1[i],wo1[i],wo2[i],zo[i],so[i],smo[i],subo[i],format=fr
endfor
free_lun,5

stop
END

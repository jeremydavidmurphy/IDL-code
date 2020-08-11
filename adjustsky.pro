PRO adjustsky, list
; this code takes the Vaccine output, uses the pefy.fits frames from
; the 2 sky nods, and generates a range of levels of sky
; subtraction. once this is complete, the output weights can be fed
; back into Vaccine to run a proper sky model and subtraction.

; Example of the list:
; jm1897_pe.fits jm1896_pefy.fits 2.0 jm1898_pefy.fits 2.0 flat_oct08_n1en.fits

; SUGGESTION! USE YOUR SCIENCE_N?.LIST.2 LISTS TO GENERATE THIS
; LIST. It helps with pairing the sky frames to the correct science
; frames.

readcol,list,f='a,a,f,a,f,a',datal,sky1l,w1l,sky2l,w2l,flatl
n1 = n_elements(datal)
wgts1 = [1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0]
wgts2 = [0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1,0.0]

for j=0,n1-1 do begin
   print,'Working on frame '+datal[j]
   print,''
   data = readfits(datal[j],hd,/silent)
;   data = double(data)
   sky1 = readfits(sky1l[j],/silent)
;   sky1 = double(sky1)
   wt = w1l[j]
   sky1 = sky1 * wt
   sky2 = readfits(sky2l[j],/silent)
;   sky2 = double(sky2)
   wt = w2l[j]
   sky2 = sky2 * wt
   flat = readfits(flatl[j],/silent)
;   flat = double(flat)
   i = where(flat eq 0)
   flat[i] = 1.0
   dataf = data/flat
   dataf[i] = 0.0
   nameout = strsplit(datal[j],'.',/extract)
   nameout = nameout[0]
   for k=0,9 do begin
      print,'Adjusting sky...'
      w1 = wgts1[k]
      w2 = wgts2[k]
      skyA = sky1*w1 + sky2*w2
      skyB = sky1*w2 + sky2*w1
      datafA = dataf - skyA
      datafB = dataf - skyB
      writefits,nameout+'_'+strn(wgts1[k])+'-'+strn(wgts2[k])+'.fits',datafA,h
      writefits,nameout+'_'+strn(wgts2[k])+'-'+strn(wgts1[k])+'.fits',datafB,h
   endfor
endfor

STOP
END

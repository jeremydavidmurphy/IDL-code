; This routine returns the average of a value, after some amount of
; sigma_clipping. It's a tweaked version of clipmean.pro

function meanWclip,input,nsigma=nsigma,niter=niter,decimals=decimals,error=error,converge=converge
  compile_opt idl2

  if keyword_set(nsigma) eq 0 then nsigma=3 else nsigma=nsigma
  if keyword_set(niter) eq 0 then niter=20 else niter=niter
  if keyword_set(error) then begin
     earr=error 
 endif else begin
     earr = findgen(n_elements(input))+1
 endelse

  arr=input
  var=1/earr
  if keyword_set(decimals) then begin
     if decimals ge 2 and decimals le 10 then begin
        len=decimals
    endif else begin
        print,' '
        print,'  ERROR: IF SET, DECIMALS KEYWORD MUST'
        print,'  BE GREATER THAN OR EQUAL TO 2 AND LESS'
        print,'  THAN OR EQUAL TO 10.'
        print,' '
        result=-99
        return,result
    endelse
endif else begin
     len=4
 endelse
  n=0
  err=0
retry:
  for i=0,niter-1 do begin
      n=n+1
      m0=strtrim(string(mean(arr)),2)
      sigma=stddev(arr)
      avg=(total(var*arr)/n_elements(arr))/(total(var)/n_elements(var))
      a=where(arr gt (avg-(nsigma*sigma)) and arr lt (avg+(nsigma*sigma)),count)
      if a[0] eq -1 or count le 1 then begin
          print,' '
          print,'  ERROR: WEIRD DISTRIBUTION!'
          print,' '
          result=-99
          return,result
      endif
      arr=arr[a]
      earr=earr[a]
      var=1/earr
      m1=strtrim(string(mean(arr)),2)
      dec0=strpos(m0,'.')
      dec1=strpos(m1,'.')
      r0=fix(strmid(m0,dec0+len+1,1))
      r1=fix(strmid(m1,dec1+len+1,1))
      if r0 ge 0 and r0 le 4 then rm0=strmid(m0,0,dec0+len+1) else rm0=strmid(strtrim(string(double(strmid(m0,0,dec0+len+1))+10.^(-len)),2),0,dec0+len+1)
      if r1 ge 0 and r1 le 4 then rm1=strmid(m1,0,dec1+len+1) else rm1=strmid(strtrim(string(double(strmid(m1,0,dec1+len+1))+10.^(-len)),2),0,dec1+len+1)
      if keyword_set(converge) then begin
          if rm0 eq rm1 then begin
              err=0
              break
          endif
          if rm0 ne rm1 and i eq (niter-1) and len ge 2 then begin
              len=len-1
              goto,retry
              err=0
          endif else begin
              err=1
          endelse
      endif else begin
          err=0
      endelse
  endfor
  if err eq 0 then begin
      result=rm1
  endif else begin
      print,'  '
     print,'  ERROR: SIGMA CLIPPED MEAN FAILED TO'
     print,'  CONVERGE TO AN ACCEPTABLY ACCURATE'
     print,'  VALUE AFTER '+strtrim(string(n),2)+' ITERATIONS.'
     print,' '
     result=-99
 endelse
  return,result

end

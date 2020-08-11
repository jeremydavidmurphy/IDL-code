;#############################################################################
;
; Copyright (C) 2008-2009, Michele Cappellari
; E-mail: cappellari_at_astro.ox.ac.uk
;
; Updated versions of the software are available from my web page
; http://purl.org/cappellari/idl
;
; If you have found this software useful for your research,
; I would appreciate an acknowledgment and a link to the website.
;
; This software is provided as is without any warranty whatsoever.
; Permission to use, for non-commercial purposes is granted.
; Permission to modify for personal or internal use is granted,
; provided this copyright and disclaimer are included unchanged
; at the beginning of the file. All other rights are reserved.
;
;#############################################################################
function ibetam, a1, b, x
compile_opt idl2
on_error, 2
;
; Incomplete beta function defined in the same 
; way as the Mathematica function Beta[x,a,b].
; This routine requires x < 1 && b > 0.
; It was tested against Mathematica.
;
; V1.0: Michele Cappellari, Oxford, 01/APR/2008
; V2.0: Use Hypergeometric function for a < 0 || b < 0.
;    From equation (6.6.8) of Abramoviz & Stegun (1965)
;    http://www.nrbook.com/abramowitz_and_stegun/page_263.htm 
;    MC, Oxford, 04/APR/2008
; V3.0: Use recurrence relation of equation (26.5.16) 
;    from Abramoviz & Stegun (1965) for a < 0 && b > 0.
;    Removed case b < 0 which is not required for JAM.
;    http://www.nrbook.com/abramowitz_and_stegun/page_944.htm
;    After suggestion by Gary Mamon. MC, Oxford, 16/APR/2009

a = a1 + 3d-7 ; Perturb to avoid singularities in gamma and ibeta
if b ge 0 then $
    if a ge 0 then $
        ib = ibeta(a,b,x)*beta(a,b) $
    else begin        
        p = ceil(abs(a))
        sum = 0d ; Summation is vectorized over x
        for j=0,p-1 do sum += gamma(a+b+j)/gamma(a+j+1d)*x^(a+j)
        ib = (sum*(1d - x)^b/gamma(b) + ibeta(a+p,b,x)) * beta(a,b)
    endelse else message, 'b < 0 not implemented'     
     
return, ib 
end
;----------------------------------------------------------------------
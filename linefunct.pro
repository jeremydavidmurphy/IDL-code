;for use in the LMFIT.PRO routine
FUNCTION polyfunct,x,a
func = a[0] + a[1]*x + a[2]*x^2

;The return is the value of the function, and the partial derivatives
;of each element, a, in the function.
RETURN, [[func],[1.0],[a[1]],[2*a[2]*x]]
END

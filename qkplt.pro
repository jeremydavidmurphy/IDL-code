;This routine is part of makebins.pro. It simply plots the input
;spectra in window,2

function qkplt, specguess
on_error,2

;this function just plots a 1-D array against its index + 1.

window,2,retain=2,xsize=650,ysize=450
;set_plot,'x'

n1 = n_elements(specguess)
x = indgen(n1)+1

plot,x,specguess,thick=1,title='The spectra (summed) of your chosen bin',$
  xrange=[100,n1-250],xstyle=1

return,specguess
end


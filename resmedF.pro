FUNCTION resmedF, file

; This function reads in a bin####.res file (the output of
; resolution.pro) and returns a file bin####.resT, which is a text
; file of three numbers, giving the median of the bin####.res file 
; The input file must be in the form
; 4046.55  4077.83  4358.33  4678.15  4799.91  4916.07  5085.82  5154.66  5460.74
;  4.8362   4.8057   4.5354   4.4667   4.4744   4.2674   4.5164   4.8774   4.6218
;  4.8186   4.7709   4.5216   4.4607   4.4585   4.2687   4.4807   5.1480   4.5776
;  4.9689   4.8800   4.4006   4.1535   4.1059   4.0006   4.0438   4.0645   3.9775

readcol,file,b1,b2,b3,g1,g2,g3,r1,r2,r3

wave = fltarr(3)
wave[0] = mean(b1[0],b2[0],b3[0])
wave[1] = mean(g1[0],g2[0],g3[0])
wave[2] = mean(r1[0],r2[0],r3[0])

readcol,file,b1,b2,b3,g1,g2,g3,r1,r2,r3,skipline=1
bmean = mean([[b1],[b2],[b3]])
bmedian = median([[b1],[b2],[b3]],even)
bvar = variance([[b1],[b2],[b3]])

END

; This function reads in the bin####.coord files and returns the mean
; of the absolute value of the RA and Dec values. 

FUNCTION avecoordF, file

readcol,file,format='i,f,f',skipline=1,fiber,ra,dec,silent=1

ira = where(ra ne -1.0)
idec = where(dec ne -1.0)
mra = mean(abs(ra[ira]))
mdec = mean(abs(dec[idec]))
rad = sqrt(mra^2 + mdec^2)
;window,0,retain=2
;loadct,4
;if (n_elements(ra) gt 1) then plot,abs(ra),abs(dec),psym=1 $
;  else  plots,abs(ra),abs(dec),psym=1;

;plots,mra,mdec,psym=4,color=150
;pause
return,[mra,mdec,rad]

END

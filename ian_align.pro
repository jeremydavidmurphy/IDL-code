pro ian_align

readcol,'wave.list',format='a',wavelist
readcol,'spectra.list',format='a',speclist

n0 = n_elements(wavelist)
if (n_elements(speclist) ne n0) then stop

wavearr = dblarr(2048,n0)
dataarr = dblarr(2048,n0)

wr1 = ''
wr2 = ''
print,'Enter a starting wavelength:'
read,wr1
wr1 = float(wr1)
print,'Enter an ending wavelength:'
read,wr2
wr2 = float(wr2)

cntr = 0
cntr1 = 0
swch = 'notdone'
loadct,0
window,0,retain=2
device,decomposed=0
set_plot,'x'

repeat begin
    w1 = readfits(wavelist[cntr1])
    w2 = readfits(wavelist[cntr1+1])
    wch2 = max(w2)
    if (wch2 lt wr1) then goto,jump1
    d1 = readfits(speclist[cntr1])
    d2 = readfits(speclist[cntr1+1])
    ti = where(w1 eq 0.0,count)
    if (count ne 0) then w1[ti] = !values.f_nan
    ti = where(w2 eq 0.0,count)
    if (count ne 0) then w2[ti] = !values.f_nan
    ti = where(d1 eq 0.0,count)
    if (count ne 0) then d1[ti] = !values.f_nan
    ti = where(d2 eq 0.0,count)
    if (count ne 0) then d2[ti] = !values.f_nan
    w1m = max(w1)
    ind = where(w2 ge w1m)
    print,w1[0]
    

    if (wch2 ge wr2) then begin
        swch = 'done'
        if (cntr eq 0) then wout = [w1,w2[ind]] else wout = [wout,w2[ind]]
        if (cntr eq 0) then dout = [d1,d2[ind]] else dout = [dout,d2[ind]]
        goto,jump1
    endif
    if (cntr eq 0) then wout = [w1,w2[ind]] else wout = [wout,w2[ind]]
    if (cntr eq 0) then dout = [d1,d2[ind]] else dout = [dout,d2[ind]]
    if (cntr eq 0) then begin
        loadct,0
        plot,w1,d1,title='The final spectra',xrange=[wr1-15,wr2+15],$
          xstyle=1,yrange=[0.0,1.2],ystyle=1
        loadct,4
        oplot,w2[ind],d2[ind],color=60
    endif else begin
        oplot,w2[ind],d2[ind],color=60+cntr*20
    endelse
    cntr = cntr + 1
    jump1:
    cntr1 = cntr1 + 1
endrep until (swch eq 'done')

;the nan's and inf's are changed back to -666 flags.    
goodindex = 0
for j=0,n_elements(wout)-1 do begin
    p = finite(wout[j])
    if (p eq 1) then goodindex = [goodindex,j]
endfor
goodindex = goodindex[1:*]

wout = wout[goodindex]
dout = dout[goodindex]
dif = fltarr(n_elements(wout))
for j=0,n_elements(dif)-2 do dif[j] = wout[j+1]-wout[j]
print,'The median dispersion is: '+strn(median(dif))
print,'The average dispersion is: '+strn(mean(dif))

disp=fltarr(1)
print,''
print,'Enter a dispersion for the interpolation:'
read,disp
print,disp

wo = fltarr(n_elements(wout))
for j=0,n_elements(wout)-1 do wo[j] = wr1+j*disp

data = interpol(dout,wout,wo)
print,'Next ENTER overplots the interpolated data...'
pause

loadct,0
oplot,wo,data

wfirst = where(wo gt wr1)
wfirst = wfirst[0]
wsecond = where(wo le wr2)
wsecond = wsecond[n_elements(wsecond)-1]

data = data[wfirst:wsecond]

name = strn(wo[wfirst])+'_'+strn(wo[wsecond])+'_'+strn(disp)+'.fits'
print,name

writefits,name,data
print,'Next ENTER deletes the plot...'
pause
wdelete
stop
end

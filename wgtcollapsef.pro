function wgtcollapsef, data, weight, maskname, aperture

; THIS VERSION REJECTS ALL -666 FLAGS AND SETS THAT COLUMN OF DATA
; EQUAL TO -666. VERSION 2 OF THIS CODE TAKES THE MEDIAN, AND
; THEREFORE CAN HANDLE UP TO 2 -666 FLAGS BEFORE REJECTING THE COLUMN
; OF DATA ENTIRELY.

; this function collapses VIRUS-P data. It requires the mask file. If
; more than 1 pixel in a given column for a given fiber is masked it
; returns a -666 flag. If one or no pixels are masked, it returns the
; following value

; DATAOUT = total (data * weight) / total(weight)

; DATA: The uncollapsed data.
; WEIGHT: The corresponding weight file
; INDEX: The mask file values - 1 (so that it's an IDL index)
; APERTURE: The width of the extraction aperture

; MODIFIED ON DEC 10, 2010: This code used to keep columns of data,
; even when a -666 flag existed. This isn't the right way to run this,
; as the output is the TOTAL of the aperture. So, setting something to
; zero doesn't work. It now just returns a -666 flag for all. (Another
; way to approach this in the future would be to take the median.

step = (aperture-1)/2.0

n1 = n_elements(data[*,0])

readcol,maskname,format='i,x',y
n3 = n_elements(y)
indexarr = y - 1.0

dataout = fltarr(n1,n3)

for jj=0,n3-1 do begin ;a loop through each fiber
    index = uint(indexarr[jj])
    rowD = data[*,index-step:index+step]
    rowW = weight[*,index-step:index+step]
    for kk=0,n1-1 do begin ;a loop through wavelength
        top = rowD[kk,*]
        bot = rowW[kk,*]
        badi = where(top eq -666,count)
        if (count ge 1) then begin
;            print,'-666!'
            dataout[kk,jj] = -666.0
        endif
        if (count eq 0) then begin
            toptop = total(top * bot)
            botbot = total(bot)
            dataout[kk,jj] = toptop / botbot
        endif
    endfor
endfor

return,dataout

stop
end

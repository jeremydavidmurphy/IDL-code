PRO Pframe, dataframe, ptowname, maskname

; A modification of Pall to plot an entire, extracted, frame. The
; fibers are plotted by numbers

n_profile = 5.
step = uint((n_profile-1)/2.)

data = readfits(dataframe)
n0 = n_elements(data[*,0]) ;wavelength
n1 = n_elements(data[0,*]) ;fiber number

data[where(data eq -666)] = !Values.F_NAN

readcol,silent=1,ptowname,format='d,d,d,d,d',c0,c1,c2,c3,c4
ptow = [[c0],[c1],[c2],[c3],[c4]];this is a 247x5 array, NOT a 5x247 array!
wave = fltarr(n0,n1)

for k=0,n_elements(c0)-1 do begin
    for j=0,n0-1 do begin
        wave[j,k] = c0[k]+c1[k]*j+c2[k]*j*j+c3[k]*j*j*j+c4[k]*j*j*j*j
    endfor
endfor

readcol,maskname,silent=1,format='i,i',mask,fibs
imask = mask-1

ans = ''

fiber = intarr(1)
window,2,retain=2
device,decomposed=0
repeat begin
    
    print,'Enter a fiber number to plot...'
    read,fiber
    fiber = fiber - 1
    onefib = data[*,imask[fiber]-step:imask[fiber]+step]
    help,onefib
    print,onefib[1600,*]
    onefib = total(onefib,2)
    
    plot,wave[*,fiber],onefib,xrange=[3600,5500],xstyle=1
   
   print,'Press ENTER to continue and "q" to quit:'
   read, ans
   
endrep until (ans eq 'q')


print,'The next ENTER deletes the plots.'
wdelete

stop
END

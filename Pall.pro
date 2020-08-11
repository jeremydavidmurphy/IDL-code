PRO Pall, file, type

;This routine is used to plot individual fits files from a list (if
;the WI and DISP keywords are used) or from a collapsed frame. It is a visual
;inspection routine to quickly look at several different spectra. The
;previous spectra is kept on screen in a different window. The
;ptow.dat file is required if you're wanting to plot a collapsed frame
;rather than a list of individual spectra.

;TYPE: enter 'single','Rframe' or 'Aframe' for a list of individual
;1-D files, a non-aligned frame, or an aligned frame.

;Modified (10-08-2009): The routine now allows for individual fibers
;in a collapsed frame to be plotted. The user will be prompted for
;this.

;Modified (02-2010): Dumped the clever sorting routine. Added 'type'
;as a keyword.

xdown = 3600
xup = 5500

if (type eq 'Rframe') then begin
    swch = 'off'
    ptowf = ''
    print,'Enter the name of the ptow.dat frame:'
    read,ptowf
    readcol,silent=1,ptowf,f='d,d,d,d,d',a,b,c,d,e
    n1 = n_elements(a)
    data = readfits(file,/silent)
    if (n_elements(data[0,*]) ne n1) then begin
        print,'The length of your ptow file and data frame are not in agreement!'
        stop
    endif
    n2 = n_elements(data[*,0])
    wave = fltarr(n2)
    cntr = 0
    ans = ''
    window,0,retain=2
    window,2,retain=2

    repeat begin
        for j=0,n2-1 do wave[j]=a[cntr]+b[cntr]*j+c[cntr]*$
          j^2+d[cntr]*j^3+e[cntr]*j^4
        wset,2
        plot,wave,data[*,cntr],xrange=[xdown,xup],xstyle=1,/ynozero,$
          title='Fiber '+strn(cntr+1),xtitle='Wavelength (A)',charsize=1.2
        if (cntr ne 0) then begin
            wset,0
            plot,wave,data[*,cntr-1],xrange=[xdown,xup],xstyle=1,/ynozero,$
              title='Fiber '+strn(cntr),xtitle='Wavelength (A)',charsize=1.2
        endif
        
        print,'Press ENTER to continue, "f" for a fiber, "p" to print, or "q" to quit:'
        read, ans
        if (ans eq 'f') then begin
            ans1 = ''
            num = ''
            jump2:
            print,'Enter a fiber number:'
            read,num
            num = uint(num)
            spec1 = data[*,num-1]
            if (swch eq 'on') then begin
                wset,0
                plot,wave,spec2,xrange=[xdown,xup],xstyle=1,/ynozero,$
                  title='Fiber '+strn(num2),xtitle='Wavelength (A)',charsize=1.2
            endif
            wset,2
            plot,wave,spec1,xrange=[xdown,xup],xstyle=1,/ynozero,$
              title='Fiber '+strn(num),xtitle='Wavelength (A)',charsize=1.2
            print,'Another fiber? (y/n)'
            read,ans1
            if (ans1 eq 'y') then begin
                spec2 = spec1
                num2 = num
                swch = 'on'
                goto, jump2
            endif else goto, jump3
        endif
        jump3:
        cntr = cntr + 1
    endrep until (ans eq 'q' or cntr eq n1-1)
endif

if (type eq 'Aframe') then begin
    data = readfits(file,/silent)
    n2 = n_elements(data[*,0])
    n1 = n_elements(data[0,*])
    wave = fltarr(n2)
    cntr = 0
    ans = ''
    window,0,retain=2
    window,2,retain=2
    wi = 0.0
    disp = 0.0
    print,'Enter an initial wavelength:'
    read,wi
    print,'Enter a dispersion:'
    read,disp
    for j=0,n2-1 do wave[j] = wi + j * disp
    repeat begin
        wset,2
        plot,wave,data[*,cntr],xrange=[xdown,xup],xstyle=1,/ynozero,$
          title='Fiber '+strn(cntr+1),xtitle='Wavelength (A)',charsize=1.2
        if (cntr ne 0) then begin
            wset,0
            plot,wave,data[*,cntr-1],xrange=[xdown,xup],xstyle=1,/ynozero,$
          title='Fiber '+strn(cntr),xtitle='Wavelength (A)',charsize=1.2
        endif
        print,'Press ENTER to continue, "f" for a fiber, or "q" to quit:'
        read, ans
        if (ans eq '') then goto, jump6
        if (ans eq 'f') then begin
            ans1 = ''
            num = ''
            jump5:
            print,'Enter a fiber number:'
            read,num
            num = uint(num)
            spec1 = data[*,num-1]
            if (swch eq 'on') then begin
                wset,0
                plot,wave,spec2,xrange=[xdown,xup],xstyle=1,/ynozero,$
                  title='Fiber '+strn(num2),xtitle='Wavelength (A)',charsize=1.2
            endif
            wset,2
            plot,wave,spec1,xrange=[xdown,xup],xstyle=1,/ynozero,$
              title='Fiber '+strn(num),xtitle='Wavelength (A)',charsize=1.2
            print,'Another fiber? (y/n)'
            read,ans1
            if (ans1 eq 'y') then begin
                spec2 = spec1
                num2 = num
                swch = 'on'
                goto, jump5
            endif else goto, jump6
        endif
        jump6:
        cntr = cntr + 1
    endrep until (ans eq 'q' or cntr eq n1-1)
endif
    
if (type eq 'single') then begin
    readcol,file,silent=1,format='a',files
    n1 = n_elements(files)
    test = readfits(files[0],/silent)
    n2 = n_elements(test)
    
    wave = fltarr(n2)
    wi = 0.0
    disp = 0.0
    print,'Enter an initial wavelength:'
    read,wi
    print,'Enter a dispersion:'
    read,disp
    for j=0,n2-1 do wave[j] = wi + j * disp

    ans = ''
    cntr = 0
    window,0,retain=2
    window,2,retain=2

    repeat begin
        
        spectra1 = readfits(files[cntr],/silent)
        i = where(spectra1 eq -666,count)
        if (count ne 0) then spectra1[i] = !values.F_NAN
        name1 = files[cntr]
        
        wset,2
        plot,wave,spectra1,xrange=[xdown,xup],xstyle=1,/ynozero,$
          title=name1,xtitle='Wavelength (A)',charsize=1.2
        
        if (cntr ne 0) then begin
            wset,0
            plot,wave,spectra2,xrange=[xdown,xup],xstyle=1,/ynozero,$
              title=name2,xtitle='Wavelength (A)',charsize=1.2
        endif
        
        spectra2 = spectra1
        name2 = name1
        if (cntr eq n1-1) then goto, jump1
        
        cntr = cntr + 1
        
        print,'Press ENTER to continue, "p" to print, or "q" to quit:'
        read, ans

        if (ans eq 'p') then begin
            set_plot,'ps'
            name = name1+'.ps'
            device,file=name,/color
            plot,wave,spectra1,xrange=[xdown,xup],xstyle=1,/ynozero,$
              title=name1,xtitle='Wavelength (A)',charsize=1.2,$
              xthick=2,ythick=2,charthick=2,thick=2
            device,/close_file
            set_plot,'x'
        endif
        
    endrep until (ans eq 'q')
endif
jump1:
    
print,'The next ENTER deletes the plots.'
pause
wdelete,0,2

END

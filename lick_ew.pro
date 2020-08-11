Pro Lick_ew, speclist, velocity_dispersion=vel_disp, resolution=res, $
             outfile=outfile, PLOT=SHOWPLOT, AIR=AIR, LINEAR=LINEAR, $
             SILENT=SILENT, error_list=error_list, err_exten=err_exten, $
             err_nslice=err_nslice, spec_exten=spec_exten, $
             spec_nslice=spec_nslice, SN=SN_FLAG, PIX=PIX
;+
; NAME:
; 		LICKEW
; VERSION:
;		4.0
; PURPOSE:
;		Measures Lick equivalent widths for a list of input spectra
; CATEGORY:
;		
; CALLING SEQUENCE:
;		LICKEW, SPECLIST [, INDEXLIST=, OUTFILE=, RES=, ERROR_LIST=
;		/PLOT, /AIR]
; INPUTS:
;	SPECLIST
;		A file containing a list of the spectra for which to
;		measure Lick indices
;       INDEXLIST
;               A file containing a list of Lick EW to measure.  The
;               file should have the starting and ending wavelengths
;               of the index bandpass (band1 & band2), as well as the
;               starting and ending wavelengths of the blue (blue1 &
;               blue2) and red (red1 & red2) pseudo-continuum
;               regions.  It should also contain the FWHM resolution
;               in Angstroms of the Lick/IDS system at the wavelength
;               of that index, and whether or not the index is
;               measured in magnitudes instead of equivalent widths (1
;               for magnitudes, 0 for EW. It should be in the
;               following format:
;
;               band1  band2  blue1  blue2  red1  red2 IndRes Mag? IndName
;                                                    
;               For example:
;               
;         4143.375 4178.375 4081.375 4118.875 4245.375 4285.375 10.3 1 CN1
;
;               If no indexlist is specified, the program will use the
;               default file '$IDL_PATH/lickindexlist.txt' which should
;               be placed in the $IDL_PATH directory.
;
;        OUTFILE
;               The output file to which the table of Lick Index
;               measurements will be written.  If no outfile is
;               specified, the default output file is index_out.txt.
;        RES    
;               The resolution of the spectrum.  RES may be an
;               array of numbers the same length as the input
;               filelist, one resolution for each 
;               input spectrum.  If RES is not set, Lickew.pro
;               assumes the spectrum has IDS variable resolution and
;               does not change the resolution to mesaure indices. 
;        ERROR_LIST
;
;        /PLOT  
;               If set, the program shows plots of each Lick Index as
;               it is measured.
;        /AIR
;               Set for input spectra in air wavelengths. Otherwise, 
;               Lickew.pro assumes the input is in vacuum wavelengths
;               and will shift the lick indices into vacuum
;               wavelengths before measuring.
;
; OUTPUTS:
;		The output file of measured Lick Indices
; COMMON BLOCKS:
;		None.
; SIDE EFFECTS:
;		None.
; RESTRICTIONS:
;		None.
; PROCEDURE:
;
; MODIFICATION HISTORY:
;		Created 16-APR-04 by Genevieve Graves
;               Fixed smoothing problem 27-APR-04 so now appears to
;               correctly smooth SDSS.
;               Smoothing was still not a good match to IRAF 'gauss'
;               12-JUN-04 Smoothing fixed!! Now matches IRAF exactly.
;               ??-JUL-05 Added /VAC as an option, so that input can
;               be in vacuum wavelengths -- lickew.pro will shift the
;               indices to vacuum wavelengths.
;               25-JUL-05 Changing input resolution to be in pixels
;               rather than in Angstroms -- this is good for SDSS
;               input where everything is binned logarithmically (so
;               that lambda/dlambda is constant over the spectrum).
;               2-AUG-05 Changed default to be input in vacuum
;               wavelengths ---> replaced keyword /VAC with keyword
;               /AIR.  If not set, program assumes vacuum
;               wavelengths.  Also changed program so that if
;               resolution input is a single number (not an array),
;               Lickew.pro will convert it into an array.
;               30-NOV-05 Added options /SILENT and /LINEAR.  /SILENT
;               allows the program to be run with minimal printing to
;               the terminal.  /LINEAR allows input of linearly-binned
;               spectra.  Code was changed to accommodate these new
;               flags and data input.
;               07-Aug-06 Fixed problem with sig_tot not getting set
;               if everything is already on the Lick/IDS system.
;-



if N_PARAMS() eq 0 then begin
    print, 'Syntax - lick_ew, speclist [, velocity_dispersion=,'
    print, '         resolution=, outfile=, spec_exten=, spec_nslice=, '
    print, '         error_list=, err_exten=, err_nslice=,'
    print, '         /PLOT, /AIR, /LINEAR, /SN, /PIX, /SILENT]'
    return
endif

d = FSC_COLOR(/AllColors,ColorStructure=ctable)    ; load color table
c = 3e5                         ; speed of light in km/s

;----------------------------------------------------------------
; Set up clean exit on error
;----------------------------------------------------------------

CATCH, error_status

if error_status ne 0 then begin
    print, !ERR_STRING
    exit_cleanly, unit_abuntab=ew_tab, unit_errtab=ewerr_tab
    return
endif 
    
;----------------------------------------------------------------
; *************** CHECK INPUT ***********************
;----------------------------------------------------------------
; Determine if the spectral resolution is on the IDS system 
; already, or constant (in angstrom space) across spectrum
;----------------------------------------------------------------

if N_ELEMENTS(res) eq 0 then begin
    flag_idsres = 1
endif else begin
    flag_idsres = 0
endelse

if N_ELEMENTS(outfile) eq 0 then begin
    outfile=speclist+'_index'
endif

;----------------------------------------------------------------
; Default spectrum and error EXTEN_NO and NSLICE to zero.
;---------------------------------------------------------------

if (N_ELEMENTS(error_list)+N_ELEMENTS(err_exten)+N_ELEMENTS(err_nslice)) gt 0 then begin
    err_flag = 1
endif else begin
    err_flag =0
endelse 


if N_ELEMENTS(spec_nslice) eq 0 then begin
    spec_nslice = 0
endif
if N_ELEMENTS(err_nslice) eq 0 then begin
    err_nslice = 0
endif else begin
    err_nslice = err_nslice - 1 
endelse 
if N_ELEMENTS(spec_exten) eq 0 then begin
    spec_exten = 0
endif
if N_ELEMENTS(err_exten) eq 0 then begin
    err_exten = 0
endif

;----------------------------------------------------------------
; Load velocity dispersion corrections
;---------------------------------------------------------------

LOAD_SIGCORRS, sigcorr_slope_arr, sigcorr_icpt_arr

;------------------------------------------------------
; ************** INPUT FILES *******************
;------------------------------------------------------
; Read in list of spectra, lick indices
;------------------------------------------------------

indexlist='$EZ_AGES_DIR/lickindexlist.txt'

if NOT(KEYWORD_SET(silent)) then begin
    print, 'READING IN LIST OF SPECTRA.'
endif

readcol, speclist, list, FORMAT='A', /silent
nspecs = (SIZE(list))[1]

if N_ELEMENTS(error_list) gt 0 then begin
    if NOT(KEYWORD_SET(silent)) then begin
        print, 'READING IN LIST OF ERROR SPECTRA.'
    endif
    readcol, error_list, errlist, FORMAT='A', /silent
    if (SIZE(errlist))[1] ne nspecs then begin
        print, 'ERROR: # of error files does not equal # of spectra.'
        return
    endif
endif 

if N_ELEMENTS(vel_disp) eq 0 then begin
    vel_disp = FLTARR(nspecs)
    print, 'Warning! Using velocity_dispersion = 0 for all input spectra.'
    print, '         Probably okay for globular clusters'
    print, '         NOT OKAY FOR GALAXIES!'
endif else if SIZE(vel_disp,/TYPE) eq 7 then begin
    READCOL,vel_disp,vel_array
    vel_disp = vel_array
endif else begin
    print, 'Using velocity dispersion = ', vel_disp, 'for all input spectra.'
    vel_array = fltarr(nspecs)
    vel_array[*] = vel_disp
    vel_disp = vel_array
endelse 

if N_ELEMENTS(res) ne 0 then begin
    if SIZE(res,/TYPE) eq 7 then begin
        if NOT(KEYWORD_SET(silent)) then begin
            print, 'Reading resolution values from input file.'
            READCOL, res, res_array,/SILENT
        endif 
    endif else begin
        if NOT(KEYWORD_SET(silent)) then begin
            print, 'Using same resolution for all input spectra:'
            print, '      Resolution = ', res
        endif 
        res_array = fltarr(nspecs)
        res_array[*] = res
    endelse 
    res = res_array
endif

if NOT(KEYWORD_SET(silent)) then begin
    print, 'READING IN LIST OF LICK INDICES.'
endif 

readcol, indexlist, band1, band2, blue1, blue2, red1, red2, indres, $
                         indmag, indexname, c1, c2, $
                         FORMAT='F,F,F,F,F,F,F,I,A,F,F', /silent
nindex = (SIZE(indexname))[1]

;--------------------------------------------------------
; If /AIR is not set, convert indices to vacuum wavelengths
;--------------------------------------------------------

if NOT(KEYWORD_SET(air)) then begin

    AIRTOVAC, band1
    AIRTOVAC, band2
    AIRTOVAC, blue1
    AIRTOVAC, blue2
    AIRTOVAC, red1
    AIRTOVAC, red2

endif

;---------------------------------------------------------
; ************** OUTPUT FILES *********************
;---------------------------------------------------------
; Open output file(s), write column titles
;---------------------------------------------------------

ew_tab = 50
openw, ew_tab, outfile
printf, ew_tab, FORMAT='($,25(A,:,"  "))',indexname
printf, ew_tab, ''

if NOT(KEYWORD_SET(silent)) then begin
    print, 'FILE ',outfile,' OPEN FOR OUTPUT.'
    print, 'Close manually using: "IDL> close, 50" if interrupt occurs.'
endif 

if err_flag then begin
    ewerr_tab = 51
    openw, ewerr_tab, outfile+'_err'
    printf, ewerr_tab, FORMAT='($,25(A,:,"  "))',indexname
    printf, ewerr_tab, ''
endif

;--------------------------------------------------------
; **************** MAIN LOOP ************************
;--------------------------------------------------------
; For each spectrum,
;--------------------------------------------------------

for i = 0,(nspecs - 1) do begin

    ew = FLTARR(nindex)          ; output arrays to store index (& errors)

    if err_flag then begin
        ew_error = FLTARR(nindex)
    endif 

    ;-----------------------------------------------------
    ; Read in spectrum and parameters, plus error spectrum
    ;     if specified
    ;-----------------------------------------------------

    if NOT(KEYWORD_SET(silent)) then begin
        print, 'READING IN SPECTRUM:',i+1
    endif 

    spectrum = READFITS(list[i],specheader,EXTEN_NO=spec_exten,/SILENT)

    if (SIZE(spectrum))[0] gt 1 then begin
        spectrum = spectrum[*,spec_nslice]
    endif

    ;---------------------------------------------------------
    ; Normalize spectrum to prevent arithmetic underflows
    ;---------------------------------------------------------
    
    max_spec = (MAX(spectrum))
    spectrum = spectrum / max_spec
    
    if err_flag then begin
        if N_ELEMENTS(error_list) eq 0 then begin
            errlist = list
        endif 
        if err_nslice gt 0 then begin
            err_spec = READFITS(errlist[i],errheader,EXTEN_NO=err_exten,NSLICE=err_nslice,$
                                /SILENT)
        endif else begin
            err_spec = READFITS(errlist[i],errheader,EXTEN_NO=err_exten,/SILENT)
        endelse 

        if KEYWORD_SET(sn_flag) then begin
            err_spec = spectrum / err_spec
        endif else begin
            err_spec = err_spec / max_spec ; Normalize error spectrum 
        endelse 
    endif 

    spec_params = READSPECINFO(specheader)  ; Get wave1,dwave from header

    specwave1 = spec_params[0]              ; set wavelength parameters
    dspecwave = spec_params[1]              ;   for spectrum

    ;----------------------------------------------
    ; Shift to rest frame, make wavelength array
    ;-------------------------------------------------

    specwave = FINDGEN((size(spectrum))[1])               ; Make wave array
    specwave = specwave*dspecwave + specwave1

    if NOT(KEYWORD_SET(linear)) then begin
        specwave = 10^specwave
    endif

    if (MAX(specwave) lt MIN(red2)) or (MIN(specwave) gt MAX(blue1)) then begin
        print, '*****No Lick Indices covered! Check linear vs.'
        print, '     log-linear wavelength scale.'
    endif 

    ;------------------------------------------------------
    ; For each Lick Index
    ;------------------------------------------------------

    for j = 0, (nindex - 1) do begin

        if ((MIN(specwave) gt blue1[j]) or (MAX(specwave) lt red2[j])) $
           then begin

            ;-----------------------------------------------
            ; If the whole Lick index isn't covered in the 
            ;   spectrum, assign NaN values for its EW and
            ;   error, skip to next index
            ;-----------------------------------------------

            ew[j] = !VALUES.F_NAN

            if err_flag then begin
                ew_error[j] = !VALUES.F_NAN
            endif 

        endif else begin   ; index is covered
        
            temparray = specwave[WHERE(specwave LE (red2[j]+50))]
            subarray = WHERE(temparray GE (blue1[j]-50))

            j_wave = specwave[subarray] ; Make subarray of the spectrum
            j_spec = spectrum[subarray] ;   in wavelength range of index

            if err_flag then begin
                j_err = err_spec[subarray]
            endif 

        ;------------------------------------------------------
        ; Smooth the subarray of the spectrum to the 
        ; Lick/IDS resolution
        ;------------------------------------------------------

            if flag_idsres then begin

                smoothspec = j_spec
                
                if err_flag then begin
                    smooth_err = j_err
                endif 

            endif else begin    ; smoothing needed

            ;-----------------------------------------------------------
            ; Find the size for the gaussian -- we want a gaussian that
            ; extends to 4.0 sigma (99.99%), but it has to have an odd
            ; number of pixels in order to be properly centered.  
            ; Force the gaussian filter size to be 8*sigma rounded to 
            ; the nearest odd number.
            ;-----------------------------------------------------------

                dlambda = (max(j_wave)-min(j_wave))/((size(j_wave))[1]-1)
                lambda = mean(j_wave)

                sig_ind = indres[j] / sqrt(8.0*alog(2.0)) / dlambda
                sig_gal = (FLOAT(vel_disp[i]) / c) * (lambda / dlambda)

                if KEYWORD_SET(pix) then begin
                    sig_res = res[i] / sqrt(8.0*alog(2.0))
                endif else begin
                    sig_res = res[i] / sqrt(8.0*alog(2.0)) / dlambda
                endelse 

                sig_tot = sqrt(sig_res^2 + sig_gal^2)
                
                if (sig_ind le sig_tot) then begin

                    smoothspec = j_spec
                    
                    if err_flag then begin
                       smooth_err = j_err
                   endif 
                
                endif else begin

                    sigma = sqrt(sig_ind^2 - sig_tot^2)
                    filtersize=(FIX(sigma*8)+(FIX(sigma*8+1) mod 2))

                    gaussfilter = PSF_GAUSSIAN(NPIXEL=filtersize, $
                                         ST_DEV=sigma, NDIMEN=1,/NORMALIZE)
            
                    smoothspec = CONVOL(j_spec, gaussfilter, /CENTER, $
                                        /EDGE_TRUNCATE)

                    ;---------------------------------------------
                    ; Smoothing reduces error -- add gaussian-
                    ;    weighted errors from surrounding pixels
                    ;    in quadrature
                    ;---------------------------------------------

                    if err_flag then begin
                        smooth_err = CONVOL((j_err)^2, gaussfilter, /CENTER, $
                                        /EDGE_TRUNCATE)
                        smooth_err = SQRT(smooth_err)
                   endif 
                    
                endelse 
            
            endelse             ; smoothing needed

        ;---------------------------------------------------------
        ; Separate out the continuum and bandpass regions
        ;---------------------------------------------------------

            blue = WHERE((j_wave gt blue1[j])and(j_wave lt blue2[j]))
            red  = WHERE((j_wave gt  red1[j])and(j_wave lt  red2[j]))
            band = WHERE((j_wave gt band1[j])and(j_wave lt band2[j]))
                
            lambda_b = (blue1[j] + blue2[j]) / 2.
            lambda_r = ( red1[j] +  red2[j]) / 2.

            theta_blue = j_wave[blue] - j_wave[blue-1]
            theta_red  = j_wave[red]  - j_wave[red -1]
            theta_band = j_wave[band] - j_wave[band-1]

            bandspec = smoothspec[band]

            ;-----------------------------------------------
            ; Include fractional pixels at end of bandpass
            ;-----------------------------------------------
            
            f_b1 = (j_wave[MIN(blue)] - blue1[j]) /               $
                   (j_wave[MIN(blue)] - j_wave[MIN(blue)-1])
            f_b2 = (blue2[j] - j_wave[MAX(blue)]) /             $
                   (j_wave[MAX(blue)+1] - j_wave[MAX(blue)])

            last_b = (SIZE(blue))[1]-1

            f_r1 = (j_wave[MIN(red)] - red1[j]) /               $
                   (j_wave[MIN(red)] - j_wave[MIN(red)-1])
            f_r2 = (red2[j] - j_wave[MAX(red)]) /             $
                   (j_wave[MAX(red)+1] - j_wave[MAX(red)])

            last_r = (SIZE(red))[1]-1

            ;-----------------------------------------------------
            ; Find "average" point in blue and red pseudocontinua
            ;-----------------------------------------------------

            s_blue = (INT_TABULATED(j_wave[blue], smoothspec[blue])      +  $
                      f_b1 * smoothspec[MIN(blue)-1] * theta_blue[0]     +  $
                      f_b2 * smoothspec[MAX(blue)] * theta_blue[last_b]) /  $
                     (blue2[j] - blue1[j])
            s_red  = (INT_TABULATED(j_wave[red] , smoothspec[red])       +  $
                       f_r1 * smoothspec[MIN(red)-1] * theta_red[0]      +  $
                       f_r2 * smoothspec[MAX(red)] * theta_red[last_r] ) /  $
                      (red2[j]  -  red1[j])

        ;--------------------------------------------------
        ; Make pseudocontinuum.
        ;--------------------------------------------------

            pseudoslope = (s_red-s_blue) / $
                          (lambda_r-lambda_b)
            pseudointercept = s_blue-pseudoslope*lambda_b

            pseudocont = pseudoslope*j_wave[band] + pseudointercept
            wholecont  = pseudoslope*j_wave + pseudointercept

        ;--------------------------------------------
        ; Show Lick Index in Plot, if option is set
        ;--------------------------------------------

            if keyword_set(showplot) then begin

                set_plot,'x'
                !P.MULTI=[0,1,1]

                window, 0
                PLOT, specwave, spectrum, title=indexname[j], $
                      xrange=[min(j_wave[blue])-50, max(j_wave[red])+50], $
                      yrange=[min(j_spec)-0.15*(max(j_spec)-min(j_spec)), $
                              max(j_spec)]
                OPLOT, [red1[j],red2[j]], $
                       [s_red,s_red], $
                       thick=2,color=ctable.green
                OPLOT, [blue1[j],blue2[j]], $
                       [s_blue,s_blue],$
                       thick=2,color=ctable.green
                OPLOT, j_wave, smoothspec, thick=2,color=ctable.red
                OPLOT, [lambda_b,lambda_r],$
                       [s_blue,s_red],$
                       linestyle=2,color=ctable.blue
                OPLOT, j_wave[band], pseudocont, thick=4,color=ctable.blue
                XYOUTS, !X.CRANGE[0]+0.1*(!X.CRANGE[1]-!X.CRANGE[0]), $
                        !Y.CRANGE[0]+0.1*(!Y.CRANGE[1]-!Y.CRANGE[0]), $
                        'Left-click for next index', charsize=1.5
                XYOUTS, !X.CRANGE[0]+0.1*(!X.CRANGE[1]-!X.CRANGE[0]), $
                        !Y.CRANGE[0]+0.05*(!Y.CRANGE[1]-!Y.CRANGE[0]), $
                        'Right-click to skip plotting', charsize=1.5
                POLYFILL,[j_wave[band[0]],j_wave[band],j_wave[REVERSE(band)]],$
                         [pseudocont[0],smoothspec[band],REVERSE(pseudocont)],$
                         /LINE_FILL, color=ctable.gold
                
                
                cursor,xcoord,ycoord,/wait ; Wait for cursor click

                if !MOUSE.button eq 4 then begin 
                    ; If right-click, do not plot again
                    showplot = 0
                endif
                
            endif               ; showplot
            
        ;--------------------------------
        ; Measure Lick Index
        ;-------------------------------

            f_band1 = (j_wave[MIN(band)] - band1[j]) /               $
                   (j_wave[MIN(band)] - j_wave[MIN(band)-1])
            f_band2 = (band2[j] - j_wave[MAX(band)]) /             $
                   (j_wave[MAX(band)+1] - j_wave[MAX(band)])

            last_band = (SIZE(band))[1]-1

            if indmag[j] EQ 0 then begin ; angstrom index
 
                ew[j] = INT_TABULATED(j_wave[band],(1-bandspec/pseudocont))
           
                ew[j] = ew[j] +                                                $
                        f_band1 * ((1-smoothspec[MIN(band)-1] /                $
                                    wholecont[MIN(band)-1]) * theta_band[0] +   $
                                   (1-smoothspec[MIN(band)]   /                $
                                    wholecont[MIN(band)]) * theta_band[0]) / 2. $ 
                        +                                                      $
                        f_band2 * ((1-smoothspec[MAX(band)+1] /                  $
                                    wholecont[MAX(band)+1]) *                  $
                                          theta_band[last_band]             +   $
                                   (1-smoothspec[MAX(band)]   /                $
                                    wholecont[MAX(band)]) *                    $
                                          theta_band[last_band]) / 2.

            endif else begin    ; magnitude index

                integral = INT_TABULATED(j_wave[band],(bandspec/pseudocont))
 
                integral = integral +                                          $
                  f_band1 * ((smoothspec[MIN(band)-1] /                        $
                                    wholecont[MIN(band)-1]) * theta_band[0] +   $
                                   (smoothspec[MIN(band)]   /                $
                                    wholecont[MIN(band)]) * theta_band[0]) / 2. $ 
                  +                                                            $
                  f_band2 * ((smoothspec[MAX(band)+1] /                          $
                                    wholecont[MAX(band)+1]) *                  $
                                          theta_band[last_band]             +   $
                                   (smoothspec[MAX(band)]   /                $
                                    wholecont[MAX(band)]) *                    $
                                          theta_band[last_band]) / 2.
         
                ew[j] = -2.5*alog10(integral/(band2[j]-band1[j]))

            endelse              ; magnitude index

        ;------------------------------------------------------
        ; If error spectra are provided, measure error for
        ;    this index using equations in Cardiel et al. 
        ;    1998, A&A, 127, 597
        ;------------------------------------------------------

            if err_flag then begin

                ; Calculate blue variance including fractional pixels 
                ;   at edge of bandpass

                var_sb = (TOTAL(smooth_err[blue]^2 * theta_blue^2)) / $
                              (blue2[j] - blue1[j])^2

                var_sb = var_sb + (                                        $
                         f_b1 * smooth_err[MIN(blue)-1]^2 * theta_blue[0]^2 +   $
                         f_b2 * smooth_err[MAX(blue)+1]^2 * theta_blue[last_b]^2  $
                         ) / (blue2[j] - blue1[j])^2                        

                ; Calculate red variance including fractional pixels 
                ;   at edge of bandpass

                var_sr = (TOTAL(smooth_err[red]^2 * theta_red^2)) / $
                              (red2[j] - red1[j])^2

                var_sr = var_sr + (                                        $
                         f_r1 * smooth_err[MIN(red)-1]^2 * theta_red[0]^2 +   $
                         f_r2 * smooth_err[MAX(red)+1]^2 * theta_red[last_r]^2  $
                         ) / (red2[j] - red1[j])^2                        

                sum = 0.

                for ii = -1, (SIZE(band))[1] do begin

                    var_c = (lambda_r - j_wave[MIN(band)+ii])^2 / $
                            (lambda_r - lambda_b)^2 * var_sb + $
                            (j_wave[MIN(band)+ii] - lambda_b)^2 / $
                            (lambda_r - lambda_b)^2 * var_sr

                    if ii eq -1 then begin
                        frac_i = f_band1
                        sub_i = 0
                    endif else if ii eq (SIZE(band))[1] then begin
                        frac_i = f_band2
                        sub_i = (SIZE(band))[1] -1
                    endif else begin
                        frac_i = 1.0
                        sub_i = ii
                    endelse 

                    sum = sum +  frac_i *                                    $
                         ( wholecont[MIN(band)+ii]^2 * smooth_err[MIN(band)+ii]^2 $
                           + smoothspec[MIN(band)+ii]^2 * var_c                  $
                          )  / wholecont[MIN(band)+ii]^4 / theta_band[sub_i]^2

                    for jj = -1, (SIZE(band))[1] do begin

                        if jj eq -1 then begin
                            frac_j = f_band1
                            sub_j = 0
                        endif else if jj eq (SIZE(band))[1] then begin
                            frac_j = f_band2
                            sub_j = (SIZE(band))[1] -1
                        endif else begin
                            frac_j = 1.0
                            sub_j = jj
                        endelse 

                        lparam1 = (lambda_r - j_wave[MIN(band)+ii]) * $
                                  (lambda_r - j_wave[MIN(band)+jj]) / $
                                  (lambda_r - lambda_b)^2

                        lparam4 = (j_wave[MIN(band)+ii] - lambda_b) * $
                                  (j_wave[MIN(band)+jj] - lambda_b) / $
                                  (lambda_r - lambda_b)^2

                        if ii ne jj then begin

                            sum = sum + frac_i * frac_j *                        $
                                  smoothspec[MIN(band)+ii]  *                    $ 
                                  smoothspec[MIN(band)+jj]  /                    $
                                  wholecont[MIN(band)+ii]^2                    / $
                                  wholecont[MIN(band)+jj]^2                    * $
                                  (lparam1 * var_sb + lparam4 * var_sr)        / $
                                  (theta_band[sub_i] * theta_band[sub_j])

                        endif    ; ii ne jj
                     
                    endfor       ;jj

                endfor            ; ii

                sigma_index = sqrt(sum)

                if indmag[j] EQ 1 then begin
                    
                    sigma_index = 2.5 * alog10(exp(1)) / 10^(-0.4*ew[j]) /  $
                                  (band2[j] - band1[j]) * sigma_index

                endif 

                ew_error[j] = sigma_index

            endif                ; measuring error

        endelse  

        ;--------------------------------------------
        ; If system has high velocity dispersion,
        ;  apply a velocity correction to the measured
        ;  index -- corrections only go to sigma=300 km/s
        ;  so cap the velocity at that
        ;--------------------------------------------

        if NOT(flag_idsres) then begin

            if (sig_ind lt sig_tot) then begin

                veldisp_eff = sqrt(sig_tot^2 - sig_ind^2) * dlambda / lambda * c

                if veldisp_eff gt 300. then begin
                    veldisp_eff = 300.0
                endif

                new_ew = SIGMA_CORRECT(veldisp_eff, ew[j], $
                          sigcorr_slope_arr[j,*], sigcorr_icpt_arr[j,*])
            
                ew[j] = new_ew
                
            endif 

        endif 

    endfor                       ; Done measuring index [j], measure next index

    printf,ew_tab,FORMAT='(A)','#'+list[i]
    printf,ew_tab,FORMAT='(25(F0,:,"    "))',[ew]

    if err_flag then begin
        printf, ewerr_tab,FORMAT='(A)','#'+list[i]
        printf, ewerr_tab,FORMAT='(25(F0,:,"    "))',[ew_error]
    endif 

endfor                           ; Done with this spectrum [i], get next

;-----------------------------------
; Close output file
;-----------------------------------

close, ew_tab

if err_flag then begin
    close, ewerr_tab
endif

end 


        

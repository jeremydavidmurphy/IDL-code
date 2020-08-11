Pro EZ_Ages, datafile, imf=imf, isochrone=isochrone, o_fe=o_fe, age_tol = age_tol, fe_tol = fe_tol, errorfile=errorfile, si_fe=si_fe, ti_fe=ti_fe, na_fe=na_fe, cr_fe=cr_fe, init_abuns=init_abuns

if N_PARAMS() eq 0 then begin
    print, 'Syntax - EZ_Ages, datafile[, imf=, '
    print, '         isochrone=, o_fe=, na_fe=, si_fe=, cr_fe='
    print, '         ti_fe=, age_tol=, fe_tol=, errorfile='
    print, '         init_abuns= ]'
    return
endif

;------------------------------------------------
; Set up clean exit on error
;------------------------------------------------

CATCH, error_status

if error_status ne 0 then begin
    print, !ERR_STRING
    exit_cleanly, unit_abuntab=unit_abuntab, unit_errtab=unit_errtab
    return
endif 

!EXCEPT=0

;------------------------------------------------
; define structures
;------------------------------------------------

@define_structs

;------------------------------------------------
; Set up deltabund directory, symbolic links
;------------------------------------------------

MAKE_LINKS

;------------------------------------------------
; Prompt for input parameters if not given
;------------------------------------------------

PROMPT_INPUT, imf=imf, isochrone=isochrone, $
  age_tol=age_tol, fe_tol=fe_tol 

;------------------------------------------------
; Set up treatment of abundances not involved
;   in fitting.  If not specified at the command
;   line, set the following defaults: 
;
;         [O/Fe] to match isochrone
;         [Cr/Fe] = 0.0 (Solar)
;         [Na/Fe], [Si/Fe], [Ti/Fe] = [Mg/Fe]
;------------------------------------------------

SET_NOFIT_ABUNS, o_fe=o_fe, na_fe=na_fe, si_fe=si_fe, $
                 cr_fe=cr_fe, ti_fe=ti_fe, mg_link=mg_link,$
                 isochrone=isochrone

;------------------------------------------------
; Set up input structure for deltabund
;------------------------------------------------

deltabund_input = deltabund_input_node
deltabund_input.x = imf
deltabund_input.iso = isochrone
deltabund_input.o = o_fe
deltabund_input.na = na_fe
deltabund_input.si = si_fe
deltabund_input.cr = cr_fe
deltabund_input.ti = ti_fe

;------------------------------------------------
; Get names of input spectra
;------------------------------------------------

READCOL, datafile, spec_names, FORMAT='A', SKIPLINE=1, /SILENT
npoints = (SIZE(spec_names))[1] / 2
spec_names = spec_names[FINDGEN(npoints)*2]

;------------------------------------------------
; Get data points to be fit with models
;------------------------------------------------

data_array = GET_DATA(datafile, npoints)

if N_ELEMENTS(errorfile) gt 0 then begin
    error_array = GET_DATA(errorfile, npoints)
endif

;----------------------------------------------------
; Open output files for plotting and abundance table
;----------------------------------------------------

OPEN_OUTPUT_FILES, datafile=datafile, errorfile=errorfile

;--------------------------------------------------------
; Find best model fit for each data point --- MAIN LOOP!
;--------------------------------------------------------

for n = 0, npoints-1 do begin       ; loop over data points

    print
    print, 'Fitting for data point #', n

    if N_ELEMENTS(errorfile) gt 0 then begin
        DEFINE_ABUN_ERRS, age_error, fe_h_error, mg_fe_error, c_fe_error, $
                          n_fe_error, ca_fe_error
    endif
    
    ;--------------------------------------------------------
    ; When fitting, use the indices selected by the 
    ;     'weights.dat' file.
    ;--------------------------------------------------------

    data_use_n = data_use_node

    if N_ELEMENTS(errorfile) gt 0 then begin

        error_use_n = data_use_node

        WEIGHT_DATA_LINES, data_array[n], data_use_n, $
                           error_set=error_array[n], $
                           error_use_n = error_use_n, $
                           weight_vectors=weight_vectors

    endif else begin

        WEIGHT_DATA_LINES, data_array[n], data_use_n, $
                          weight_vectors=weight_vectors

    endelse 

    ;------------------------------------------------
    ; Zero the model
    ;------------------------------------------------

    if KEYWORD_SET(init_abuns) then begin
        SET_INIT_ABUNS, init_abuns, deltabund_input
    endif else begin
        SET_TO_SOLAR, deltabund_input
    endelse

    SET_NOFIT_ABUNS, o_fe=o_fe, na_fe=na_fe, si_fe=si_fe, $
                     cr_fe=cr_fe, ti_fe=ti_fe, mg_link=mg_link,$
                     isochrone=isochrone

    ;------------------------------------------------------
    ; Spawn a new model and weight to match the data
    ;    lines
    ;------------------------------------------------------

    model_array = MAKE_NEW_MODEL(model_use_array, deltabund_input, mg_link, $
                                 weight_vectors)

    ;--------------------------------------------------
    ; Find the best-fitting model, then loop to make
    ;   sure the fiducial "true" values haven't changed -- as
    ;   long as the abundance has changed during the 
    ;   loop, repeat the loop.
    ;--------------------------------------------------

    change_flag = 1
    count_iter = 0

    while ((change_flag eq 1) and FINITE(data_use_n.balmer) $
                              and FINITE(data_use_n.fe))   do begin
    
        change_flag = 0   ; if we get to the end, and it's still
                          ;  zero, we're done!

        count_iter = count_iter + 1    
        
        ;-----------------------------------------------------
        ; ***
        ; *** Find fiducial "true" age & [Fe/H] from Fe vs. Balmer  ***
        ; ***
        ;-----------------------------------------------------

        print, 'Getting fiducial values for age, [Fe/H]'
        
        if N_ELEMENTS(errorfile) gt 0 then begin
       
            FIND_AGE_FE, data_use_n, model_use_array, age_true, fe_h_true, $
              error_use_n=error_use_n, age_error=age_error, $
              fe_h_error=fe_h_error

        endif else begin
        
            FIND_AGE_FE, data_use_n, model_use_array, age_true, fe_h_true

        endelse

        if FINITE(age_true) then begin 

            ;---------------------------------------------------
            ; ***
            ; *** Next check Mg vs. Balmer fit, adjust [Mg/Fe]
            ; ***
            ;---------------------------------------------------

            if FINITE(data_use_n.mg) then begin

                FIND_MG_FE, data_array, data_use_n, model_use_array, $
                  age_true, fe_h_true, fe_tol, deltabund_input, $
                  mg_link, change_flag, weight_vectors

            endif
            
            ;---------------------------------------------------
            ; ***
            ; *** Next check C fit, adjust [C/Fe]
            ; ***
            ; *** NOTE: If using C4668, combine with Balmer line
            ; ***       If using G4300, combine with Fe line
            ; ***
            ;---------------------------------------------------
            
            if FINITE(data_use_n.c) then begin

                if data_use_n.g4300_flag then begin

                    FIND_C_FE_G4300, data_array, data_use_n, model_use_array, $
                      age_true, age_tol, fe_h_true, deltabund_input, mg_link, $
                      change_flag, c_fit_found, weight_vectors

                endif else begin

                    FIND_C_FE_C4668, data_array, data_use_n, model_use_array, $
                      age_true, fe_h_true, fe_tol, deltabund_input, mg_link, $
                      change_flag, c_fit_found, weight_vectors

                endelse 

            ;---------------------------------------------------
            ; ***
            ; *** Next check CN vs. Balmer fit, adjust [N/Fe]
            ; ***
            ;---------------------------------------------------

                if (FINITE(data_use_n.n) and c_fit_found) then begin

                    FIND_N_FE, data_array, data_use_n, model_use_array, $
                      age_true, fe_h_true, fe_tol, deltabund_input, mg_link, $
                      change_flag, n_fit_found, weight_vectors

            ;---------------------------------------------------
            ; ***
            ; *** Next check Ca vs. Balmer fit, adjust [Ca/Fe]
            ; ***
            ;---------------------------------------------------

                    if (FINITE(data_use_n.ca) and n_fit_found) then begin

                        FIND_CA_FE, data_array, data_use_n, model_use_array, $
                          age_true, fe_h_true, fe_tol, deltabund_input, $
                          mg_link, change_flag, weight_vectors

                    endif       ; We have Ca line
                endif           ; We have CN line
            endif               ; We have C line
           
        endif                   ; Data is on Fe vs. Balmer plot

        ;----------------------------------------
        ; If we've already iterated 4 times, 
        ;  assume we can't do better and move on!
        ;----------------------------------------

        if count_iter eq 4 then begin
            print, '4 iterations -- move to next point!'
            change_flag = 0
        endif

        ;---------------------------------------
        ;  **** Done fitting *****
        ;    Loop again to double-check fit,
        ;    unless nothing has changed.
        ;---------------------------------------
 
    endwhile      ; check to see if model changed
      
    ;------------------------------------------------
    ; If index errors are available, calculate
    ;    errors in the derived ages and abundances
    ;------------------------------------------------

    if N_ELEMENTS(errorfile) gt 0 then begin

        ESTIMATE_ERRORS, data_array, data_use_n, model_use_array, $
                     age_true, fe_h_true, error_use_n, weight_vectors, $
                     age_tol, fe_tol, deltabund_input, mg_link, $
                     age_error=age_error, fe_h_error=fe_h_error, $
                     mg_fe_error=mg_fe_error, c_fe_error=c_fe_error, $
                     n_fe_error=n_fe_error, ca_fe_error=ca_fe_error

    endif                       ; Measure errors

    ;-------------------------------------------
    ; Now that we've found the best fit, make
    ;   plots so we can verify visually.
    ;
    ;   ---> Remake the best fitting model if it 
    ;        has been modified to estimate erros.
    ;-------------------------------------------

    if N_ELEMENTS(error_use_n) eq 0 then begin

        PLOT_GRIDS, data_array[n], data_use_n, isochrone, spec_names[n]

    endif else begin
        
        WRITE_INPUTFILE, deltabund_input, mg_link
        SPAWN, '$EZ_AGES_DIR/deltabund_lector_batch'

        PLOT_GRIDS_ERR, data_array[n], data_use_n, isochrone, $
          spec_names[n], error_set_n=error_array[n], $
          error_use_n=error_use_n

    endelse

    ;-------------------------------------------------------------
    ; Get ages from all three Balmer lines (Hb, HgF, HdF), 
    ;     if available, with or w/o errors
    ;-------------------------------------------------------------

    if N_ELEMENTS(error_use_n) eq 0 then begin

        GET_BALMER_AGES, data_array[n], data_use_n, $
                         model_use_array, model_array, $
                         age_true, age_hb, age_hg, age_hd,$
                         deltabund_input, mg_link, weight_vectors

    endif else begin  

        GET_BALMER_AGES, data_array[n], data_use_n, $
                         model_use_array, model_array, $
                         age_true, age_hb, age_hg, age_hd, $
                         deltabund_input, mg_link, weight_vectors, $
                         error_use_n=error_use_n, $
                         error_array=error_array[n], $
                         age_error=age_error, age_hb_error=age_hb_error, $
                         age_hg_error=age_hg_error, $
                         age_hd_error=age_hd_error, $
                         c_fe_error=c_fe_error

    endelse   ; Get Balmer ages if we're measuring errors

    ;--------------------------------------------------
    ; Write to output file: Age, [Fe/H], all abundaces, 
    ;   Age(Hb), Age(Hg), Age(Hd), isochrone, x(imf)
    ; 
    ; If we're measuring errors, also write an output
    ;   error file.
    ;--------------------------------------------------

    WRITE_OUTPUT, age_true, fe_h_true, deltabund_input, mg_link, $
      isochrone, imf, age_hb, age_hg, age_hd, spec_names[n], $
      errorfile=errorfile, age_error=age_error, fe_h_error=fe_h_error, $
      mg_fe_error=mg_fe_error, c_fe_error=c_fe_error, $
      n_fe_error=n_fe_error, ca_fe_error=ca_fe_error, $
      age_hb_error=age_hb_error, $
      age_hg_error=age_hg_error, age_hd_error=age_hd_error

endfor        ; loop over data points

;-------------------------------------------------
; Close output files
;-------------------------------------------------

CLOSE_OUTPUT_FILES, errorfile=errorfile

;-------------------------------------------------
; Clean up Fortran-generated files
;-------------------------------------------------

spawn,'rm -f fort.98'
spawn,'rm -f fort.92'
spawn,'rm -f fort.91'

end           ; program!




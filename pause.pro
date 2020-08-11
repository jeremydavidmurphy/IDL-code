pro PAUSE
    ON_ERROR, 2

    prompt_save = !prompt
    rp = ''
    read,rp,prompt='Press ENTER to continue...'
    rp = strupcase(rp)
    !prompt=prompt_save
    if rp eq 'Q' then begin
         close,/all
         message,'USER-BREAK'
    endif
end

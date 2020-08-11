;This routine is used to plot the rms of various levels of sky
;subtraction. The "run" file is the actual run file used in Hoku. The
;rms file is the saved screen output. This is now piped into an output
;named rmsall.out

;This routine makes no assumptions about how many files in a row you
;want to consider- it can be iterative. The minimum rms value for each
;user defined block is written into a file that has the same name as
;the rms file, only with ".min" appended to the end.

;THIS CODE NOW READS THE FRAME NUMBERS AND BREAKS UP THE RMS ANALYSIS
;ON THIS NUMBER (the #### in a frame named jm####bcTTfc.fits)

;------------------------------------------------------
PRO plotrmsky, runlist, rmslist
;------------------------------------------------------

;THIS CODE CAN NOT HANDLE RUNLISTS WITH DUPLICATE ####'S IN THE SAME
;FILE, MEANING, YOU MUST RUN YOUR DIFFERENT SPECTRAL REGIONS SEPARATELY.

iora = ''
print,'Run the routine interactively or automatically?'
print,'Enter "i" or "a"'
read,iora

readcol,runlist,format='x,a,x,x,x,x,x,a',files,regions
readcol,rmslist,format='d,x,x,x',rmsvalues

t0 = n_elements(files)
t1 = n_elements(rmsvalues)

if (t0 ne t1) then begin
    print,'The RMS file and the RUN file are not the same length!'
    STOP
endif

;------------------------------------------------------------------
;The correction for the different levels of sky-subtraction is made
;------------------------------------------------------------------

idarr = strarr(t0)
corarr = fltarr(t0)
namearr = strarr(2,t0)
for j=0,t0-1 do begin
    p = strsplit(files[j],'jmvp0123456789Tf',/extract)
    pp = strsplit(files[j],'jmvpabcxyzTf',/extract)
    idarr[j] = p[0]
    namearr[0,j] = pp
    if (n_elements(p) gt 1) then namearr[1,j] = p[1]
endfor

i1 = where(idarr eq 'cc')
i2 = [where(idarr eq 'cb'), where(idarr eq 'bc')] 
i3 = [where(idarr eq 'ca'), where(idarr eq 'ac'), where(idarr eq 'bb')]
i4 = [where(idarr eq 'ba'), where(idarr eq 'ab')]
i5 = where(idarr eq 'aa')
i6 = where(idarr eq 'zz')
i7 = [where(idarr eq 'zy'), where(idarr eq 'yz')]
i8 = [where(idarr eq 'zx'), where(idarr eq 'xz'), where(idarr eq 'yy')]
i9 = [where(idarr eq 'yx'), where(idarr eq 'xy')]
i10 = where(idarr eq 'xx')

;*****************************************************************
stop
; CHECK THIS CODE BEFORE USING IT AGAIN! It looks like the rmsvalues
; array is a different size than the corarr array. (a t0 vs t1)
;*****************************************************************

corarr[i1] = 4/(2.15+2.15)
corarr[i2] = 4/(2.00+2.15)
corarr[i3] = 4/(2.00+2.00)
corarr[i4] = 4/(2.00+1.85)
corarr[i5] = 4/(1.85+1.85)
corarr[i6] = 4/(1.70+1.70)
corarr[i7] = 4/(1.70+1.60)
corarr[i8] = 4/(1.60+1.60)
corarr[i9] = 4/(1.60+1.50)
corarr[i10] = 4/(1.50+1.50)

rmsvalues = rmsvalues*corarr
;------------------------------------------------------------------

ans = ''
print,'Plot the output?'
read,ans

if (iora eq 'i') then begin
   c1 = 0
   c2 = 0
   print,'Enter the maximum length of a section:'
   read,c2
   step = c2
   done = 0
   repeat begin
      for j=c1,c2 do begin
         print,j,'   ',files[j]
      endfor
      print,'Enter a LAST element:'
      read,c2
      print,files[c1:c2]
      n1 = uint(c2-c1)+1
      xaxis = findgen(n1)+1
      
;---------------------------------------------------------------
;An addition that outputs the SORTED file names and rms values.
;--------------------------------------------------------------- 
      piece = namearr[*,c1:c2]
      rmsp = rmsvalues[c1:c2]
      indy = bsort(rmsp)
      stdrms = rmsp[indy]
      filesp = files[c1:c2]
      stdnames = filesp[indy]
      tmane = strsplit(stdnames[0],'abcxyz',/extract)
      if (namearr[1,c1] eq 'c') then back = 'rmsc.min'
      if (namearr[1,c1] eq '') then back = 'rms.min'
      rr = regions[c1]
      tmane = tmane[0]+'_'+rr+back
      n0 = n_elements(filesp)
      print,tmane
      openw,5,tmane
      for j=0,n0-1 do printf,5,stdnames[j],stdrms[j]
      free_lun,5
;---------------------------------------------------------------
       
      winner = where(rmsvalues[c1:c2] eq min(rmsvalues[c1:c2]))
      mins = piece[winner]
      
      if (ans eq 'y') then begin
         yup = max(rmsvalues[c1:c2])+0.002
         ydown = min(rmsvalues[c1:c2])-0.003
         xup = max(xaxis)+1
         xdown = min(xaxis)-1
         window,0,retain=2
         device,decomposed=0
         loadct,0
         plot,xaxis,rmsvalues[c1:c2],psym=-2,xrange=[xdown,xup],$
              yrange=[ydown,yup],xstyle=1,ystyle=1,title='Sky Subtraction RMS Values',$
              ytitle='RMS',charsize=1.5
         xyouts,xaxis,ydown,files[c1:c2],orientation=90,charsize=1.5
         loadct,4
         xyouts,xaxis[winner],ydown,piece[winner],color=150,charthick=2,$
                orientation=90,charsize=1.5
      endif
      
      c1 = c2+1
      c2 = c1 + step -1
      if (c2 ge t1-1) then begin
         c2 = t1-1 
         done = done+1
      endif
      
   endrep until (done eq 2)
   
endif

;-------------------------------------------------------------
;THE AUTOMATIC SECTION
;-------------------------------------------------------------
i1 = 0
if (iora eq 'a') then begin
    repeat begin
        n1 = namearr[0,i1]
        lindx = where(namearr[0,*] eq n1)
;        print,''
;        print,'the lindx is:'
;        print,lindx
; jm0083ccTTfc     0.037766507

;        pause
;        l1 = lindx[0]
;        cntr = 0
;         repeat begin
;            l1 = lindx[cntr]
;            l2 = lindx[cntr+1]
;            eval = l2-l1
;            cntr = cntr + 1
;        endrep until (eval gt 1)
;        print,lindx[0:cntr-1] & pause
;        sindx = lindx[0:cntr-1]
        n = n_elements(lindx) - 1
        i1 = lindx[0]
        i2 = lindx[n]
        print,'The files considered in the RMS minimization are:'
        print,files[i1:i2]
        piece = files[i1:i2]
        rmsp = rmsvalues[i1:i2]
        indy = bsort(rmsp)
        stdrms = rmsp[indy]
        filesp = files[i1:i2]
        stdnames = filesp[indy]
        tmane = strsplit(stdnames[0],'abcxyz',/extract)
        if (namearr[1,i1] eq 'c') then back = '_rmsc.min'
        if (namearr[1,i1] eq '') then back = '_rms.min'
        rr = regions[i1]
        tmane = tmane[0]+'_'+rr+back
        n0 = n_elements(filesp)
        openw,5,tmane
        for j=0,n0-1 do printf,5,stdnames[j],stdrms[j]
        free_lun,5
        
        winner = where(rmsvalues[i1:i2] eq min(rmsvalues[i1:i2]))
        mins = piece[winner]
        
        nn1 = uint(i2-i1)+1
        xaxis = findgen(nn1)+1
        
        if (ans eq 'y') then begin
            yup = max(rmsvalues[i1:i2])+0.002
            ydown = min(rmsvalues[i1:i2])-0.003
            xup = max(xaxis)+1
            xdown = min(xaxis)-1
            window,0,retain=2
            plot,xaxis,rmsvalues[i1:i2],psym=-2,xrange=[xdown,xup],$
              yrange=[ydown,yup],xstyle=1,ystyle=1,title='Sky Subtraction RMS Values',$
              ytitle='RMS',charsize=1.5
            xyouts,xaxis,ydown,files[i1:i2],orientation=90,charsize=1.5
            loadct,4
            xyouts,xaxis[winner],ydown,piece[winner],color=150,charthick=2,$
              orientation=90,charsize=1.5
            pause
        endif
        
        i1 = i2+1

    endrep until (i1 eq t0)
wdelete,0
endif

stop
end

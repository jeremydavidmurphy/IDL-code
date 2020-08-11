;THIS VERSION NOW ACCEPTS UNCOLLAPSED DATA, AND THE WEIGHT
;FILES, THE BIN LISTS, AND RETURNS THE FINAL BINNED DATA AND THE
;CORRESPONDING SIGNAL TO NOISE CALCULATION.

;This is the second half of the pipeline and assumes you've run
;PIPE1.PRO in the same directory you're running PIPE2.PRO
;from. This is necessary as PIPE1.PRO generates the name.R files
;required by PIPE2.PRO. This procedure works from lists with the
;ultimate goal of returning the final, aligned, and combined file for
;a given pointing. Pipe2 is a reworking of alignall.pro, realign.pro
;and the combiwgt.pro procedures. 

;*****************************************************************
pro pipe2, LIST=listoflists, BIN=yourbinlist, WGTYN=wgtyn
;*****************************************************************

;***************************************************************************
;The initial wavelength and dispersion values for the linear interpolation
wi = 3545.0 ;(M87:3545.0)
wdisp = 1.125 ;(M87:1.125)
;wdisp = 0.4
;The size of the fiber profile, in pixels.
;nprofile = 7
nprofile = 5
step = (nprofile-1)/2
plotswitch = 'off'
;plotswitch = 'on'
;***************************************************************************

;The listoflists is just that. A list of different data lists. The
;only purpose of this is to allow different mixes of sky-subtraction
;to be run sequentially and without supervision.
;Each of the datalists must have the form below. The listoflists will
;look something like
;  M87cc.list
;  M87bc.list
;  etc.

;where each line is the name of a text file. Each of these text files
;must have the structure shown below.

;If the keyword call isn't used, the routine searches for a
;list called list.list 


;The datalist:
;THE FORMAT FOR THE DATALIST (i.e. M87cc.list)
;The list is composed of 4 columns. The first is the data
;file, sans any pefsm.fits. The second indicates the pointing,
;corresponding to the binlist letters. The third column is the
;corresponding ptow.dat AND mask.dat files. This will be the name you
;gave VACCINE which gets used in the naming convention for the data
;for that run. PIPE2.PRO will append the necessary 'mask' and 'ptow'
;as needed. The 4th column are the "name.R" files that are output from 
;PIPE1.PRO. 

;This is a list of the form:
;  jm0123cb_   a   _jan08_n2.dat    jm0123.R
;  jm0125cd_   a   _jan08_n2.dat    jm0125.R
;  etc.

;This list must include ALL THE DATA for a given galaxy, INCLUDING ALL THE
;POINTINGS, that you want to go into the final data. This is 
;necessary due to the overlapping pointings.

;YOURBINLIST:
;This is just a list of the bin lists. So, the bin.list should look
;like,
;  bin1601
;  bin1602
;  etc.
;With bin1601 being a text file of the form 
;  203_a
;  218_a
;  1_c
;  2_c
;  etc.

;IF THE BIN KEYWORD IS NOT USED THE CODE WILL LOOK FOR A FILE NAMED BIN.LIST

;STEPS OF THE REDUCTION:
; All the data and weight arrrays are read in. The weight arrays
; contain values outside of 0 and 1. These are generally data of very
; poor quality, and so the weights are set equal to zero.

; The data is run through realign.pro, which makes the linearization
; of the wavelength. This runs on both the data and weight array.

; All the data and weights for a given bin are sorted out. Currently
; the code is hardwired to accept 5 pointings, maximum. Increasing
; this is a simple manner of copying over a complete chunk of the
; a,b,c,etc. code.

; WGTYN: This sets whether the weighted biweight or the straight
; biweight is run. Enter 'wgtn' for no weights. The default is to NOT use
; the weighted biweight.

ans = ''

if (n_elements(listoflists) eq 0) then listoflists = 'list.list'
readcol,silent=1,listoflists,format='A',allists

;The bins are read in
if (n_elements(yourbinlist) eq 0) then readcol,silent=1,'bin.list',format='a',binlist $
else readcol,silent=1,yourbinlist,format='a',binlist
n0 = n_elements(binlist)

if (n_elements(wgtyn) eq 0) then wgtyn = 'wgtn'

for p=0,n_elements(allists)-1 do begin ;a loop through each of the lists.
    onelist = allists[p]
    print,''
    print,'Working on list '+onelist+'...'
    print,''
    
;the parameter file is read in, with all the fits files to be
;combined. The size of the data array is determined.
    readcol,silent=1,onelist,format='a,a,a,a',files,pointing,dat,dotR
    test = readfits(files[0]+'pefsm.fits',/silent)
    n1 = n_elements(test[*,0]) ;N1 IS WAVELENGTH
    n2 = n_elements(test[0,*]) ;N2 IS HEIGHT OF ARRAY
    n3 = n_elements(files)
    
;this needs to be float rather than double. As all the data
;necessarily gets read in at once, galaxies with several pointings run
;into the memory limit. 
    if (p eq 0) then darr = fltarr(n1,n2,n3)
    if (p eq 0) then woffset = strarr(n3)
    if (p eq 0) then warr = darr
    
;the size of the dotR array is determined. this contains which fibers
;to mask and the dRA and dDec offsets.
    readcol,silent=1,dotR[0],format='i,x,x,x,x,x',temp
    n4 = n_elements(temp)
    fibmaskarr = intarr(n3,n4)
    fibradarr = fltarr(n3,n4)
    dRAarr = fltarr(n3,n4)
    dDecarr = fltarr(n3,n4)
    lumarr = fltarr(n3,n4)
;The data, and weight files are read in.
    for j=0,n3-1 do begin
        print,'Reading in file ',j+1
        darr[*,*,j] = readfits(files[j]+'pefsm.fits',/silent,h)
        hn = n_elements(h)
        cntr = 0
        out = 'lost'
        repeat begin
            hh = strsplit(h[cntr],' ',/extract)
            if (hh[0] eq 'WAVEZP') then out = 'found' else cntr = cntr + 1
        endrep until (out eq 'found')
        woffset[j] = float(hh[2])
        print,'The wavelength offset is '+hh[2]
        warr[*,*,j] = readfits(files[j]+'pefw.fits',/silent)
    endfor
    
;the dotR array is read in. it makes it into three different
;arrays. the fibmaskarr array, the dRAarr array and dDecarr array.
    for j=0,n3-1 do begin
        readcol,silent=1,dotR[j],format='i,f,f,f,f,x',fiber,fibrad,dra,ddec,lum
        bfi = where(fiber eq -1) ;The bad fibers from PIPE1 are flagged.
        readcol,silent=1,'mask'+dat[j],format='i,x',y
        y = y-1  ;it's now an index
        ny = n_elements(y)
        if (ny ne n4) then begin ;this is just done to get the files to the same length
            dif = n4-ny
            attach = intarr(dif)
            y = [y,attach]
            fibrad = [fibrad,attach]
            dra = [dra,attach]
            ddec = [ddec,attach]
            lum = [lum,attach]
        endif
        fibmaskarr[j,*] = y
        dRAarr[j,*] = dra
        dDecarr[j,*] = ddec
        lumarr[j,*] = lum
        fibradarr[j,*] = fibrad

;FIBMASKARR is the INDEX of the CENTER of the fiber.
;The fibers get the same mask values given to them by Vaccine. The
;step value is used to mask the data for all the rows for a given
;fiber.
        if (bfi[0] ne -1) then begin
            for k=0,n_elements(bfi)-1 do begin
                darr[*,y[bfi[k]]-step:y[bfi[k]]+step,j] = -666
                warr[*,y[bfi[k]]-step:y[bfi[k]]+step,j] = 0.0
            endfor
        endif
    endfor
   
;the wavelength linearization is made on both the data and the weight
;array. this now handles the issue with interpolating over the -666
;flags

;a labor intensive test. toss after a go round
;i = where(darr le -666,count)
;print,'The percentage of pixels tossed in the data set is: 'strn(float(count/n1*n2*n3)*100)

    for k=0,n3-1 do begin
        print,'Realigning file ',k+1
        mask = 'mask'+dat[k]
        ptow = 'ptow'+dat[k]
        darr[*,*,k] = realignf(darr[*,*,k],ptow,mask,wi,wdisp,step,woffset[k])
        warr[*,*,k] = realignf(warr[*,*,k],ptow,mask,wi,wdisp,step,woffset[k])
    endfor
    
; The weights in the weight array are not masked to zero. This is done here.
    i = where(darr le -666,count)
    if (count ne 0) then begin
;        print,'The number of -666 flags in the data is '+strn(count)
        darr[i] = -666
        warr[i] = 0.0
    endif
    i = where(warr lt 0.0,count)
    if (count ne 0) then begin
;        print,'The number of -0.0 weights in the weight array is '+strn(count)
        darr[i] = -666
        warr[i] = 0.0
    endif
    i = where(warr gt 1.0,count)
    if (count ne 0) then begin
;        print,'The number of +1.0 weights in the weight array is '+strn(count)
        darr[i] = -666
        warr[i] = 0.0
    endif

;-----------------------------------------------------------------------------
;This section is hard-coded for a maximum of 5 different pointings on
;one galaxy. This will need to be rewritten if more than 5 exist.
;-----------------------------------------------------------------------------

;Description of the indexes...
;ipa, ipb,... The index in the data array (dimension 3) for a given
;pointing (a,b,c,etc.) 
;iba, ibb,... The index which dictates which elements in the bin
;arrays correspond to a given letter. These are then used to get the
;actual fiber index, the ifa's
;ifa, ifb,... The index for the fiber for a specific pointing to
;include. Again, this is an intermediate step and is used in the
;fibmaskarr to get the actual index for the data array.
;ima, imb,... The fiber array numbers. Feed these into the data or
;weight arrays and you are at the center of the fiber for a given pointing.

;the indexs for the pointings are set outside the loop as these don't change.
    ipa = where(pointing eq 'a')
    ipb = where(pointing eq 'b')
    ipc = where(pointing eq 'c')
    ipd = where(pointing eq 'd')
    ipe = where(pointing eq 'e')
    
    for j=0,n0-1 do begin
        binname = binlist[j]
        print,''
        print,'Working on '+binname+'...'
        readcol,silent=1,binname,format='a',fibs
        bin = strarr(2,n_elements(fibs))
        for k=0,n_elements(fibs)-1 do bin[*,k] = strsplit(fibs[k],'_',/extract)
        
;all indexes and data arrays are reset for each new bin
        iba = -1 & ifa = -1 & ima = -1
        ibb = -1 & ifb = -1 & imb = -1
        ibc = -1 & ifc = -1 & imc = -1
        ibd = -1 & ifd = -1 & imd = -1 
        ibe = -1 & ife = -1 & ime = -1 
        d1=0 & d2=0 & d3=0 & d4=0 & d5=0
        
        alldata = fltarr(n1)
        allwgts = fltarr(n1)
        allwgtsn = fltarr(n1)
        alloffsets = [0.0,0.0,0.0,0.0,0.0]
        fibnorm = fltarr(1)
        
;-------------------------------------------------------------------------------
        
        iba = where(bin[1,*] eq 'a')
        ibb = where(bin[1,*] eq 'b')
        ibc = where(bin[1,*] eq 'c')
        ibd = where(bin[1,*] eq 'd')
        ibe = where(bin[1,*] eq 'e')

;if the pointing is included in the bin, then the ifa array is set,
;giving the coordinates for the CENTERS of each of the fibers as per
;the mask file.
        if (iba[0] ne -1) then begin
            print,'Num of fibers from pointing A: ',n_elements(iba)
            ifa = bin[0,iba]-1
        endif
        if (ibb[0] ne -1) then begin
            print,'Num of fibers from pointing B: ',n_elements(ibb)
            ifb = bin[0,ibb]-1
        endif
        if (ibc[0] ne -1) then begin
            print,'Num of fibers from pointing C: ',n_elements(ibc)
            ifc = bin[0,ibc]-1
        endif
        if (ibd[0] ne -1) then begin
            print,'Num of fibers from pointing D: ',n_elements(ibd)
            ifd = bin[0,ibd]-1
        endif
        if (ibe[0] ne -1) then begin
            print,'Num of fibers from pointing E: ',n_elements(ibe)
            ife = bin[0,ibe]-1
        endif
      
;if the pointing exists AND the bin requires that pointing, the data
;is placed into an array called DA, DB, etc. and the weights are
;placed into arrays called WA, WB, etc.
      
        if (ipa[0] ne -1) and (iba[0] ne -1) then begin
            for k=0,n_elements(ipa)-1 do begin ;a loop over each pointing
                ima = fibmaskarr[ipa[k],ifa] ;the index for THE CENTER of all fibers in this pointing
                for l=0,n_elements(ima)-1 do begin ;a loop over each fiber
                    if (k eq 0 and l eq 0) then begin
                        da = darr[*,ima[l]-step:ima[l]+step,ipa[k]]
                        wa = warr[*,ima[l]-step:ima[l]+step,ipa[k]]
                        wan = wa/max(wa) ;normalized to compete with the fiber weights
                        fibnorm = [fibnorm,median(da)]
                    endif else begin
                        dpa = darr[*,ima[l]-step:ima[l]+step,ipa[k]]
                        da = [[da],[dpa]]
                        wpa = warr[*,ima[l]-step:ima[l]+step,ipa[k]]
                        wpan = wpa/max(wpa) ;normalized to compete with the fiber weights
                        wa = [[wa],[wpa]]
                        wan = [[wan],[wpan]]
                        fibnorm = [fibnorm,median(dpa)]
                    endelse
                endfor
            endfor
            d1 = n_elements(da[0,*])
                                ;the data and weights get concatenated with the other pointings
            alldata = [[alldata],[da]]
            allwgts = [[allwgts],[wa]]
            allwgtsn = [[allwgtsn],[wan]]

            for k=0,n_elements(ifa)-1 do begin
                piece = [ifa[k]+1,dRAarr[ipa[0],ifa[k]],dDecarr[ipa[0],ifa[k]],$
                         fibradarr[ipa[0],ifa[k]],lumarr[ipa[0],ifa[k]]]
                alloffsets = [[alloffsets],[piece]]
            endfor
        endif
        
        if (ipb[0] ne -1) and (ibb[0] ne -1) then begin
            for k=0,n_elements(ipb)-1 do begin
                imb = fibmaskarr[ipb[k],ifb]
                for l=0,n_elements(imb)-1 do begin
                    if (k eq 0 and l eq 0) then begin
                        db = darr[*,imb[l]-step:imb[l]+step,ipb[k]]
                        wb = warr[*,imb[l]-step:imb[l]+step,ipb[k]]
                        wbn = wb/max(wb)
                        fibnorm = [fibnorm,median(db)]
                    endif else begin
                        dpb = darr[*,imb[l]-step:imb[l]+step,ipb[k]]
                        db = [[db],[dpb]]
                        wpb = warr[*,imb[l]-step:imb[l]+step,ipb[k]]
                        wpbn =wpb/max(wpb)
                        wb = [[wb],[wpb]]
                        wbn = [[wbn],[wpbn]]
                        fibnorm = [fibnorm,median(dpb)]
                    endelse
                endfor
            endfor
            d2 = n_elements(db[0,*])
            alldata = [[alldata],[db]]
            allwgts = [[allwgts],[wb]]
            allwgtsn = [[allwgtsn],[wbn]]
            for k=0,n_elements(ifb)-1 do begin
                piece = [ifb[k]+1,dRAarr[ipb[0],ifb[k]],dDecarr[ipb[0],ifb[k]],$
                         fibradarr[ipb[0],ifb[k]],lumarr[ipb[0],ifb[k]]]
                alloffsets = [[alloffsets],[piece]]
            endfor
        endif
        
        if (ipc[0] ne -1) and (ibc[0] ne -1) then begin
            for k=0,n_elements(ipc)-1 do begin
                imc = fibmaskarr[ipc[k],ifc]
                for l=0,n_elements(imc)-1 do begin
                    if (k eq 0 and l eq 0) then begin
                        dc = darr[*,imc[l]-step:imc[l]+step,ipc[k]]
                        wc = warr[*,imc[l]-step:imc[l]+step,ipc[k]]
                        wcn = wc/max(wc)
                        fibnorm = [fibnorm,median(dc)]
                    endif else begin
                        dpc = darr[*,imc[l]-step:imc[l]+step,ipc[k]]
                        dc = [[dc],[dpc]]
                        wpc = warr[*,imc[l]-step:imc[l]+step,ipc[k]]
                        wpcn = wpc/max(wpc)
                        wc = [[wc],[wpc]]
                        wcn = [[wcn],[wpcn]]
                        fibnorm = [fibnorm,median(dpc)]
                    endelse
                endfor
            endfor
            d3 = n_elements(dc[0,*])
            alldata = [[alldata],[dc]]
            allwgts = [[allwgts],[wc]]
            allwgtsn = [[allwgtsn],[wcn]]
            for k=0,n_elements(ifc)-1 do begin
                piece = [ifc[k]+1,dRAarr[ipc[0],ifc[k]],dDecarr[ipc[0],ifc[k]],$
                         fibradarr[ipc[0],ifc[k]],lumarr[ipc[0],ifc[k]]]
                alloffsets = [[alloffsets],[piece]]
            endfor
        endif
        
        if (ipd[0] ne -1) and (ibd[0] ne -1) then begin
            for k=0,n_elements(ipd)-1 do begin
                imd = fibmaskarr[ipd[k],ifd]
                for l=0,n_elements(imd)-1 do begin ;a loop over each fiber
                    if (k eq 0 and l eq 0) then begin
                        dd = darr[*,imd[l]-step:imd[l]+step,ipd[k]]
                        wd = warr[*,imd[l]-step:imd[l]+step,ipd[k]]
                        wdn = wd/max(wd)
                        fibnorm = [fibnorm,median(dd)]
                    endif else begin
                        dpd = darr[*,imd[l]-step:imd[l]+step,ipd[k]]
                        dd = [[dd],[dpd]]
                        wpd = warr[*,imd[l]-step:imd[l]+step,ipd[k]]
                        wpdn = wpd/max(wpd)
                        wd = [[wd],[wpd]]
                        wdn = [[wdn],[wpdn]]
                        fibnorm = [fibnorm,median(dpd)]
                    endelse
                endfor
            endfor
            d4 = n_elements(dd[0,*])
            alldata = [[alldata],[dd]]
            allwgts = [[allwgts],[wd]]
            allwgtsn = [[allwgtsn],[wdn]]
            for k=0,n_elements(ifd)-1 do begin
                piece = [ifd[k]+1,dRAarr[ipd[0],ifd[k]],dDecarr[ipd[0],ifd[k]],$
                         fibradarr[ipd[0],ifd[k]],lumarr[ipd[0],ifd[k]]]
                alloffsets = [[alloffsets],[piece]]
            endfor
        endif
        
        if (ipe[0] ne -1) and (ibe[0] ne -1) then begin
            for k=0,n_elements(ipe)-1 do begin ;a loop over each pointing
                ime = fibmaskarr[ipe[k],ife]
                for l=0,n_elements(ime)-1 do begin ;a loop over each fiber
                    if (k eq 0 and l eq 0) then begin
                        de = darr[*,ime[l]-step:ime[l]+step,ipe[k]]
                        we = warr[*,ime[l]-step:ime[l]+step,ipe[k]]
                        wen = we/max(we)
                        fibnorm = [fibnorm,median(de)]
                    endif else begin
                        dpe = darr[*,ime[l]-step:ime[l]+step,ipe[k]]
                        de = [[de],[dpe]]
                        wpe = warr[*,ime[l]-step:ime[l]+step,ipe[k]]
                        wpen = wpe/max(wpe)
                        we = [[we],[wpe]]
                        wen = [[wen],[wpen]]
                        fibnorm = [fibnorm,median(dpe)]
                    endelse
                endfor
            endfor
            d5 = n_elements(de[0,*])
            alldata = [[alldata],[de]]
            allwgts = [[allwgts],[we]]
            allwgtsn = [[allwgtsn],[wen]]
            for k=0,n_elements(ife)-1 do begin
                piece = [ife[k]+1,dRAarr[ipe[0],ife[k]],dDecarr[ipe[0],ife[k]],$
                         fibradarr[ipe[0],ife[k]],lumarr[ipe[0],ife[k]]]
                alloffsets = [[alloffsets],[piece]]
            endfor
        endif
        
;now the initial zeros in the arrays are removed.
;if there are no fibers that belong in a given bin (as when you are
;rerunning a single pointing) then this skips over to the next bin.

        if (n_elements(alldata[0,*]) ge 5) then begin
            alldata = alldata[*,1:*]
            allwgts = allwgts[*,1:*]
            allwgtsn = allwgtsn[*,1:*]
            alloffsets = alloffsets[*,1:*]
            fibnorm = fibnorm[1:*]
        endif else begin
            print,'There is no data for '+binname+' in the list '+onelist
            print,'Skipping to the next bin...'
            goto, jump1
        endelse

;-------------------------------------------------------------------------------
; The bad fibers still exist in the arrays. This is a necessary step
; as the fiber indexing done above gets confused if you REMOVE the bad
; fibers from the data. But, now they are just baggage and are tossed
; from ALL THE ARRAYS (this includes the dRA and dDec arrays)
;-------------------------------------------------------------------------------

        flags = median(alldata,dimension=1)
        index = where(flags ne -666,count)
        print,''
        print,'The number of total rows in '+binname+' is: '+strn(n_elements(alldata[0,*]))
        print,''
        print,'The number of rows kept, after rejecting bad fibers is: '+strn(count)
        print,''
        
        dt = count

        alldata = alldata[*,index]
        allwgts = allwgts[*,index]
        allwgtsn = allwgtsn[*,index]
        fibnorm = fibnorm[index]
        
        fibnorm = fibnorm/max(fibnorm) ;the weights get normalized to 1.
        
;the fibnorm array is stretched to match up with the data
        fibnorml = fltarr(1)
        for k=0,n_elements(fibnorm)-1 do $
          fibnorml = [fibnorml,replicate(fibnorm[k],nprofile)]
        fibnorml = fibnorml[1:*]
        
;---------------------------------------------------------------------------------
;The normalization for the wbiweight routine is made.
;---------------------------------------------------------------------------------
;I've elected here to give the entire fiber a weight, then let the
;profile weight add in quadrature to that. This constitutes the total
;weight that's fed into the weighted biweight (wbiweight.pro).

        alldatanorm = fltarr(n1,dt)
        flags = where(alldata eq -666,count)
        for k=0,dt-1 do alldatanorm[*,k] = alldata[*,k]/fibnorml[k]
        if (count ne 0) then begin
            alldatanorm[flags] = -666
            allwgtsn[flags] = 0.0
            allwgts[flags] = 0.0
        endif

; A plotting routine is done for each bin as a visual inspection. This
; generates two plots. The first is every row of data BEFORE
; NORMALIZATION BY THE FIBER WEIGHTS. The second are the normalized
; data. There is a switch that sends either the raw or normalized data
; into the wbiweight. (At least there will be at some point.)
        
        if (plotswitch eq 'on') then begin
            wave = fltarr(n1)
            yup = max(alldata)+25
            ydown = min(alldata)-25
            yupn = max(alldatanorm)+0.1
            ydownn = min(alldatanorm)-0.1
            for k=0,n1-1 do wave[k] = wi + wdisp*k
            
            repeat begin
                window,0,retain=2
                device,decomposed=0
                loadct,0
                plot,wave,alldata[*,0],xtitle='Wavelength',ytitle='Flux',$
                  yrange=[ydown,yup],xrange=[min(wave),max(wave)],ystyle=1,$
                  xstyle=1,title='Data for '+binname,charsize=1.5,/nodata
                loadct,27
                for k=0,dt-1 do oplot,wave,alldata[*,k],color=10+k
                
                window,1,retain=2
                device,decomposed=0
                loadct,0
                plot,wave,alldatanorm[*,0],xtitle='Wavelength',ytitle='Flux (normalized)',$
                  yrange=[ydownn,yupn],xrange=[min(wave),max(wave)],ystyle=1,$
                  xstyle=1,title='Normalized data for '+binname,charsize=1.5,/nodata
                loadct,27
                for k=0,dt-1 do oplot,wave,alldatanorm[*,k],color=10+k
                
                print,'Change the plotting range for the unnormalized data? ("y" or "n")'
                read,ans
                if (ans eq 'y') then begin
                    print,'Enter a new upper limit:'
                    read,yup
                    print,'Enter a new lower limit:'
                    read,ydown
                endif
            endrep until (ans eq 'n')
            
            wdelete,0,1
        endif
        
;***************************************************************************
;the allwgts array needs to be scaled up so that when added in
;quadrature with the fiber weights, they have some degree of
;influence. as the weight of the fiber is already taken into account,
;each fiber profile weight is scaled (done above, during the process
;of creating them) so as to be normalized to 1.
;***************************************************************************
        if (wgtyn eq 'wgty') then begin
            wgtarr = fltarr(n1,dt)
            for k=0,dt-1 do begin
                for l=0,n1-1 do wgtarr[l,k] = sqrt(allwgtsn[l,k]^2 + fibnorml[k]^2)
            endfor
;THESE WEIGHTS CAN NOW RANGE BETWEEN SQRT(2) AND ZERO. To fix this,
;they get another rescaling... how best to do this?
        
            wgtarr = wgtarr/max(wgtarr) ;the temporary solution...
        
            print,''
            print,'The max weight value is '+strn(max(wgtarr))
            print,'The min weight value is '+strn(min(wgtarr))
        
                                ;another visual check
        
            if(plotswitch eq 'on') then begin
                temp = wgtarr[bsort(wgtarr)]
                window,0,retain=2
                plot,temp,title='THE WEIGHTS, SORTED',charsize=1.2
                
                wdelete,0
            endif
        endif

;------------------------------------------------------------------------------
;THE dRA AND dDEC CALCULATION IS MADE
;------------------------------------------------------------------------------

; MODIFIED ON 04/11/09: The old method doesn't work as bins commonly
; combine fibers across the major axis. The new routine now simply
; returns the dRA and dDec values for each bin. 

; A calculation of the radius is going to be included in the near
; future...
;xxx
;This is still broken. if you want to include a luminosity weighted radius
;then you will need to get the weights from either the luminosity weighting from the 
;dotR files or (probably better) from the actual data in ALLDATA.

;USE THE FIBNORM ARRAY.


        if (p eq 0) then begin
            openw,5,binname+'.coord'
            printf,5,'F #','dRA','dDec',format='(a4,2x,11x,a3,10x,a4)'
            for k=0,n_elements(alloffsets[0,*])-1 $
              do printf,5,alloffsets[*,k],format='(i4,2x,f14.5,f14.5)'
            free_lun,5
        endif
        
;-----------------------------------------------------------------------------
;THE SIGNAL-TO-NOISE CALCULATION
        
        onelisttrim = strsplit(onelist,'.',/extract)
        outnamed = binname+'_'+onelisttrim[0]+'.fits'
        outnamesn = binname+'_'+onelisttrim[0]+'_s2n.fits'

        flags = where(alldata eq -666,count)
        if (count ne 0) then alldata[flags] = 0.0
        
        signal = dblarr(n1)
        noise = dblarr(n1)
        for s=0,n1-1 do signal[s] = total(alldata[s,*])
        for s=0,n1-1 do begin
            ind = where(allwgts[s,*] ne 0.0)
            noise[s] = sqrt(total(1/allwgts[s,ind]))
        endfor

        S2N = signal / noise
        writefits,outnamesn,S2N
        
;------------------------------------------------------------------------------
;THE WEIGHTED BIWEIGHT
;------------------------------------------------------------------------------
        if (wgtyn eq 'wgty') then outdata = wbiweightf(alldatanorm,wgtarr)
        if (wgtyn eq 'wgtn') then outdata = biweightf(alldatanorm)
        writefits,outnamed,outdata
;------------------------------------------------------------------------------
        jump1: ;this jump occurs when no data exists for a given bin
    endfor ;The end of the bin loop.
endfor ;The end of the list loop

print,''
print,'PIPE2 has finished successfully.'
print,''

end

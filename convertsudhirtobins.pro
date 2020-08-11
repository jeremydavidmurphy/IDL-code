PRO convertsudhirtobins, file, binend, angleend

;this rountine takes Sudhir's bin format and converts it into
;the dithers and fiber numbers my code expects. He gives me a list
;with a row showing the radial bins (in IDL index) followed by a row
;of the fiber numbers (with D2 and D3 denoted as fiber#+246 and
;fiber#+(246 * 2) respectively). The fiber numbers are also in IDL
;index! This codes just generates the output bin text files (to create
;your bin.list). 

;FILE: The name of Sudhir's text file. (ex: CoAddBinInfo.txt)
                                ;note: you've got to modify his
                                ;text files to cut the top line
                                ;(start things with the first bin) and
                                ;get rid of all the empty lines.
;BINEND: The last radial bin (in real values, not IDL indices)
;ANGLEEND: The last angular bin

free_lun,5
openr,5,file
finish = 'nope'

bin = 0
angle = 0
number = 0

repeat begin
   count = strarr(1)
   readf,5,count
   reads,count,format='(I0,I0,I0)',bin,angle,number
   bin = bin + 1 ;to adjust for idl index
   angle = angle + 1 ;to adjust for idl index
   print,bin,angle,number
   fibers = intarr(number)
   readf,5,fibers
   fibers = fibers + 1 ;to adjust for IDL index
   print,fibers
   openw,6,'bin'+strn(bin)+strn(angle)+'.txt'
   i1 = where(fibers le 246)
   if i1[0] ne -1 then for j=0,n_elements(i1)-1 do printf,6,strn(fibers[i1[j]])+'_D1'
   i2 = where(fibers gt 246 and fibers le 492)
   if i2[0] ne -1 then begin
      fibers[i2] = fibers[i2]-246 
      for j=0,n_elements(i2)-1 do printf,6,strn(fibers[i2[j]])+'_D2'
   endif
   i3 = where(fibers gt 492)
   if i3[0] ne -1 then begin
      fibers[i3] = fibers[i3]-492 
      for j=0,n_elements(i3)-1 do printf,6,strn(fibers[i3[j]])+'_D3'
   endif
   free_lun,6
   if bin eq binend and angle eq angleend then finish = 'yep'
endrep until (finish eq 'yep')

stop
END

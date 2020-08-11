pro Pmcfit2, list

readcol,list,silent=1,format='a',mcfit2s
n0 = n_elements(mcfit2s)

window,0,retain=2
device,decomposed=0
loadct,0
ans = ''

for l=0,n0-1 do begin
    file = mcfit2s[l]
    name = strsplit(file,'_.',/extract)
    binname = name[0]
    region = name[1]
    readcol,silent=1,file,format='f,f,f,f',vel,int,down,up

    plot,vel,int,title='LOSVD for '+binname+' Region: '+region,xrange=[-1500,1500],xstyle=1,$
      xtitle='Velocity (km/sec)',charsize=1.2
    oplot,vel,down,linestyle=2
    oplot,vel,up,linestyle=2
    print,'Plot as a .ps file? (y/n)'
    print,'Type "q" to quit...'
    read,ans
    if ans eq 'q' then goto,jumpend
    if ans eq 'y' then begin
       set_plot,'ps'
       device,file=binname+'_mcfitLOSVD.ps'
       plot,vel,int,title='LOSVD for '+binname+' Region: '+region,xrange=[-1500,1500],xstyle=1,$
            xtitle='Velocity (km/sec)',charsize=1.2
       oplot,vel,down,linestyle=2
       oplot,vel,up,linestyle=2
       device,/close_file
       set_plot,'x'
    endif

endfor

jumpend:
print,'Next ENTER closes the plotting window...'
pause
wdelete

stop
end

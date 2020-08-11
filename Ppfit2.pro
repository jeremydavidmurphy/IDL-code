PRO plotpfit

readcol,'bin.radius.N2832',f='a,f',bins,rad
readcol,'pfit.list',f='a',lists
n0 = n_elements(lists)

readcol,lists[0],f='x,f,f,x,f,f,x,x,x',v,d,h3,h4,silent=1
window,0,retain=2,xsize=1200
device,decomposed=0

colors = [10,50,90,130,180,220]

plot,rad,v,psym=1,/nodata,xrange=[-300,300],yrange=[-100,100]

for j=0,n0-1 do begin
    readcol,lists[j],f='x,f,f,x,f,f,x,x,x',v,d,h3,h4,silent=1
    oplot,rad,v,psym=1,color=colors[j]
    pause
endfor

stop
END

PRO sky_eval, a, b

a = 'vp'+strn(a)+'_p.fits'
b = 'vp'+strn(b)+'_p.fits'
aa = readfits(a,/silent)
bb = readfits(b,/silent)

print,a+'  '+strn(median(aa))
print,b+'  '+strn(median(bb))

END

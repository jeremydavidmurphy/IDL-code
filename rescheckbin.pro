PRO rescheckbin, binlist, datalist

;The difference between this code and rescheck is that this does not
;generate the plots for each night. Rather it reads in a list of bins
;and calculates the statistics on all the data going into that bin.

;Compile: rescheckF
readcol,binlist,format='a',list
n0 = n_elements(list)

readcol,datalist,format='x,a,a,x',pointing,dat
test = ['a','b','c','d','e','f']

for j=0,n0-1 do begin
    onelist = list[j]
    readcol,onelist,format='a',binlist
    n1 = n_elements(binlist)
    fibers = intarr(n1)
    letters = strarr(n1)
    for l=0,n1-1 do begin
        temp = strsplit(binlist[l],'_',/extract)
        fibers[l] = uint(temp[0])
        letters[l] = temp[1]
        print,fibers
        print,letters
    endfor
    i = 0
    for k=0,n_elements(letters)-1 do begin
        let = letters[k]
        letold = [letold,let]
        i = [i,where(pointing eq let)]
    endfor
    print,'Working on '+onelist+'...'
    print,n_elements(i)*n1
    onedat = dat[i]
    n2 = n_elements(onedat)
    Bkeep1 = fltarr(3,1)
    Gkeep1 = fltarr(3,1)
    Rkeep1 = fltarr(3,1)
    Bkeep2 = fltarr(3,1)
    Gkeep2 = fltarr(3,1)
    Rkeep2 = fltarr(3,1)
    for k=0,n2-1 do begin
        out = rescheckF(onedat[k])
        Bkeep1 = [[Bkeep1],[out[0:2,ifiber,0]]]
        Bkeep2 = [[Bkeep2],[out[0:2,ifiber,1]]]
        Gkeep1 = [[Gkeep1],[out[3:5,ifiber,0]]]
        Gkeep2 = [[Gkeep2],[out[3:5,ifiber,1]]]
        Rkeep1 = [[Rkeep1],[out[6:8,ifiber,0]]]
        Rkeep2 = [[Rkeep2],[out[6:8,ifiber,1]]]
    endfor
    Bkeep1 = Bkeep1[*,1:*]
    Bkeep2 = Bkeep2[*,1:*]
    Gkeep1 = Gkeep1[*,1:*]
    Gkeep2 = Gkeep2[*,1:*]
    Rkeep1 = Rkeep1[*,1:*]
    Rkeep2 = Rkeep2[*,1:*]

    out = fltarr(6,3)

    out[0,0] = median(Bkeep1)
    out[1,0] = mean(Bkeep1)
    out[2,0] = variance(Bkeep1)
    out[3,0] = median(Bkeep2)
    out[4,0] = mean(Bkeep2)
    out[5,0] = variance(Bkeep2)

    out[0,1] = median(Gkeep1)
    out[1,1] = mean(Gkeep1)
    out[2,1] = variance(Gkeep1)
    out[3,1] = median(Gkeep2)
    out[4,1] = mean(Gkeep2)
    out[5,1] = variance(Gkeep2)

    out[0,2] = median(Rkeep1)
    out[1,2] = mean(Rkeep1)
    out[2,2] = variance(Rkeep1)
    out[3,2] = median(Rkeep2)
    out[4,2] = mean(Rkeep2)
    out[5,2] = variance(Rkeep2)

    openw,5,onelist+'.RESstats'
    for m=0,2 do printf,5,out[*,m]
    free_lun,5
endfor

end

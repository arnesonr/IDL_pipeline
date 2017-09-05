pro bpm

;read in the dark that is already bias subtracted and in ADU/s
dark = readfits('~/clip_median_dark_20.fits',hdr)

n_cols = sxpar(hdr, 'NAXIS1')
n_rows = sxpar(hdr, 'NAXIS2')
;create the bpm
bpm = intarr(n_cols,n_rows)
bpm = bpm*0

djs_iterstat, dark, mean =mean, sigma =sigma
;loop through the dark and set anything > 5sigma to 1 (bad pixel)
for i=0L, n_cols-1, 1L do begin
    for j=0L, n_rows-1, 1L do begin
        if dark[i,j] gt mean + 5.*sigma then bpm[i,j]=1
    endfor
endfor
;exclude everything within 20 pixels of edge
bpm[0:19,*]=1
bpm[n_cols-21:n_cols-1,*]=1
bpm[*,0:19]=1
bpm[*,n_rows-21:n_rows-1]=1
;write the bpm to fits
fxaddpar,hdr,'BZERO',0.0,''
writefits, 'bpm.fits',bpm, hdr
end


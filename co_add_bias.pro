pro co_add_bias

readcol, 'bias.txt', filename, format = ('A')

filename = '/uufs/astro.utah.edu/common/astro_data/kdawson/ASTR5015/data/Sep2/'+filename
;read in all the biases
bias_0 = readfits(filename[0],hdr)
bias_1 = readfits(filename[1],hdr)
bias_2 = readfits(filename[2],hdr)
bias_3 = readfits(filename[3],hdr)
bias_4 = readfits(filename[4],hdr)
bias_5 = readfits(filename[5],hdr)
bias_6 = readfits(filename[6],hdr)
bias_7 = readfits(filename[7],hdr)
bias_8 = readfits(filename[8],hdr)
bias_9 = readfits(filename[9],hdr)
bias_10 = readfits(filename[10],hdr)
bias_11 = readfits(filename[11],hdr)
bias_12 = readfits(filename[12],hdr)
bias_13 = readfits(filename[13],hdr)
bias_14 = readfits(filename[14],hdr)
bias_15 = readfits(filename[15],hdr)
bias_16 = readfits(filename[16],hdr)
bias_17 = readfits(filename[17],hdr)
bias_18 = readfits(filename[18],hdr)
bias_19 = readfits(filename[19],hdr)

n_cols = sxpar(hdr, 'NAXIS1')
n_rows = sxpar(hdr, 'NAXIS2')

;=============================================================================
;MEAN

;create the stacked array to put values into
mean_bias_20 = fltarr(n_cols, n_rows)
;calculate the stacked mean
mean_bias_20 = (bias_0 + bias_1 + bias_2 + bias_3 + bias_4 + bias_5 + bias_6 + bias_7 + bias_8 + bias_9 + bias_10 + bias_11 + bias_12 + bias_13 + bias_14 + bias_15 + bias_16 + bias_17 + bias_18 + bias_19) / 20L

;create a stack of only 10 biases
mean_bias_10 = fltarr(n_cols, n_rows)
;calculate the stacked mean
mean_bias_10 = (bias_0 + bias_2 + bias_4 + bias_5 + bias_9 + bias_11 + bias_13 + bias_16 + bias_18 + bias_19) / 10L

;create a stack of only 4 biases
mean_bias_4 = fltarr(n_cols, n_rows)
;calculate the stacked mean
mean_bias_4 = (bias_2 + bias_5 + bias_11 + bias_18) / 4L

fxaddpar,hdr,'BZERO',0.0,''

;write fits files
writefits, 'mean_bias_20.fits', mean_bias_20, hdr
writefits, 'mean_bias_10.fits', mean_bias_10, hdr
writefits, 'mean_bias_4.fits', mean_bias_4, hdr

;=============================================================================
;MEDIAN

;create the median arrays
median_bias_20 = fltarr(n_cols, n_rows)
median_bias_10 = fltarr(n_cols, n_rows)
median_bias_4 = fltarr(n_cols, n_rows)
stacked_im = fltarr(n_cols, n_rows, 20)
;stacking 20 images
for i=0L, 19L, 1L do begin
    im = readfits(filename[i])
    stacked_im[*,*, i]=im
endfor

;iterate through each pixel and calcualte the median of the stacked im
for i=0L, n_cols-1, 1L do begin
    for j=0L, n_rows-1, 1L do begin
        median_bias_20[i,j] = median(stacked_im[i,j,*])
    endfor
endfor

stacked_im = fltarr(n_cols, n_rows, 10)
;stacking 10 images
for i=0L, 9L, 1L do begin
    im = readfits(filename[i])
    stacked_im[*,*, i]=im
endfor

;iterate through each pixel and calcualte the median of the stacked im
for i=0L, n_cols-1, 1L do begin
    for j=0L, n_rows-1, 1L do begin
        median_bias_10[i,j] = median(stacked_im[i,j,*])
    endfor
endfor

stacked_im = fltarr(n_cols, n_rows, 4)
;stacking 4 images
for i=0L, 3L, 1L do begin
    im = readfits(filename[i])
    stacked_im[*,*, i]=im
endfor

;iterate through each pixel and calcualte the median of the stacked im
for i=0L, n_cols-1, 1L do begin
    for j=0L, n_rows-1, 1L do begin
        median_bias_4[i,j] = median(stacked_im[i,j,*])
    endfor
endfor

fxaddpar,hdr,'BZERO',0.0,''

writefits, 'median_bias_20.fits', median_bias_20, hdr
writefits, 'median_bias_10.fits', median_bias_10, hdr
writefits, 'median_bias_4.fits', median_bias_4, hdr
;==============================================================================
;CLIPPED MEAN

;create the clipped median arrays
clip_mean_bias_20 = fltarr(n_cols, n_rows)
clip_mean_bias_10 = fltarr(n_cols, n_rows)
clip_mean_bias_4 = fltarr(n_cols, n_rows)
stacked_im = fltarr(n_cols, n_rows, 20)
;stacking 20 images
for i=0L, 19L, 1L do begin
    im = readfits(filename[i])
    stacked_im[*,*, i]=im
endfor

;iterate through each pixel and calcualte the median of the stacked im
for i=0L, n_cols-1, 1L do begin
    for j=0L, n_rows-1, 1L do begin
        djs_iterstat, stacked_im[i,j,*], mean=m
        clip_mean_bias_20[i,j] = m
    endfor
endfor

stacked_im = fltarr(n_cols, n_rows, 10)
;stacking 10 images
for i=0L, 9L, 1L do begin
    im = readfits(filename[i])
    stacked_im[*,*, i]=im
endfor

;iterate through each pixel and calcualte the median of the stacked im
for i=0L, n_cols-1, 1L do begin
    for j=0L, n_rows-1, 1L do begin
        djs_iterstat, stacked_im[i,j,*], mean=m
        clip_mean_bias_10[i,j] = m
    endfor
endfor

stacked_im = fltarr(n_cols, n_rows, 4)
;stacking 4 images
for i=0L, 3L, 1L do begin
    im = readfits(filename[i])
    stacked_im[*,*, i]=im
endfor

;iterate through each pixel and calcualte the median of the stacked im
for i=0L, n_cols-1, 1L do begin
    for j=0L, n_rows-1, 1L do begin
        djs_iterstat, stacked_im[i,j,*], mean=m
        clip_mean_bias_4[i,j] = m
    endfor
endfor

fxaddpar,hdr,'BZERO',0.0,''

writefits, 'clip_mean_bias_20.fits', clip_mean_bias_20, hdr
writefits, 'clip_mean_bias_10.fits', clip_mean_bias_10, hdr
writefits, 'clip_mean_bias_4.fits', clip_mean_bias_4, hdr

;==============================================================================
;CLIPPED MEDIAN
;create the clipped median arrays
clip_median_bias_20 = fltarr(n_cols, n_rows)
clip_median_bias_10 = fltarr(n_cols, n_rows)
clip_median_bias_4 = fltarr(n_cols, n_rows)
stacked_im = fltarr(n_cols, n_rows, 20)
;stacking 20 images
for i=0L, 19L, 1L do begin
    im = readfits(filename[i])
    stacked_im[*,*, i]=im
endfor


;iterate through each pixel and calcualte the median of the stacked im
for i=0L, n_cols-1, 1L do begin
    for j=0L, n_rows-1, 1L do begin
        djs_iterstat, stacked_im[i,j,*], median=m
        clip_median_bias_20[i,j] = m
    endfor
endfor

stacked_im = fltarr(n_cols, n_rows, 10)
;stacking 10 images
for i=0L, 9L, 1L do begin
    im = readfits(filename[i])
    stacked_im[*,*, i]=im
endfor

;iterate through each pixel and calcualte the median of the stacked im
for i=0L, n_cols-1, 1L do begin
    for j=0L, n_rows-1, 1L do begin
        djs_iterstat, stacked_im[i,j,*], median=m
        clip_median_bias_10[i,j] = m
    endfor
endfor

stacked_im = fltarr(n_cols, n_rows, 4)
;stacking 4 images
for i=0L, 3L, 1L do begin
    im = readfits(filename[i])
    stacked_im[*,*, i]=im
endfor

;iterate through each pixel and calcualte the median of the stacked im
for i=0L, n_cols-1, 1L do begin
    for j=0L, n_rows-1, 1L do begin
        djs_iterstat, stacked_im[i,j,*], median=m
        clip_median_bias_4[i,j] = m
    endfor
endfor

fxaddpar,hdr,'BZERO',0.0,''

writefits, 'clip_median_bias_20.fits', clip_median_bias_20, hdr
writefits, 'clip_median_bias_10.fits', clip_median_bias_10, hdr
writefits, 'clip_median_bias_4.fits', clip_median_bias_4, hdr
;==============================================================================
end


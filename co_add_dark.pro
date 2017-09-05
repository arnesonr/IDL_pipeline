pro co_add_dark

readcol, 'dark.txt', filename, format = ('A')

filename = '/uufs/astro.utah.edu/common/astro_data/kdawson/ASTR5015/data/55441/'+filename
mean_bias_20 = readfits('/uufs/astro.utah.edu/common/home/u0450159/mean_bias_20.fits')
;read in all the darks
dark_0 = readfits(filename[0],hdr)
exp_time = sxpar(hdr, 'EXPTIME')
dark_0 = (dark_0*1.-mean_bias_20)/exp_time
dark_1 = readfits(filename[1],hdr)
exp_time = sxpar(hdr, 'EXPTIME')
dark_1 = (dark_1*1.-mean_bias_20)/exp_time
dark_2 = readfits(filename[2],hdr)
exp_time = sxpar(hdr, 'EXPTIME')
dark_2 = (dark_2*1.-mean_bias_20)/exp_time
dark_3 = readfits(filename[3],hdr)
exp_time = sxpar(hdr, 'EXPTIME')
dark_3 = (dark_3*1.-mean_bias_20)/exp_time
dark_4 = readfits(filename[4],hdr)
exp_time = sxpar(hdr, 'EXPTIME')
dark_4 = (dark_4*1.-mean_bias_20)/exp_time
dark_5 = readfits(filename[5],hdr)
exp_time = sxpar(hdr, 'EXPTIME')
dark_5 = (dark_5*1.-mean_bias_20)/exp_time
dark_6 = readfits(filename[6],hdr)
exp_time = sxpar(hdr, 'EXPTIME')
dark_6 = (dark_6*1.-mean_bias_20)/exp_time
dark_7 = readfits(filename[7],hdr)
exp_time = sxpar(hdr, 'EXPTIME')
dark_7 = (dark_7*1.-mean_bias_20)/exp_time
dark_8 = readfits(filename[8],hdr)
exp_time = sxpar(hdr, 'EXPTIME')
dark_8 = (dark_8*1.-mean_bias_20)/exp_time
dark_9 = readfits(filename[9],hdr)
exp_time = sxpar(hdr, 'EXPTIME')
dark_9 = (dark_9*1.-mean_bias_20)/exp_time
dark_10 = readfits(filename[10],hdr)
exp_time = sxpar(hdr, 'EXPTIME')
dark_10 = (dark_10*1.-mean_bias_20)/exp_time
dark_11 = readfits(filename[11],hdr)
exp_time = sxpar(hdr, 'EXPTIME')
dark_11 = (dark_11*1.-mean_bias_20)/exp_time
dark_12 = readfits(filename[12],hdr)
exp_time = sxpar(hdr, 'EXPTIME')
dark_12 = (dark_12*1.-mean_bias_20)/exp_time
dark_13 = readfits(filename[13],hdr)
exp_time = sxpar(hdr, 'EXPTIME')
dark_13 = (dark_13*1.-mean_bias_20)/exp_time
dark_14 = readfits(filename[14],hdr)
exp_time = sxpar(hdr, 'EXPTIME')
dark_14 = (dark_14*1.-mean_bias_20)/exp_time
dark_15 = readfits(filename[15],hdr)
exp_time = sxpar(hdr, 'EXPTIME')
dark_15 = (dark_15*1.-mean_bias_20)/exp_time
dark_16 = readfits(filename[16],hdr)
exp_time = sxpar(hdr, 'EXPTIME')
dark_16 = (dark_16*1.-mean_bias_20)/exp_time
dark_17 = readfits(filename[17],hdr)
exp_time = sxpar(hdr, 'EXPTIME')
dark_17 = (dark_17*1.-mean_bias_20)/exp_time
dark_18 = readfits(filename[18],hdr)
exp_time = sxpar(hdr, 'EXPTIME')
dark_18 = (dark_18*1.-mean_bias_20)/exp_time
dark_19 = readfits(filename[19],hdr)
exp_time = sxpar(hdr, 'EXPTIME')
dark_19 = (dark_19*1.-mean_bias_20)/exp_time

n_cols = sxpar(hdr, 'NAXIS1')
n_rows = sxpar(hdr, 'NAXIS2')
exp_time = sxpar(hdr, 'EXPTIME')
;=============================================================================
;MEAN

;create the stacked array to put values into
mean_dark_20 = intarr(n_cols, n_rows)
;calculate the stacked mean
mean_dark_20 = (dark_0 + dark_1 + dark_2 + dark_3 + dark_4 + dark_5 + dark_6 + dark_7 + dark_8 + dark_9 + dark_10 + dark_11 + dark_12 + dark_13 + dark_14 + dark_15 + dark_16 + dark_17 + dark_18 + dark_19) / 20L

;create a stack of only 10 darks
mean_dark_10 = intarr(n_cols, n_rows)
;calculate the stacked mean
mean_dark_10 = (dark_0 + dark_2 + dark_4 + dark_5 + dark_9 + dark_11 + dark_13 + dark_16 + dark_18 + dark_19) / 10L

;create a stack of only 4 darks
mean_dark_4 = intarr(n_cols, n_rows)
;calculate the stacked mean
mean_dark_4 = (dark_2 +dark_5 + dark_11 + dark_18) / 4L

fxaddpar,hdr,'BZERO',0.0,''
;write fits files
writefits, 'mean_dark_20.fits', mean_dark_20, hdr
writefits, 'mean_dark_10.fits', mean_dark_10, hdr
writefits, 'mean_dark_4.fits', mean_dark_4, hdr

;=============================================================================
;MEDIAN

;create the median arrays
median_dark_20 = intarr(n_cols, n_rows)
median_dark_10 = intarr(n_cols, n_rows)
median_dark_4 = intarr(n_cols, n_rows)
stacked_im = fltarr(n_cols, n_rows, 20)
median_bias_20 = readfits('/uufs/astro.utah.edu/common/home/u0450159/median_bias_20.fits')
;stacking 20 images
for i=0L, 19L, 1L do begin
    im = readfits(filename[i])
    im = im*1.
    stacked_im[*,*, i]=im
endfor

;iterate through each pixel and calcualte the median of the stacked im
for i=0L, n_cols-1, 1L do begin
    for j=0L, n_rows-1, 1L do begin
        median_dark_20[i,j] = median(stacked_im[i,j,*])
    endfor
endfor
stacked_im = fltarr(n_cols, n_rows, 10)
median_bias_10 = readfits('/uufs/astro.utah.edu/common/home/u0450159/median_bias_10.fits')
;stacking 10 images
for i=0L, 9L, 1L do begin
    im = readfits(filename[i])
    im = im*1.
    stacked_im[*,*, i]=im
endfor

;iterate through each pixel and calcualte the median of the stacked im
for i=0L, n_cols-1, 1L do begin
    for j=0L, n_rows-1, 1L do begin
        median_dark_10[i,j] = median(stacked_im[i,j,*])
    endfor
endfor

stacked_im = fltarr(n_cols, n_rows, 4)
median_bias_4 = readfits('/uufs/astro.utah.edu/common/home/u0450159/median_bias_4.fits')
;stacking 4 images
for i=0L, 3L, 1L do begin
    im = readfits(filename[i])
    im = im*1.
    stacked_im[*,*, i]=im
endfor

;iterate through each pixel and calcualte the median of the stacked im
for i=0L, n_cols-1, 1L do begin
    for j=0L, n_rows-1, 1L do begin
        median_dark_4[i,j] = median(stacked_im[i,j,*])
    endfor
endfor
median_dark_20 = (median_dark_20-median_bias_20)/exp_time
median_dark_10 = (median_dark_10-median_bias_10)/exp_time
median_dark_4 = (median_dark_4-median_bias_4)/exp_time

fxaddpar,hdr,'BZERO',0.0,''
writefits, 'median_dark_20.fits', median_dark_20, hdr
writefits, 'median_dark_10.fits', median_dark_10, hdr
writefits, 'median_dark_4.fits', median_dark_4, hdr
;==============================================================================
;CLIPPED MEAN

;create the clipped median arrays
clip_mean_dark_20 = intarr(n_cols, n_rows)
clip_mean_dark_10 = intarr(n_cols, n_rows)
clip_mean_dark_4 = intarr(n_cols, n_rows)
stacked_im = fltarr(n_cols, n_rows, 20)
clip_mean_bias_20 = readfits('/uufs/astro.utah.edu/common/home/u0450159/clip_mean_bias_20.fits')
;stacking 20 images
for i=0L, 19L, 1L do begin
    im = readfits(filename[i])
    im = im*1.
    stacked_im[*,*, i]=im
endfor

;iterate through each pixel and calcualte the median of the stacked im
for i=0L, n_cols-1, 1L do begin
    for j=0L, n_rows-1, 1L do begin
        djs_iterstat, stacked_im[i,j,*], mean=m
        clip_mean_dark_20[i,j] = m
    endfor
endfor

stacked_im = fltarr(n_cols, n_rows, 10)
clip_mean_bias_10 = readfits('/uufs/astro.utah.edu/common/home/u0450159/clip_mean_bias_10.fits')
;stacking 10 images
for i=0L, 9L, 1L do begin
    im = readfits(filename[i])
    im = im*1.
    stacked_im[*,*, i]=im
endfor

;iterate through each pixel and calcualte the median of the stacked im
for i=0L, n_cols-1, 1L do begin
    for j=0L, n_rows-1, 1L do begin
        djs_iterstat, stacked_im[i,j,*], mean=m
        clip_mean_dark_10[i,j] = m
    endfor
endfor

stacked_im = fltarr(n_cols, n_rows, 4)
clip_mean_bias_4 = readfits('/uufs/astro.utah.edu/common/home/u0450159/clip_mean_bias_4.fits')
;stacking 4 images
for i=0L, 3L, 1L do begin
    im = readfits(filename[i])
    im = im*1.
    stacked_im[*,*, i]=im
endfor

;iterate through each pixel and calcualte the median of the stacked im
for i=0L, n_cols-1, 1L do begin
    for j=0L, n_rows-1, 1L do begin
        djs_iterstat, stacked_im[i,j,*], mean=m
        clip_mean_dark_4[i,j] = m
    endfor
endfor
clip_mean_dark_20 = (clip_mean_dark_20-clip_mean_bias_20)/exp_time
clip_mean_dark_10 = (clip_mean_dark_10-clip_mean_bias_10)/exp_time
clip_mean_dark_4 = (clip_mean_dark_4-clip_mean_bias_4)/exp_time

fxaddpar,hdr,'BZERO',0.0,''
writefits, 'clip_mean_dark_20.fits', clip_mean_dark_20, hdr
writefits, 'clip_mean_dark_10.fits', clip_mean_dark_10, hdr
writefits, 'clip_mean_dark_4.fits', clip_mean_dark_4, hdr

;==============================================================================
;CLIPPED MEDIAN
;create the clipped median arrays
clip_median_dark_20 = intarr(n_cols, n_rows)
clip_median_dark_10 = intarr(n_cols, n_rows)
clip_median_dark_4 = intarr(n_cols, n_rows)
stacked_im = fltarr(n_cols, n_rows, 20)
clip_median_bias_20 = readfits('/uufs/astro.utah.edu/common/home/u0450159/clip_median_bias_20.fits')
;stacking 20 images
for i=0L, 19L, 1L do begin
    im = readfits(filename[i])
    im = im*1.
    stacked_im[*,*, i]=im
endfor

;iterate through each pixel and calcualte the median of the stacked im
for i=0L, n_cols-1, 1L do begin
    for j=0L, n_rows-1, 1L do begin
        djs_iterstat, stacked_im[i,j,*], median=m
        clip_median_dark_20[i,j] = m
    endfor
endfor

stacked_im = fltarr(n_cols, n_rows, 10)
clip_median_bias_10 = readfits('/uufs/astro.utah.edu/common/home/u0450159/clip_median_bias_10.fits')
;stacking 10 images
for i=0L, 9L, 1L do begin
    im = readfits(filename[i])
    im = im*1.
    stacked_im[*,*, i]=im
endfor

;iterate through each pixel and calcualte the median of the stacked im
for i=0L, n_cols-1, 1L do begin
    for j=0L, n_rows-1, 1L do begin
        djs_iterstat, stacked_im[i,j,*], median=m
        clip_median_dark_10[i,j] = m
    endfor
endfor

stacked_im = fltarr(n_cols, n_rows, 4)
clip_median_bias_4 = readfits('/uufs/astro.utah.edu/common/home/u0450159/clip_median_bias_4.fits')
;stacking 4 images
for i=0L, 3L, 1L do begin
    im = readfits(filename[i])
    im = im*1.
    stacked_im[*,*, i]=im
endfor

;iterate through each pixel and calcualte the median of the stacked im
for i=0L, n_cols-1, 1L do begin
    for j=0L, n_rows-1, 1L do begin
        djs_iterstat, stacked_im[i,j,*], median=m
        clip_median_dark_4[i,j] = m
    endfor
endfor

clip_median_dark_20 = (clip_median_dark_20-clip_median_bias_20)/exp_time
clip_median_dark_10 = (clip_median_dark_10-clip_median_bias_10)/exp_time
clip_median_dark_4 = (clip_median_dark_4-clip_median_bias_4)/exp_time

fxaddpar,hdr,'BZERO',0.0,''
writefits, 'clip_median_dark_20.fits', clip_median_dark_20, hdr
writefits, 'clip_median_dark_10.fits', clip_median_dark_10, hdr
writefits, 'clip_median_dark_4.fits', clip_median_dark_4, hdr
;==============================================================================
end




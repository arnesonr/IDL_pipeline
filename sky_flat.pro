pro sky_flat

readcol, 'SBIG_g.txt', filename, format = ('A')

filename = '/uufs/astro.utah.edu/common/astro_data/kdawson/ASTR5015/data/55441/'+filename
how_many = n_elements(filename)
bpm = readfits('~/bpm.fits',hdr)*1.
n_cols = sxpar(hdr, 'NAXIS1')
n_rows = sxpar(hdr, 'NAXIS2')
bias = readfits('~/clip_median_bias_20.fits')
dark = readfits('~/clip_median_dark_20.fits')
fxaddpar,hdr,'BZERO',0.0,''

;stack the SBIGs in a 3D array
stacked_im = fltarr(n_cols, n_rows, how_many)

for i=0L, how_many-1, 1L do begin
    im = readfits(filename[i],hdr)
    exp_time = sxpar(hdr, 'EXPTIME')
    im = (im * 1. - bias)/exp_time - dark
    stacked_im[*,*,i]=im
endfor
;create 3D array of object masks (good pixels = 0)
obj_mask = intarr(n_cols,n_rows,how_many)*0
;use bpm to calculate mean of each SBIG and flag any +/- 5 sigma
for i=0L, how_many-1L, 1L do begin
    trigger = 1.
    while trigger gt 0.01 do begin
        obj = stacked_im[*,*,i]
        a = where(obj_mask[*,*,i] eq 0 AND bpm eq 0)
        good_pixels = obj[a]
        djs_iterstat, good_pixels, sigma = sigma, mean=mean
        ;flag the bad pixels as 1 and grow it by 5 pixels
        for j=0L, n_cols-1, 1L do begin
            for k=0L, n_rows-1, 1L do begin
                if bpm[j,k] lt 1 then begin
                    if obj[j,k] gt (mean+5.*sigma) then obj_mask[j-5:j+5,k-5:k+5,i]=1
                    ;if obj[j,k] lt (mean-5.*sigma) then obj_mask[j-5:j+5,k-5:k+5,i]=1
                    vector = [obj[j,k], obj[j,k-1], obj[j,k+1], obj[j-1,k], obj[j+1,k]]
                    if min(vector) gt (mean+2.*sigma) then obj_mask[j-5:j+5,k-5:k+5,i]=1
                endif
            endfor
        endfor
        b = where(obj_mask[*,*,i] eq 0 AND bpm eq 0)
        good_pixels_2 = obj[b]
        djs_iterstat, good_pixels_2, sigma = sigma_2
        trigger = (sigma-sigma_2)/sigma
    endwhile
endfor
;put bpm through all of obj_mask
for i=0L, how_many-1L, 1L do begin
    for j=0L, n_cols-1, 1L do begin
        for k=0L, n_rows-1, 1L do begin
            if bpm[j,k] eq 1 AND obj_mask[j,k,i] eq 0 then obj_mask[j,k,i]=1
        endfor
    endfor
endfor

;write the masks to fits files to inspect
fxaddpar,hdr,'BZERO',0.0,''
writefits, 'g_band_mask.fits', obj_mask, hdr

;compute the denominator of M (weighted mean)
denom = fltarr(how_many)
for i =0L, how_many-1, 1L do begin
    c = where(obj_mask[*,*,i] eq 0)
    obj = stacked_im[*,*,i]
    array = obj[c]
    djs_iterstat, array, sigma=sig, mean=m
    stacked_im[*,*,i]=obj/m
    ;obj_mask[*,*,i]=obj_mask[*,*,i]/m
    ;sig=sig/m
    denom[i] = 1./sig^2
endfor
denominator = total(denom)
flat = fltarr(n_cols,n_rows)
;compute the numerator of M (weighted mean)
for j=0L, n_cols-1, 1L do begin
    for k=0L, n_rows-1, 1L do begin
        vector=stacked_im[j,k,*]
        a=where(obj_mask[j,k,*] eq 0)
        if min(a) ge 0 then begin
            vector_good = vector(a)
            djs_iterstat, vector_good, median=median
            flat[j,k]=median;/total(denom[a])
            ;flat[j,k]=total(vector[a]*denom[a])/total(denom[a])
        endif
; set the global bad pixels around the border equal to 0.0
        if min(a) lt 0 then flat[j,k]= 0.0
    endfor
endfor
;divide the flat by its mean
a = where(flat gt 0)
djs_iterstat, flat[a], mean=M, sigma=sigma
flat = flat/M

print, M, sigma

;write the super_flat
fxaddpar,hdr,'BZERO',0.0,''
writefits,'G_super_flat.fits',flat, hdr

end



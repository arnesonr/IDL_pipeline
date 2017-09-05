
;SDSS1:
spAll.fits
 
;BOSS:
spAll-v5_4_14.fits
 
;Arguments to readspec:
 
;SDSS1 data:
sdss1_topdir='/uufs/chpc.utah.edu/common/home/bolton_data2/SDSS1/spectro'
sdss1_run1d=''
sdss1_run2d=''
 
;BOSS data:
boss_topdir='/uufs/chpc.utah.edu/common/home/bolton_data2/SDSS3/BOSS/bossredux_lbl'
boss_run1d='v5_4_14'
boss_run2d='v5_4_14'

;select only the LRG sample
sdss1_spall = sdss1_topdir + "/spAll.fits"
columns = ['PLATE','FIBERID','MJD','OBJTYPE','PLUG_RA','PLUG_DEC', $
      'PROGNAME','SPECPRIMARY','PRIMTARGET','SECTARGET','CLASS','Z','ZWARNING']
spall = hogg_mrdfits(sdss1_spall,1,columns=columns)
ilrg = where(spall.specprimary EQ 1 $
      AND strmatch(spall.progname,'main*') $
      AND (spall.primtarget AND 2L^5 + 2L^26) NE 0 $
      AND strmatch(spall.class,'GALAXY*') $
      AND spall.zwarning EQ 0, nlrg)
;SDSS1 has 173,748 LRGs
;
;
;keep only 1,000  spectra randomly selected
nkeep=1000l
srt=sort(randomu(rseed,nlrg))
ilrg=ilrg[srt[0:nkeep-1]]
readspec, spall[ilrg].plate,spall[ilrg].fiberid,mjd=spall[ilrg].mjd,$
flux=flux, wave=wave, invvar=invvar, zans=zans, synflux=synflux, /align,$
topdir=sdss1_topdir, run1d=sdss1_run1d, run2d=sdss1_run2d
;plotting a spectra, template solution, and residual spectra respectively
splot, wave, flux[*,500]
soplot, wave, synflux[*,500], color=2
soplot, wave, (flux-synflux)[*,500], color=3
;define the residual
rflux = flux - synflux
;form a guassian kernel
pixbase = findgen(15)
pixbase = pixbase - mean(pixbase)
sigma = 2.4
kernel = exp(-0.5 * pixbase^2 / sigma^2)
kernel = kernel / total(kernel)
;plotting the gaussian kernel
splot, pixbase, kernel
splot, pixbase, kernel, ps=10
;define the residual spectra
rflux_test = rflux[*,500]
;define the invvar
invvar_test = invvar[*,500]
;convolve the invvar and kernel
conv = convolve(invvar_test, kernel^2)
;trick to keep only noise greater than 0
lnoise = (conv gt 0.) / sqrt(conv + (conv le 0.))
;plot the noise
splot, wave, lnoise
;need to find a solution based on wavelength not pixel
;pixel to wavelength related by:
print, alog10(wave[300:350])

;convolve the kernel with all of the invvars
conv = convolve(invvar, kernel^2)
;find the noise of all the spectra
lnoise = (conv gt 0.) / sqrt(conv + (conv le 0.))
;stack all the lnoise spectra
dims = size(lnoise,/dimensions)
ave_lnoise = fltarr(nkeep)
for i=0L, nkeep-1, 1L do ave_lnoise[i] = mean(lnoise[*,i])

srt = sort(ave_lnoise)
srt_lnoise = lnoise[*,srt]
;plot the 50th, 84th, and 16th percentile average lnoise spectra
splot, wave, srt_lnoise[*,499]
soplot, wave, srt_lnoise[*,839],color=2
soplot, wave, srt_lnoise[*,159],color=3

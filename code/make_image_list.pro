pro make_image_list

; make list of images and filenames for downloading cutouts
; this code has been expanded to make a catalog including more
; variables for use in analysis

acs=mrdfits("~/data/cosmos/code/lensing18_20110914.fits",1)
group=mrdfits("~/data/cosmos/code/group5_20110914.fits",1) ; for halo masses
ext=".list"
listFile="../images/acs_mem"
openw,u,listFile+ext,/get_lun,width=1000

; select only members of interest
sel=where(acs.kevin_mstar GT 9.4 $
          AND acs.p_mem_best GT 0.5 $
         )

acs=acs[sel]

; get other properties from morph catalog
;morph=mrdfits("~/data/cosmos/catalogs/griffith/cosmos_i_public_catalog_V1.0.fits.gz",1)
morph=mrdfits("~/data/cosmos/catalogs/griffith/ACS-GC_published_catalogs/cosmos_i_public_catalog_V1.0.fits.gz",1)
close_match_radec,acs.alpha_j2000,acs.delta_j2000,morph.ra,morph.dec,m1,m2,1./3600,1,miss1
acsmiss=acs[miss1]
acs=acs[m1]
morph=morph[m2]

; get group masses and orientation angles for member galaxies
halomass=fltarr(n_elements(acs))
theta=fltarr(n_elements(acs))
for gg=0,n_elements(group)-1 do begin
   mem=where(acs.group_id_best EQ group[gg].id,nmem)
   if(nmem GT 0) then begin
      halomass[mem]=group[gg].lensing_m200
      
      ; get angle between satellite and central (to rotate images in plot)
      ; define theta: angle between central, satellite, and east (left)
      ; rotate the (N up) image clockwise by theta in order the have the
      ; central-facing side of the satellite to point left on a plot
      dx=group[gg].alpha_mmgg_scale - acs[mem].alpha_j2000
      dy=-group[gg].delta_mmgg_scale + acs[mem].delta_j2000
      theta[mem]=atan(dy,dx)*!radeg ; theta is CW rotation angle in degrees
   endif
endfor

; this is where the image cutout filename is recorded
; (following IRSA's naming convention)
printf,u,"|          id |   flag |             ra |            dec |           z |          sm |       mhalo |           r |  ztype | zbulge | nomatch | rgflag |      rgsize |    rgsersic |        rgba |        rgpa |       color |        ssfr |         ebv |       theta |                                              filename |"

; convert size from pixels to kpc
pixScale=0.05 ; arcsec/pixel

for ii=0,n_elements(acs)-1 do begin
   printf,u,$
          acs[ii].ident,$
          acs[ii].group_flag_best,$
          acs[ii].alpha_j2000,$
          acs[ii].delta_j2000,$
          acs[ii].photoz_non_comb,$
          acs[ii].kevin_mstar,$
          halomass[ii],$
          acs[ii].dist_bcg_r200,$
          acs[ii].zest_type,$
          acs[ii].zest_bulge,$
          0, $ ; nomatch flag
          morph[ii].flag_galfit_hi,$
          as2kpc(morph[ii].re_galfit_hi * pixScale, acs[ii].photoz_non_comb),$ ; kpc
          morph[ii].n_galfit_hi,$
          morph[ii].ba_galfit_hi,$
          morph[ii].pa_galfit_hi,$
          acs[ii].mnuv_mr,$
          acs[ii].ssfr_med,$
          acs[ii].ebv,$
          theta[ii],$
          "   "+string(acs[ii].alpha_j2000,format='(F9.5)')+"_"+string(acs[ii].delta_j2000,format='(F9.7)')+"_acs_I_mosaic_30mas_sci.fits"
endfor

; add in objects for which there is no match in morph catalog
; first get halo masses and angles for missed matches
halomassmiss=fltarr(n_elements(acsmiss))
thetamiss=fltarr(n_elements(acsmiss))
for gg=0,n_elements(group)-1 do begin
   mem=where(acsmiss.group_id_best EQ group[gg].id,nmem)
   if(nmem GT 0) then begin
      halomassmiss[mem]=group[gg].lensing_m200
      
      ; get angle between satellite and central (to rotate images in plot)
      ; define theta: angle between central, satellite, and east (left)
      ; rotate the (N up) image clockwise by theta in order the have the
      ; central-facing side of the satellite to point left on a plot
      dx=group[gg].alpha_mmgg_scale - acsmiss[mem].alpha_j2000
      dy=-group[gg].delta_mmgg_scale + acsmiss[mem].delta_j2000
      thetamiss[mem]=atan(dy,dx)*!radeg ; theta is CW rotation angle in degrees
   endif
endfor

for ii=0,n_elements(acsmiss)-1 do begin
   printf,u,$
          acsmiss[ii].ident,$
          acsmiss[ii].group_flag_best,$
          acsmiss[ii].alpha_j2000,$
          acsmiss[ii].delta_j2000,$
          acsmiss[ii].photoz_non_comb,$
          acsmiss[ii].kevin_mstar,$
          halomassmiss[ii],$
          acsmiss[ii].dist_bcg_r200,$
          acsmiss[ii].zest_type,$
          acsmiss[ii].zest_bulge,$
          1, $ ; nomatch flag
          -99, $ ; flag_galfit_hi
          -99., $ ; re_galfit_hi
          -99., $ ; n_galfit_hi
          -99., $ ; ba_galfit_hi
          -99., $ ; pa_galfit_hi
          acsmiss[ii].mnuv_mr,$
          acsmiss[ii].ssfr_med,$
          acsmiss[ii].ebv,$
          thetamiss[ii],$
          "   "+string(acsmiss[ii].alpha_j2000,format='(F9.5)')+"_"+string(acsmiss[ii].delta_j2000,format='(F9.7)')+"_acs_I_mosaic_30mas_sci.fits"
endfor

close,u
free_lun,u

; split into multiple lists if needed (IRSA can only take 1000 at a time)
if(n_elements(acs) GT 1000) then begin
   nLists=ceil(float(n_elements(acs))/1000)
   for ii=0,nLists-1 do begin
      spawn,'head -1 '+listFile+ext+" > "+listFile+"_"+string(ii,format='(I0)')+ext
      if(ii LT nLists-1) then $
         spawn,'head -'+string(1+(ii+1)*1000,format='(I0)')+" "+listFile+ext+" | tail -1000 >> "+listFile+"_"+string(ii,format='(I0)')+ext $
      else spawn,'head -'+string(1+(ii+1)*1000,format='(I0)')+" "+listFile+ext+" | tail -"+string(n_elements(acs) MOD 1000,format='(I0)')+" >> "+listFile+"_"+string(ii,format='(I0)')+ext

   endfor
   print,nLists,' sub-lists printed'
endif

stop

end

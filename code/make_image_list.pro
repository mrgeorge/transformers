pro make_image_list

; make list of images and filenames for downloading cutouts

acs=mrdfits("../../code/lensing18_20110914.fits",1)
listFile="../images/acs_mem.list"
openw,u,listFile,/get_lun,width=1000

sel=where(acs.kevin_mstar GT 9.8 $
          AND acs.p_mem_best GT 0.5 $
          AND acs.photoz_non_comb GT 0.2 $
          AND acs.photoz_non_comb LT 0.5)

; this is where the image cutout filename is recorded
; (following IRSA's naming convention)
printf,u,"|          id |             ra |            dec |           z |          sm |           r |   type |  bulge |       color |                                              filename |"
for ii=0,n_elements(sel)-1 do $
   printf,u,$
          acs[sel[ii]].ident,$
          acs[sel[ii]].alpha_j2000,$
          acs[sel[ii]].delta_j2000,$
          acs[sel[ii]].photoz_non_comb,$
          acs[sel[ii]].kevin_mstar,$
          acs[sel[ii]].dist_bcg_r200,$
          acs[sel[ii]].zest_type,$
          acs[sel[ii]].zest_bulge,$
          acs[sel[ii]].mnuv_mr,$
          "   "+string(ii+1,format='(I04)')+"_"+string(acs[sel[ii]].alpha_j2000,format='(F9.5)')+"_"+string(acs[sel[ii]].delta_j2000,format='(F9.7)')+"_acs_I_mosaic_30mas_sci.fits"

close,u
free_lun,u

end

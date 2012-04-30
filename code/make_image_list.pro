pro make_image_list

; make list of images and filenames for downloading cutouts

acs=mrdfits("../../code/lensing18_20110914.fits",1)
group=mrdfits("~/data/cosmos/code/group5_20110914.fits",1) ; for halo masses
ext=".list"
listFile="../images/acs_mem"
openw,u,listFile+ext,/get_lun,width=1000

sel=where(acs.kevin_mstar GT 9.4 $
          AND acs.p_mem_best GT 0.5 $
          ,nSel)

theta=fltarr(nSel)

; get group masses for member galaxies
halomass=fltarr(n_elements(acs))
for gg=0,n_elements(group)-1 do begin
   mem=where(acs.group_id_best EQ group[gg].id,nmem)
   if(nmem GT 0) then halomass[mem]=group[gg].lensing_m200
endfor


; this is where the image cutout filename is recorded
; (following IRSA's naming convention)
printf,u,"|          id |   flag |             ra |            dec |           z |          sm |       mhalo |           r |   type |  bulge |       color |       theta |                                              filename |"
for ii=0,nSel-1 do begin

   ; get angle between satellite and central (to rotate images in plot)
   ; define theta: angle between central, satellite, and east (left)
   ; rotate the (N up) image clockwise by theta in order the have the
   ; central-facing side of the satellite to point left on a plot
   if(acs[sel[ii]].mmgg_scale EQ 1) then theta[ii]=0. $
   else begin
      cen=where(acs.mmgg_scale EQ 1 $
                AND acs.group_id_best EQ acs[sel[ii]].group_id_best $
                ,nCen)
      if(nCen NE 1) then begin
         print,'Error finding central'
         stop
      endif
      dy=acs[sel[ii]].delta_j2000-acs[cen].delta_j2000 ; north is up
      dx=acs[cen].alpha_j2000-acs[sel[ii]].alpha_j2000 ; east (higher RA) is left
      theta[ii]=atan(dy,dx)*!radeg ; theta is CW rotation angle in degrees
   endelse


   printf,u,$
          acs[sel[ii]].ident,$
          acs[sel[ii]].group_flag_best,$
          acs[sel[ii]].alpha_j2000,$
          acs[sel[ii]].delta_j2000,$
          acs[sel[ii]].photoz_non_comb,$
          acs[sel[ii]].kevin_mstar,$
          halomass[sel[ii]],$
          acs[sel[ii]].dist_bcg_r200,$
          acs[sel[ii]].zest_type,$
          acs[sel[ii]].zest_bulge,$
          acs[sel[ii]].mnuv_mr,$
          theta[ii],$
          "   "+string(acs[sel[ii]].alpha_j2000,format='(F9.5)')+"_"+string(acs[sel[ii]].delta_j2000,format='(F9.7)')+"_acs_I_mosaic_30mas_sci.fits"

endfor

close,u
free_lun,u

; split into multiple lists if needed (IRSA can only take 1000 at a time)
if(nSel GT 1000) then begin
   nLists=ceil(float(nSel)/1000)
   for ii=0,nLists-1 do begin
      spawn,'head -1 '+listFile+ext+" > "+listFile+"_"+string(ii,format='(I0)')+ext
      if(ii LT nLists-1) then $
         spawn,'head -'+string(1+(ii+1)*1000,format='(I0)')+" "+listFile+ext+" | tail -1000 >> "+listFile+"_"+string(ii,format='(I0)')+ext $
      else spawn,'head -'+string(1+(ii+1)*1000,format='(I0)')+" "+listFile+ext+" | tail -"+string(nSel MOD 1000,format='(I0)')+" >> "+listFile+"_"+string(ii,format='(I0)')+ext

   endfor
   print,nLists,' sub-lists printed'
endif

end

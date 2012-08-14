pro print_result,fileName,magval,distval,nPhotMemMockMem,nPhotMem,nMockMem,nPhotMem_w,nMockMem_w
; print values used for individual plot to a file
fileExists=file_test(fileName)
if(NOT fileExists) then begin
   openw,u,fileName,/get_lun,width=1000
   printf,u,'# magvals distvals nPhotMemMockMem   nPhotMem   nMockMem   nPhotMem_w   nMockMem_w'
endif else openu,u,fileName,/get_lun,/append,width=1000
printf,u,magval,distval,nPhotMemMockMem,nPhotMem,nMockMem,nPhotMem_w,nMockMem_w
close,u
free_lun,u
end

pro memstat_separability

; The purity of group membership assignment depends on a number of
; variables, most importantly magnitude (for photoz quality) and
; distance from the group center (the group population is falling, so
; of the ones selected more are field galaxies in the outskirts).
; This code should test the separability of these two factors; 
; does P(m,R) = P(m) * P(R)?

; Based on code from ~/data/cosmos/code/mocks/memstat_mock.pro and
; ~/data/cosmos/code/analysis/plot_memstat_all.pro for a more
; specialized analysis

mockDir="../../mocks/"
nMocks=10
suffix='' ; could change for SZ, Chandra, etc.
acsFileArr=mockDir+"acs_mock4."+string(indgen(nMocks),format='(I0)')+suffix+".fits"
groupFileArr=mockDir+"group_mock2."+string(indgen(nMocks),format='(I0)')+suffix+".fits"
outDir="../outfiles/"
spawn,'rm '+outDir+'magdist.?.dat'

minZ=0.2
memThresh=0.5

;------------------------------------
; Magnitude and Distance from center
;------------------------------------
;   magBins=[16,19,20,21,22.5,24.2]
magBins=[17,20,21,22.5,24.2]
nMagBins=n_elements(magBins)

;   distBins=findgen(6)/5
distBins=findgen(5)/4
distBins[0]=0.01                ; exclude centrals
nDistBins=n_elements(distBins)

magVals=fltarr(nMagBins-1,nDistBins-1)
distVals=fltarr(nMagBins-1,nDistBins-1)

totPhotNonMockMem=fltarr(nMagBins-1,nDistBins-1)
totPhotMemMockNon=fltarr(nMagBins-1,nDistBins-1)
totPhotMemMockMem=fltarr(nMagBins-1,nDistBins-1)

for nm=0,nMocks-1 do begin
   nmStr=string(nm,format='(I0)') ; for use with filenames later

   ; open catalogs
   group=mrdfits(groupFileArr[nm],1)
   acs=mrdfits(acsFileArr[nm],1)

   ; exclude galaxies fainter than 24.2 and stellar mass < 8.5 and z < 0.2
   acs=acs[where(acs.mag_auto LT 24.2 AND acs.kevin_mstar GT 8.5 AND acs.z GT minZ)]

;;groups in cat have now been cut to z range and LX sensitivity limit
;;in find_bcg_mock.pro
;   ; exclude groups outside 0.2 < z < 1 and M < 13
;   group=group[where(group.zphot GT minZ AND group.zphot LT 1.0 AND group.lensing_m200 GT 13.)]

   nGroups=n_elements(group)
   nAcs=n_elements(acs)

; exclude groups for which id_mmgg_scale != id_mmgg_scale_specz
;bad=where(group.flag_include EQ 1 AND group.id_mmgg_scale NE group.id_mmgg_scale_specz,nbad)
;group[bad].flag_include=0
;for i=0,nbad-1 do begin
;   tmp=where(acs.group_id_best EQ group[bad[i]].id,ntmp)
;   if(ntmp GT 0) then acs[tmp].group_flag_best=0
;endfor

   ; keep track of galaxies that belong in good halos that are NOT in the
   ; catalog due to edge effects, etc.
   haloInd=intarr(nAcs)
   groupInd=lonarr(nAcs)
   for ii=0L,nAcs-1 do begin
      haloInd[ii]=(where(group.id EQ acs[ii].haloid))[0] ; -1 if halo is missing
      groupInd[ii]=where(group.id EQ acs[ii].group_id_best)
   endfor
   bad=where(haloInd EQ -1,complement=good)
   d_halo_r200=fltarr(nAcs)
   d_halo_r200[bad]=-1.
   d_halo_r200[good]=distance(acs[good].alpha_j2000,acs[good].delta_j2000,group[haloInd[good]].alpha_j2000,group[haloInd[good]].delta_j2000)*3600./group[haloInd[good]].lensing_r200_as

   for ii=0,nMagBins-2 do begin
      for jj=0,ndistbins-2 do begin
         photNonMockMem=where(acs.p_mem_best LT memThresh $
                              AND acs.mag_auto GT magbins[ii] $
                              AND acs.mag_auto LE magbins[ii+1] $
                              AND haloInd GE 0 $
                              AND d_halo_r200 LT 1. $
                              AND d_halo_r200 GE distbins[jj] $
                              AND d_halo_r200 LT distbins[jj+1] $
                              , nPhotNonMockMem)
         photMemMockNon=where(acs.p_mem_best GE memThresh $
                              AND acs.mag_auto GT magbins[ii] $
                              AND acs.mag_auto LE magbins[ii+1] $
                              AND haloInd EQ -1 $
                              AND acs.dist_bcg_r200 GE distbins[jj] $
                              AND acs.dist_bcg_r200 LT distbins[jj+1] $
                              , nPhotMemMockNon)
         photMemMockMem=where(acs.p_mem_best GE memThresh $
                              AND acs.mag_auto GT magbins[ii] $
                              AND acs.mag_auto LE magbins[ii+1] $
                              AND acs.group_id_best EQ acs.haloid $
                              AND d_halo_r200 LT 1. $
                              AND acs.dist_bcg_r200 GE distbins[jj] $
                              AND acs.dist_bcg_r200 LT distbins[jj+1] $
                              , nPhotMemMockMem)
         
         totPhotNonMockMem[ii,jj]+=nPhotNonMockMem
         totPhotMemMockNon[ii,jj]+=nPhotMemMockNon
         totPhotMemMockMem[ii,jj]+=nPhotMemMockMem
         
                                ; alt method - used sum of membership
                                ;              probabilties vs. total
                                ;              number of true members
                                ;              to get purity and completeness
         photMem_w=where(acs.p_mem_best GT 0 $
                         AND acs.mag_auto GT magbins[ii] $
                         AND acs.mag_auto LE magbins[ii+1] $
                         AND acs.dist_bcg_r200 GE distbins[jj] $
                         AND acs.dist_bcg_r200 LT distbins[jj+1] $
                         , nPhotMem_w)
         if(nPhotMem_w GT 0) then nPhotMem_w=total(acs[photMem_w].p_mem_best)
         mockMem_w=where(acs.mag_auto GT magbins[ii] $
                         AND acs.mag_auto LE magbins[ii+1] $
                         AND acs.dist_bcg_r200 GE distbins[jj] $
                         AND acs.dist_bcg_r200 LT distbins[jj+1] $
                         AND haloInd GE 0 $
                         AND d_halo_r200 LT 1. $
                         , nMockMem_w)
         
         magVals[ii,jj]+=mean(acs[photMem_w].mag_auto)
         distVals[ii,jj]+=mean(acs[photMem_w].dist_bcg_r200)
         
;         print_result,outDir+'magdist.'+nmStr+'.dat',magVals[ii,jj],distVals[ii,jj],nPhotMemMockMem,nPhotMemMockNon+nPhotMemMockMem,nPhotNonMockMem+nPhotMemMockMem,nPhotMem_w,nMockMem_w
      endfor
   endfor
endfor   ; end loop over lightcones

magVals=magVals/nMocks
distVals=distVals/nMocks

totPhotMem=totPhotMemMockNon + totPhotMemMockMem
totMockMem=totPhotNonMockMem + totPhotMemMockMem

purity=float(totPhotMemMockMem)/totPhotMem
contamination=1.-purity
completeness=float(totPhotMemMockMem)/totMockMem


!p.multi=[0,2,1]
!p.font=0
!p.charthick=1.2
!p.charsize=1.2
!p.thick=3
!x.thick=3
!y.thick=3

simpctable

; P vs R
plot,/nodata,[0,1],[0,1],xtitle=textoidl('R/R_{200c}'),ytitle='Contamination (%)',color=!black,background=!white
for ii=0,nMagBins-2 do begin
   oplot,distVals[ii,*],contamination[ii,*],linestyle=ii,color=ii
endfor

; P vs mag
plot,/nodata,[19,24],[0,1],xtitle=textoidl('F814W'),ytitle='Contamination (%)',color=!black,background=!white
for ii=0,nDistBins-2 do begin
   oplot,magVals[*,ii],contamination[*,ii],linestyle=ii,color=ii
endfor


stop
end

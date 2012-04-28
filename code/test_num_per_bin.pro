pro test_num_per_bin

acs=mrdfits("../../code/lensing18_20110914.fits",1)
group=mrdfits("../../code/group5_20110914.fits",1)
group=group[where(group.flag_include EQ 1,nGroups)]

;smbins=[9.4,9.8,10.3,10.7,11.5]
smbins=[9.4,11.5]
;rbins=[0.01,0.33,0.66,1.0]
rbins=[0.01,1.0]
zbins=[0.2,0.5,0.85,1.0]

minGroupMass=13.6
maxGroupMass=14.0

; get group masses for member galaxies
halomass=fltarr(n_elements(acs))
for gg=0,ngroups-1 do begin
   mem=where(acs.group_id_best EQ group[gg].id,nmem)
   halomass[mem]=group[gg].lensing_m200
endfor


morph=acs.zest_type + 0.1*acs.zest_bulge
morph[where(acs.zest_type EQ 1)]=1.0
morph[where(acs.zest_type EQ 2 AND acs.zest_bulge EQ 0)]=2.0
morph[where(acs.zest_type EQ 2 AND acs.zest_bulge EQ 1)]=2.3
morph[where(acs.zest_type EQ 2 AND acs.zest_bulge EQ 2)]=2.6
morph[where(acs.zest_type EQ 2 AND acs.zest_bulge EQ 3)]=2.9
morph[where(acs.zest_type EQ 3)]=3.3
morph[where(acs.zest_type EQ 0)]=4.0
morph[where(acs.zest_type EQ 9)]=4.3

color=acs.mnuv_mr
simpctable

for zz=0,n_elements(zbins)-2 do begin
   for sm=0,n_elements(smbins)-2 do begin
      for rr=0,n_elements(rbins)-2 do begin
         sel=where(acs.p_mem_best GT 0.5 $
                   AND halomass GE minGroupMass $
                   AND halomass LT maxGroupMass $
                   AND acs.kevin_mstar GE smbins[sm] $
                   AND acs.kevin_mstar LT smbins[sm+1] $
                   AND acs.dist_bcg_r200 GE rbins[rr] $
                   AND acs.dist_bcg_r200 LT rbins[rr+1] $
                   AND acs.photoz_non_comb GE zbins[zz] $
                   AND acs.photoz_non_comb LT zbins[zz+1] $
                  ,nSel)
         
         print,zbins[zz],smbins[sm],rbins[rr],nSel
         
         plotcircle
;         plot,color[sel],morph[sel],psym=8,xtit='M(NUV)-M(R)',ytit='Morph',yr=[0,5],xr=[-2,7]
;         oplot,[3.5,3.5],[0,5]
;         oplot,[-2,7],[2.1,2.1]
         plot,/nodata,[0,1],[9,12],xtit='R',ytit='SM',yst=1,xst=1
         be=where(color[sel] LT 3.5 AND morph[sel] LE 2.0, nbe)
         re=where(color[sel] GT 3.5 AND morph[sel] LE 2.0, nre)
         bd=where(color[sel] LT 3.5 AND morph[sel] GT 2.0, nbd)
         rd=where(color[sel] GT 3.5 AND morph[sel] GT 2.0, nrd)
         if(nbe GT 0) then oplot,acs[sel[be]].dist_bcg_r200,acs[sel[be]].kevin_mstar,color=!blue,ps=8
         if(nre GT 0) then oplot,acs[sel[re]].dist_bcg_r200,acs[sel[re]].kevin_mstar,color=!red,ps=8
         if(nbd GT 0) then oplot,acs[sel[bd]].dist_bcg_r200,acs[sel[bd]].kevin_mstar,color=!blue,ps=4
         if(nrd GT 0) then oplot,acs[sel[rd]].dist_bcg_r200,acs[sel[rd]].kevin_mstar,color=!red,ps=4

         stop
      endfor
   endfor
endfor

stop
end

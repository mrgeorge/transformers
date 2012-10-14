pro emulate_kb_quench

; create a quenched catalog similar to Bundy et al. 2010 color-color
; cut for selecting quenched galaxies, excluding dusty star-formers.

acsFile="~/data/cosmos/code/lensing18_20110914.fits"
quenchOutFile="~/data/cosmos/catalogs/quench_mrg.fits"

acs=mrdfits(acsFile,1)

quench=intarr(n_elements(acs))
; -2 = bad photometry (some abs mag < -30)
; -1 = not selected, low SM or outside z range
;  0 = selected, but not quenched
;  1 = selected and quenched
quench[*]=-1

; only consider galaxies that go into the transformer analysis, make 3
; rough cuts for different z bins

sel=where(acs.kevin_mstar GT 9.8 $
          AND acs.photoz_non_comb GE 0.2 $
          AND acs.photoz_non_comb LE 1.0 $
          AND acs.type EQ 0)
quench[sel]=0

lowz=where(acs[sel].photoz_non_comb LT 0.5)
lowzred=where(acs[sel[lowz]].mnuv_mr GT 3.5)
midz=where(acs[sel].photoz_non_comb GE 0.5 AND acs.photoz_non_comb LT 0.8)
midzred=where(acs[sel[midz]].mnuv_mr GT 3.5)
highz=where(acs[sel].photoz_non_comb GE 0.8)
highzred=where(acs[sel[highz]].mnuv_mr GT 3.5)


xr=[-0.5,2.]
yr=[0,6]
simpctable
plotcircle,0.5
!p.multi=[0,3,1]

;low z
plot,acs[sel[lowz]].mr-acs[sel[lowz]].mj,acs[sel[lowz]].mnuv-acs[sel[lowz]].mr,ps=3,xr=xr,yr=yr,xst=1,yst=1,xtit='R-J',ytit='NUV-R'
oplot,acs[sel[lowz[lowzred]]].mr-acs[sel[lowz[lowzred]]].mj,acs[sel[lowz[lowzred]]].mnuv-acs[sel[lowz[lowzred]]].mr,ps=3,color=!red

bluesph=where(acs[sel[lowz]].kevin_mstar GT 10.7 AND $
              (acs[sel[lowz]].zest_type EQ 1 OR (acs[sel[lowz]].zest_type EQ 2 AND acs[sel[lowz]].zest_bulge EQ 0))); $
;              AND acs[sel[lowz]].mnuv_mr LT 3.5)
oplot,acs[sel[lowz[bluesph]]].mr-acs[sel[lowz[bluesph]]].mj,acs[sel[lowz[bluesph]]].mnuv-acs[sel[lowz[bluesph]]].mr,ps=8,color=!blue

nuvrCut=4.
rjSlope=3.5
rjIntercept=1.5

oplot,!x.crange,[1,1]*nuvrCut
oplot,!x.crange,rjSlope*!x.crange+rjIntercept

lowzquench=where(acs[sel[lowz]].mnuv - acs[sel[lowz]].mr GE nuvrCut $
                 AND acs[sel[lowz]].mnuv - acs[sel[lowz]].mr GE rjSlope*(acs[sel[lowz]].mr - acs[sel[lowz]].mj) + rjIntercept)
quench[sel[lowz[lowzquench]]]=1


;mid z
plot,acs[sel[midz]].mr-acs[sel[midz]].mj,acs[sel[midz]].mnuv-acs[sel[midz]].mr,ps=3,xr=xr,yr=yr,xst=1,yst=1,xtit='R-J',ytit='NUV-R'
oplot,acs[sel[midz[midzred]]].mr-acs[sel[midz[midzred]]].mj,acs[sel[midz[midzred]]].mnuv-acs[sel[midz[midzred]]].mr,ps=3,color=!red

bluesph=where(acs[sel[midz]].kevin_mstar GT 10.7 AND $
              (acs[sel[midz]].zest_type EQ 1 OR (acs[sel[midz]].zest_type EQ 2 AND acs[sel[midz]].zest_bulge EQ 0))); $
;              AND acs[sel[midz]].mnuv_mr LT 3.5)
oplot,acs[sel[midz[bluesph]]].mr-acs[sel[midz[bluesph]]].mj,acs[sel[midz[bluesph]]].mnuv-acs[sel[midz[bluesph]]].mr,ps=8,color=!blue

nuvrCut=3.7

oplot,!x.crange,[1,1]*nuvrCut
oplot,!x.crange,rjSlope*!x.crange+rjIntercept

midzquench=where(acs[sel[midz]].mnuv - acs[sel[midz]].mr GE nuvrCut $
                 AND acs[sel[midz]].mnuv - acs[sel[midz]].mr GE rjSlope*(acs[sel[midz]].mr - acs[sel[midz]].mj) + rjIntercept)
quench[sel[midz[midzquench]]]=1


;high z
plot,acs[sel[highz]].mr-acs[sel[highz]].mj,acs[sel[highz]].mnuv-acs[sel[highz]].mr,ps=3,xr=xr,yr=yr,xst=1,yst=1,xtit='R-J',ytit='NUV-R'
oplot,acs[sel[highz[highzred]]].mr-acs[sel[highz[highzred]]].mj,acs[sel[highz[highzred]]].mnuv-acs[sel[highz[highzred]]].mr,ps=3,color=!red

bluesph=where(acs[sel[highz]].kevin_mstar GT 10.7 AND $
              (acs[sel[highz]].zest_type EQ 1 OR (acs[sel[highz]].zest_type EQ 2 AND acs[sel[highz]].zest_bulge EQ 0))); $
;              AND acs[sel[highz]].mnuv_mr LT 3.5)
oplot,acs[sel[highz[bluesph]]].mr-acs[sel[highz[bluesph]]].mj,acs[sel[highz[bluesph]]].mnuv-acs[sel[highz[bluesph]]].mr,ps=8,color=!blue

nuvrCut=3.7

oplot,!x.crange,[1,1]*nuvrCut
oplot,!x.crange,rjSlope*!x.crange+rjIntercept

highzquench=where(acs[sel[highz]].mnuv - acs[sel[highz]].mr GE nuvrCut $
                 AND acs[sel[highz]].mnuv - acs[sel[highz]].mr GE rjSlope*(acs[sel[highz]].mr - acs[sel[highz]].mj) + rjIntercept)
quench[sel[highz[highzquench]]]=1





bad=where(acs.mnuv LT -30 $
          OR acs.mr LT -30 $
          OR acs.mj LT -30)
quench[bad]=-2

quench_cat=create_struct('IDENT',0L,'RA',0.D,'DEC',0.D,'QUENCHED',0)
quench_cat=replicate(quench_cat,n_elements(acs))
quench_cat.ident=acs.ident
quench_cat.ra=acs.alpha_j2000
quench_cat.dec=acs.delta_j2000
quench_cat.quenched=quench

;mwrfits,quench_cat,quenchOutFile,/create

end

pro emulate_kb_quench

; create a quenched catalog similar to Bundy et al. 2010 color-color
; cut for selecting quenched galaxies, excluding dusty star-formers.

acsFile="~/data/cosmos/code/lensing18_20110914.fits"

acs=mrdfits(acsFile,1)

; only consider galaxies that go into the transformer analysis, make 3
; rough cuts for different z bins

sel=where(acs.kevin_mstar GT 9.8 $
          AND acs.photoz_non_comb GE 0.2 $
          AND acs.photoz_non_comb LE 1.0 $
          AND acs.type EQ 0)

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
              (acs[sel[lowz]].zest_type EQ 1 OR (acs[sel[lowz]].zest_type EQ 2 AND acs[sel[lowz]].zest_bulge EQ 0)) $
              AND acs[sel[lowz]].mnuv_mr LT 3.5)
oplot,acs[sel[lowz[bluesph]]].mr-acs[sel[lowz[bluesph]]].mj,acs[sel[lowz[bluesph]]].mnuv-acs[sel[lowz[bluesph]]].mr,ps=8,color=!blue


oplot,!x.crange,[1,1]*4.
oplot,!x.crange,3.5*!x.crange+1.5

;mid z
plot,acs[sel[midz]].mr-acs[sel[midz]].mj,acs[sel[midz]].mnuv-acs[sel[midz]].mr,ps=3,xr=xr,yr=yr,xst=1,yst=1,xtit='R-J',ytit='NUV-R'
oplot,acs[sel[midz[midzred]]].mr-acs[sel[midz[midzred]]].mj,acs[sel[midz[midzred]]].mnuv-acs[sel[midz[midzred]]].mr,ps=3,color=!red

oplot,!x.crange,[1,1]*3.7
oplot,!x.crange,3.5*!x.crange+1.5


;high z
plot,acs[sel[highz]].mr-acs[sel[highz]].mj,acs[sel[highz]].mnuv-acs[sel[highz]].mr,ps=3,xr=xr,yr=yr,xst=1,yst=1,xtit='R-J',ytit='NUV-R'
oplot,acs[sel[highz[highzred]]].mr-acs[sel[highz[highzred]]].mj,acs[sel[highz[highzred]]].mnuv-acs[sel[highz[highzred]]].mr,ps=3,color=!red

oplot,!x.crange,[1,1]*3.7
oplot,!x.crange,3.5*!x.crange+1.5


stop
end

pro add_morph

acsInFile="~/data/cosmos/code/lensing18_20110914.fits"
acsOutFile="~/data/cosmos/code/lensing18_20110914_morph.fits"

morphFile="~/data/cosmos/catalogs/cosmos_morph_tasca_1.1.tbl"
quenchFile="~/data/cosmos/catalogs/quenched_cat_nocuts.fits"

acs=mrdfits(acsInFile,1)
sel=where(acs.mag_auto LT 24.2 $
          AND acs.kevin_mstar GT 9 $
          AND acs.photoz_non_comb GT 0 $
          AND acs.photoz_non_comb LT 1)
; keep the full acs file intact, but only match on these good sources

readcol,morphFile,id,ra,dec,mag,ab,rh,gini,conc,asym,m1,m2,m3
; 1=early type, 2=spirals, 3=irregulars (0=unclassified?)
; m1=class_int - Tasca preferred
; m2=class_linee - Abraham
; m3=class_svn - Huertas-Company

close_match_radec,acs[sel].alpha_j2000,acs[sel].delta_j2000,ra,dec,match1,match2,1./3600,1

; create arrays to be added to ACS catalog (0 will be default for missing or unclassified)
nAcs=n_elements(acs)
morph_tasca1=intarr(nAcs)
morph_tasca2=intarr(nAcs)
morph_tasca3=intarr(nAcs)

morph_tasca1[sel[match1]]=m1[match2]
morph_tasca2[sel[match1]]=m2[match2]
morph_tasca3[sel[match1]]=m3[match2]

acs=mrg_addcol(acs,'MORPH_TASCA1',morph_tasca1)
acs=mrg_addcol(acs,'MORPH_TASCA2',morph_tasca2)
acs=mrg_addcol(acs,'MORPH_TASCA3',morph_tasca3)


; now match with quenched catalog
quench_cat=mrdfits(quenchFile,1)

close_match_radec,acs[sel].alpha_j2000,acs[sel].delta_j2000,quench_cat.ra,quench_cat.dec,match1,match2,1./3600,1,miss

; create array for ACS catalog, -2 will be default for objects w/o matches (quench cat has values of -1, 0, and 1)
kb_quench=intarr(nAcs)
kb_quench[*]=-2
kb_quench[sel[match1]]=quench_cat[match2].quenched

acs=mrg_addcol(acs,"KB_QUENCH",kb_quench)

mwrfits,acs,acsOutFile,/create

end

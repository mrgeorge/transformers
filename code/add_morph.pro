pro add_morph

acsInFile="~/data/cosmos/code/lensing18_20110914.fits"
acsOutFile="~/data/cosmos/code/lensing18_20110914_morph.fits"

tascaFile="~/data/cosmos/catalogs/cosmos_morph_tasca_1.1.tbl"
cassataFile="~/data/cosmos/catalogs/cosmos_morph_cassata_1.1.tbl"
morph2005File="~/data/cosmos/catalogs/cosmos_morphology_2005.tbl"
quenchFile="~/data/cosmos/catalogs/quench_mrg.fits"

acs=mrdfits(acsInFile,1)
sel=where(acs.mag_auto LT 24.2 $
          AND acs.kevin_mstar GT 9 $
          AND acs.photoz_non_comb GT 0 $
          AND acs.photoz_non_comb LT 1)
; keep the full acs file intact, but only match on these good sources

; create arrays to be added to ACS catalog (0 will be default for missing or unclassified)
nAcs=n_elements(acs)
morph_tasca1=intarr(nAcs)
morph_tasca2=intarr(nAcs)
morph_tasca3=intarr(nAcs)
morph_cassata=intarr(nAcs)
morph2005_gini=fltarr(nAcs)
morph2005_conc=fltarr(nAcs)


; TASCA MORPHOLOGY
readcol,tascaFile,id,ra,dec,mag,ab,rh,gini,conc,asym,m1,m2,m3
; 1=early type, 2=spirals, 3=irregulars (0=unclassified?)
; m1=class_int - Tasca preferred
; m2=class_linee - Abraham
; m3=class_svn - Huertas-Company

close_match_radec,acs[sel].alpha_j2000,acs[sel].delta_j2000,ra,dec,match1,match2,1./3600,1

morph_tasca1[sel[match1]]=m1[match2]
morph_tasca2[sel[match1]]=m2[match2]
morph_tasca3[sel[match1]]=m3[match2]

acs=mrg_addcol(acs,'MORPH_TASCA1',morph_tasca1)
acs=mrg_addcol(acs,'MORPH_TASCA2',morph_tasca2)
acs=mrg_addcol(acs,'MORPH_TASCA3',morph_tasca3)

undefine,id
undefine,ra
undefine,dec
undefine,mag
undefine,ab
undefine,rh
undefine,gini
undefine,conc
undefine,asym
undefine,match1
undefine,match2
undefine,morph_tasca1
undefine,morph_tasca2
undefine,morph_tasca3

; CASSATA MORPHOLOGY
readcol,cassataFile,id,ra,dec,mag,rpetro,rhalf,conc,asym,gini,m20,ab,class,weight

close_match_radec,acs[sel].alpha_j2000,acs[sel].delta_j2000,ra,dec,match1,match2,1./3600,1

morph_cassata[sel[match1]]=class[match2]

acs=mrg_addcol(acs,'MORPH_CASSATA',morph_cassata)

undefine,id
undefine,ra
undefine,dec
undefine,mag
undefine,rpetro
undefine,rhalf
undefine,conc
undefine,asym
undefine,gini
undefine,m20
undefine,ab
undefine,class
undefine,weight
undefine,match1
undefine,match2
undefine,moprh_cassata

; 2005 MORPHOLOGY
morph2005=read_ipac_table(morph2005File)

close_match_radec,acs[sel].alpha_j2000,acs[sel].delta_j2000,morph2005.ra,morph2005.dec,match1,match2,1./3600,1

morph2005_gini[sel[match1]]=(morph2005.gini)[match2] ; NOTE RHS indexing is *much* faster than morph2005[match2].gini
morph2005_conc[sel[match1]]=(morph2005.con)[match2]

acs=mrg_addcol(acs,'MORPH2005_GINI',morph2005_gini)
acs=mrg_addcol(acs,'MORPH2005_CONC',morph2005_conc)

undefine,morph2005
undefine,match1
undefine,match2
undefine,morph2005_gini
undefine,morph2005_conc

; QUENCHED CATALOG (from emulate_kb_quench)
quench_cat=mrdfits(quenchFile,1)

; quench_cat is already matched line by line to acs
acs=mrg_addcol(acs,"QUENCH_MRG",quench_cat.quenched)

undefine,quench_cat

mwrfits,acs,acsOutFile,/create

end

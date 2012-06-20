pro bootstrap, arr, average, sigma
  nSample=1000
  size=n_elements(arr)
  if(size LT nSample^(1./size)) then begin
     print, 'Too few elements in array for bootstrapping'
     stop
  endif

  sampleInd=round(randomu(seed,size,nSample)*(size-1))
  avgArr=total(arr[sampleInd],1)/size
  average=mean(avgArr)
  sigma=stddev(avgArr)
end

pro plot_quenchr, loz=loz,midz=midz,hiz=hiz,lom=lom,red=red,blue=blue

; like plot_quenchz but with radial bins at 0.2<z<0.8, 13.6<M<14.

acs=mrdfits("~/data/cosmos/code/lensing18_20110914.fits",1)
group=mrdfits("~/data/cosmos/code/group5_20110914.fits",1)
;quench_cat=mrdfits("/Users/mgeorge/data/cosmos/catalogs/quenched_cat_nocuts.fits",1)

plotDir="~/data/cosmos/transformers/plots/"

smbins=[9.4,9.8,10.3,10.7,11.5]
;rbins=[0.0,0.25,0.5,1.0]
;rbins=[0.01,0.25,0.5,0.75,1.0]
rbins=[0.01,0.33,0.66,1.0]
sun=sunsymbol()
;rNames=textoidl(['0<R/R_{200}<0.25','0.25<R/R_{200}<0.5','0.5<R/R_{200}<1'])
smNames=['[9.4,9.8)','[9.8,10.3)','[10.3,10.7)','[10.7,11.5)']



nSMbins=n_elements(smbins)-1
nrbins=n_elements(rbins)-1

if(keyword_set(loz)) then begin
   minz=0.2
   maxz=0.8
   fillArr=[0,1,1,1]
   minGroupMass=13.6
   maxGroupMass=14.0
   plot_suffix='_loz'
endif else if(keyword_set(midz)) then begin
   minz=0.5
   maxz=0.85
   fillArr=[0,1,1,1]
   minGroupMass=13.6
   maxGroupMass=14.0
   plot_suffix='_midz'
endif else if(keyword_set(hiz)) then begin
   minz=0.80
   maxz=1.0
   fillArr=[0,0,1,1]
   minGroupMass=13.6
   maxGroupMass=14.0
   plot_suffix='_hiz'
endif else if(keyword_set(lom)) then begin
   minz=0.2
   maxz=0.5
   fillArr=[1,1,1,1]
   minGroupMass=13.3
   maxGroupMass=14.6
   plot_suffix='_lom'
endif

; trim catalogs to relevant cuts
sel=where(acs.kevin_mstar GT 2 $
          AND acs.photoz_non_comb GE minz $
          AND acs.photoz_non_comb LT maxz $
          AND acs.mag_auto LT 24.2 $
          AND acs.good_zphot_lens EQ 1 $
          AND acs.mu_class EQ 1 $
         )
acs=acs[sel]
sel=where(group.flag_include EQ 1 $
          AND group.zphot GT minz $
          AND group.zphot LE maxz $
          ,ngroups)
group=group[sel]

; associate quenched catalog (older photometric cat) with new acs file
;close_match_radec,acs.alpha_j2000,acs.delta_j2000,quench_cat.ra,quench_cat.dec,asub,qsub,1./3600,1,miss

; array of quench flags, one entry for each object in acs. 1=quenched
;quench=quench_cat[qsub].quenched

; eliminate misses from acs catalog
;acs=acs[asub]

; eliminate objects with quench flagged as -1
;noflag=where(quench GE 0)
;quench=quench[noflag]
;acs=acs[noflag]

; use Ilbert's color classification
quench=(acs.mnuv_mr GT 3.5)

; use ZEST type 1,2.0
nMorph=3 ; 3 morphological categories
morph=intarr(nMorph,n_elements(acs))
morph[0,*]=(acs.zest_type EQ 1 $       ; ellipticals and spheroidals
     OR (acs.zest_type EQ 2 AND acs.zest_bulge EQ 0))
morph[1,*]=(acs.zest_type EQ 2 AND acs.zest_bulge EQ 1) ; intermediate disks
morph[2,*]=(acs.zest_type EQ 2 AND (acs.zest_bulge EQ 2 OR acs.zest_bulge EQ 3)) ; bulgeless disks

if(keyword_set(red)) then begin
   color_suffix='_red'
   morph *= rebin(transpose(quench),nMorph,n_elements(quench))
   quenchSel=quench
endif else if(keyword_set(blue)) then begin
   color_suffix='_blue'
   morph *= rebin(~(transpose(quench)),nMorph,n_elements(quench))
   quenchSel=~quench
endif else begin
   color_suffix=''
   quenchSel=replicate(1,n_elements(quench))
endelse

; get group masses for member galaxies
halomass=fltarr(n_elements(acs))
for gg=0,ngroups-1 do begin
   mem=where(acs.group_id_best EQ group[gg].id,nmem)
   halomass[mem]=group[gg].lensing_m200
endfor

; group fractions by bin (i.e. lump all galaxies in groups together)
quenchFracMem=fltarr(nMorph,nSMbins,nrbins)
quenchFracMemErr=fltarr(nMorph,nSMbins,nrbins)
quenchFracMemBoot=fltarr(nMorph,nSMbins,nrbins)
quenchFracMemErrBoot=fltarr(nMorph,nSMbins,nrbins)
morphFracMem=fltarr(nMorph,nSMbins,nrbins)
morphFracMemErr=fltarr(nMorph,nSMbins,nrbins)
morphFracMemBoot=fltarr(nMorph,nSMbins,nrbins)
morphFracMemErrBoot=fltarr(nMorph,nSMbins,nrbins)
memSM=fltarr(nSMbins,nrbins)
memSMLo=fltarr(nSMbins,nrbins)
memSMHi=fltarr(nSMbins,nrbins)
memR=fltarr(nSMbins,nrbins)
memRLo=fltarr(nSMbins,nrbins)
memRHi=fltarr(nSMbins,nrbins)

; centrals
quenchFracCen=fltarr(nMorph,nSMbins)
quenchFracCenBoot=fltarr(nMorph,nSMbins)
quenchFracCenErrBoot=fltarr(nMorph,nSMbins)
morphFracCen=fltarr(nMorph,nSMbins)
morphFracCenBoot=fltarr(nMorph,nSMbins)
morphFracCenErrBoot=fltarr(nMorph,nSMbins)

; field fractions (lumping all field galaxies together)
quenchFracField=fltarr(nMorph,nSMbins)
quenchFracFieldErr=fltarr(nMorph,nSMbins)
morphFracField=fltarr(nMorph,nSMbins)
morphFracFieldErr=fltarr(nMorph,nSMbins)
fieldSM=fltarr(nSMbins)
fieldSMLo=fltarr(nSMbins)
fieldSMHi=fltarr(nSMbins)

for sm=0,nSMbins-1 do begin
   ; centrals
   cen=where(acs.mmgg_scale EQ 1 $
             AND haloMass GE minGroupMass $
             AND haloMass LT maxGroupMass $
             AND acs.kevin_mstar GE smbins[sm] $
             AND acs.kevin_mstar LT smbins[sm+1] $
             AND acs.photoz_non_comb GE minz $
             AND acs.photoz_non_comb LT maxz $
             AND quenchSel $
             , ncen)
   if(ncen GT 5) then begin
      for mm=0,nMorph-1 do begin
         tmp=where(morph[mm,cen] EQ 1,nMorphCen)
         morphFracCen[mm,sm]=float(nMorphCen)/ncen
         bootstrap,morph[mm,cen],mean,sigma
         morphFracCenBoot[mm,sm]=mean
         morphFracCenErrBoot[mm,sm]=sigma

         if(nMorphCen GT 5) then begin
            tmp2=where(quench[cen[tmp]] EQ 1,nQuenchCen)
            quenchFracCen[mm,sm]=float(nQuenchCen)/nMorphCen
            bootstrap,quench[cen[tmp]],mean,sigma
            quenchFracCenBoot[mm,sm]=mean
            quenchFracCenErrBoot[mm,sm]=sigma
         endif else begin
            quenchFracCen[mm,sm]=!values.f_nan
            quenchFracCenBoot[mm,sm]=!values.f_nan
            quenchFracCenErrBoot[mm,sm]=!values.f_nan
         endelse
      endfor
   endif else begin             ; not enough objects, don't plot
      morphFracCen[*,sm]=!values.f_nan
      morphFracCenBoot[*,sm]=!values.f_nan
      morphFracCenErrBoot[*,sm]=!values.f_nan

      quenchFracCen[*,sm]=!values.f_nan
      quenchFracCenBoot[*,sm]=!values.f_nan
      quenchFracCenErrBoot[*,sm]=!values.f_nan
   endelse

   ; get info for field in this SM bin
   field=where(acs.p_mem_best EQ 0 $
               AND acs.kevin_mstar GE smbins[sm] $
               AND acs.kevin_mstar LT smbins[sm+1] $
               AND acs.photoz_non_comb GE minz $
               AND acs.photoz_non_comb LT maxz $
               AND quenchSel $
               ,nfield)
   if(nfield GT 10) then begin
      for mm=0,nMorph-1 do begin
         tmp=where(morph[mm,field] EQ 1,nMorphField)
         morphFracField[mm,sm]=float(nMorphField)/nfield
         ; if used, should be changed to binomial or bootstrap - morphFracFieldErr[mm,sm,rr]=morphFracField[mm,sm,rr]*sqrt(1./nMorphField+1./nfield)

         if(nMorphField GT 10) then begin
            tmp2=where(quench[field[tmp]] EQ 1,nQuenchField)
            quenchFracField[mm,sm]=float(nQuenchField)/nMorphField
            ; if used, should be changed to binomial or bootstrap - quenchFracFieldErr[sm]=quenchFracField[sm]*sqrt(1./nQuenchField+1./nfield)
         endif else begin
            quenchFracField[mm,sm]=!values.f_nan
         endelse
      endfor

      fieldSM[sm]=alog10(mean(10.^(acs[field].kevin_mstar)))
      fieldSMLo[sm]=fieldSM[sm]-smbins[sm]
      fieldSMHi[sm]=smbins[sm+1]-fieldSM[sm]
   endif else begin
      morphFracField[*,sm]=!values.f_nan
      quenchFracField[*,sm]=!values.f_nan
   endelse

   ; Satellites in radial bins
   for rr=0,nrbins-1 do begin
      ; get info for galaxies in groups in this r,SM bin
      mem=where(acs.group_flag_best EQ 1 $
                AND acs.p_mem_best GT 0.5 $
                AND acs.dist_bcg_r200 GE rbins[rr] $
                AND acs.dist_bcg_r200 LT rbins[rr+1] $
                AND haloMass GE minGroupMass $
                AND haloMass LT maxGroupMass $
                AND acs.kevin_mstar GE smbins[sm] $
                AND acs.kevin_mstar LT smbins[sm+1] $
                AND acs.photoz_non_comb GE minz $
                AND acs.photoz_non_comb LT maxz $
                AND quenchSel $
                ,nmem)
      if(nmem GT 5) then begin
         for mm=0,nMorph-1 do begin
            tmp=where(morph[mm,mem] EQ 1,nMorphMem)
            morphFracMem[mm,sm,rr]=float(nMorphMem)/nmem
;         ellFracMemErr[sm,rr]=ellFracMem[sm,rr]*sqrt(1./nEllMem+1./nmem)
            morphFracMemErr[mm,sm,rr]=sqrt(1.*(nMorphMem+1)*(nmem-nMorphMem+1)/((nmem+3)*(nmem+2)^2)) ; binomial stddev - see http://www.roma1.infn.it/~dagos/proportions/node3.html
            bootstrap,morph[mm,mem],mean,sigma
            morphFracMemBoot[mm,sm,rr]=mean
            morphFracMemErrBoot[mm,sm,rr]=sigma

            if(nMorphMem GT 5) then begin
               tmp2=where(quench[mem[tmp]] EQ 1,nQuenchMem)
               quenchFracMem[mm,sm,rr]=float(nQuenchMem)/nMorphMem
     ;         quenchFracMemErr[sm,rr]=quenchFracMem[sm,rr]*sqrt(1./nQuenchMem+1./nmem)
               quenchFracMemErr[mm,sm,rr]=sqrt(1.*(nQuenchMem+1)*(nmem-nQuenchMem+1)/((nmem+3)*(nmem+2)^2)) ; binomial stddev - see http://www.roma1.infn.it/~dagos/proportions/node3.html
               bootstrap,quench[mem[tmp]],mean,sigma
               quenchFracMemBoot[mm,sm,rr]=mean
               quenchFracMemErrBoot[mm,sm,rr]=sigma
            endif else begin
               quenchFracMem[mm,sm,rr]=!values.f_nan
               quenchFracMemBoot[mm,sm,rr]=!values.f_nan
               quenchFracMemErrBoot[mm,sm,rr]=!values.f_nan
            endelse
         endfor

         memSM[sm,rr]=alog10(mean(10.^(acs[mem].kevin_mstar)))
         memSMLo[sm,rr]=memSM[sm,rr]-smbins[sm]
         memSMHi[sm,rr]=smbins[sm+1]-memSM[sm,rr]
         memR[sm,rr]=mean(acs[mem].dist_bcg_r200)
         memRLo[sm,rr]=memR[sm,rr]-rbins[rr]
         memRHi[sm,rr]=rbins[rr+1]-memR[sm,rr]
      endif else begin
         morphFracMem[*,sm,rr]=!values.f_nan
         morphFracMemBoot[*,sm,rr]=!values.f_nan
         morphFracMemErrBoot[*,sm,rr]=!values.f_nan

         quenchFracMem[*,sm,rr]=!values.f_nan
         quenchFracMemBoot[*,sm,rr]=!values.f_nan
         quenchFracMemErrBoot[*,sm,rr]=!values.f_nan
      endelse
   endfor
endfor

set_plot,'ps'
simpctable
!p.font=0
!p.thick=3
!x.thick=3
!y.thick=3
!p.charthick=1.2
!p.charsize=1.1

;xr=[9.25,11.65]
xr=[-0.03,1.03]
yr=[-0.03,1.03]
xticknames=['0','0.2','0.4','0.6','0.8','1']
xtickv=[0,0.2,0.4,0.6,0.8,1]
xticks=n_elements(xticknames)-1
; for redshift labels
;xtext=9.35
;ytext=0.93
xtext=0.1
ytext=0.93
; for legend
yspace=0.07
xsymoffset=0.05
ysymoffset=0.02

groupsym=[0,8,4] ; circle, square, triangle
smColors=replicate(!black,nsmbins)
morphColors=[!darkred,!darkgreen,!purple,!orange]
morphOffsets=[-1,0,1]*0.03

; QUENCHING PLOT
device,filename=plotDir+'quenchr'+plot_suffix+color_suffix+'.eps',/color,/helvetica,xsize=10,ysize=4,/inches,/encapsul
multiplot,/default
multiplot,[nsmbins,1],mxtitle=textoidl('Group-centric Distance (R/R_{200c})'),mxtitsize=1.4,mxtitoffset=0.5,mytitle='Quenched Fraction',mytitsize=1.4,mytitoffset=0,/rowmajor

for sm=0,nsmbins-1 do begin
   ; top row - binned
   plot,xr,yr,xstyle=1,ystyle=1,/nodata,ytitle=ytitle,xticks=xticks,xtickv=xtickv,xtickname=xticknames,xminor=4
   xyouts,xtext,ytext,smNames[sm],/data,color=smColors[sm]

   for mm=0,nMorph-1 do begin
      ; plot centrals
      plotsym,groupsym[mm],thick=3,/fill
      if(finite(quenchFracCenErrBoot[mm,sm])) then oploterror,[0],[quenchFracCen[mm,sm]],[quenchFracCenErrBoot[mm,sm]],psym=8,color=morphColors[mm],errcolor=morphColors[mm],errstyle=0,errthick=3

      ; satellites
      if NOT(fillArr[sm]) then plotsym,groupsym[mm],thick=3 ; unfilled
;      oploterror,memR[sm,*],quenchFracMem[sm,*],memRLo[sm,*],quenchFracMemErrBoot[sm,*],psym=8,color=smColors[sm],linestyle=0,errcolor=smColors[sm],errstyle=0,errthick=3,/lobar
;      oploterror,memR[sm,*],quenchFracMem[sm,*],memRHi[sm,*],quenchFracMemErrBoot[sm,*],psym=8,color=smColors[sm],linestyle=0,errcolor=smColors[sm],errstyle=0,errthick=3,/hibar
      if(total(finite(quenchFracMemErrBoot[mm,sm,*])) GT 0) then $
         oploterror,memR[sm,*]+morphOffsets[mm],quenchFracMem[mm,sm,*],quenchFracMemErrBoot[mm,sm,*],psym=-8,color=morphColors[mm],linestyle=0,errcolor=morphColors[mm],errstyle=0,errthick=3

      ; field
      if(finite(quenchFracField[mm,sm])) then $
         oplot,!x.crange,replicate(quenchFracField[mm,sm],2),linestyle=1,color=morphColors[mm]
   endfor

   multiplot
endfor
multiplot,/reset
device,/close


; MORPHOLOGY PLOT
device,filename=plotDir+'morphr'+plot_suffix+color_suffix+'.eps',/color,/helvetica,xsize=10,ysize=4,/inches,/encapsul
multiplot,/default
multiplot,[nsmbins,1],mxtitle=textoidl('Group-centric Distance (R/R_{200c})'),mxtitsize=1.4,mxtitoffset=0.5,mytitle='Morphological Fraction',mytitsize=1.4,mytitoffset=0,/rowmajor

for sm=0,nsmbins-1 do begin
   ; top row - binned
   plot,xr,yr,xstyle=1,ystyle=1,/nodata,ytitle=ytitle,xticks=xticks,xtickv=xtickv,xtickname=xticknames,xminor=4
   xyouts,xtext,ytext,smNames[sm],/data,color=smColors[sm]

   for mm=0,nMorph-1 do begin
      ; plot centrals
      plotsym,groupsym[mm],thick=3,/fill
      if(finite(morphFracCenErrBoot[mm,sm])) then oploterror,[0],[morphFracCen[mm,sm]],[morphFracCenErrBoot[mm,sm]],psym=8,color=morphColors[mm],errcolor=morphColors[mm],errstyle=0,errthick=3

      ; satellites
      if NOT(fillArr[sm]) then plotsym,groupsym[mm],thick=3 ; unfilled
;      oploterror,memR[sm,*],morphFracMem[mm,sm,*],memRLo[sm,*],morphFracMemErrBoot[mm,sm,*],psym=8,color=smColors[sm],linestyle=0,errcolor=smColors[sm],errstyle=0,errthick=3,/lobar
;      oploterror,memR[sm,*],morphFracMem[mm,sm,*],memRHi[sm,*],morphFracMemErrBoot[mm,sm,*],psym=8,color=smColors[sm],linestyle=0,errcolor=smColors[sm],errstyle=0,errthick=3,/hibar
      if(total(finite(morphFracMem[mm,sm,*])) GT 0) then $
         oploterror,memR[sm,*]+morphOffsets[mm],morphFracMem[mm,sm,*],morphFracMemErrBoot[mm,sm,*],psym=-8,color=morphColors[mm],linestyle=0,errcolor=morphColors[mm],errstyle=0,errthick=3
      
      ; field
      if(finite(morphFracField[mm,sm])) then $
         oplot,!x.crange,replicate(morphFracField[mm,sm],2),linestyle=1,color=morphColors[mm]
   endfor

   multiplot
endfor
multiplot,/reset
device,/close

!p.multi=0

end

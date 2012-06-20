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

pro plot_quenchr, loz=loz,midz=midz,hiz=hiz,lom=lom

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
   maxz=0.5
   fillArr=[1,1,1,1]
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
   minz=0.85
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
ell=(acs.zest_type EQ 1 $       ; ellipticals and spheroidals
     OR (acs.zest_type EQ 2 AND acs.zest_bulge EQ 0))
; other ZEST disk categories
disk21=(acs.zest_type EQ 2 AND acs.zest_bulge EQ 1)
disk22=(acs.zest_type EQ 2 AND acs.zest_bulge EQ 2)
disk23=(acs.zest_type EQ 2 AND acs.zest_bulge EQ 3)

; get group masses for member galaxies
halomass=fltarr(n_elements(acs))
for gg=0,ngroups-1 do begin
   mem=where(acs.group_id_best EQ group[gg].id,nmem)
   halomass[mem]=group[gg].lensing_m200
endfor

; group fractions by bin (i.e. lump all galaxies in groups together)
quenchFracMem=fltarr(nSMbins,nrbins)
quenchFracMemErr=fltarr(nSMbins,nrbins)
quenchFracMemBoot=fltarr(nSMbins,nrbins)
quenchFracMemErrBoot=fltarr(nSMbins,nrbins)
ellFracMem=fltarr(nSMbins,nrbins)
ellFracMemErr=fltarr(nSMbins,nrbins)
ellFracMemBoot=fltarr(nSMbins,nrbins)
ellFracMemErrBoot=fltarr(nSMbins,nrbins)
memSM=fltarr(nSMbins,nrbins)
memSMLo=fltarr(nSMbins,nrbins)
memSMHi=fltarr(nSMbins,nrbins)
memR=fltarr(nSMbins,nrbins)
memRLo=fltarr(nSMbins,nrbins)
memRHi=fltarr(nSMbins,nrbins)

; centrals
quenchFracCen=fltarr(nSMbins)
quenchFracCenBoot=fltarr(nSMbins)
quenchFracCenErrBoot=fltarr(nSMbins)
ellFracCen=fltarr(nSMbins)
ellFracCenBoot=fltarr(nSMbins)
ellFracCenErrBoot=fltarr(nSMbins)

; field fractions (lumping all field galaxies together)
quenchFracField=fltarr(nSMbins,nrbins)
quenchFracFieldErr=fltarr(nSMbins,nrbins)
ellFracField=fltarr(nSMbins,nrbins)
ellFracFieldErr=fltarr(nSMbins,nrbins)
fieldSM=fltarr(nSMbins,nrbins)
fieldSMLo=fltarr(nSMbins,nrbins)
fieldSMHi=fltarr(nSMbins,nrbins)

for rr=0,nrbins-1 do begin
   for sm=0,nSMbins-1 do begin
      ; centrals
      if(rr EQ 0) then begin ; only need to do this once
         cen=where(acs.mmgg_scale EQ 1 $
                   AND haloMass GE minGroupMass $
                   AND haloMass LT maxGroupMass $
                   AND acs.kevin_mstar GE smbins[sm] $
                   AND acs.kevin_mstar LT smbins[sm+1] $
                   AND acs.photoz_non_comb GE minz $
                   AND acs.photoz_non_comb LT maxz $
                   , ncen)
         if(ncen GT 5) then begin
            tmp=where(quench[cen] EQ 1,nQuenchCen)
            tmp=where(ell[cen] EQ 1,nEllCen)
            quenchFracCen[sm]=float(nQuenchCen)/ncen
            bootstrap,quench[cen],mean,sigma
            quenchFracCenBoot[sm]=mean
            quenchFracCenErrBoot[sm]=sigma
            ellFracCen[sm]=float(nEllCen)/ncen
            bootstrap,ell[cen],mean,sigma
            ellFracCenBoot[sm]=mean
            ellFracCenErrBoot[sm]=sigma
         endif else begin ; not enough objects, don't plot
            quenchFracCen[sm]=!values.f_nan
            quenchFracCenBoot[sm]=!values.f_nan
            quenchFracCenErrBoot[sm]=!values.f_nan
            ellFracCen[sm]=!values.f_nan
            ellFracCenBoot[sm]=!values.f_nan
            ellFracCenErrBoot[sm]=!values.f_nan
         endelse
      endif

      ; get info for field in this r,SM bin
      field=where(acs.p_mem_best EQ 0 $
                  AND acs.kevin_mstar GE smbins[sm] $
                  AND acs.kevin_mstar LT smbins[sm+1] $
                  AND acs.photoz_non_comb GE minz $
                  AND acs.photoz_non_comb LT maxz $
                  ,nfield)
      if(nfield GT 0) then begin
         tmp=where(quench[field] EQ 1,nQuenchField)
         tmp=where(ell[field] EQ 1,nEllField)
         quenchFracField[sm,rr]=float(nQuenchField)/nfield
         quenchFracFieldErr[sm,rr]=quenchFracField[sm,rr]*sqrt(1./nQuenchField+1./nfield)
         ellFracField[sm,rr]=float(nEllField)/nfield
         ellFracFieldErr[sm,rr]=ellFracField[sm,rr]*sqrt(1./nEllField+1./nfield)
         fieldSM[sm,rr]=alog10(mean(10.^(acs[field].kevin_mstar)))
         fieldSMLo[sm,rr]=fieldSM[sm,rr]-smbins[sm]
         fieldSMHi[sm,rr]=smbins[sm+1]-fieldSM[sm,rr]
      endif

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
                ,nmem)
      if(nmem GT 0) then begin
         tmp=where(quench[mem] EQ 1,nQuenchMem)
         tmp=where(ell[mem] EQ 1,nEllMem)
         quenchFracMem[sm,rr]=float(nQuenchMem)/nmem
;         quenchFracMemErr[sm,rr]=quenchFracMem[sm,rr]*sqrt(1./nQuenchMem+1./nmem)
         quenchFracMemErr[sm,rr]=sqrt(1.*(nQuenchMem+1)*(nmem-nQuenchMem+1)/((nmem+3)*(nmem+2)^2)) ; binomial stddev - see http://www.roma1.infn.it/~dagos/proportions/node3.html
         bootstrap,quench[mem],mean,sigma
         quenchFracMemBoot[sm,rr]=mean
         quenchFracMemErrBoot[sm,rr]=sigma
         ellFracMem[sm,rr]=float(nEllMem)/nmem
;         ellFracMemErr[sm,rr]=ellFracMem[sm,rr]*sqrt(1./nEllMem+1./nmem)
         ellFracMemErr[sm,rr]=sqrt(1.*(nEllMem+1)*(nmem-nEllMem+1)/((nmem+3)*(nmem+2)^2)) ; binomial stddev - see http://www.roma1.infn.it/~dagos/proportions/node3.html
         bootstrap,ell[mem],mean,sigma
         ellFracMemBoot[sm,rr]=mean
         ellFracMemErrBoot[sm,rr]=sigma
         memSM[sm,rr]=alog10(mean(10.^(acs[mem].kevin_mstar)))
         memSMLo[sm,rr]=memSM[sm,rr]-smbins[sm]
         memSMHi[sm,rr]=smbins[sm+1]-memSM[sm,rr]
         memR[sm,rr]=mean(acs[mem].dist_bcg_r200)
         memRLo[sm,rr]=memR[sm,rr]-rbins[rr]
         memRHi[sm,rr]=rbins[rr+1]-memR[sm,rr]
      endif

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

groupsym=4
fieldsym=8
rColors=[!magenta,!darkorange,!darkgreen,!purple]
smColors=[!darkred,!darkgreen,!purple,!orange]
fieldColor=!black

; QUENCHING PLOT
device,filename=plotDir+'quenchr'+plot_suffix+'.eps',/color,/helvetica,xsize=10,ysize=4,/inches,/encapsul
multiplot,/default
multiplot,[nsmbins,1],mxtitle=textoidl('Group-centric Distance (R/R_{200c})'),mxtitsize=1.4,mxtitoffset=0.5,mytitle='Red Fraction',mytitsize=1.4,mytitoffset=0,/rowmajor

for sm=0,nsmbins-1 do begin
   ; top row - binned
   plot,xr,yr,xstyle=1,ystyle=1,/nodata,ytitle=ytitle,xticks=xticks,xtickv=xtickv,xtickname=xticknames,xminor=4
   xyouts,xtext,ytext,smNames[sm],/data,color=smColors[sm]

   plotsym,groupsym,thick=3
   oploterror,memR[sm,*],quenchFracMem[sm,*],memRLo[sm,*],quenchFracMemErrBoot[sm,*],psym=8,color=smColors[sm],linestyle=0,errcolor=smColors[sm],errstyle=0,errthick=3,/lobar
   oploterror,memR[sm,*],quenchFracMem[sm,*],memRHi[sm,*],quenchFracMemErrBoot[sm,*],psym=8,color=smColors[sm],linestyle=0,errcolor=smColors[sm],errstyle=0,errthick=3,/hibar
   ; plot centrals
   plotsym,groupsym,thick=3,/fill
   if(finite(quenchFracCenErrBoot[sm])) then oploterror,[0],[quenchFracCen[sm]],[quenchFracCenErrBoot[sm]],psym=8,color=smColors[sm],errcolor=smColors[sm],errstyle=0,errthick=3
      
   oplot,!x.crange,replicate(quenchFracField[sm,0],2),linestyle=1,color=smColors[sm]

   ; fill in points with SM completeness
   for rr=0,nrbins-1 do begin
      plotsym,groupsym,thick=3,fill=fillArr[sm]
      oplot,[memR[sm,rr]],[quenchFracMem[sm,rr]],psym=8,color=smColors[sm]
   endfor
   
;   if(rr EQ 0) then begin
;      xlegend=!x.crange[0]+0.88*(!x.crange[1]-!x.crange[0])
;      ylegend=!y.crange[0]+0.05*(!y.crange[1]-!y.crange[0])
;      xyouts,xlegend-xsymoffset,ylegend+1*yspace,'Field',alignment=1,charsize=1
;      xyouts,xlegend-xsymoffset,ylegend+0*yspace,'Group Members',alignment=1,charsize=1
;      plotsym,fieldsym,/fill
;      oplot,[xlegend+xsymoffset],[ylegend+1*yspace+ysymoffset],psym=8,color=fieldColor
;      plotsym,groupsym,/fill
;      oplot,[xlegend+xsymoffset],[ylegend+0*yspace+ysymoffset],psym=8,color=smColors[sm]
;   endif

   multiplot
endfor
multiplot,/reset
device,/close


; MORPHOLOGY PLOT
device,filename=plotDir+'morphr'+plot_suffix+'.eps',/color,/helvetica,xsize=10,ysize=4,/inches,/encapsul
multiplot,/default
multiplot,[nsmbins,1],mxtitle=textoidl('Group-centric Distance (R/R_{200c})'),mxtitsize=1.4,mxtitoffset=0.5,mytitle='Morphological Fraction',mytitsize=1.4,mytitoffset=0,/rowmajor

for sm=0,nsmbins-1 do begin
   ; top row - binned
   plot,xr,yr,xstyle=1,ystyle=1,/nodata,ytitle=ytitle,xticks=xticks,xtickv=xtickv,xtickname=xticknames,xminor=4
   xyouts,xtext,ytext,smNames[sm],/data,color=smColors[sm]

   plotsym,groupsym,thick=3
   oploterror,memR[sm,*],ellFracMem[sm,*],memRLo[sm,*],ellFracMemErrBoot[sm,*],psym=8,color=smColors[sm],linestyle=0,errcolor=smColors[sm],errstyle=0,errthick=3,/lobar
   oploterror,memR[sm,*],ellFracMem[sm,*],memRHi[sm,*],ellFracMemErrBoot[sm,*],psym=8,color=smColors[sm],linestyle=0,errcolor=smColors[sm],errstyle=0,errthick=3,/hibar
   ; plot centrals
   plotsym,groupsym,thick=3,/fill
   if(finite(ellFracCenErrBoot[sm])) then oploterror,[0],[ellFracCen[sm]],[ellFracCenErrBoot[sm]],psym=8,color=smColors[sm],errcolor=smColors[sm],errstyle=0,errthick=3
      
   oplot,!x.crange,replicate(ellFracField[sm,0],2),linestyle=1,color=smColors[sm]

   ; fill in points with SM completeness
   for rr=0,nrbins-1 do begin
      plotsym,groupsym,thick=3,fill=fillArr[sm]
      oplot,[memR[sm,rr]],[ellFracMem[sm,rr]],psym=8,color=smColors[sm]
   endfor
   
;   if(rr EQ 0) then begin
;      xlegend=!x.crange[0]+0.88*(!x.crange[1]-!x.crange[0])
;      ylegend=!y.crange[0]+0.05*(!y.crange[1]-!y.crange[0])
;      xyouts,xlegend-xsymoffset,ylegend+1*yspace,'Field',alignment=1,charsize=1
;      xyouts,xlegend-xsymoffset,ylegend+0*yspace,'Group Members',alignment=1,charsize=1
;      plotsym,fieldsym,/fill
;      oplot,[xlegend+xsymoffset],[ylegend+1*yspace+ysymoffset],psym=8,color=fieldColor
;      plotsym,groupsym,/fill
;      oplot,[xlegend+xsymoffset],[ylegend+0*yspace+ysymoffset],psym=8,color=smColors[sm]
;   endif

   multiplot
endfor
multiplot,/reset
device,/close

!p.multi=0

stop

end

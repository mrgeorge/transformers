pro pcorr

; compute partial correlation matrix of a bunch of galaxy properties

inFile="../images/acs_mem.list"
readcol,inFile,skipline=1,id,flag,ra,dec,z,sm,mhalo,r,type,bulge,color,theta,imageFile,format='L,I,D,D,F,F,F,F,I,I,F,F,A'


data=transpose([ $
     [z],$
     [sm],$
     [mhalo],$
     [r],$
     [color],$
     [theta]])

minZ=0.1
maxZ=1.0
minMh=13.5
maxMh=14.0
sel=where(flag EQ 1 $
          AND z GT minZ $
          AND z LT maxZ $
          AND mhalo GT minMh $
          AND mhalo LT maxMh)

data=data[*,sel]

nVar=(size(data,/dim))[0]
nGal=(size(data,/dim))[1]
corrMatrix=fltarr(nVar,nVar)

; compute ranks for each property
rankData=fltarr(nVar,nGal)
for ii=0,nVar-1 do begin
   rankData[ii,*]=mrg_rank(data[ii,*])
endfor

for ii=0,nVar-1 do begin
   for jj=0,nVar-1 do begin
      indepVar=ii
      depVar=jj
      controlVar=where(findgen(nVar) NE indepVar AND findgen(nVar) NE depVar)
      corrMatrix[ii,jj]=p_correlate(reform(rankData[depVar,*]),reform(rankData[indepVar,*]),rankData[controlVar,*])
   endfor
endfor

print, corrMatrix

stop
end

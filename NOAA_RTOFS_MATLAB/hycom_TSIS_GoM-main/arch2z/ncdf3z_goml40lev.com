#!/bin/csh
#
#PBS -N goml004_40lev2
#PBS -o goml004_40lev2.log
#PBS -A ONRDC10855122
##PBS -l mppwidth="1"
#PBS -l walltime="12:00:00"
#PBS -q transfer
#
set echo
#
# --- interpolate to 3-d z-levels from a single HYCOM archive file.
# --- z-levels, via linear interpolation, at Levitus depths.
#
# --- output is netCDF.
# --- this is an example, customize it for your datafile needs.
#
#
# --- optional title and institution.
#
setenv CDF_TITLE	"HYCOM GOMl0.04 FSU COAPS"
setenv CDF_INST 	"40 levels nested within 1/12 Atl."
#
# --- D,y,d select the archive files.
#
setenv X 04.2
setenv E `echo ${X} | awk '{printf("%3.3i",$1*10)}'`
setenv R GOMl0.04
setenv ALL /u/home/wallcraf/hycom/ALL/archive/src
setenv yr 097
setenv S /scr/${user}/hycom/${R}/expt_${X}
setenv D ${S}/netcdf2
setenv DR /u/home/ddukhovs/hycom/${R}/topo
setenv O ${S}/extract${yr}
setenv T 72
setenv Y `echo ${yr} | awk '{printf("%04d",$1+1900)}'`
setenv DRO /u/home/${user}/DATA
setenv hs 0   # output hour stamp

cd ${DRO}

if (! -e ${D}) then
  mkdir -p ${D}
endif
#
if (! -e ${S}/topo) then
  mkdir -p ${S}/topo
endif

if (! -e ${S}/ncout) then
  mkdir -p ${S}/ncout
endif

# Clean netcdf directory
/bin/rm -f ${D}/*gz
#touch regional.depth.a regional.depth.b
#if (-z regional.depth.a) then
if (-e ${S}/topo/depth_${R}_${T}.a && ! -z ${S}/topo/depth_${R}_${T}.a) then
  /bin/rm regional.depth.a
  /bin/ln -sf ${S}/topo/depth_${R}_${T}.a regional.depth.a
  /bin/rm regional.depth.b
  /bin/ln -sf ${S}/topo/depth_${R}_${T}.b regional.depth.b
else
  /usr/bin/rcp newton:${DR}/depth_${R}_${T}.a ${S}/topo/. &
  wait
  /usr/bin/rcp newton:${DR}/depth_${R}_${T}.b ${S}/topo/. &
  wait
  /bin/rm regional.depth.a
  /bin/ln -sf ${S}/topo/depth_${R}_${T}.a regional.depth.a
  /bin/rm regional.depth.b
  /bin/ln -sf ${S}/topo/depth_${R}_${T}.b regional.depth.b
endif
#
touch regional.grid.a regional.grid.b
#if (-z regional.grid.a) then
if (-e ${S}/topo/regional.grid.a && ! -z ${S}/topo/regional.grid.a) then
   /bin/rm regional.grid.a
  /bin/ln -sf ${S}/topo/regional.grid.a regional.grid.a
  /bin/rm regional.grid.b
  /bin/ln -sf ${S}/topo/regional.grid.b regional.grid.b
else
  /usr/bin/rcp newton:${DR}/regional.grid.a ${S}/topo/. &
  wait
  /usr/bin/rcp newton:${DR}/regional.grid.b ${S}/topo/. &
  wait
  /bin/rm regional.grid.a
  /bin/ln -sf ${S}/topo/regional.grid.a regional.grid.a
endif
if (-z regional.grid.b) then
  /bin/rm regional.grid.b
  /bin/ln -sf ${S}/topo/regional.grid.b regional.grid.b
endif


set ds=0
set de=366
@ icnt = 0
set fnmb=0


#foreach d ( a b c d e f g h i j k l )
while ($ds <= $de)
#    setenv O ${D}/tar_${Y}${d}
    setenv d `echo ${ds} | awk '{printf("%03d",$1)}'`
     @ ds++
    setenv hh `echo ${hs} | awk '{printf("%02d",$1)}'`
#    setenv CDF030  ${O}/${E}_archv.${yr}_${d}_00.nc
#  ls ${O}/${E}_archc_*tar.gz >& /dev/null
#  if ( $status == 0) then
#    mv ${O}/${E}_archc_*tar.gz ${O}/../ncout/.
#  endif
#
#  foreach f (${O}/${E}_archv.*.a)
#  setenv FL `basename ${f} .a | sed 's|\(.[0-9]*_[0-9]*\)_\([0-9]*\)|\1|'` 
    set f = ${O}/${E}_archv.${Y}_${d}_${hh}.a
    set FL = ${E}GOMl0.04-${Y}_${d}_${hh}_
    setenv CDF031  ${D}/${FL}u.nc
    setenv CDF032  ${D}/${FL}v.nc
    setenv CDF033  ${D}/${FL}t.nc
    setenv CDF034  ${D}/${FL}s.nc
   echo $CDF031
   setenv FLR ${f} 
   echo ${FLR} 
   ls $CDF031 >& /dev/null
  if ($status == 0) then
    /bin/rm -f  $CDF031  $CDF032  $CDF033  $CDF034    
  endif
#/u/home/wallcraf/hycom/ALL/archive/src/archv2ncdf3z <<E-o-D
${ALL}/archv2ncdf3z <<E-o-D
${FLR}
netCDF
  32 	'iexpt ' = experiment number x10 (000=from archive file)
  3	'yrflag' = days in year flag (0=360J16,1=366J16,2=366J01,3-actual)
 541 	'idm   ' = longitudinal array size
 385	'jdm   ' = latitudinal  array size
  40	'kdm   ' = number of layers
  25.0	'thbase' = reference density (sigma units)
  0	'smooth' = smooth the layered fields (0=F,1=T)
  1  	'iorign' = i-origin of plotted subregion
  1 	'jorign' = j-origin of plotted subregion
  0 	'idmp  ' = i-extent of plotted subregion (<=idm; 0 implies idm)
  0 	'jdmp  ' = j-extent of plotted subregion (<=jdm; 0 implies jdm)
  1	'itype ' = interpolation type (0=sample,1=linear)
  33	'kz    ' = number of depths to sample
  5.0	'z     ' = sample depth  1
  10.0	'z     ' = sample depth  2
  20.0	'z     ' = sample depth  3
  30.0	'z     ' = sample depth  4
  50.0	'z     ' = sample depth  5
  75.0	'z     ' = sample depth  6
 100.0	'z     ' = sample depth  7
 125.0	'z     ' = sample depth  8
 150.0	'z     ' = sample depth  9
 200.0	'z     ' = sample depth 10
 250.0	'z     ' = sample depth 11
 300.0	'z     ' = sample depth 12
 400.0	'z     ' = sample depth 13
 500.0	'z     ' = sample depth 14
 600.0	'z     ' = sample depth 15
 700.0	'z     ' = sample depth 16
 800.0	'z     ' = sample depth 17
 900.0	'z     ' = sample depth 18
1000.0	'z     ' = sample depth 19
1100.0	'z     ' = sample depth 20
1200.0	'z     ' = sample depth 21
1300.0	'z     ' = sample depth 22
1400.0	'z     ' = sample depth 23
1500.0	'z     ' = sample depth 24
1750.0	'z     ' = sample depth 25
2000.0	'z     ' = sample depth 26
2500.0	'z     ' = sample depth 27
3000.0	'z     ' = sample depth 28
3500.0	'z     ' = sample depth 29
4000.0	'z     ' = sample depth 30
4500.0	'z     ' = sample depth 31
5000.0	'z     ' = sample depth 32
5500.0	'z     ' = sample depth 33
  0 	'botio ' = bathymetry  I/O unit (0 no I/O)
  0	'mltio ' = mix.l.thk.  I/O unit (0 no I/O)
  0     'tempml' = temperature jump across mixed-layer (degC,  0 no I/O)
  0     'densml' =   density jump across mixed-layer (kg/m3, 0 no I/O)
  0	'infio ' = intf. depth I/O unit (0 no I/O, <0 label with layer #)
  0	'wviio ' = intf. veloc I/O unit (0 no I/O)
  0	'wvlio ' = w-velocity  I/O unit (0 no I/O)
  31 	'uvlio ' = u-velocity  I/O unit (0 no I/O)
  32 	'vvlio ' = v-velocity  I/O unit (0 no I/O)
  0	'splio ' = speed       I/O unit (0 no I/O)
  33	'temio ' = temperature I/O unit (0 no I/O)
  34	'salio ' = salinity    I/O unit (0 no I/O)
  0	'tthio ' = density     I/O unit (0 no I/O)
  0     'keio  ' = kinetic egy I/O unit (0 no I/O)
E-o-D

  @ icnt++
  if ( ${icnt} == 30 || ${ds} >= ${de} ) then
#    unset $icnt
    @ icnt = 0
    cd ${D}
    @ fnmb++
    set FTR = z3d${E}_${Y}_${fnmb}
    tar -cvf ${FTR}.tar ${E}GOMl0.04-${Y}_*.nc
    ls -l *.tar
    gzip ${FTR}.tar
    wait

    if (! $status) then
      /bin/mv *.gz ${S}/ncout/.
      /bin/rm ${E}*${Y}_*.nc
    endif
    cd ${DRO}
  endif

end    # while  

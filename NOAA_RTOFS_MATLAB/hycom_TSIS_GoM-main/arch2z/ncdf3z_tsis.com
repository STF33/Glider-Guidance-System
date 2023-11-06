#!/bin/csh
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
setenv CDF_TITLE	"HYCOM TSIS IAS0.03  Tendral/FSU COAPS"
setenv CDF_INST 	"30 levels North Atlantic Data Assimilated"
#
# --- D,y,d select the archive files.
#
setenv X 04.2
setenv E `echo ${X} | awk '{printf("%3.3i",$1*10)}'`
setenv R IAS
setenv ALL /Net/ocean/ddmitry/HYCOM/hycom/ALL4/archive/src_2.2.35
#setenv ALL /Net/yucatan/abozec/PACIFIC/HYCOM/hycom/ALL4/archive/src_2.2.35_core
setenv yr 110
setenv Y `echo ${yr} | awk '{printf("%04d",$1+1900)}'`
setenv S /Net/gleam/dmitry/hycom/TSIS/IASx0.03/output
setenv Dtopo /Net/kronos/ddmitry/hycom/TSIS/topo_IAS0.03
setenv Dab   /Net/gleam/dmitry/hycom/TSIS/IASx0.03/output/${Y}
setenv D ${S}/nc3z/${Y}
#setenv DRO ${D}
setenv hs 0   # output hour stamp

mkdir -pv ${D}
#
cd ${S}

if (! -e topo) then
   ln -sf ${Dtopo} topo
endif


# Clean netcdf directory
/bin/rm -f ${D}/*.nc
touch regional.depth.a regional.depth.b
/bin/rm regional.depth.a
/bin/ln -sf ${S}/topo/regional.depth.a
/bin/rm regional.depth.b
/bin/ln -sf ${S}/topo/regional.depth.b

touch regional.grid.a regional.grid.b
/bin/ln -sf ${Dtopo}/regional.grid.a
/bin/ln -sf ${Dtopo}/regional.grid.b

set ds=0
set de=365
@ icnt = 0
set fnmb=0


#foreach d ( a b c d e f g h i j k l )
while ($ds <= $de)
    setenv d `echo ${ds} | awk '{printf("%03d",$1)}'`
     @ ds++
    setenv hh `echo ${hs} | awk '{printf("%02d",$1)}'`

    set FLR = ${Dab}/archv.${Y}_${d}_${hh}.a
    set FL = TSIS_IAS0.03-${Y}_${d}_${hh}_
    setenv CDF031  ${D}/${FL}u.nc
    setenv CDF032  ${D}/${FL}v.nc
    setenv CDF033  ${D}/${FL}t.nc
    setenv CDF034  ${D}/${FL}s.nc
   echo $CDF031
   echo ${FLR} 
   ls $CDF031 >& /dev/null
  if ($status == 0) then
    /bin/rm -f  $CDF031  $CDF032  $CDF033  $CDF034    
  endif
#/u/home/wallcraf/hycom/ALL/archive/src/archv2ncdf3z <<E-o-D
${ALL}/archv2ncdf3z <<E-o-D
${FLR}
netCDF
  22 	'iexpt ' = experiment number x10 (000=from archive file)
  3	'yrflag' = days in year flag (0=360J16,1=366J16,2=366J01,3-actual)
 1401 	'idm   ' = longitudinal array size
 891	'jdm   ' = latitudinal  array size
  30	'kdm   ' = number of layers
  34.0	'thbase' = reference density (sigma units)
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

#  /bin/mv ${FL}

end    # while  

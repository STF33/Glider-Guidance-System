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
# --- optional title and institution.
#
cd /Net/gleam/dmitry/hycom/TSIS/IASx0.03/output
touch regional.depth.a regional.depth.b
if (-z regional.depth.a) then
  /bin/rm regional.depth.a
  /bin/ln -s /Net/kronos/ddmitry/hycom/TSIS/topo_IAS0.03/regional.depth.a .
endif
if (-z regional.depth.b) then
  /bin/rm regional.depth.b
  /bin/ln -s /Net/kronos/ddmitry/hycom/TSIS/topo_IAS0.03/regional.depth.b .
endif
#
touch regional.grid.a regional.grid.b
if (-z regional.grid.a) then
  /bin/rm regional.grid.a
  /bin/ln -s /Net/kronos/ddmitry/hycom/TSIS/topo_IAS0.03/regional.grid.a regional.grid.a
endif
if (-z regional.grid.b) then
  /bin/rm regional.grid.b
  /bin/ln -s /Net/kronos/ddmitry/hycom/TSIS/topo_IAS0.03/regional.grid.b regional.grid.b
endif
#
# --- D,y,d select the archive files.
setenv CDF_TITLE        "HYCOM ARCTIC"
setenv CDF_INST         "Naval Research Laboratory"
#
setenv D /Net/gleam/dmitry/hycom/TSIS/IASx0.03/output/2009
#setenv O /Net/gleam/dmitry/hycom/TSIS/IASx0.03/output/nc3z/2009
setenv O /Net/ocean/ddmitry
setenv E 500 
setenv Y 2009
setenv m 1
#1948 setenv d 10
#
    foreach d ( 003 )
    setenv d1 `echo $d | awk '{printf("%03d", $1)}'`
   setenv CDF031  ${O}/archv.${Y}_${d1}_temp.nc
   setenv CDF032  ${O}/archv.${Y}_${d1}_sal.nc
   setenv CDF033  ${O}/archv.${Y}_${d1}_dens.nc
   /bin/rm $CDF031  $CDF032 $CDF033

/Net/yucatan/abozec/PACIFIC/HYCOM/hycom/ALL4/archive/src/archv2ncdf3z <<E-o-D
${D}/archv.${Y}_${d}_00.a
netCDF
 000	'iexpt ' = experiment number x10 (000=from archive file)
   0	'yrflag' = days in year flag (0=360J16,1=366J16,2=366J01,3-actual)
 1401    'idm   ' = longitudinal array size
 891	'jdm   ' = latitudinal  array size
  30	'kdm   ' = number of layers
  34.0	'thbase' = reference density (sigma units)
   0	'smooth' = smooth the layered fields (0=F,1=T)
   1	'iorign' = i-origin of plotted subregion
   1	'jorign' = j-origin of plotted subregion
   0	'idmp  ' = i-extent of plotted subregion (<=idm; 0 implies idm)
   0	'jdmp  ' = j-extent of plotted subregion (<=jdm; 0 implies jdm)
   1	'itype ' = interpolation type (0=sample,1=linear)
  33    'kz    ' = number of depths to sample
  5.0   'z     ' = sample depth  1
  10.0  'z     ' = sample depth  2
  20.0  'z     ' = sample depth  3
  30.0  'z     ' = sample depth  4
  50.0  'z     ' = sample depth  5
  75.0  'z     ' = sample depth  6
 100.0  'z     ' = sample depth  7
 125.0  'z     ' = sample depth  8
 150.0  'z     ' = sample depth  9
 200.0  'z     ' = sample depth 10
 250.0  'z     ' = sample depth 11
 300.0  'z     ' = sample depth 12
 400.0  'z     ' = sample depth 13
 500.0  'z     ' = sample depth 14
 600.0  'z     ' = sample depth 15
 700.0  'z     ' = sample depth 16
 800.0  'z     ' = sample depth 17
 900.0  'z     ' = sample depth 18
1000.0  'z     ' = sample depth 19
1100.0  'z     ' = sample depth 20
1200.0  'z     ' = sample depth 21
1300.0  'z     ' = sample depth 22
1400.0  'z     ' = sample depth 23
1500.0  'z     ' = sample depth 24
1750.0  'z     ' = sample depth 25
2000.0  'z     ' = sample depth 26
2500.0  'z     ' = sample depth 27
3000.0  'z     ' = sample depth 28
3500.0  'z     ' = sample depth 29
4000.0  'z     ' = sample depth 30
4500.0  'z     ' = sample depth 31
5000.0  'z     ' = sample depth 32
5500.0  'z     ' = sample depth 33
  0	'botio ' = bathymetry  I/O unit (0 no I/O)
  0	'mltio ' = mix.l.thk.  I/O unit (0 no I/O)
  0     'tempml' = temperature jump across mixed-layer (degC,  0 no I/O)
  0     'densml' =   density jump across mixed-layer (kg/m3, 0 no I/O)
  0	'infio ' = intf. depth I/O unit (0 no I/O, <0 label with layer #)
  0	'wviio ' = intf. veloc I/O unit (0 no I/O)
  0	'wvlio ' = w-velocity  I/O unit (0 no I/O)
  0	'uvlio ' = u-velocity  I/O unit (0 no I/O)
  0	'vvlio ' = v-velocity  I/O unit (0 no I/O)
  0	'splio ' = speed       I/O unit (0 no I/O)
  31	'temio ' = temperature I/O unit (0 no I/O)
  32	'salio ' = salinity    I/O unit (0 no I/O)
  33	'tthio ' = density     I/O unit (0 no I/O)
  0	'keio  ' = kinetic egy I/O unit (0 no I/O)
E-o-D

echo Salut
    set m=`expr $m + 1`
end


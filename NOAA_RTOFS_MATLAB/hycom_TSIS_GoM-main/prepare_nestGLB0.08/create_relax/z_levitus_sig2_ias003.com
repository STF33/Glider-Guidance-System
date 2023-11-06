#!/bin/csh -x
# This script prepares 
# PHC climatology T&S& dens fields in Z levels
# interpolated into HYCOM grid
#>& /dev/null
if (! $status) then 
  set path = (${path} /Net/ocean/ddmitry/HYCOM/hycom/ALL4/bin)  
  set path = (${path} /Net/ocean/ddmitry/HYCOM/hycom/ALL4/bin/C)
endif
setenv pget cp
setenv pput cp
#
#--- P is primary path,
#--- S is scratch directory,
#--- D is scalar permanent directory,
#--- L is Levitus directory
#--- T topography version: 07, 09, or 11
#setenv P /Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.04/relax/PHC3.0
#setenv DCL /Net/data/PHC_3.0

#
setenv P /Net/kronos/ddmitry/hycom/TSIS/relax_41lrs
setenv S $P/SCRATCH
setenv D $P
setenv DT /Net/kronos/ddmitry/hycom/TSIS/topo_IAS0.03
setenv ALL /nexsan/people/ddmitry/Net_ocean/HYCOM/hycom/ALL4
setenv L /Net/data/PHC_3.0

#
mkdir -pv $S
cd       $S
#
touch   z_levitus
/bin/rm z_levitus
#
#--- 12 months
#
foreach MM ( 01 02 03 04 05 06 07 08 09 10 11 12 )
#
#--- Input.
#
touch      fort.71 fort.73 fort.72
/bin/rm -f fort.71 fort.73 fort.72
${pget} $L/r_m${MM}.d fort.71 &
${pget} $L/t_m${MM}.d fort.72 &
${pget} $L/s_m${MM}.d fort.73 &
#
#touch regional.grid.a regional.grid.b
touch regional.grid.a regional.grid.b
if (-z regional.grid.b) then
    ${pget} ${DT}/regional.grid.b regional.grid.b &
endif
if (-z regional.grid.a) then
    ${pget} ${DT}/regional.grid.a regional.grid.a &
endif
#
touch z_levitus
if (-z z_levitus) then
  ${pget} ${ALL}/relax/src_2.2.35/z_levitus . &
endif
wait
chmod a+rx z_levitus
#
/bin/rm -f core
touch core

setenv FOR010A fort.10A
setenv FOR011A fort.11A
setenv FOR012A fort.12A
/bin/rm -f fort.10 fort.10A fort.11 fort.11A fort.12 fort.12A
./z_levitus <<E-o-D
 &AFTITL
  CTITLE = '1234567890123456789012345678901234567890',
  CTITLE = 'Levitus monthly',
 /
 &AFFLAG
  ICTYPE =   3,
  KSIGMA =   2,
  INTERP =   0,
  ITEST  =   7,
  JTEST  =  29,
  MONTH  = $MM,
 /
E-o-D
#
#--- Required Output, potential density and temperature.
#
${pput} fort.10  ${D}/temp_sig2_m${MM}.b
${pput} fort.10A ${D}/temp_sig2_m${MM}.a
${pput} fort.12  ${D}/dens_sig2_m${MM}.b
${pput} fort.12A ${D}/dens_sig2_m${MM}.a
#
#--- Optional Output.
#
${pput} fort.11  ${D}/saln_sig2_m${MM}.b
${pput} fort.11A ${D}/saln_sig2_m${MM}.a
#
#--- end of month foreach loop
#
/bin/rm fort.7[123]
end
#
#--- Delete all files.
#
#/bin/rm -f *

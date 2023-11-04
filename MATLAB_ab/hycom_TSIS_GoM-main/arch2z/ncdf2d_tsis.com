#!/bin/csh
#
#PBS -N GOMl_2D_20lev
#PBS -o GOMl_2D_20lev.log
#PBS -A ONRDC10855122
##PBS -l mppwidth="1"
#PBS -l walltime="12:00:00"
#PBS -q transfer
#
set echo
#
# --- extract 2-d fields from a single HYCOM archive file.
#
# --- output can be formatted, unformatted (BINARY), .[ab] (HYCOM).
# --- or use ${E}archv2ncdf2d for netCDF output.
#
# --- output is HYCOM .a files, converted to "raw" .A files.
# --- this is an example, customize it for your datafile needs.
#
# --- D,y,d select the archive files.
setenv CDF_TITLE	"HYCOM TSIS IAS0.03  Tendral/FSU COAPS"
setenv CDF_INST 	"30 levels North Atlantic Data Assimilated"
#
#setenv y 109
setenv X 04.2
setenv E `echo ${X} | awk '{printf("%3.3i",$1*10)}'`
setenv R IAS
setenv ALL /Net/ocean/ddmitry/HYCOM/hycom/ALL4/archive/src_2.2.35
setenv yr 110
setenv Y `echo ${yr} | awk '{printf("%04d",$1+1900)}'`
setenv S /Net/gleam/dmitry/hycom/TSIS/IASx0.03/output
setenv Dtopo /Net/kronos/ddmitry/hycom/TSIS/topo_IAS0.03
setenv Dab   /Net/gleam/dmitry/hycom/TSIS/IASx0.03/output/${Y}
setenv D ${S}/nc2d/${Y}
#setenv DR /u/home/${user}/hycom/${R}/topo
#setenv O ${S}/extract${y}
#setenv T 72 
#setenv Y `echo ${y} | awk '{printf("%04d",$1+1900)}'`
#setenv DR0 /p/home/${user}/data_gom
setenv hs 0   # output hour stamp

mkdir -pv ${D}

#
cd ${S}

if (! -e topo) then
   ln -sf ${Dtopo} topo
endif

# Clean netcdf directory
/bin/rm -f ${D}/*gz
touch regional.depth.a regional.depth.b
/bin/rm regional.depth.a
/bin/ln -sf ${S}/topo/regional.depth.a .
/bin/rm regional.depth.b
/bin/ln -sf ${S}/topo/regional.depth.b .

touch regional.grid.a regional.grid.b
/bin/ln -sf ${Dtopo}/regional.grid.a .
/bin/ln -sf ${Dtopo}/regional.grid.b .


set ds=0
set de=365
set icnt=0
set fnmb=0
#set f1= `echo ${ds} | awk '{printf("%03d",$1+1)}'`

#
#foreach d ( a b c d e f g h i j k l )
while ( ${ds} <= ${de} )
  @ ds++
  setenv d `echo ${ds} | awk '{printf("%03d",$1)}'`
  setenv hh `echo ${hs} | awk '{printf("%02d",$1)}'`

  set FLR = ${Dab}/archv.${Y}_${d}_${hh}.a
  set FL = TSIS_IAS0.03-${Y}_${d}_${hh}_
  setenv CDF021 ${D}/${FL}ssh.nc
#  setenv CDF022 ${D}/${FL}mlu.nc
#  setenv CDF023 ${D}/${FL}mlv.nc
#  setenv CDF024 ${D}/${FL}mlt.nc
#  setenv CDF025 ${D}/${FL}mls.nc
  echo $CDF021 
  echo ${FLR}
  ls $CDF021 >& /dev/null
  if ($status == 0) then
    /bin/rm $CDF021  
  endif

  if ( -e ${FLR}) then

${ALL}/archv2ncdf2d << E-o-D
${FLR}
netCDF
 23	'iexpt ' = experiment number x10 (000=from archive file)
  3	'yrflag' = days in year flag (0=360J16,1=366J16,2=366J01,3-actual)
 1401	'idm   ' = longitudinal array size
 891	'jdm   ' = latitudinal  array size
 30	'kdm   ' = number of layers
 34.0	'thbase' = reference density (sigma units)
  0	'smooth' = smooth fields before plotting (0=F,1=T)
  0     'mthin ' = mask thin layers from plots (0=F, 1=T)
  1	'iorign' = i-origin of plotted subregion
  1	'jorign' = j-origin of plotted subregion
  0	'idmp  ' = i-extent of plotted subregion (<=idm; 0 implies idm)
  0	'jdmp  ' = j-extent of plotted subregion (<=jdm; 0 implies jdm)
  0	'botio ' = bathymetry       I/O unit (0 no I/O)
  0	'flxio ' = surf. heat flux  I/O unit (0 no I/O)
  0	'empio ' = surf. evap-pcip  I/O unit (0 no I/O)
  0	'ttrio ' = surf. temp trend I/O unit (0 no I/O)
  0	'strio ' = surf. saln trend I/O unit (0 no I/O)
  0	'icvio ' = ice coverage     I/O unit (0 no I/O)
  0	'ithio ' = ice thickness    I/O unit (0 no I/O)
  0	'ictio ' = ice temperature  I/O unit (0 no I/O)
 21	'sshio ' = sea surf. height I/O unit (0 no I/O)
  0	'bsfio ' = baro. strmfn.    I/O unit (0 no I/O)
  0 	'uvmio ' = mix. lay. u-vel. I/O unit (0 no I/O)
  0	'vvmio ' = mix. lay. v-vel. I/O unit (0 no I/O)
  0	'spmio ' = mix. lay. speed  I/O unit (0 no I/O)
  0	'bltio ' = bnd. lay. thick. I/O unit (0 no I/O)
  0	'mltio ' = mix. lay. thick. I/O unit (0 no I/O)
  0	'sstio ' = mix. lay. temp.  I/O unit (0 no I/O)
  0 	'sssio ' = mix. lay. saln.  I/O unit (0 no I/O)
  0	'ssdio ' = mix. lay. dens.  I/O unit (0 no I/O)
 -1	'kf    ' = first output layer (=0 end output; <0 label with layer #)
  1	'kl    ' = last  output layer
  0	'uvlio ' = layer k   u-vel. I/O unit (0 no I/O)
  0	'vvlio ' = layer k   v-vel. I/O unit (0 no I/O)
  0	'splio ' = layer k   speed. I/O unit (0 no I/O)
  0	'wvlio ' = layer k   i-vel. I/O unit (0 no I/O)
  0	'infio ' = layer k   i.dep. I/O unit (0 no I/O)
  0	'thkio ' = layer k   thick. I/O unit (0 no I/O)
  0	'temio ' = layer k   temp   I/O unit (0 no I/O)
  0	'salio ' = layer k   saln.  I/O unit (0 no I/O)
  0	'tthio ' = layer k   dens,  I/O unit (0 no I/O)
  0	'sfnio ' = layer k  strmfn. I/O unit (0 no I/O)
  0	'kf    ' = first output layer (=0 end output; <0 label with layer #)
E-o-D

  @ icnt++
  endif # if exist $FLR
end     # while ds


exit 0

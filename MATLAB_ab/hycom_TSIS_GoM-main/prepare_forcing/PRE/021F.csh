#!/bin/csh -x
#PBS -N XXX
#PBS -j oe
#PBS -o XXX.log
#PBS -W umask=027
# single node
#PBS -l select=1:ncpus=24
#PBS -l walltime=04:00:00
#PBS -l application=hycom
#PBS -A NRLSS03755018
#PBS -q standard
#
set echo
set time = 1
set timestamp
date +"START  %c"
#
# --- Create model interpolated fluxes for HYCOM.
# --- From NWP NRL .nc files.
#
# --- prebuild this script similar to a model script
# --- awk -f 284.awk y01=098 ab=a 284F.csh > 284f098a.csh
#
# --- EX is experiment directory
#
setenv EX /home/ddmitry/codes/anls_mtlb_utils/hycom_TSIS/prepare_forcing 
#
# --- E is expt number.
# --- P is primary path.
# --- D is permanent directory.
# --- S is scratch   directory, must not be the permanent directory.
# --- N is data-set name, e.g. cfsr.
# --- W is permanent NWP NRL .nc files directory.
#
source ${EX}/EXPT.src
#
# --- HUM is vapmix or spchum
#
setenv HUM vapmix
#setenv HUM spchum
#
# --- System Type
#
switch ($OS)
case 'SunOS':
#   assumes /usr/5bin is before /bin and /usr/bin in PATH.
    setenv BINRUN   ""
    breaksw
case 'Linux':
case 'OSF1':
case 'IRIX64':
case 'AIX':
    setenv BINRUN   ""
    breaksw
case 'HPE':
case 'HPEI':
#   load modules
    unset echo
    source ~/HYCOM-tools/Make_ncdf.src
    set echo
    setenv BINRUN   ""
    breaksw
case 'XC30':
case 'XC40':
# --- ~/HYCOM-tools/bin is for CNL (ftn)
    setenv BINRUN   "aprun -n 1"
default:
    echo 'Unknown Operating System: ' $OS
    exit (1)
endsw
#
# --- pget, pput "copy" files between scratch and permanent storage.
# --- Can both be cp if the permanent filesystem is mounted locally.
#
switch ($OS)
case 'SunOS':
case 'Linux':
case 'HPE':
case 'HPEI':
case 'XT4':
case 'OSF1':
case 'AIX':
case 'unicos':
case 'unicosmk':
    if (-e        ~${user}/bin/pget) then
      setenv pget ~${user}/bin/pget
      setenv pput ~${user}/bin/pput
    else
      setenv pget /bin/cp
      setenv pput /bin/cp
    endif
    breaksw
default:
      setenv pget /bin/cp
      setenv pput /bin/cp
endsw
#
switch ($OS)
case 'XT4':
case 'XC30':
case 'XC40':
    mkdir -p      $S
    lfs setstripe $S -S 1048576 -i -1 -c 8
    breaksw
endsw
#
mkdir -p $S/flux
cd       $S/flux
#
# --- For whole year runs.
# ---   ymx number of years per model run.
# ---   Y01 initial model year of this run.
# ---   YXX is the last model year of this run, and the first of the next run.
# ---   A and B are identical, typically blank.
# --- For part year runs.
# ---   A is this part of the year, B is next part of the year.
# ---   Y01 initial model year of this run.
# ---   YXX year at end of this part year run.
# ---   ymx is 1.
# --- Note that these variables and the .awk generating script must
# ---  be consistant to get a good run script.
#
# --- One year spin-up run.
#
@ ymx =  1
#
setenv A ""
setenv B ""
setenv Y01 "101"
#
switch ("${B}")
case "${A}":
    setenv YXX `echo $Y01 $ymx | awk '{printf("%03d", $1+$2)}'`
    breaksw
case "a":
    setenv YXX `echo $Y01 | awk '{printf("%03d", $1+1)}'`
    breaksw
default:
    setenv YXX $Y01
endsw
#
echo "Y01 =" $Y01 "YXX = " $YXX  "A =" ${A} "B =" ${B}
#
# --- time limits.
#
if (-e ${D}/../${E}y${Y01}${A}.limits) then
  setenv TS `sed -e "s/-/ /g" ${D}/../${E}y${Y01}${A}.limits | awk '{print $1}'`
  setenv TM `cat              ${D}/../${E}y${Y01}${A}.limits | awk '{print $2}'`
else
# use "LIMITI" when starting a run after day zero.
# use "LIMITS9" (say) for a 9-day run.
  setenv TS `echo "LIMITS" | awk -f ${S}/021.awk y01=${Y01} ab=${A} | awk '{print $1}'`
  setenv TM `echo "LIMITS" | awk -f ${S}/021.awk y01=${Y01} ab=${A} | awk '{print $2}'`
endif
#
echo "TS =" $TS "TM =" $TM
#
# --- input files from file server.
#
touch  regional.grid.a regional.grid.b
if (-z regional.grid.b) then
  ${pget} ${DS}/topo/regional.grid.b regional.grid.b &
endif
if (-z regional.grid.a) then
  ${pget} ${DS}/topo/regional.grid.a regional.grid.a &
endif
#
setenv YC1 `echo $Y01 | awk '{printf("%04d\n",$1+1900)}'`
setenv YCX `echo $YXX | awk '{printf("%04d\n",$1+1900)}'`
setenv CDF071 ${N}-sea_${YC1}_01hr_TaqaQrQp.nc
setenv CDF072 ${N}-sea_${YCX}_01hr_TaqaQrQp.nc
#
#touch  $CDF071
#if (-z $CDF071) then
	#  ${pget} ${W}/$CDF071 . &
ln -s ${W}/$CDF071 . 
  #endif
#
#if ($CDF071 != $CDF072) then
	#  touch  $CDF072
	#  if (-z $CDF072) then
       #${pget} ${W}/$CDF072 . &
ln -s ${W}/$CDF072 .
    #  endif
    #endif
#
# --- executable
#
/bin/cp /nexsan/people/abozec/HYCOM-tools/force/src/ap_nc . &
#/bin/cp /nexsan/people/ddmitry/Net_ocean/HYCOM/HYCOM-tools/force/src/ap_nc . &
wait
chmod ug+rx ap_nc
ls -laFq
#
# --- NAMELIST input.
#
touch   fort.05i
/bin/rm fort.05i
cat <<E-o-D  > fort.05i
 &AFTITL
  CTITLE = '1234567890123456789012345678901234567890',
  CTITLE = '${N}-sea, 1hrly, MKS',
 /
 &AFTIME
  HMKS   =   1.0,          !kg/kg             to kg/kg
  RMKS   =   1.0,          !W/m**2 into ocean to W/m**2 into ocean
  PMKS   =   1.1574074E-5, !m/day  into ocean to m/s    into ocean
  BIASPC =   0.0,
  PCMEAN =   0.0,
  BIASRD =   0.0,
  RDMEAN =   0.0,
  FSTART = ${TS},
  TSTART = ${TS},
  TMAX   = ${TM},
 /
 &AFFLAG
  IFFILE =     5,  !3:monthly-climo; 5:actual-day;
  IFTYPE =     4,  !5:Ta-Ha-Qr-Qp-Pc; 4:Ta-Ha-Qr-Qp; 3:Ta-Ha; 2:Qr; 1:Pc;
  INTERP =     0,  !0:bilinear; 1:cubic spline; 3:piecewise bi-cubic;
  INTMSK =     0,  !0:no mask; 1:land/sea=0/1; 2:land/sea=1/0;
 /
E-o-D
#
# --- run the flux interpolation.
#
date +"lwdflx %c"
touch      fort.10     fort.10a
/bin/rm -f fort.1[0-4] fort.1[0-4]a
#
setenv FOR010A fort.10a
setenv FOR011A fort.11a
setenv FOR012A fort.12a
setenv FOR013A fort.13a
setenv FOR014A fort.14a
#
/bin/rm -f core
touch      core
${BINRUN} ./ap_nc < fort.05i
#
# --- Output.
#
/bin/mv fort.10  airtmp_${Y01}${A}.b
/bin/mv fort.10a airtmp_${Y01}${A}.a
if ($HUM == "vapmix") then
  /bin/mv fort.11  vapmix_${Y01}${A}.b
  /bin/mv fort.11a vapmix_${Y01}${A}.a
else
  sed -e "s/vapmix/spchum/g" fort.11 >! spchum_${Y01}${A}.b
  /bin/mv fort.11a                      spchum_${Y01}${A}.a
endif
touch     fort.12
if (! -z  fort.12) then
  /bin/mv fort.12  radflx_${Y01}${A}.b
  /bin/mv fort.12a radflx_${Y01}${A}.a
  /bin/mv fort.13  shwflx_${Y01}${A}.b
  /bin/mv fort.13a shwflx_${Y01}${A}.a
endif
touch     fort.14
if (! -z  fort.14) then
  /bin/mv fort.14  precip_${Y01}${A}.b
  /bin/mv fort.14a precip_${Y01}${A}.a
endif
touch      fort.10     fort.10a
/bin/rm -f fort.1[0-4] fort.1[0-4]a
#
if (-e ./SAVE) then
  ln ??????_${Y01}${A}.a ./SAVE
  ln ??????_${Y01}${A}.b ./SAVE
endif
if (-e ../cice) then
  date +"CICE   %c"
  setenv IDM  `grep idm regional.grid.b | awk '{print $1}'`
  setenv JDM  `grep jdm regional.grid.b | awk '{print $1}'`
  setenv MAP  `grep map regional.grid.b | awk '{print $1}'`
  if ($MAP < 10) then
    setenv JDA $JDM
  else
#   global tripole region
    setenv JDA `expr $JDM - 1`
  endif
  cd ../cice
  /bin/rm -f vapmix_${Y01}${A}.r airtmp_${Y01}${A}.a airtmp_${Y01}${A}.r netrad_${Y01}${A}.r netQlw_${Y01}${A}.a
# airtmp is in degK, used immediately so not run in background
  ${BINRUN} ~/HYCOM-tools/bin/hycom_expr ../flux/airtmp_${Y01}${A}.a ONE ${IDM} ${JDM} 1.0 273.15 airtmp_${Y01}${A}.a >! airtmp_${Y01}${A}.b
#
  ${BINRUN} ~/HYCOM-tools/bin/hycom2raw8         airtmp_${Y01}${A}.a ${IDM} ${JDM} 1 1 ${IDM} ${JDA} airtmp_${Y01}${A}.r >! airtmp_${Y01}${A}.B &
  ${BINRUN} ~/HYCOM-tools/bin/hycom2raw8 ../flux/${HUM}_${Y01}${A}.a ${IDM} ${JDM} 1 1 ${IDM} ${JDA} vapmix_${Y01}${A}.r >! vapmix_${Y01}${A}.B &
  if (! -e ../lrad) then
#   netrad may not be needed, create it anyway
    ${BINRUN} ~/HYCOM-tools/bin/hycom2raw8 ../flux/shwflx_${Y01}${A}.a ${IDM} ${JDM} 1 1 ${IDM} ${JDA} netrad_${Y01}${A}.r >! netrad_${Y01}${A}.B &
#   incoming Qlw from netQlw and Tsur, just calculate Qlw here
    ${BINRUN} ~/HYCOM-tools/bin/hycom_expr ../flux/radflx_${Y01}${A}.a ../flux/shwflx_${Y01}${A}.a ${IDM} ${JDM} 1.0 -1.0 netQlw_${Y01}${A}.a >! netQlw_${Y01}${A}.b &
  endif
  wait
  if (-e ./SAVE) then
    foreach f ( vapmix_${Y01}${A}.[rB] airtmp_${Y01}${A}.[rB] )
      ln ${f} ./SAVE/${f}
    end
  endif
#cice
endif
#
#  --- END OF JUST IN TIME FLUXES GENERATION SCRIPT.
#
date +"END    %c"

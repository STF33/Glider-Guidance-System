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
# --- Create model interpolated 10m winds for HYCOM. Use bi-cubic interpolation.
# --- From NWP NRL .nc files.
#
# --- prebuild this script similar to a model script
# --- awk -f 284.awk y01=098 ab=a 284V.csh > 284v098a.csh
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
mkdir -p $S/wvel
cd       $S/wvel
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
touch      fort.44 fort.44a fort.45 fort.45a
#/bin/rm -f fort.44 fort.44a fort.45 fort.45a
if (-z fort.44) then
  ${pget} ${DS}/force_otf/offset/tauewd_zero.b fort.44  &
endif
if (-z fort.44a) then
  ${pget} ${DS}/force_otf/offset/tauewd_zero.a fort.44a &
endif
if (-z fort.45) then
  ${pget} ${DS}/force_otf/offset/taunwd_zero.b fort.45  &
endif
if (-z fort.45a) then
  ${pget} ${DS}/force_otf/offset/taunwd_zero.a fort.45a &
endif
#
setenv YC1 `echo $Y01 | awk '{printf("%04d\n",$1+1900)}'`
setenv YCX `echo $YXX | awk '{printf("%04d\n",$1+1900)}'`
setenv CDF071 ${N}-sec2_${YC1}_01hr_uv-10m.nc
setenv CDF072 ${N}-sec2_${YCX}_01hr_uv-10m.nc
#
#touch  $CDF071
#if (-z $CDF071) then
#  ${pget} ${W}/$CDF071 . &
#endif
#C
#if ($CDF071 != $CDF072) then
#  touch  $CDF072
#  if (-z $CDF072) then
#    ${pget} ${W}/$CDF072 . &
#  endif
#endif
ln -s ${W}/$CDF071 .
ln -s ${W}/$CDF072 .
#
# --- executable
#
/bin/cp /nexsan/people/abozec/HYCOM-tools/force/src/kp_nc . &
#/bin/cp /nexsan/people/ddmitry/Net_ocean/HYCOM/HYCOM-tools/force/src/wi_nc . &
wait
chmod ug+rx wi_nc
ls -laFq
#
# --- NAMELIST input.
#
touch   fort.05i
/bin/rm fort.05i
cat <<E-o-D  > fort.05i
 &WWTITL
  CTITLE = '123456789012345678901234567890123456789012345678901234567890',
  CTITLE = '${N}, 1-hrly, Scat Corrected, MKS',
  CNAME  = 'wndewd', 'wndnwd',
 /
 &WWTIME
  SPDMIN =   0.0,  !minimum allowed wind speed
  WSCALE =   1.0,  !scale factor to mks
  WSTART = ${TS},
  TSTART = ${TS},
  TMAX   = ${TM},
 /
 &WWFLAG
  IGRID  =   2,  !0=p; 1=u&v; 2=p
  ISPEED =   0,  !0:none; 1:const; 2:kara; 3:coare; -1,-2:input wind velocity
  INTERP =   3,  !0:bilinear; 1:cubic spline; 2:piecewise bessel; 3:piecewise bi-cubic;
  INTMSK =   0,  !0:no mask; 1:land/sea=0/1; 2:l/s=1/0;
  IFILL  =   3,  !0,1:tx&ty; 2,3:magnitude; 1,3:smooth; (intmsk>0 only)
  IOFILE =   0,  !0:single offset; 1:multiple offsets; 2:multi-off no .b check
  IWFILE =   4,  !1:ann/mon; 2:multi-file; 4:actual wind day
 /
E-o-D
#
# --- run the wnd10m interpolation.
#
date +"wnd10m %c"
touch      fort.10     fort.10a
/bin/rm -f fort.1[012] fort.1[012]a
#
setenv FOR010A fort.10a
setenv FOR011A fort.11a
setenv FOR012A fort.12a
#
setenv FOR044A fort.44a
setenv FOR045A fort.45a
#
/bin/rm -f core
touch      core
${BINRUN} ./wi_nc < fort.05i
#
# --- Output.
#
/bin/mv fort.10  wndewd_${Y01}${A}.b
/bin/mv fort.10a wndewd_${Y01}${A}.a
/bin/mv fort.11  wndnwd_${Y01}${A}.b
/bin/mv fort.11a wndnwd_${Y01}${A}.a
if (-e fort.12) then
  /bin/mv fort.12  wndspd_${Y01}${A}.b
  /bin/mv fort.12a wndspd_${Y01}${A}.a
endif
#
if (-e ./SAVE) then
  ln wnd???_${Y01}${A}.a ./SAVE
  ln wnd???_${Y01}${A}.b ./SAVE
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
  /bin/rm wndewd_${Y01}${A}.r wndnwd_${Y01}${A}.r
  ${BINRUN} ~/HYCOM-tools/bin/hycom2raw8 ../wvel/wndewd_${Y01}${A}.a ${IDM} ${JDM} 1 1 ${IDM} ${JDA} wndewd_${Y01}${A}.r >! wndewd_${Y01}${A}.B &
  ${BINRUN} ~/HYCOM-tools/bin/hycom2raw8 ../wvel/wndnwd_${Y01}${A}.a ${IDM} ${JDM} 1 1 ${IDM} ${JDA} wndnwd_${Y01}${A}.r >! wndnwd_${Y01}${A}.B &
  wait
  if (-e ./SAVE) then
    foreach f ( wndewd_${Y01}${A}.[rB] wndnwd_${Y01}${A}.[rB] )
      ln ${f} ./SAVE/${f}
    end
  endif
#cice
endif
#
#  --- END OF JUST IN TIME WND10M GENERATION SCRIPT.
#
date +"END    %c"

#! /bin/csh -x
#PBS -N Rnest_GOMb0.08
#PBS -j oe
#PBS -o Rnest_GOMb0.08.log
#PBS -W umask=027 
#PBS -l select=1:ncpus=48
#PBS -l place=scatter:excl
#PBS -l walltime=1:30:00
#PBS -A ONRDC10855122
#PBS -q standard
#
module list
set echo
set time = 1
hostname
date
#
## --- form interpolated subregion daily mean archive files, GLBb0.08 to IAS0.03
#
# --- R  is the original region
# --- U  is the target   region
# --- E  is the experiment number
# --- Y  is year
# --- P  is month
# --- S  is the location to run from
# --- D  is the location of the outer domain   archive files
# --- N  is the location of the subregion (nested)  archive files
# --- TR is the location of the outer domain topo 
# --- TU is the location of the subregion topo    files
# --- VR is the original  depth version (${TR}/depth_${R}_${VR}.[ab])
# --- VU is the subregion depth version (${TU}/depth_${U}_${VU}.[ab])
# 
setenv R GLBb0.08
setenv U IAS0.03
setenv E 737
setenv Y 119

setenv YY `echo ${Y} | awk '{printf("%04d", $1+1900)}'`
setenv EX `echo ${E} | awk '{printf("%04.1f", $1*0.1)}'`

setenv S  /home/ddmitry/codes/anls_mtlb_utils/hycom_TSIS/prepare_nestGLB0.08_py/subregion
setenv TR /nexsan/people/ddmitry/hycom/ARCc0.04_022/misc
setenv TU /Net/kronos/ddmitry/hycom/TSIS/topo_IAS0.03 
setenv VR 09m11
setenv VU 09m11_ias0.03
setenv D /nexsan/people/ddmitry/hycom/ARCc0.04_022/misc 
setenv N /Net/kronos/ddmitry/hycom/TSIS/nest_files_gofs3.1

setenv ALL /nexsan/people/ddmitry/Net_ocean/HYCOM/hycom/ALL4
setenv BINRUN ""
mkdir -pv $N
#
#mkdir -pv $S
cd        $S
#
/bin/rm   regional.depth.a regional.depth.b
touch     regional.depth.a regional.depth.b

if (-z    regional.depth.a) then
  /bin/rm regional.depth.a
  /bin/ln -s ${TR}/depth_${R}_${VR}.a regional.depth.a
endif
if (-z    regional.depth.b) then
  /bin/rm regional.depth.b
  /bin/ln -s ${TR}/depth_${R}_${VR}.b regional.depth.b
endif
#
touch     regional.grid.a regional.grid.b
if (-z    regional.grid.a) then
  /bin/rm regional.grid.a
  /bin/ln -s ${TR}/regional.grid.a .
endif
if (-z    regional.grid.b) then
  /bin/rm regional.grid.b
  /bin/ln -s ${TR}/regional.grid.b .
endif
#
## determine start and end days, one day overlap for mean files
#
#
foreach iday ( 001 )
	setenv FLG ${E}_archm.${YY}_${iday}_12
#   no leading expt number for nesting files
	setenv FL archm.topo09m11ias003_${YY}_${iday}   # output file with interp archm IAS 0.03
	/bin/rm -f ${N}/${FL}.[ab]
	${ALL}/subregion/src/isubaregion <<E-o-D
${TU}/regional.grid.a
${TU}/regional.gmapi_${R}.a
${TU}/depth_${VU}.a
regional.depth.a
${D}/${FLG}.a
${N}/${FL}.a
${R} interpolated to ${U}
1401   'idm   ' = target longitudinal array size
 891   'jdm   ' = target latitudinal  array size
	 0   'iceflg' = ice in output archive flag (0=none,1=energy loan model)
	 0   'smooth' = smooth interface depths    (0=F,1=T)
E-o-D

	touch  ${N}/${FL}.b
	if (-z ${N}/${FL}.b) then
		echo "missing archive file: " ${N}/${FL}.b
		exit (2)
	endif

	cd ${S}
  date
end

\rm region*






#!/bin/csh -x
#
set echo
set time = 1
# --- Convert z-level interpolated files with 
# NEMO+GLORYS interpolated onto 75 Z-levels NEMO
# for the IAS HYCOM-TSIS domain on
# 0.03-defree HYCOM-TSIS grid
# Considered as ("climatology files")
# Use Alan W's code to interpolate into 
# HYCOM-TSIS vertical layers
setenv pget /bin/cp
setenv pput /bin/cp

setenv mm 07   # month
setenv SG sig2
setenv NL 30 # vert layers

# P - work directory for creating realx files
#     from NEMO-GLORYS Z-level interpolated files
# S - temporary scratch dir
setenv P /Net/kronos/ddmitry/hycom/TSIS/relax
setenv S $P/SCRATCH
setenv D /Net/kronos/ddmitry/hycom/TSIS/nemo_glorys_interp
setenv DT /Net/kronos/ddmitry/hycom/TSIS/topo_IAS0.03 
setenv ALL /nexsan/people/ddmitry/Net_ocean/HYCOM/hycom/ALL4 
setenv DUMP ${P}/DUMP
setenv DOUT $P

mkdir -pv $DUMP
mkdir -pv $S
#mkdir -pv $DOUT
cd $S
#
touch   relax_${SG} blkdat.input fort.51 fort.51A
/bin/rm relax_${SG} blkdat.input fort.51 fort.51A

#
# Input fields
#
touch      fort.73 fort.73A fort.72 fort.72A
/bin/rm -f fort.73 fort.73A fort.72 fort.72A

setenv dstamp 09072011
setenv flnm zlev_hycom_sig2_${dstamp}
${pget} ${D}/temp_${flnm}.b fort.72 &
${pget} ${D}/temp_${flnm}.a fort.72A &
${pget} ${D}/saln_${flnm}.b fort.73 &
${pget} ${D}/saln_${flnm}.a fort.73A &

#
# Grid and topo
touch fort.51 fort.51A
if (-z fort.51) then
			/bin/cp $DT/regional.depth.b fort.51 &
endif
if (-z fort.51A) then
			/bin/cp $DT/regional.depth.a fort.51A &
endif

touch regional.grid.a regional.grid.b
if (-z regional.grid.b) then
		${pget} ${DT}/regional.grid.b regional.grid.b &
endif
if (-z regional.grid.a) then
		${pget} ${DT}/regional.grid.a regional.grid.a &
endif
#    C
touch blkdat.input
if (-z blkdat.input) then
  ${pget} ${D}/../blkdat.input_clim blkdat.input &
endif

setenv HEX relaxi
touch ${HEX}
if (-z ${HEX}) then
  /bin/cp ${ALL}/relax/src_2.2.35/${HEX} . &
endif
wait

chmod a+rx ${HEX}
sed -e "s/^[  0-9]*'month ' =/  ${mm}   'month ' =/" blkdat.input >! fort.99

/bin/rm -f core
touch core

setenv FOR010A fort.10A
setenv FOR011A fort.11A
setenv FOR012A fort.12A
setenv FOR021A fort.21A
setenv FOR051A fort.51A
setenv FOR072A fort.72A
setenv FOR073A fort.73A
/bin/rm -f fort.10  fort.11  fort.12  fort.21
/bin/rm -f fort.10A fort.11A fort.12A fort.21A

echo "Submitting ${HEX}"

./${HEX}

mv fort.10  ${DOUT}/relax_tem_m${dstamp}.b
mv fort.10A ${DOUT}/relax_tem_m${dstamp}.a
mv fort.11  ${DOUT}/relax_sal_m${dstamp}.b
mv fort.11A ${DOUT}/relax_sal_m${dstamp}.a
mv fort.12  ${DOUT}/relax_int_m${dstamp}.b
mv fort.12A ${DOUT}/relax_int_m${dstamp}.a

mv fort.21  ${DOUT}/relax.m${dstamp}.b
mv fort.21A ${DOUT}/relax.m${dstamp}.a

/bin/rm fort.7[12]
/bin/rm fort.7[12]A

# Move used climat files
#cd ${D}
#foreach ff (temp saln)
#  /bin/mv ${ff}_${flnm}.b ${ff}_${flnm}.b-bkp
#  /bin/mv ${ff}_${flnm}.a ${ff}_${flnm}.a-bkp
#end

exit 0











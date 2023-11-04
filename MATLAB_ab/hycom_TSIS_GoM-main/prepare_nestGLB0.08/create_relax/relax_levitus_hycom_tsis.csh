#!/bin/csh -x
#
set echo
set time = 1
# --- Convert z-level Levitus climatology
# Use Alan W's code to interpolate into 
# HYCOM-TSIS vertical layers
setenv pget /bin/cp
setenv pput /bin/cp

#setenv mm 01   # month
setenv SG sig2
setenv NL 41 # vert layers

# P - work directory for creating realx files
# S - temporary scratch dir
setenv P /Net/kronos/ddmitry/hycom/TSIS/relax_41lrs
setenv S $P/SCRATCH
setenv D $P
setenv DT /Net/kronos/ddmitry/hycom/TSIS/topo_IAS0.03 
setenv ALL /nexsan/people/ddmitry/Net_ocean/HYCOM/hycom/ALL4 
setenv DUMP ${P}/DUMP
setenv DOUT $P

mkdir -pv $DUMP
#mkdir -pv $S
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

foreach mm( 01 02 03 04 05 06 07 08 09 10 11 12 )
  /bin/rm -f fort.7*

	${pget} ${D}/temp_${SG}_m${mm}.b fort.72 &
	${pget} ${D}/temp_${SG}_m${mm}.a fort.72A &
	${pget} ${D}/saln_${SG}_m${mm}.b fort.73 &
	${pget} ${D}/saln_${SG}_m${mm}.a fort.73A &

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
		${pget} ${P}/../blkdat.input_clim_41lrs blkdat.input &
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

	mv fort.10  ${DOUT}/relax_tem_m${mm}.b
	mv fort.10A ${DOUT}/relax_tem_m${mm}.a
	mv fort.11  ${DOUT}/relax_sal_m${mm}.b
	mv fort.11A ${DOUT}/relax_sal_m${mm}.a
	mv fort.12  ${DOUT}/relax_int_m${mm}.b
	mv fort.12A ${DOUT}/relax_int_m${mm}.a

	mv fort.21  ${DOUT}/relax.m${mm}.b
	mv fort.21A ${DOUT}/relax.m${mm}.a

	/bin/rm fort.7[12]
	/bin/rm fort.7[12]A
end  # month loop

#
##
# --- Merge or not "climatology" files into one file.
# #
setenv MRG yes
if (${MRG} == 'yes') then
  cd ${DOUT}
  /bin/cp -f relax_int_m01.b relax_int.b
  /bin/cp -f relax_sal_m01.b relax_sal.b
  /bin/cp -f relax_tem_m01.b relax_tem.b
  ls -l relax_*.b

  foreach MM ( 02 03 04 05 06 07 08 09 10 11 12 )
  tail --lines=+6 relax_int_m${MM}.b >> relax_int.b
  tail --lines=+6 relax_sal_m${MM}.b >> relax_sal.b
  tail --lines=+6 relax_tem_m${MM}.b >> relax_tem.b
  end

  /bin/cp -f relax_int_m01.a relax_int.a
  /bin/cp -f relax_sal_m01.a relax_sal.a
  /bin/cp -f relax_tem_m01.a relax_tem.a

  foreach MM ( 02 03 04 05 06 07 08 09 10 11 12 )
  cat relax_int_m${MM}.a >> relax_int.a
  cat relax_sal_m${MM}.a >> relax_sal.a
  cat relax_tem_m${MM}.a >> relax_tem.a
  end
#
## --- delete the monthly files
#
  /bin/rm relax_int_m??.[ab]
  /bin/rm relax_sal_m??.[ab]
  /bin/rm relax_tem_m??.[ab]
else
  /bin/rm relax_*_[12]*.[ab]
endif



exit 0











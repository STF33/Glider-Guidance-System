#! /bin/csh -x
#
## --- form subregion grid array index map file, GLBb0.08 --> IAS 0.03
#
# ALL - directory with hycom util codes
# R - directory with subregion topo
# NGLB - outer domain that is being subset
# NSB  - name of the subset domain
#
setenv ALL /nexsan/people/ddmitry/Net_ocean/HYCOM/hycom/ALL4
setenv NSB IAS0.03
setenv R /Net/kronos/ddmitry/hycom/TSIS/topo_IAS0.03
#
setenv NGLB GLBb0.08 
setenv RGLB /nexsan/people/ddmitry/hycom/ARCc0.04_022/misc
ln -s ${RGLB}/regional.grid.? .
#
/bin/rm -f ${R}/regional.gmapi_${NGLB}.[ab]
${ALL}/subregion/src/isuba_gmapi <<E-o-D
${R}/regional.grid.a
${R}/regional.gmapi_${NGLB}.a
${NGLB} (4500x3298) to ${NSB} (maxinc=9)
 1401     'idm   ' = target longitudinal array size
  891     'jdm   ' = target latitudinal  array size
    9     'maxinc' = maximum input array index jump on target grid
E-o-D



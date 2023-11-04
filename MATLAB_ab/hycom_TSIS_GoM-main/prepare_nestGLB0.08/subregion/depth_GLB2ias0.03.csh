#! /bin/csh -x
#
## --- form subregion bathymetry GLBb0.08 --> IAS 0.03
#
# ALL - directory with hycom util codes
# R - directory with subregion topo
# NGLB - outer domain that is being subset
# NSB  - name of the subset domain
# TGLB - topo for outer domain
#
setenv ALL /nexsan/people/ddmitry/Net_ocean/HYCOM/hycom/ALL4
setenv NSB IAS0.03
setenv R /Net/kronos/ddmitry/hycom/TSIS/topo_IAS0.03
#
setenv NGLB GLBb0.08 
setenv RGLB /nexsan/people/ddmitry/hycom/ARCc0.04_022/misc
setenv TGLB 09m11
ln -s ${RGLB}/regional.grid.? .
#
/bin/rm -f ${R}/depth_GLBb0.08${TGLB}_ias0.03.[ab]
${ALL}/subregion/src/isuba_topog <<E-o-D
${R}/regional.gmapi_${NGLB}.a
${RGLB}/depth_GLBb0.08_09m11.b
${R}/depth_09m11_ias0.03.b
depth_GLBb0.08${TGLB} subregioned to ${NSB} via isuba_topog
 1401     'idm   ' = target longitudinal array size
  891     'jdm   ' = target latitudinal  array size
E-o-D



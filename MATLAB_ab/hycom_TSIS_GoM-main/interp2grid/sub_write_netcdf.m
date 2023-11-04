function sub_write_netcdf(flcdf,nlon,nlat,Znas,sshi,Tzx,Szx,Uzx,Vzx,dnmb,huge);
% Create new netcdf file
% write fields
% Note fields will be rotated for netcdf 
%
%
d0=datenum(1900,12,31);
nday=dnmb-d0;
[mnas,nnas]=size(nlon);
knas=length(Znas);

fprintf('Writing %s\n',flcdf);

% Open file      
ncid=netcdf.create(flcdf,'NETCDF4');

% Define dimensions
dimidt = netcdf.defDim(ncid,'time',1);
dimidy = netcdf.defDim(ncid,'north',mnas);
dimidx = netcdf.defDim(ncid,'east',nnas);
dimidz = netcdf.defDim(ncid,'depth',knas);

% Define variables
dayID  = netcdf.defVar(ncid,'time','double',1);
lonID  = netcdf.defVar(ncid,'longitude','double',[dimidx, dimidy]);
latID  = netcdf.defVar(ncid,'latitude','double',[dimidx, dimidy]);
dpthID = netcdf.defVar(ncid,'depth','double',[dimidz]);
sshID  = netcdf.defVar(ncid,'srfhgt','double',[dimidx dimidy]);
tempID = netcdf.defVar(ncid,'temp','double',[dimidx dimidy dimidz]);
salnID = netcdf.defVar(ncid,'saln','double',[dimidx dimidy dimidz]);
uvelID = netcdf.defVar(ncid,'u-vel','double',[dimidx dimidy dimidz]);
vvelID = netcdf.defVar(ncid,'v-vel','double',[dimidx dimidy dimidz]);

netcdf.defVarFill(ncid,sshID,false,huge);
netcdf.defVarFill(ncid,tempID,false,huge);
netcdf.defVarFill(ncid,salnID,false,huge);
netcdf.defVarFill(ncid,uvelID,false,huge);
netcdf.defVarFill(ncid,vvelID,false,huge);

% Global Attributes:
varid = netcdf.getConstant('GLOBAL');
netcdf.putAtt(ncid,varid,'info','HYCOM-TSIS GoM 0.03, COAPS FSU');
netcdf.putAtt(ncid,varid,'creation_date',datestr(now));


netcdf.endDef(ncid);

% Time
%netcdf.putVar(ncid,dayID,nday);
ncwrite(flcdf,'time',nday)
ncwriteatt(flcdf, 'time', 'standard_name', 'Days');
ncwriteatt(flcdf, 'time', 'units','days since 1900-12-31 00:00:00');

% Put grid
netcdf.putVar(ncid,lonID,nlon');
ncwriteatt(flcdf, 'longitude', 'standard_name', 'Longitude');
ncwriteatt(flcdf, 'longitude', 'units','degrees_east');
ncwriteatt(flcdf, 'longitude', 'axis','X');

netcdf.putVar(ncid,latID,nlat');
ncwriteatt(flcdf, 'latitude', 'standard_name', 'Latitude');
ncwriteatt(flcdf, 'latitude', 'units','degrees_north');
ncwriteatt(flcdf, 'latitude', 'axis','Y');

netcdf.putVar(ncid,dpthID,Znas);
ncwriteatt(flcdf, 'depth', 'standard_name', 'Depths of fixed Z levels');
ncwriteatt(flcdf, 'depth', 'units','m');
ncwriteatt(flcdf, 'depth', 'axis','Z');

% 2D SSH
netcdf.putVar(ncid,sshID,sshi');
ncwriteatt(flcdf, 'srfhgt', 'standard_name', 'sea surface height');
ncwriteatt(flcdf, 'srfhgt', 'units','m');
%netcdf.reDef(ncid);
%netcdf.defVarFill(ncid,sshID,false,huge);
%%netcdf.endDef(ncid);

% 3D Temp
%fprintf('    Writing temperature\n');
% rotate fields for netcdf writing  
Fout=[];
for kk=1:knas
  A=squeeze(Tzx(kk,:,:));
  Fout(:,:,kk)=A';
end
netcdf.putVar(ncid,tempID,Fout);
ncwriteatt(flcdf, 'temp', 'standard_name', 'potential temperature');
ncwriteatt(flcdf, 'temp', 'units','degrees Celsius');

% 3D Salin
Fout=[];
for kk=1:knas
  A=squeeze(Szx(kk,:,:));
  Fout(:,:,kk)=A';
end
netcdf.putVar(ncid,salnID,Fout);
ncwriteatt(flcdf, 'saln', 'standard_name', 'salinity');
ncwriteatt(flcdf, 'saln', 'units','salinity units');

% 3D U-comp
Fout=[];
for kk=1:knas
  A=squeeze(Uzx(kk,:,:));
  Fout(:,:,kk)=A';
end
netcdf.putVar(ncid,uvelID,Fout);
ncwriteatt(flcdf, 'u-vel', 'standard_name', 'u-component of ocean velocity vector');
ncwriteatt(flcdf, 'u-vel', 'units','m/s');
ncwriteatt(flcdf, 'u-vel', 'positive dir','eastward');


% 3D V-comp
Fout=[];
for kk=1:knas
  A=squeeze(Vzx(kk,:,:));
  Fout(:,:,kk)=A';
end
netcdf.putVar(ncid,vvelID,Fout);
ncwriteatt(flcdf, 'v-vel', 'standard_name', 'v-component of ocean velocity vector');
ncwriteatt(flcdf, 'v-vel', 'units','m/s');
ncwriteatt(flcdf, 'v-vel', 'positive dir','northward');

netcdf.close(ncid);

%fprintf(' File created \n\n');

return
% from a struct array with POOL(Nfgr).Time(Ntime).WK() 
% get monthly/weekly data from MHD struct array to POOL
% Nfgr - # of forecast groups
% Ntime - # of time levels
%
function POOL = sub_pool_month_week_OSE(POOL,dmm,ifc,WK);

%dmm=MHD(ifc).MHD;
[a1,a2] = size(dmm);

if isempty(dmm); 
  fprintf(' WARN: sub_pool_month_week - empty input array dmm\n\n');
  return; 
end;

% Pool all runs for the same forecast group and time period
pm1=POOL(ifc).mo(1).pm;
pm2=POOL(ifc).mo(2).pm;
pm3=POOL(ifc).mo(3).pm;

A=dmm(1:30,:);
[b1,b2]=size(A);
A=reshape(A,[b1*b2,1]);
pm1=[pm1;A];

A=dmm(31:60,:);
[b1,b2]=size(A);
A=reshape(A,[b1*b2,1]);
pm2=[pm2;A];

A=dmm(61:91,:);
[b1,b2]=size(A);
A=reshape(A,[b1*b2,1]);
pm3=[pm3;A];

POOL(ifc).mo(1).pm=pm1;
POOL(ifc).mo(2).pm=pm2;
POOL(ifc).mo(3).pm=pm3;

% Weekly pools
nwk=length(WK)-1;
for iwk=1:nwk
	id1=WK(iwk);
	id2=WK(iwk+1)-1;
	pw=POOL(ifc).WK(iwk).pw;
	A=dmm(id1:id2,:);
	[b1,b2]=size(A);
	A=reshape(A,[b1*b2,1]);
	pw=[pw;A];
	POOL(ifc).WK(iwk).pw=pw;
end

return

 

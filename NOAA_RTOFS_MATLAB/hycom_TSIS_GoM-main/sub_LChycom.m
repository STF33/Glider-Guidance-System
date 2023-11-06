function LCH = sub_LChycom(fina,finb,LON,LAT,Bisol,INH);
% Extract LC contour
% from HYCOM simulations
rg=9806;  % convert pressure to depth, m
huge=1e20;

fprintf('Reading %s\n',fina);
fld = 'srfhgt';
[F,nn,mm,ll] = read_hycom(fina,finb,fld);
F(F>huge)=nan;
ssh=squeeze(F)./(1e-3*rg);  % ssh m
%
% Subtract anomaly:
dmm=ssh;
dmm(INH==0)=nan;
%  dmm(HH>-200)=nan;
sshM=nanmean(nanmean(dmm));
ssh=ssh-sshM;

%
% Derive LC contour:
% 
dmm=ssh;
dmm(INH==0)=nan;
LCH = identify_LC(LON,LAT,dmm,Bisol);


return

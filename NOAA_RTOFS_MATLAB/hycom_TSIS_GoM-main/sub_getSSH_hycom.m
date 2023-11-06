% Read SSH hycom fields
% specified demeaned or not: dmean
% INH only needed if dmean is not 0
%
function ssh = sub_getSSH_hycom(fina,finb,INH,dmean);

rg=9806;
huge=1e20;

ie = exist(fina,'file');

if ~ie
		fprintf('sub_getSSH_hycom:  ERR ** Missing HYCOM file: %s\n',fina);
		keyboard
end

fprintf('Reading %s\n',fina);
fld = 'srfhgt';
[F,nn,mm,ll] = read_hycom(fina,finb,fld);
F(F>huge)=nan;
ssh=squeeze(F)./(1e-3*rg);  % ssh m
%
% Subtract anomaly:
if dmean>0
		dmm=ssh;
%		dmm(INH==0)=nan;
		%  dmm(HH>-200)=nan;
%		sshM=nanmean(nanmean(dmm));
  I=find(INH==1);
  sshM=nanmean(ssh(I));
		ssh=ssh-sshM;
end

return

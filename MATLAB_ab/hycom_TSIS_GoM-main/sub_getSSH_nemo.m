% Extract SSH field
% demean if dmean>0
% from NEMO 1/100 simulation 
% used for OSSE hindcasts
%
function ssh = sub_getSSH_nemo(DV,LONN,LATN,INN,dmean);

fpwd='https://intercambio:5w9SHURN@cigom.cicese.mx/thredds/dodsC/intercambio/GOLFO108-NAS00-S/';


yr=DV(1);
mo=DV(2);
dm=DV(3);
dnmb1=datenum(yr,mo,1);
dnmb2=dnmb1+32;
v2=datevec(dnmb2);
dnmb2=datenum(v2(1),v2(2),1);
d2=dnmb2-datenum(yr,mo,1);



fnemo=sprintf('GOLFO108-NAS00_1d_%4.4i%2.2i%2.2i_%4.4i%2.2i%2.2i_grid_T.nc',yr,mo,1,yr,mo,d2);
fprintf('NEMO requested date: %i/%i/%i\n',DV(1:3));
fprintf('Reading NEMO: %s\n',fnemo);

fin = sprintf('%s/%i/%s',fpwd,yr,fnemo);
%keyboard

enm = squeeze(ncread(fin,'ssh',[1 1 dm],[Inf Inf 1]));
enm = enm';
I=find(enm==0);
enm(I)=nan;

% Subtract spatial mean ssh
if dmean>0
%		dmm=enm;
%		dmm(INN==0)=nan;
%		sshM=nanmean(nanmean(dmm));
  I=find(INN==1);
  sshM = nanmean(enm(I));
		enm = enm-sshM;
end

ssh=enm;


return

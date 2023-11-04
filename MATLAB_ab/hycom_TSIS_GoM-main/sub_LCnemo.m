function LCN = sub_LCnemo(yr,mo,d2,fpwd,LONN,LATN,Bisol,INN,dm);
% Extract LC contour
% from NEMO 1/100 simulation 
% used for OSSE hindcasts

fnemo=sprintf('GOLFO108-NAS00_1d_%4.4i%2.2i%2.2i_%4.4i%2.2i%2.2i_grid_T.nc',yr,mo,1,yr,mo,d2);
fprintf('Reading NEMO: %s\n',fnemo);

fin = sprintf('%s/%i/%s',fpwd,yr,fnemo);

if isempty(LONN)
		fmesh=sprintf('%smesh_mask.nc',fpwd);
		dmm = ncread(fmesh,'nav_lon');
		LONN = dmm';
		dmm = squeeze(ncread(fmesh,'nav_lat'));
		LATN = dmm';

		[mm,nn] = size(LONN);

		[XM,YM]=meshgrid([1:nn],[1:mm]);
		INN = inpolygon(XM,YM,GOMN(:,1),GOMN(:,2));
		clear XM YM

end


enm = squeeze(ncread(fin,'ssh',[1 1 dm],[Inf Inf 1]));
enm = enm';
I=find(enm==0);
enm(I)=nan;

% Subtract spatial mean ssh
dmm=enm;
dmm(INN==0)=nan;
sshM=nanmean(nanmean(dmm));
enm = enm-sshM;

Bisol=0.17;
dmm=enm;
dmm(INN==0)=nan;
LCN = identify_LC(LONN,LATN,dmm,Bisol);



return

function mntg1 = sub_montg1_regr2(ssh,fmat);
% Use regression to estimate montg1
% from the SSH field
% see: fit_regr_montgomery.m
%
% Make sure that land mask of montg matches the SSH
% otherwise HYCOM will generate huge velocities
% pbavg=SSH-montg1
%
fprintf('Calculating montgomery ...\n');
gg     = 9.806;
%pthmat = sprintf('/Net/ocean/ddmitry/HYCOM/GoM/GOMl0.04/080/data_anls/');
%pthmat='/Net/mars/ddmitry/hycom/hycom_TSIS/data_mat/';
%fmat=sprintf('%sregression_montg_ssh_IAS003.mat',pthmat);

load(fmat);

hg=2^100;
%[F,n,m,l] = read_hycom(finaOLD,finbOLD,'srfhgt');
%F=F./gg;  % SSH pressure -> m
%I=find(F>1e20);
%F(F>1e20)=nan;
%ssh=squeeze(F);
I = find(isnan(ssh));

B0=RGR.intercept;
B1=RGR.slope;

mntg1 = B0 + B1.*ssh;

% Filter montgomery field
% to avoid high-freq noise near the OB
nij=2; % Filter is (2*nij+1)x(2*nij+1)
eMf=sub_fltr_Gauss(nij,mntg1);
mntg1=eMf;
mntg1(I)=hg;

%keyboard

return

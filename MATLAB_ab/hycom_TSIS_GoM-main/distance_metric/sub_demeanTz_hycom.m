% Demean Tz interpolated onto z level
% dTmn - demeaned T inside GoM domain
% tMh - spatially mean T inside GoM
%
function [dTmn,tMh] = sub_demeanTz_hycom(dmm);

%
% Subtract spatial mean:
%      dmm(INH==0)=nan;
%      tMh=nanmean(nanmean(dmm));
dmm(1:121,312:end)=nan; % Caribbean 
dmm(:,504:end)=nan;

tMh = nanmean(nanmean(dmm));
dTmn = dmm-tMh;


return

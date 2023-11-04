% Find best LC/LCE combination given contorol contour (e.g., NEMO Nature Run)
% For NEMO - try different LC/LCE combinations as well
% Xc,Yc - control contours
% Use MHD - min MHD is the best contour combination
% for specified day ihc - f/cast day #
%
function [Xhc,Yhc,Xc,Yc,mhd_fcst] = sub_best_contour_mhd(LCE,LCXY,xn,yn,ihc,LCEN,irc,lonW);

xh1 = LCXY.XY(ihc).X;
yh1 = LCXY.XY(ihc).Y;

Nlc1 = sub_numbLCE(LCEN,irc);  % # of LCEs 
[ECombN, icombN] = sub_lce_comb(Nlc1);  % Possible LC/LCE comb for NEMO

Nlc2 = sub_numbLCE(LCE,ihc); % # of LCE on this date 
[EComb, icomb] = sub_lce_comb(Nlc2);  % all possible comb. of LC/LCE


% NEMO NR - all LCEs 
% Try all LC/LCE combinations for NEMO and HYCOM OSSE 
% to get best MHD
MHD_LCLCE = zeros(icombN,1)*nan;
MHD_Prs = MHD_LCLCE;
flg_lce1 = [];
for ilce=1:icombN   % LC/LCE combination for NEMO NR
  [Xc,Yc,flg] = sub_combineLCLCE_EComb(xn,yn,LCEN,irc,ECombN,ilce,'lonW',lonW);

  if flg < 0
    MHD_LCLCE(ilce) = 1.e30;  % no LCEs comb - west of lonW
    continue
  end

% For HYCOM forecast - try all possible combinations of LC/LCEs and find the best
% that matches NEMO 
  [Xhc,Yhc,mhd_osse] = sub_best_contour_mhd(LCE,LCXY,Xc,Yc,ihc);
  MHD_LCLCE(ilce) = mhd_osse;
end

mhd_min = min(MHD_LCLCE);
iNbest = find(MHD_LCLCE==mhd_min,1);
[Xc,Yc,flg] = sub_combineLCLCE_EComb(xn,yn,LCEN,irc,ECombN,iNbest,'lonW',lonW);
[Xhc,Yhc,mhd_fcst] = sub_best_contour_mhd(LCE,LCXY,Xc,Yc,ihc);




return

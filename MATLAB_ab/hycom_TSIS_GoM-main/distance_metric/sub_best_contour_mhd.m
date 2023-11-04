% Find best LC/LCE combination given contorol contour (e.g., NEMO Nature Run)
% Xc,Yc - control contours
% Use MHD - min MHD is the best contour combination
% for specified day ihc - f/cast day #
%
function [Xhc,Yhc,mhd_fcst] = sub_best_contour_mhd(LCE,LCXY,Xc,Yc,ihc);

xh1 = LCXY.XY(ihc).X;
yh1 = LCXY.XY(ihc).Y;

Nlc2 = sub_numbLCE(LCE,ihc); % # of LCE on this date 
[EComb, icomb] = sub_lce_comb(Nlc2);  % all possible comb. of LC/LCE

% For HYCOM forecast - try all possible combinations of LC/LCEs and find the best
% that matches the NEMO 
for jlce=1:icomb   % adding eddies to the f/cast
  [Xhc,Yhc,flg] = sub_combineLCLCE_EComb(xh1,yh1,LCE,ihc,EComb,jlce);
  if jlce>1 & flg==0 % no LCE 
    MHD_LCLCE(jlce,1)=1.e20;
    continue;
  end

  P = [Xc,Yc];    % NEMO NR
  Q = [Xhc,Yhc];  % HYCOM fcst LC 
  mhd = modified_hausdorff_distance(P,Q,'geo');
  MHD_LCLCE(jlce,1) = mhd;
end

mhd_fcst = min(min(MHD_LCLCE));
imhd0 = find(MHD_LCLCE==mhd_fcst,1);
% Combine with the contour that produces best MHD
[Xhc,Yhc,flg] = sub_combineLCLCE_EComb(xh1,yh1,LCE,ihc,EComb,imhd0);

%keyboard

return

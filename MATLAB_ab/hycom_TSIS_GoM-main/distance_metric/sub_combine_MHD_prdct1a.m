% Predictability 1a forecasts - initialized from NEMO+GLORYS 
% interpolated onto HYCOM horiz/vert grid 
% 7 forecasts ran with 1 week shift starting May 1
% combine all 7 runs as 1 forecast group to compare 
% with OSSE forecast groups (AVISO, PIES, etc.)
% mhd - see mhd_LCLCEcntrPrdct_nemo_fcsthycom.m
% Grab only runs 1 e.g. MHD_LCLCE_nemo_persist_OSSEfcst16-0201.mat
% 0202 is IC with 2day shit, 0203 -1 day, 0204 +1, 0205 +2 - not needed
function MHD = sub_combine_MHD_prdct1a(MHD,pthmat);
ifc = length(MHD);
ifc = ifc+1;
MHD(ifc).Fcst_OSSE = 'Predictability 1a NEMO+GLORYS expts10-16 run1';
MHD(ifc).Fcst_nmb = '10-16';
for itime = 1:2
  irr = 0;
  for iFcst = 10:16
    FCST  = sub_fcstPrdct_info(iFcst);
    Nhnd  = FCST.Nhind;  % initial cond from the  hindcast #
    nmexp = FCST.Hindcast_Name;
    pthHcst = FCST.Hindcast_path;    % hindcast path with initial fields

    irun = 1;  % pick only prdct 1a runs
    nmexp = sprintf('fcst%2.2i-%2.2i%2.2i',Nhnd,itime,irun);
%    fmat1 = sprintf('%sMHD_LCLCE_nemo_persist_hycom%s.mat',...
%                pthmat,nmexp);
    fmat1 = sprintf('%sMHD_LCLCE_nemo_persist_OSSE%s.mat',...
                  pthmat,nmexp);

    fprintf('Loading %s\n',fmat1);
    A=load(fmat1);  % MHD for f/cast and persistence
    DV=datevec(A.TM);

    irr = irr+1;
    fcst_name = sprintf('prdct1a%2.2i-%2.2i%2.2i',Nhnd,itime,irun); 
    MHD(ifc).Time(itime).FcstName = fcst_name;
    MHD(ifc).Time(itime).MHD(:,irr) = A.MHD(:,1);
    MHD(ifc).Time(itime).MHD_prst(:,irr) = A.MHD(:,2); % persist HYCOM on day 1 modified code
%    MHD(ifc).Time(itime).MHD_prst(:,irr) = A.MHD(:,3); % persist HYCOM on day 1

%    MHD(ifc).nemo0 = A.MHD(:,2); % NEMO true contour on day 1
%    MHD(ifc).hycom0= A.MHD(:,3); % HYCOM contour on day 1
  end
end

%keyboard

return

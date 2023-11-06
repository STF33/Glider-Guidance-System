% Prepare mat file with hindcast information
% Epxeriments <10 - first round of forecasts experiments
%                   Initial state - from IAS HYCOM-TSIS data assimilative hindcast
%                   NEMO synthetic observations are assimilated
%                   mimicking real obs (AVISO, SST, PIES T/S profiles)
%
%             10  - second set of predictability experiments
%                   Initial state - from 3D interpolated NEMO+GLORYS fields
%                   mapped onto HYCOM horiz and vertical grid
ii=0;
EXPT = struct;
% Experiments:
% freerun
ii=ii+1;
EXPT(ii).Name = 'FreeRun';
EXPT(ii).Name_short = 'FreeRun';
EXPT(ii).path = '/Net/gleam/abozec/HYCOM/TSIS/IASx0.03/NEMO/freerun/';

% No SSS or SST analysis:
% full field SSH + no pies  
ii=ii+1;
EXPT(ii).Name = '2DSSH noPIES noSSS noSST';
EXPT(ii).Name_short = '2DSSH';
EXPT(ii).path = '/Net/gleam/abozec/HYCOM/TSIS/IASx0.03/NEMO/hindcast_full_nopies_newtsis/';
% AVISO tracks SSH + no pies 
ii=ii+1;
EXPT(ii).Name = 'AVISOSwathsSSH noPIES noSSS noSST';
EXPT(ii).Name_short = 'AVISO';
EXPT(ii).path = '/Net/gleam/abozec/HYCOM/TSIS/IASx0.03/NEMO/hindcast_aviso_nopies_newtsis/';
% “one track” SSH + no pies
ii=ii+1;
EXPT(ii).Name = '1swathSSH noPIES noSSS noSST';
EXPT(ii).Name_short = 'AVISO 1 swath';
EXPT(ii).path = '/Net/gleam/abozec/HYCOM/TSIS/IASx0.03/NEMO/hindcast_onetrack_nopies_newtsis/';
% AVISO tracks SSH + ugos pies (no SSS and no SST)
ii=ii+1;
EXPT(ii).Name = 'AVISOSwathSSH ugosPIES noSSS noSST';
EXPT(ii).Name_short = 'AVISO PIES noSST';
EXPT(ii).path = '/Net/gleam/abozec/HYCOM/TSIS/IASx0.03/NEMO/hindcast_aviso_ugos_nosss_nosst_newtsis/';

% No SSS analysis but with SST analysis:
%No SSH track + pies distributed all over the GOM domain (1/30 points)
ii=ii+1;
EXPT(ii).Name = 'noSSH allGoMPIES noSSS SST';
EXPT(ii).Name_short = 'TSprof dx=30';
EXPT(ii).path = '/Net/gleam/abozec/HYCOM/TSIS/IASx0.03/NEMO/hindcast_nosla_gompies_nosss_newtsis/';
%AVISO tracks SSH + ugos pies (small area distribution)
ii=ii+1;
EXPT(ii).Name = 'AVISOSwathSSH ugosPIES noSSS SST';
EXPT(ii).Name_short = 'AVISO PIES';
EXPT(ii).path = '/Net/gleam/abozec/HYCOM/TSIS/IASx0.03/NEMO/hindcast_aviso_ugos_nosss_newtsis/';
% AVISO tracks SSH + extd pies (bigger area distribution)
ii=ii+1;
EXPT(ii).Name = 'AVISOSwathSSH extendedPIES noSSS SST';
EXPT(ii).Name_short = 'AVISO extdPIES';
EXPT(ii).path = '/Net/gleam/abozec/HYCOM/TSIS/IASx0.03/NEMO/hindcast_aviso_extd_nosss_newtsis/';
% No SSH track + pies distributed all over the GOM domain  (1/60 points)
ii=ii+1;
EXPT(ii).Name = 'NoSSH allGoMPIES60 noSSS SST';
EXPT(ii).Name_short = 'TSprof dx=60';
EXPT(ii).path = '/Net/gleam/abozec/HYCOM/TSIS/IASx0.03/NEMO/hindcast_nosla_gompies60_nosss_newtsis/';
%
%
%  PREDICTABILITY experiments:
%  Same time periods as previous f/casts (May 2011 and Jan 2012)
%  but initial conditions for HYCOM are interpolated NEMO+GLORYS fields
%  (not assimilative hindcasts)
for ii=10:16;
		EXPT(ii).Name = 'Predictability NEMO+GLORYS 3D interpolated';
		EXPT(ii).Name_short = sprintf('Predictab %i',ii-9);
  EXPT(ii).path = '/nexsan/people/ddmitry/hycom/TSIS/IASx0.03/output_predictability/';		
%	EXPT(ii).path = '/Net/kronos/ddmitry/hycom/TSIS/nemo_glorys_interp/'; %initial fields
%EPXT(ii).path = '/Net/gleam/dmitry/hycom/TSIS/IASx0.03/output_predictability/'; 
end

fhnd = 'hycom_tsis_expts.mat';
fprintf('Saving %s\n',fhnd);
save(fhnd,'EXPT');



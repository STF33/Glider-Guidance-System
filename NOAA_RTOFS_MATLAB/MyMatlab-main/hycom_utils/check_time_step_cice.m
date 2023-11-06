function check_time_step_cice(tbcl,tbtr,tice,tcpl,dsurfq);
%function check_time_step(tbcl,tbtr,tice,tcpl,dsurfq);
% dsurfq - check blkdat.input, # of days btw model diagnstc. at surface
%          can be 0
% HYCOM: barocl & barotr time steps have 
% CICe: tice - time step for CICE, sec. (600)
% tcpl: sea ice coupling, seconds (see blkdat.input cplifq, : 
%             tcpl=3600*24*cplifq)
% to be related and satisfy several requirements
%tbcl=120;
%tbtr=4;
%dsurfq=0.125;

c0=86400/tbcl;
if dsurfq>0
  c1=86400/tbcl*dsurfq;
end
c2=tbcl/2*1/tbtr;

fprintf('T baroclinic = %8.3f\n',tbcl);
fprintf('T barotropic = %8.3f\n',tbtr);
fprintf('Number of days - surface diagn = %8.3f\n',dsurfq);

fprintf('T barocl - integer devisor of 86400\n');
if round(c0)~=c0
  fprintf('T barocl is not an integer devisor of 86400\n');
  t1=86400/round(c0);
  t2=86400/round(c0+1);
  fprintf('Suggested T bcl= %9.3f, %9.3f\n',t1,t2);
  error('CHange T barocl');
else
  fprintf('Passed\n');
end

fprintf('\nT barocl and T barotrop\n')
if round(c2)~=c2
  fprintf('Barocl. & barotr. time steps do not match\n');
  t1=2*round(c2)*tbtr;
  t2=2*round(c2+1)*tbtr;
  fprintf('Suggested T bcl= %9.3f, %9.3f\n',t1,t2);
  error('Adjust T barocl. or T barotrop');
end

% Check coupling time
% hycom_cice.F:
%      ocn_cpl_day = OCN_nts_day/OCN_nts_cpl
%      ice_nts_cpl = ICE_nts_day/ocn_cpl_day
 
ocn_nts_cpl = round(tcpl/tbcl); % should be integer
ocn_nts_day=86400/tbcl;
ocn_cpl_day=round(ocn_nts_day/ocn_nts_cpl);
if ocn_nts_day ~= ocn_cpl_day*ocn_nts_cpl
  fprintf(' /n  Need to make tbcl to match ice coupling freq\n');
  fprintf(' ocn_nts_day ~- ocn_cpl_day*ocn_nts_cpl\n');
  error('Coupling and time steps do not match');
else
  fprintf(' \n Ice coupling matches T barocl, OK\n');
end

if dsurfq>0
fprintf('\nT barocl and dsurfq\n')
if round(c1)~=c1
  fprintf('Barocl. time step does not match dsurfq\n');
  t1=86400/round(c1)*dsurfq;
  t2=86400/round(c1+1)*dsurfq;
  fprintf('Suggested T bcl= %9.3f, %9.3f\n',t1,t2);
  error('Adjust T barocl. or dsurfq');
else
  fprintf('OK\n');
end
end


fprintf('Time steps are ok\n');


function check_time_step(tbcl,tbtr,dsurfq,tcouple);
%function check_time_step(tbcl,tbtr,dsurfq,tcouple);
% HYCOM: barocl & barotr time steps have 
% to be related and satisfy several requirements
%tbcl=120;
%tbtr=4;
%dsurfq=0.125;
%if coupling time is setup (ocea-cice, e.g.) then
% the time step should be a divisor of coupling time (e.g., 3600 sec)
% in blkdat.input: cplifq = .04166667 --> ~3600 sec 
% 

c0=86400/tbcl;
c1=86400/tbcl*dsurfq;
c2=tbcl/2*1/tbtr;

if tcouple>0 & ~isempty(tcouple)
  c3=tcouple/tbcl;
else 
  c3=[];
end

fprintf('T baroclinic = %8.3f\n',tbcl);
fprintf('T barotropic = %8.3f\n',tbtr);
fprintf('Number of days - surface diagn = %8.3f\n',dsurfq);
fprintf('Coupling time = %8.3f\n',tcouple);

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

if ~isempty(c3)
  fprintf('\n If Coupling time - Tbcl should be a divisor of  Coupling time: \n');
  if round(c3)~=c3
    fprintf('!!!  Tbcl is not an integer divisor of Coupling time, ERROR !!!\n');

    t1=(round(tbcl*10))/10;
    c3=tcouple/t1;
    dd=abs(round(c3)-c3);
    while dd>0
      t1=t1-0.1;
      c3=tcouple/t1;
%      fprintf('t1=%8.3f, c3=%8.3f\n',t1,c3);
      dd=abs(round(c3)-c3);
      if t1<=0;
        fprintf(' Could not find tbcl for tcouple=%8.3f\n',tcouple);
        break;
      end
      if dd<1e-12; break; end;
    end

    t2=(round(tbcl*10))/10;
    c3=tcouple/tbcl;
    dd=abs(round(c3)-c3);
    while dd>0
      t2=t2+0.1;
      c3=tcouple/t2;
      dd=abs(round(c3)-c3);
      if t2>tcouple;
        fprintf(' Could not find tbcl for tcouple=%8.3f\n',tcouple);
        break;
      end
      if dd<1e-12; break; end;
    end

    fprintf('Suggested T bcl= %9.3f, %9.3f\n',t1,t2);

  else
    fprintf(' OK \n');
  end
end

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

fprintf('\nT barocl and T barotrop\n')
if round(c2)~=c2
  fprintf('Barocl. & barotr. time steps do not match\n');
  t1=2*round(c2)*tbtr;
  t2=2*round(c2+1)*tbtr;
  fprintf('Suggested T bcl= %10.7d, %10.7d\n',t1,t2);
  error('Adjust T barocl. or T barotrop');
end

fprintf('Time steps are ok\n');


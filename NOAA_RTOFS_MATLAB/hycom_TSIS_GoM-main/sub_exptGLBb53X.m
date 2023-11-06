function expt = sub_exptGLBb53X(dnmb);
% HYCOM+NCODA GLBb0.08 analysis GOFS3.1
% return GLBb experiment # for
% selected date

i=1;
DE(i,1) = datenum(2007,01,1);
DE(i,2) = datenum(2007,07,01);
DE(i,3) = 534;
i=i+1;
DE(i,1) = datenum(2008,09,19);
DE(i,2) = datenum(2009,06,30);
DE(i,3) = 535;
i=i+1;
DE(i,1) = datenum(2009,07,01);
DE(i,2) = datenum(2011,7,1);
DE(i,3) = 536;
i=i+1;
DE(i,1) = datenum(2011,01,03);
DE(i,2) = datenum(2013,7,1);
DE(i,3) = 537;
i=i+1;
DE(i,1) = datenum(2013,8,20);
DE(i,2) = datenum(2014,12,31);
DE(i,3) = 538;
i=i+1;
DE(i,1) = datenum(2015,1,1);
DE(i,2) = datenum(2016,4,18);
DE(i,3) = 539;

if dnmb<DE(1,1);
  fprintf('Selected date is out of range: %s\n',datestr(dnmb));
  fprintf('  First date of analysis: %s\n',datestr(DE(1,1)));
  expt = 0;
  return
end

ii = max(find(DE(:,1)<=dnmb));
expt = DE(ii,3);

return
function expt = sub_exptGLBb(dnmb);
% HYCOM+NCODA GLBb0.08 analysis GOFS3.1
% return GLBb experiment # for
% selected date

i=1;
DE(i,1) = datenum(2008,09,19);
DE(i,2) = datenum(2009,05,06);
DE(i,3) = 906;
i=i+1;
DE(i,1) = datenum(2009,05,07);
DE(i,2) = datenum(2011,1,2);
DE(i,3) = 908;
i=i+1;
DE(i,1) = datenum(2011,01,03);
DE(i,2) = datenum(2013,8,20);
DE(i,3) = 909;
i=i+1;
DE(i,1) = datenum(2013,8,20);
DE(i,2) = datenum(2014,4,4);
DE(i,3) = 910;
i=i+1;
DE(i,1) = datenum(2014,4,5);
DE(i,2) = datenum(2016,4,18);
DE(i,3) = 911;
i=i+1;
DE(i,1) = datenum(2016,4,18);
DE(i,2) = datenum(2017,5,31);
DE(i,3) = 912;

if dnmb<DE(1,1);
  fprintf('Selected date is out of range: %s\n',datestr(dnmb));
  fprintf('  First date of analysis: %s\n',datestr(DE(1,1)));
  expt = 0;
  return
end

ii = max(find(DE(:,1)<=dnmb));
expt = DE(ii,3);

return
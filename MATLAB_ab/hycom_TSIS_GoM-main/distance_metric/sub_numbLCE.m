% Find actual number of LCEs at given day/record
% ihc - day number, current date
% LCE - LCE struct array
% LCE(Numb_LCE).XY(Nrecords).X, Y

function Nlc = sub_numbLCE(LCE,ihc);

Nlc = length(LCE);  % max # of LCEs during the forecast at any given day

Irc=[];
for inlc=1:Nlc
  Irc(inlc) = length(LCE(inlc).XY);
end

ii = find(Irc >= ihc);
Nlc = length(ii); 

ilg = [];
for inlc=1:Nlc
  ilg(inlc) = length(LCE(inlc).XY(ihc).X);
end

ii = find(ilg > 0);
Nlc = length(ii);

return
  



function [F,IFL]=sub_fill_land(F,IFL);
% Find indices of closest ocean pnts
% To Fill land in horizontal plane
% Simply grab the closest value
% Needed for interpolation
[ma,na]=size(F);

if isempty(IFL),
  IFL=struct;
  In=find(isnan(F));
  [JJ,II]=find(~isnan(F));

  fprintf('Deriving Fill Land indices ...\n');
  for ik=1:length(In);
    if mod(ik,10000)==0
      fprintf('  Done %5.2f%% ...\n',ik/length(In)*100);
    end
    
    [j0,i0]=ind2sub([ma,na],In(ik));
    d=sqrt((JJ-j0).^2+(II-i0).^2);
    kk=find(d==min(d),1);
    jm=JJ(kk);
    im=II(kk);

    IFL.Land_Indx(ik)=In(ik);
    IFL.Land_IJ(ik)=i0;
    IFL.Land_IJ(ik)=j0;
    IFL.Ocean_Indx(ik)=sub2ind([ma,na],jm,im);
    IFL.Ocean_IJ(ik,1)=im;
    IFL.Ocean_IJ(ik,2)=jm;
  end
end
%keyboard
%
% Fill land:
In=IFL.Land_Indx;
Io=IFL.Ocean_Indx;
F(In)=F(Io);

return
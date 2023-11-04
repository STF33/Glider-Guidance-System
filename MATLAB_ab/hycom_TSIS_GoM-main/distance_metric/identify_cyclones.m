      function CYCL = identify_cyclones(elon,alat,ssh,cisol);
%
%  Identify cyclonic eddies, ssh - demeaned field
%         cisol - ssh contour to identify cyclonic eddies
%          
%   Start search cyclones from ssh = -0.25 m and 
%       going up, pick up weaker cyclones
%
%  This eliminates the problem related to identification of cyclones  
%      that are  inside large-scale SSH depression and are enclosed by a contour 
%    of smaller (in magnitude) SSH
%
%  Dmitry Dukhovskoy 2017 - created for DeepStar projects

cS = -0.25;
dS = 0.05;
cE = -0.1;

f_plot=0;  %=0 - keep figure visible off

LCssh=0;

if f_plot==0
  ff=figure('Visible','off');
end

cisol = cS-dS;
LCLCE = struct;
CYCL(1).Label = 'Cyclonic Eddy';
CYCL(1).ssh_contour=[];
CYCL(1).xx=[];
CYCL(1).yy=[];
ccT=0;

%keyboard

while cisol < cE
  cisol = cisol+dS;

  [cc1,cc2]=contour(elon,alat,ssh,[cisol cisol],'k');
% Identify cyclone
  np=size(cc1,2);
  kk=1;

  if ~isempty(cc1); 
    while kk<np
%    disp(kk);
      nrd=cc1(2,kk);
      iiS=kk+1;
      iiE=kk+nrd;
%    disp(nrd);

      if nrd>20
        xx=cc1(1,iiS:iiE)';
        yy=cc1(2,iiS:iiE)';

%  Check if the contour is closed:
        dD=distance_spheric_coord(yy(1),xx(1),yy(2),xx(2));
        Dend=distance_spheric_coord(yy(1),xx(1),yy(end),xx(end));

% Check if this contour track a cyclone and not an a/cyclone (center SSH has to be <0):
        xxo=mean(xx);
        yyo=mean(yy);
        D = distance_spheric_coord(alat,elon,yyo,xxo);
        [jjo,iio]=find(D == min(min(D)));
        ssh0 = ssh(jjo,iio);
        
        if Dend<0.01*dD & ssh0<cisol, 
%
%  Check if this eddy has been already tracked:
          IN=0;
          for jj=1:ccT
            xxo=CYCL(jj).xx;
            yyo=CYCL(jj).yy;
            IN = inpolygon(mean(xxo),mean(yyo),xx,yy);
            if IN>0; break; end;
          end;
          if IN==0
            ccT=ccT+1;
            CYCL(ccT).Label='Cyclonic Eddy';
            CYCL(ccT).ssh_contour=cisol;
            CYCL(ccT).xx=xx;
            CYCL(ccT).yy=yy;
          end
        end;  % Dend
      end;  % if nrd>

      kk=iiE+1;
    end    % while
  end;    % ~isempty

end;      % while cisol


if f_plot==0
 close(ff);
end;

return


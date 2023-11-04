    function draw_ellipse(lon1,lat1,U1,V1,sclv,x1,y1,Uvct,unts);
%
% Calculate and draw STD ellipses and mean vector
% from U,V time series at location lon1, lat1
% sclv - scaling to draw mean and std ellipse
% If x1,y1 are not empty, than
% draw a vector = Uvct units/s (same units as U1,V1) at x1,y1 
% unts - 'cm/s' or 'm/s'

XYd(:,1)=U1-mean(U1);
XYd(:,2)=V1-mean(V1);
Cv=1/length(XYd)*XYd'*XYd;

% Find eigenvectors (Ev) and eigenvalues (El):
[Ev El]=eig(Cv);
e1=Ev(:,1)';
e2=Ev(:,2)';
lmb1=El(1,1);
lmb2=El(2,2);

% Major axis (Applied multivariate statistical analysis, p. 463):
mj=(lmb1)^0.5*e1;
mn=(lmb2)^0.5*e2;

ii=lon1;
jj=lat1;
mux=mean(U1);
muy=mean(V1);

%  sclv=0.01;   % scaling=# of degrees per 1cm/s of the velocity, ellipse
  cf=0.3;
  beta=15;
  col=[0 0 1];
  lwd=1.8;
  draw_arrowF(ii,ii+mux*sclv,jj,jj+muy*sclv,cf,beta,col,lwd);

%      plot([ii ii+mux],[jj jj+muy],'r');
  pp8=plot(ii,jj,'r.');
  set(pp8,'Color',col,'MarkerSize',15);
%      hold on
  col_ell=[1 0 0];   % ellipse color
  lwde=1.8;
  plot([ii+(mux-mj(1))*sclv ii+(mux+mj(1))*sclv],...
       [jj+(muy-mj(2))*sclv jj+(muy+mj(2))*sclv],...
       'Color',col_ell,'LineWidth',lwde);
  plot([ii+(mux-mn(1))*sclv ii+(mux+mn(1))*sclv],...
       [jj+(muy-mn(2))*sclv jj+(muy+mn(2))*sclv],...
       'Color',col_ell,'LineWidth',lwde);
% Plot ellipse:
  tt=(0:0.1:2*pi+0.1);
  aa=sqrt((mj*sclv)*(mj*sclv)');
  bb=sqrt((mn*sclv)*(mn*sclv)');
  A=[aa*cos(tt);bb*sin(tt)];
  alf=atan2(mj(2),mj(1));
  Rot=[cos(-alf), sin(-alf); -sin(-alf), cos(-alf)];
  A2=Rot*A;
  Ell=[A2(1,:)+ii+mux*sclv;A2(2,:)+jj+muy*sclv];
  plot(Ell(1,:),Ell(2,:),'Color',col_ell,'LineWidth',lwde)


if ~isempty(x1)
  x2=x1+Uvct*sclv;
  draw_arrowF(x1,x2,y1,y1,cf,beta,col,lwd);
  if (round(Uvct)==Uvct)
    ustr=sprintf('%i %s',Uvct,unts);
  else
    ustr=sprintf('%3.1f %s',Uvct,unts);
  end
  text(x1,y1+0.02,ustr,'Fontsize',12);
end





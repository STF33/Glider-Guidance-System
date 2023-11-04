    function sub_plot_section(A,dH,fld,SCT,sttl0,f_layer,offst);
% Plot section, specified in SCT - str. array
% if it is empty - default: OBs
% Data - T, or S field and dH - hybrid vertical grid
pthtopo = '/Net/ocean/ddmitry/HYCOM/GoM/topo/GOMl0.04/';

rg = 9806;
Tv= 72; % topography version

if isempty(offst),
  offst=0;
end


%f_layer=1;  % plot isopyc.
plr=15;
%fld = 'salin';
%fld = 'temp';
%fld = 'u-vel.';
%fld = 'v-vel.';

ftopo = sprintf('%s/depth_GOMl0.04_%2.2i.nc',pthtopo,Tv); % 
HH  = nc_varget(ftopo,'Bathymetry');
%LON = nc_varget(ftopo,'Longitude');
%LAT = nc_varget(ftopo,'Latitude');
[mm,nn]=size(HH);

% Sections:
if isempty(SCT)
  SCT = sub_sections('ob');
end;
Nsc=length(SCT);

switch(lower(fld))
 case('salin');
  c1=34.2;
  c2=37;
  
  CMP = colormap_A(c1,c2);
  cmp=CMP.colormap;
  cnt=CMP.intervals;
  nint=length(cmp);
 case('temp');
  c1=2;
  c2=30;

  CMP = colormap_A(c1,c2);
  cmp=CMP.colormap;
  cnt=CMP.intervals;
  nint=length(cmp);
 case {'u','v'};
  c1=-0.5;
  c2=0.5;
  cl1=flipud(colormap_cold(50));
  cl2=colormap_red(50);
  cmp=[cl1;cl2];
  cmp=smooth_colormap(cmp,5,3);
  nint=length(cmp);
  cnt=[c1:(c2-c1)/nint:c2];
 otherwise
  fprintf('UNKNOWN input field %s',fld);
end

%CMP = colormap_A(c1,c2);
%cmp=CMP.colormap;
%cnt=CMP.intervals;
%nint=length(cmp);

dH(dH<0.01)=0;

[l,mm,nn]=size(dH);

for ii=1:Nsc
  IND=SCT(ii).IND;
  Dsec=squeeze(dH(:,IND));
  Dsec(Dsec==0)=nan;
% Create Depth array of interface depths:
% Note these are BOTTOM interfaces 
% So Layer 1 is between interfaces 0m and ZZ(1)
  clear ZZb
  Dsec(isnan(Dsec))=0;
  ZZb(1,:)=-Dsec(1,:);
  for kk=2:l
    ZZb(kk,:)=ZZb(kk-1,:)-Dsec(kk,:);
  end
  I0=find(Dsec(1,:)==0);
  for ka=1:length(I0)
    i0=I0(ka);
    ZZb(:,i0)=nan;
  end
  
% For plotting need to have #of layers +1 interfaces
% add surface layer, otherwise all values
% will be shifted 1 layer down
  [nl,npb]=size(ZZb);
  ZZ=zeros(nl+1,npb);
  ZZ(2:nl+1,:)=ZZb;
  hb=HH(IND);

  
%  x=SCT(ii).X;
  x=[1:npb];
  [XX,dmm]=meshgrid(x,[1:l]);
  rgnm=SCT(ii).Name;
  
  figure(ii+offst);
  clf
  axes('position',[0.11 0.2 0.8 0.7]);

  Asec=squeeze(A(:,IND));
  Asec=[Asec(1,:);Asec(:,:)];
  XX=[XX(1,:);XX(:,:)];

  In=find(~isnan(Asec));
  zb=min(ZZ(In));

  
  pcolor(XX,ZZ,Asec); shading flat;
  caxis([c1 c2]);
  hold on;
  colormap(cmp);

  if f_layer>0
    nl=size(ZZ,1);
    for k=1:nl
      z=ZZ(k,:);
      x=XX(k,:);
      plot(x,z,'k-','linewidth',1);
      if (k==plr+1), % bottom interface of the layer
	plot(x,z,'w-','linewidth',1.8);
      end
    end
  end;

  set(gca,'Color',[0 0 0],'tickdir','out');
  set(gca,'xlim',[min(min(XX)) max(max(XX))],...
	  'ylim',[1.12*zb 0],...
	  'xminortick','on')
%	  'ytick',[-4500:500:0]);
  set(gca,'fontsize',16);

  sttl=sprintf('%s %s %s',sttl0,fld,rgnm);
  title(sttl,'Fontsize',14);
  %  lt0=mean(LAT(j1,i1:i2));



  hght=[];
  lngth=[];
  mint=20;
  mbx=20;
  fsz=14;
  bxc='k';
  posc=[0.08 0.08 0.8 0.05];
  aend=1;
  [az,axc]  = colorbar_horiz(cmp,cnt,hght,lngth,mint,fsz,bxc,posc,mbx,aend);

end


return





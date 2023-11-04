    function SCT = sub_sections(sname);
% Prepare Structured array with grid indices
% for specified sections
pthtopo = '/Net/ocean/ddmitry/HYCOM/GoM/topo/GOMl0.04/';
Tv=72;
ftopo = sprintf('%s/depth_GOMl0.04_%2.2i.nc',pthtopo,Tv); % 
HH  = nc_varget(ftopo,'Bathymetry');
[mm,nn]=size(HH);

switch(lower(sname));
 case('ob');  % section along open boundaries
% Sections:
  SCT(1).I=[250,nn];
  SCT(1).J=[2,2];
  SCT(2).I=[nn-1,nn-1];
  SCT(2).J=[2,mm];
  SCT(3).I=[430,nn];
  SCT(3).J=[mm-1,mm-1];

  %figure(1); clf;
  %contour(HH,[0 0],'k');
  %hold on;

  Nsc=length(SCT);
  for ii=1:Nsc
    i1=SCT(ii).I(1);
    i2=SCT(ii).I(2);
    j1=SCT(ii).J(1);
    j2=SCT(ii).J(2);
    I=[i1:i2]';
    J=[j1:j2]';

    switch(ii),
     case(1)
      SCT(ii).Name='South OB';
     case(2)
      SCT(ii).Name='East OB';
     case(3)
      SCT(ii).Name='North OB';
    end

    if j2==j1
      di=i2-i1+1;
      J=ones(di,1)*j1;
      SCT(ii).X=I;
    end
    if i2==i1
      dj=j2-j1+1;
      I=ones(dj,1)*i1;
      SCT(ii).X=J;
    end

    Il=sub2ind([mm,nn],J,I);
    SCT(ii).IND=Il;

  %  plot(I,J,'r.');
  end
 case('ob2');  % section along open boundaries
% Sections:  step in the domain N pnts away from relaxation
  N=2;
  SCT(1).I=[250,nn];
  SCT(1).J=[10,10];
  SCT(2).I=[nn-N,nn-N];
  SCT(2).J=[2,mm];
  SCT(3).I=[430,nn];
  SCT(3).J=[mm-N,mm-N];

  %figure(1); clf;
  %contour(HH,[0 0],'k');
  %hold on;

  Nsc=length(SCT);
  for ii=1:Nsc
    i1=SCT(ii).I(1);
    i2=SCT(ii).I(2);
    j1=SCT(ii).J(1);
    j2=SCT(ii).J(2);
    I=[i1:i2]';
    J=[j1:j2]';

    switch(ii),
     case(1)
      SCT(ii).Name='South OB';
     case(2)
      SCT(ii).Name='East OB';
     case(3)
      SCT(ii).Name='North OB';
    end

    if j2==j1
      di=i2-i1+1;
      J=ones(di,1)*j1;
      SCT(ii).X=I;
    end
    if i2==i1
      dj=j2-j1+1;
      I=ones(dj,1)*i1;
      SCT(ii).X=J;
    end

    Il=sub2ind([mm,nn],J,I);
    SCT(ii).IND=Il;

  %  plot(I,J,'r.');
  end

 case('centr_gom');
  SCT(1).I=[309,309];
  SCT(1).J=[2,335];
  SCT(2).I=[7,539];
  SCT(2).J=[147,147];

  %figure(1); clf;
  %contour(HH,[0 0],'k');
  %hold on;

  Nsc=length(SCT);
  for ii=1:Nsc
    i1=SCT(ii).I(1);
    i2=SCT(ii).I(2);
    j1=SCT(ii).J(1);
    j2=SCT(ii).J(2);
    I=[i1:i2]';
    J=[j1:j2]';

    switch(ii),
     case(1)
      SCT(ii).Name='Carib-Yuc-PCity';
     case(2)
      SCT(ii).Name='Mexico-Campeche-FlStrait-Bahamas';
    end

    if j2==j1
      di=i2-i1+1;
      J=ones(di,1)*j1;
      SCT(ii).X=I;
    end
    if i2==i1
      dj=j2-j1+1;
      I=ones(dj,1)*i1;
      SCT(ii).X=J;
    end

    Il=sub2ind([mm,nn],J,I);
    SCT(ii).IND=Il;

  %  plot(I,J,'r.');
  end
  
end

return
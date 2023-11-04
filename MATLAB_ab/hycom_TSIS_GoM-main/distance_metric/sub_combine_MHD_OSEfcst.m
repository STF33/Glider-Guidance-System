% Pool together all OSE f/casts with PIES/noPIES IC 
function MHD = sub_combine_MHD_OSEfcst;
pthmat  = '/Net/mars/ddmitry/hycom/hycom_TSIS/data_mat/';
ixx=0;
for isim=1:2
  if isim==1
    esim = 'PIES';
  else
    esim = 'noPIES';
  end

  for YR=2009:2010
    im1=5;
    if YR==2010; im1=1; end;
    for imo=im1:12
      fmat = sprintf('%sMHD_LCLCE_501GOMu_persist_OSEfcst%s_%4.4i%2.2i.mat',...
                    pthmat,esim,YR,imo);

      fprintf('Loading %s\n',fmat);
      A = load(fmat);

      ixx=ixx+1;
      DV=datevec(A.TM);
      MHD(ixx).mhd = A.MHD(:,1);
      MHD(ixx).Name = sprintf('Fcst %s %4.4i%2.2i',esim,YR,imo);
      MHD(ixx).Date_str = sprintf('%2.2i/%2.2i/%4.4i',DV(1,3:-1:1));
      MHD(ixx).esim = esim;
      MHD(ixx).Time = A.TM;

% Persistence
      ixx=ixx+1;
      MHD(ixx).mhd = A.MHD(:,2);
      MHD(ixx).Name = sprintf('Persist %s %4.4i%2.2i',esim,YR,imo);
      MHD(ixx).Date_str = sprintf('%2.2i/%2.2i/%4.4i',DV(1,3:-1:1));
      MHD(ixx).esim = sprintf('Persist %s',esim);
      MHD(ixx).Time = A.TM;
    end
  end
end



return 

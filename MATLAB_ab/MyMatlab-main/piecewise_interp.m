function Pns = piecewise_interp(Xi,Yi,dn,Xnew);
% Piecewise continuous polynomial
% Use Newton interpolating polynomial
% Given mesh point Xi and values Yi
% fit a piecewise polynomial of degree dn
% and interpolate original data
% onto a higher-resolution mesh Xnew
% Output Pns - function at Xnew points
%
% Dmitry Dukhovskoy, FSU, 2018-2020
%
fprintf('picewise_interpolation: Newton Polynomial, degree=%i\n',dn);
np = length(Xi);
%dn = 5; % set the degree of the piecewise polynomial
is = dn+1; % # of points in the segment
           % note that last and 1st pnts are overlapped
	   % between adjacent segments
nsgm = floor((np-1)/(is-1)); % # of full segments
np_res = np-(nsgm*(is-1)+1); % left points

res_sgm = 0;
dn_res = 0;
if np_res>0
  res_sgm = 1; % need to fit a final piece
  is_res = np_res; % degree of last polynom (+1 pnt from previous segm)
  dn_res = is_res; % degree of the "residual" polynomial
end

fprintf('# of full segments: %i + %i points\n',nsgm,np_res);
% Create points where to evaluate the polynomial:
%dx = 0.01;
%xx = [x0:dx:xN];
xx=Xnew;

Pns = [];  % piecewise polynomial evaluated at xx points
xns = [];
for isgm =1:nsgm
%  fprintf('\nisgm=%i\n',isgm);
  
  i1 = (isgm-1)*(is-1)+1; % 1st pnt of segment
  i2 = i1+is-1;           % last pnt of segment
  Xis = Xi(i1:i2);
  Yis = Yi(i1:i2);

  Fdf = sub_newton_bern(Xis,Yis);

% Identify points inside the segment
% and evaluate polynom at xxs
  [gs,xxs] = sub_evaluate_polynom(Xis,Yis,Fdf,xx,dn);
%  keyboard
% Note: 1st pnt should be identical to last from previous segment
  if isgm>1
    gs = gs(2:end);
    xxs = xxs(2:end);
  end
  
  Pns = [Pns;gs];
  xns = [xns;xxs]; % should be identical to xx after all pieces put together
                   % if all xx pnts are within the [a,b] interval
%  xNs = 
end

% Do the last segment if needed
if np_res>0
  i1 = nsgm*(is-1)+1; % 1st pnt in the residual segment
  i2 = np;            % last pnt = end pnt
  di = i2-i1;
% Sanity check # of points in the resid. segm =
% degree of the resid. polyn + 1
  if di ~= dn_res,
    error('Check residual polynom, last segm, #pnts and degree');
  end
  
  Xis = Xi(i1:i2);
  Yis = Yi(i1:i2);

  Fdf = sub_newton_bern(Xis,Yis);
% Identify points inside the segment
% and evaluate polynom at xxs
  [gs,xxs] = sub_evaluate_polynom(Xis,Yis,Fdf,xx,dn_res);
% Chop off 1st pnt
  gs = gs(2:end);
  xxs = xxs(2:end);
  Pns = [Pns;gs];
  xns = [xns;xxs]; % should be identical to xx after all pieces put together
  
end


f_chck=0;
if f_chck==1
  figure(1); clf;
  plot(Xi,Yi,'o-');
  hold on;
  plot(xns,Pns,'r');
  stt = sprintf('Piecwise Newton Polyn.(red), degr=%i, last segm dgr=%i',...
	dn,dn_res);
  title(stt);
end

return


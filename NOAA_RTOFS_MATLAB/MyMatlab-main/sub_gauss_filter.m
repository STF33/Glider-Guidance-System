function Flf = sub_gauss_filter(Fld,sgmx,sgmy,npnts,ilm1,ilm2,jlm1,jlm2);
% Apply 2D gauss Filter
% sgm1 and smg2 defines filtering scale
% fill land with nans
% npnts - # of points in the filter in 1 direction
%         usually, npnts ~=3*sgmx
% Construct filter:
[mm,nn]=size(Fld);

dx = round(npnts/2);
dy = dx;
x=[-dx:dx];
y=[-dy:dy];
nx=length(x);
ny=length(y);
x0=mean(x);
y0=mean(y);

[X,Y]=meshgrid(x,y);

fprintf('Filtering 2D Gaussian: sgmx=%6.4f, sgmy=%6.4f, Npnts=%3.1f\n',...
	sgmx,sgmy,nx);


aa = (X-x0).^2/(2*sgmx^2);
bb = (Y-y0).^2/(2*sgmy^2);

F = 1/(2*pi*sgmx*sgmy).*exp(-(aa+bb));
%keyboard

% Check weight
% should be 1, adjust if needed
w=nansum(nansum(F));
F=F/w;
%keyboard

ntt = ((ilm2-ilm1)+1)*((jlm2-jlm1)+1);
nnmd= round(ntt/10);
cc = 0;
Flf = Fld*nan;
for ii=ilm1:ilm2
  if ii<=dx,    cc=cc+(jlm2-jlm1+1); continue; end;
  if ii>nn-dx, cc=cc+(jlm2-jlm1+1); continue; end;
  for jj=jlm1:jlm2
    cc=cc+1;
    if mod(cc,nnmd)==0,
      fprintf('... done %6.4f%%\n', cc/ntt*100);
    end
    
    if jj<=dy, continue; end;
    if jj>mm-dy; continue; end;
    if isnan(Fld(jj,ii)), continue; end;
    A=Fld(jj-dy:jj+dy,ii-dx:ii+dx);
    I=find(~isnan(A));
    IN=find(isnan(A));
    if length(I)<0.01*nx*ny; continue; end;
    b=nanmean(nanmean(A));
    A(IN)=b;
    Af = nansum(nansum(A.*F));
    Flf(jj,ii)=Af;
  end
end

      


return

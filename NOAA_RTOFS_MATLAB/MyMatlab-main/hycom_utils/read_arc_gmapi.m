function [xmap,ymap] = read_arc_gmapi(IDM,JDM,fina);
% The gmapi file is a copy of GLBa0.0X to ARCc0.0X, 
% since the actual grid does not change between GLBa0.0X
% and GLBb0.0X (and ARCc0.0X and ARCb0.0X)."
% These are remapping indices GLBa-> ARCc
%
% Input IDM, JDM - X, Y dimensions of the 
% ARCc0.0X grid

IJDM = IDM*JDM;
npad=4096-mod(IJDM,4096);

Lrec=(IJDM+npad)*4; % 1 record length, bytes

fid  = fopen(fina,'r','ieee-be');

xmap = fread(fid,IJDM,'float32');
xmap = reshape(xmap,[IDM JDM])';
fseek(fid,4*(npad+IJDM),-1);
ymap = fread(fid,IJDM,'float32');
ymap = reshape(ymap,[IDM JDM])';

fclose(fid);

mnx = min(min(xmap));
mxx = max(max(xmap));
mny = min(min(ymap));
mxy = max(max(ymap));

fprintf('Reading gmapi %s\n',fina);
fprintf('xmap:  min,max = %8.4f  %8.4f\n',mnx,mxx);
fprintf('ymap:  min,max = %8.4f  %8.4f\n',mny,mxy);


return


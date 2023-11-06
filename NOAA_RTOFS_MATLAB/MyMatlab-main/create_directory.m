  function create_directory(pthfig);
%
%  function create_directory(pthfig);
%  
% Creates all missing directories from root
% analogus to tcsh: mkdir -p ...

S=strread(deblank(pthfig(2:end)),'%s','delimiter','/');
pp='/';
%keyboard
for k=1:length(S);
  pp=[pp,char(S(k)),'/'];
  if ~exist(pp,'dir');
    fprintf('Creating new directory: %s\n',pp);
%    REX = sprintf('!mkdir %s',pp);
%    eval(REX);
    REX = sprintf('mkdir %s',pp);
    system(REX);
  end;
end;




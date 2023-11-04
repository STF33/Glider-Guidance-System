function sub_add_fields2archv(fldnm,Fdmm,fida_NEW,fidb_NEW,aa,toto);
% Add missing fields in archive files 
% Both *a and *b files are already open for writing
% and are at the position where the record is needed
% to match HYCOM reading protocol
% for the nest files
% Mean or instanteneous archive files
% used for nesting have different 
% fields in *a and and *b 
% file identifiers fida_* fidb_* - open files for writing
% Fdmm is a 2D field (JD x ID)
hg=2^100;
sz=size(Fdmm);
if length(sz)~=2,
  error('sub_add_fields2archv: Array Fdmm has to be 2D');
end
JD=sz(1);
ID=sz(2);
IJDM=ID*JD;
Fdmm=Fdmm';
dmm=reshape(Fdmm,IJDM,1);
I=find(isnan(dmm));
dmm(I)=hg;
Inh=find(dmm<hg/10);
minv=min(dmm(Inh));
maxv=max(dmm(Inh));

fwrite(fida_NEW,dmm,'float32','ieee-be');
fwrite(fida_NEW,toto,'float32','ieee-be');

% Parse string from old nest *b:
I=strfind(aa,'=');
aa(1:I-1)=' ';
aa(1:length(fldnm))=fldnm;
aa2 = sub_parse_string_B(aa,[],[],minv,maxv);
fprintf('<=  NO FIELD -----  \n');
fprintf('=> %s\n',aa2);
fprintf(fidb_NEW,[aa2,'\n']);

return
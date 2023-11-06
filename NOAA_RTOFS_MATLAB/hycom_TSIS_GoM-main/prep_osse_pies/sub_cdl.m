function sub_cdl(fnmb,pthout,ntime);
% Cread cdl file for writing netcdf
sp='                                        ';
%
fcdl=sprintf('%s%s.cdl',pthout,fnmb);
fid = fopen(fcdl,'wt');
fprintf(fid,'%s \n',sprintf('netcdf %s {',fnmb));
fprintf(fid,'%s \n','dimensions:');
fprintf(fid,'%sdepth = 75 ;\n',sp(1:8));
fprintf(fid,'%stime = UNLIMITED ; // (%i currently)\n',...
          sp(1:8),ntime);
fprintf(fid,'%s \n','variables:');
fprintf(fid,'%sfloat depth(depth) ;\n',sp(1:8));
%fprintf(fid,'\n');
% Variable attributes:
str=sprintf('\"Depth below Surface of Measurements (data) or ADCP bin\" ;');
fprintf(fid,'%sdepth:long_name = %s\n',sp(1:16),str);
fprintf(fid,'%sdepth:units = \"m\" ;\n',sp(1:16));
fprintf(fid,'%sdepth:positive = \"down\" ;\n',sp(1:16));
fprintf(fid,'%sdepth:axis = \"Z\" ;\n',sp(1:16));
fprintf(fid,'%sdouble time(time) ;\n',sp(1:8));
fprintf(fid,'%stime:standard_name = \"time\" ;\n',sp(1:16));
fprintf(fid,'%stime:long_name = \"Time\" ;\n',sp(1:16));
str=sprintf('\"minutes since 2009-04-21 12:00:00\" ;');
fprintf(fid,'%stime:units = %s\n',sp(1:16),str);
fprintf(fid,'%stime:calendar = \"standard\" ;\n',sp(1:16));
fprintf(fid,'%sfloat T_var(time, depth) ;\n',sp(1:8));
fprintf(fid,'%sT_var:long_name = \"TEMPERATURE (WATER)\" ;\n',sp(1:16));
fprintf(fid,'%sT_var:units = \"degC\" ; \n',sp(1:16));
fprintf(fid,'%sT_var:_FillValue = 999.f ;\n',sp(1:16));
fprintf(fid,'%sT_var:missing_value = 999.f ;\n',sp(1:16));
fprintf(fid,'%sfloat S_var(time, depth) ;\n',sp(1:8));
fprintf(fid,'%sS_var:long_name = \"SALINITY\" ;\n',sp(1:16));
fprintf(fid,'%sS_var:units = \"PSU\" ; \n',sp(1:16));
fprintf(fid,'%sS_var:_FillValue = 999.f ;\n',sp(1:16));
fprintf(fid,'%sS_var:missing_value = 999.f ;\n',sp(1:16));
fprintf(fid,' \n');
fprintf(fid,'// global attributes:\n');
fprintf(fid,'%s:Info1 = \"Synthetic PIES T/S profiles from 1km NEMO\" ;\n',...
              sp(1:16));
fprintf(fid,'%s:Info2 = \"NEMO simulations: CICESE \";\n',sp(1:16));
fprintf(fid,'%s:Info3 = \"Synthetic profiles prepared by NCSU \";\n',sp(1:16));
fprintf(fid,'%s:history = \"%s: created from ugos_mooring_2011.mat\" ;\n',...
       sp(1:16),date);
fprintf(fid,'%s:Time_Zone = \"GMT\" ;\n',sp(1:16));
fprintf(fid,'%s:Date_Time_Format = \"yyyy-mm-dd hh:mm:ss\" ;\n',sp(1:16));
fprintf(fid,'%s:Start_Time = \"2011-01-01 12:00:00\" ;\n',sp(1:16));
fprintf(fid,'%s:Stop_Time = \"2011-12-31 23:59:00\" ;\n',sp(1:16));
fprintf(fid,'%s:Type = \"Scalar\" ;\n',sp(1:16));
fprintf(fid,'%s:Reference =  \"COAPS FSU\" ;\n',sp(1:16));
fprintf(fid,'%s:PROJECT = \"NAS Loop Current Predictability\" ;\n',sp(1:16));
fprintf(fid,'}');
fclose(fid);


stt=sprintf('system(''ncgen -b %s'')',fcdl);
eval(stt);
stt=sprintf('system(''mv %s.nc %s.'')',fnmb,pthout);
eval(stt);
stt=sprintf('system(''rm -f %s'')',fcdl);
eval(stt);

return

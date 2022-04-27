% /* --------------------------------------------------------------------------------------
%  * File:    MigDescriptor2movit.m
%  * Date:    01/11/2018
%  * Author:  David Pastor Escuredo, research@dpastoresc.org
%  * Version: 0.2
%  * License: BSD
%  * --------------------------------------------------------------------------------------
%  Copyright (c) 2015-2019, David Pastor Escuredo

function writeHeaderMigration_v12(fid, Tave,Xave, maxp, minp, bins)

%WRiting the metada in a file
%TODO: convert this format to JSON

fprintf(fid, '%s', ['Tave ' num2str(Tave)]);
fprintf(fid, '\n');
fprintf(fid, '%s', ['Xave ' num2str(Xave)]);
fprintf(fid, '\n');
fprintf(fid, '%s', ['cid;cnum;tstep;speed_m;speed_x;speed_y;speed_z;speedave_m;speedave_x;speedave_y;speedave_z;speedr_m;speedr_x;speedr_y;speedr_z;sl;meanslx;meansly;meanslz;meanslm;neighs;neighgain;celldiameter;convergence']);
fprintf(fid, '\n');

fprintf(fid, '%s', ['cellid;cellnum;timestep;maxspeed;maxspeedx;maxspeedy;maxspeedz;maxspeedave;maxspeedavex;maxspeedavey;maxspeedavez;' ...
    'maxspeedr;maxspeedrx;maxspeedry;maxspeedrz;maxsl;maxmeanslx;maxmeansly;maxmeanslz;maxmeanslm;maxeneigh;maxneighgain;maxcelldiam;maxconverg;']);
fprintf(fid, '\n');
fprintf(fid, '%s',['-1;-1;-1;' num2str(maxp(1)) ';-1;-1;-1;' num2str(maxp(2)) ';-1;-1;-1;' num2str(maxp(3)) ';-1;-1;-1;' ...
    num2str(maxp(4)) ';' num2str(maxp(5)) ';' num2str(maxp(6)) ';' num2str(maxp(7)) ';' num2str(maxp(8)) ...
    ';' num2str(maxp(9)) ';' num2str(maxp(10)) ';' num2str(maxp(11)) ';' num2str(maxp(12))]);
fprintf(fid, '\n');
fprintf(fid, '%s', ['cellid;cellnum;timestep;minspeed;minspeedx;minspeedy;minspeedz;minspeedave;minspeedavex;minspeedavey;minspeedavez;' ...
    'minspeedr;minspeedrx;minspeedry;minspeedrz;minsl;minmeanslx;minmeansly;minmeanslz;minmeanslm;mineneigh;minneighgain;mincelldiam;minconverg;']);
fprintf(fid, '\n');
fprintf(fid, '%s',['-1;-1;-1;' num2str(minp(1)) ';-1;-1;-1;' num2str(minp(2)) ';-1;-1;-1;' num2str(minp(3)) ';-1;-1;-1;' ...
    num2str(minp(4)) ';' num2str(minp(5)) ';' num2str(minp(6)) ';' num2str(minp(7)) ';' num2str(minp(8)) ...
    ';' num2str(minp(9)) ';' num2str(minp(10)) ';' num2str(minp(11)) ';' num2str(minp(12))]);
fprintf(fid, '\n');
fprintf(fid, '%s', 'nodim;nodim;secs;nodim;nodim;nodim;umsec;nodim;nodim;nodim;umsec;nodim;nodim;nodim;nodim;nodim;nodim;nodim;nodim;nodim;nodim;nodim;nodim;nodim');
fprintf(fid, '\n');
fprintf(fid, '%s', 'numberOfBins');
fprintf(fid, '\n');
fprintf(fid, '%s', num2str(bins));
fprintf(fid, '\n');
fclose(fid);
        
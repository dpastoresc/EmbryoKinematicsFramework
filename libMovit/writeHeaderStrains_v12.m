function writeHeaderStrains_v12(fid, Tave,Xave, maxp, minp, bins)
% David Pastor Escuredo. 2012/2015 BIT-UPM
% Tracking Kinematics Framework
% (C) All rights reserved

%Writing the metada in a file
%TODO: convert this format to JSON

fprintf(fid, '%s', ['Tave ' num2str(Tave)]);
fprintf(fid, '\n');
fprintf(fid, '%s', ['Xave ' num2str(Xave)]);
fprintf(fid, '\n');
fprintf(fid, '%s', ['cid;cnum;tstep;Qs;e1_m;e1_x;e1_y;e1_z;e2_m;e2_x;e2_y;e2_z;e3_m;e3_x;e3_y;e3_z;Qd;DistortionRatio;MaxShearAngle;de1_m;de1_x;de1_y;de1_z;de2_m;de2_x;de2_y;de2_z;de3_m;de3_x;de3_y;de3_z']);
fprintf(fid, '\n');

fprintf(fid, '%s', ['cellid;cellnum;timestep;maxqs;maxe1;maxe1x;maxe1y;maxe1z;maxe2;maxe2x;maxe2y;maxe2z;maxe3;maxe3x;maxe3y;maxe3z;maxqd;maxDR;maxAngle;maxe1;maxe1x;maxe1y;maxe1z;maxe2;maxe2x;maxe2y;maxe2z;maxe3;maxe3x;maxe3y;maxe3z']);
fprintf(fid, '\n');
fprintf(fid, '%s',['-1;-1;-1;' num2str(maxp(1)) ';' num2str(maxp(2)) ';-1;-1;-1;' num2str(maxp(3)) ';-1;-1;-1;' ...
    num2str(maxp(4)) ';-1;-1;-1;' num2str(maxp(5)) ';' num2str(maxp(6)) ';' num2str(maxp(7)) ';' ...
    num2str(maxp(8)) ';-1;-1;-1;' num2str(maxp(9)) ';-1;-1;-1;' num2str(maxp(10)) ';-1;-1;-1']);
fprintf(fid, '\n');
fprintf(fid, '%s', ['cellid;cellnum;timestep;minqs;mine1;mine1x;mine1y;mine1z;mine2;mine2x;mine2y;mine2z;mine3;mine3x;mine3y;mine3z;minqd;minDR;minAngle;minde1;minde1x;min1y;minde1z;minde2;minde2x;minde2y;minde2z;minde3;minde3x;minde3y;minde3z']);
fprintf(fid, '\n');
fprintf(fid, '%s',['-1;-1;-1;' num2str(minp(1)) ';' num2str(minp(2)) ';-1;-1;-1;' num2str(minp(3)) ';-1;-1;-1;' ...
    num2str(minp(4)) ';-1;-1;-1;' num2str(minp(5)) ';' num2str(minp(6)) ';' num2str(minp(7)) ';' ...
    num2str(minp(8)) ';-1;-1;-1;' num2str(minp(9)) ';-1;-1;-1;' num2str(minp(10)) ';-1;-1;-1']);
fprintf(fid, '\n');
fprintf(fid, '%s', ['nodim;nodim;secs;nodim;um2;nodim;nodim;nodim;um2;nodim;nodim;nodim;um2;nodim;nodim;nodim;' ...
    'nodim;nodim;nodim;um2;nodim;nodim;nodim;um2;nodim;nodim;nodim;um2;nodim;nodim;nodim']);
fprintf(fid, '\n');
fprintf(fid, '%s', 'numberOfBins');
fprintf(fid, '\n');
fprintf(fid, '%s', num2str(bins));
fprintf(fid, '\n');
fprintf(fid, '%s', 'W');
fprintf(fid, '\n');
fprintf(fid, '%s', num2str(1));
fprintf(fid, '\n');
fclose(fid);

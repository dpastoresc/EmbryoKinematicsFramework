function writeHeaderTopo_v12(fid, Tave,Xave, maxp, minp, bins)

% David Pastor Escuredo. 2012/2015 BIT-UPM
% Tracking Kinematics Framework
% (C) All rights reserved

%Writing the metada in a file

ms=255;
mis=0;%labels for estabilization 50-100-150-200
fprintf(fid, '%s', ['Tave ' num2str(Tave)]);
fprintf(fid, '\n');
fprintf(fid, '%s', ['Xave ' num2str(Xave)]);
fprintf(fid, '\n');
%We create the header with max and min global values for all timesteps
fprintf(fid, '%s', ['cid;cnum;tstep;Speed;AveragedSpeed_m;AveragedSpeed_x;AveragedSpeed_y;AveragedSpeed_z;Topology;P;Q;RotationDiscriminant;PlanarExpansion_m;PlanarExpansion_x;PlanarExpansion_y;PlanarExpansion_z;PlanarCompression_m;PlanarCompression_x;PlanarCompression_y;PlanarCompression_z;Rotation_m;Rotation_x;Rotation_y;Rotation_z']);
fprintf(fid, '\n');
%We don't need to put the max and min of the coordinates of the vector field!!
fprintf(fid, '%s', ['cellid;cellnum;timestep;maxspeed;maxspeedave;maxspeedx;maxspeedy;maxspeedz;maxS;maxP;maxQ;maxD;maxexp;maxexpx;maxexpy;maxexpz;maxcomp;maxcompx;maxcompy;maxcompz;maxrot;maxrotx;maxroty;maxrotz']);
fprintf(fid, '\n');
fprintf(fid, '%s',['-1;-1;-1;' num2str(maxp(1)) ';' num2str(maxp(2)) ';-1;-1;-1;' num2str(ms) ';' num2str(maxp(3)) ';' num2str(maxp(4)) ';' num2str(maxp(5)) ';' ...
    num2str(maxp(6)) ';-1;-1;-1;'  num2str(maxp(7)) ';-1;-1;-1;' num2str(maxp(8)) ';-1;-1;-1']);
%fprintf(fid, '%s',['-1;-1;-1;' num2str(ms) ';' num2str(maxp(1)) ';' num2str(maxp(2)) ';' num2str(maxp(3)) ';' num2str(maxp(4)) ';' num2str(maxp(5)) ';-1;-1;-1;'  num2str(maxp(6)) ';-1;-1;-1;' num2str(maxp(7)) ';-1;-1;-1;' num2str(maxp(8))]);
fprintf(fid, '\n');
fprintf(fid, '%s', ['cellid;cellnum;timestep;minspeed;maxspeedave;maxspeedx;maxspeedy;maxspeedz;minS;minP;minQ;minD;minexp;minexpx;minexpy;minexpz;mincomp;mincompx;mincompy;mincompz;minrot;maxrotx;minroty;minrotz']);
fprintf(fid, '\n');
fprintf(fid, '%s',['-1;-1;-1;' num2str(minp(1)) ';' num2str(minp(2)) ';-1;-1;-1;' num2str(mis) ';' num2str(minp(3)) ';' num2str(minp(4)) ';' ...
    num2str(minp(5)) ';' num2str(minp(6)) ';-1;-1;-1;'  num2str(minp(7)) ';-1;-1;-1;' num2str(minp(8)) ';-1;-1;-1' ]);  %fprintf(fid, '%s',['-1;-1;-1;' num2str(mis) ';' num2str(minp(4)) ';' num2str(mipn) ';' num2str(miqp) ';' num2str(miqn) ';' num2str(mie) ';-1;-1;-1;'  num2str(mic) ';-1;-1;-1;' num2str(mir) ';-1;-1;-1;' num2str(mirr) ]);
fprintf(fid, '\n');
fprintf(fid, '%s', 'nodim;nodim;secs;um/sec;um/sec;nodim;nodim;nodim;nodim;nodim;nodim;nodim;nodim;nodim;nodim;nodim;nodim;nodim;nodim;nodim;degrees;nodim;nodim;nodim');
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

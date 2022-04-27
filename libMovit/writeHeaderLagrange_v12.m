function writeHeaderLagrange_v12(fid, Tave,Xave, maxp, minp, bins, T)

% David Pastor Escuredo. 2012/2015 BIT-UPM
% Tracking Kinematics Framework
% (C) All rights reserved

%Writing the metada in a file
%TODO: convert this format to JSON along with the data so it can describe
%itself

        ms=255;
        mis=0;%labels for estabilization 50-100-150-200
        fprintf(fid, '%s', ['Tave ' num2str(Tave)]);
        fprintf(fid, '\n');
        fprintf(fid, '%s', ['Xave ' num2str(Xave)]);
        fprintf(fid, '\n');
        fprintf(fid, '%s', ['cid;cnum;tstep;FTLE;isocoricFTLE;Velocity;Topology;J;MIC1;MIC2;Rotation_m;rotation_x;rotation_y;rotation_z;' ...
            'c1_m;c1_x;c1_y;c1_z;c2_m;c2_x;c2_y;c2_z;c3_m;c3_x;c3_y;c3_z;MaxShear;IntermediateShear']);
        fprintf(fid, '\n');      
        fprintf(fid, '%s', ['cid;cnum;tstep;maxFTLE;maxisoFTLE;maxAveVelocity;maxs;maxJ;maxMIC1;maxMIC2;maxrotation_m;maxrotation_x;maxrotation_y;maxrotation_z;' ...
                'maxe1_m;maxe1_x;maxe1_y;maxe1_z;maxe2_m;maxe2_x;maxe2_y;maxe2_z;maxe3_m;maxe3_x;maxe3_y;maxe3_z;maxMaxShear;maxIntermediateShear']);
        fprintf(fid, '\n');
        maxp
        minp
        fprintf(fid, '%s',['-1;-1;-1;' num2str(maxp(1)) ';' num2str(maxp(2)) ';' num2str(maxp(3)) ';' num2str(ms) ';' num2str(maxp(5)) ';' num2str(maxp(6)) ';' num2str(maxp(7)) ';' num2str(maxp(8)) ...
            ';-1;-1;-1;' num2str(maxp(9)) ';-1;-1;-1;' num2str(maxp(10)) ';-1;-1;-1;' num2str(maxp(11)) ';-1;-1;-1;' num2str(maxp(12)) ';' num2str(maxp(13))]);
        fprintf(fid, '\n');
        fprintf(fid, '%s', ['cid;cnum;tstep;minFTLE;minisoFTLE;minAveVelocity;mins;minJ;minMIC1;minMIC2;minrotation_m;minrotation_x;minrotation_y;minrotation_z;' ...
                'mine1_m;mine1_x;mine1_y;mine1_z;mine2_m;mine2_x;mine2_y;mine2_z;mine3_m;mine3_x;mine3_y;mine3_z;minMaxShear;minIntermediateShear']);
        fprintf(fid, '\n');
        fprintf(fid, '%s',['-1;-1;-1;' num2str(minp(1)) ';' num2str(minp(2)) ';' num2str(minp(3)) ';' num2str(mis) ';' num2str(minp(5)) ';' num2str(minp(6)) ';' num2str(minp(7)) ';' num2str(minp(8)) ...
            ';-1;-1;-1;' num2str(minp(9)) ';-1;-1;-1;' num2str(minp(10)) ';-1;-1;-1;' num2str(minp(11)) ';-1;-1;-1;' num2str(minp(12)) ';' num2str(minp(13))]);
        fprintf(fid, '\n');
        fprintf(fid, '%s', 'nodim;nodim;secs;nodim;nodim;micronsxmin;nodim;nodim;nodim;nodim;nodim;nodim;nodim;nodim;nodim;nodim;nodim;nodim;nodim;nodim;nodim;nodim;nodim;nodim;nodim;nodim;nodim');
        fprintf(fid, '\n');
        fprintf(fid, '%s', 'numberOfBins');
        fprintf(fid, '\n');
        fprintf(fid, '%s', num2str(bins));
        fprintf(fid, '\n');
        fprintf(fid, '%s', 'W');
        fprintf(fid, '\n');
        fprintf(fid, '%s', num2str(T));
        fprintf(fid, '\n');
        fclose(fid);  

function writeHeaderTopo_v12(fid, Tave,Xave, maxp, minp, bins)

% David Pastor Escuredo. 2012/2015 BIT-UPM
% Tracking Kinematics Framework
% (C) All rights reserved

%WRiting the metada in a file
%TODO: convert this format to JSON

        ms=255;
        mis=0;%labels for estabilization 50-100-150-200
        fprintf(fid, '%s', ['Tave ' num2str(Tave)]);
        fprintf(fid, '\n');
        fprintf(fid, '%s', ['Xave ' num2str(Xave)]);
        fprintf(fid, '\n');
        %We create the header with max and min global values for all timesteps
        fprintf(fid, '%s', ['cid;cnum;tstep;mode1;mode2;mode3;mode4;mode5;mode6']);
        fprintf(fid, '\n');
        %We don't need to put the max and min of the coordinates of the vector field!!
        fprintf(fid, '%s', ['cellid;cellnum;timestep;maxmode1;maxmode2;maxmode3;maxmode4;maxmode5;maxmode6']);
        fprintf(fid, '\n');
        fprintf(fid, '%s',['-1;-1;-1;' num2str(maxp(1)) ';' num2str(maxp(2)) ';' num2str(maxp(3)) ';' num2str(maxp(4)) ';' num2str(maxp(5)) ';' ...
               num2str(maxp(6))]);
        %fprintf(fid, '%s',['-1;-1;-1;' num2str(ms) ';' num2str(maxp(1)) ';' num2str(maxp(2)) ';' num2str(maxp(3)) ';' num2str(maxp(4)) ';' num2str(maxp(5)) ';-1;-1;-1;'  num2str(maxp(6)) ';-1;-1;-1;' num2str(maxp(7)) ';-1;-1;-1;' num2str(maxp(8))]);
        fprintf(fid, '\n');
        fprintf(fid, '%s', ['cellid;cellnum;timestep;minmode1;minmode2;minmode3;minmode4;minmode5;minmode6']);
        fprintf(fid, '\n');
        fprintf(fid, '%s',['-1;-1;-1;' num2str(minp(1)) ';' num2str(minp(2)) ';-1;-1;-1;' num2str(minp(3)) ';' num2str(minp(4)) ';' ...
                num2str(minp(5)) ';' num2str(minp(6)) ]);  %fprintf(fid, '%s',['-1;-1;-1;' num2str(mis) ';' num2str(minp(4)) ';' num2str(mipn) ';' num2str(miqp) ';' num2str(miqn) ';' num2str(mie) ';-1;-1;-1;'  num2str(mic) ';-1;-1;-1;' num2str(mir) ';-1;-1;-1;' num2str(mirr) ]);
        fprintf(fid, '\n');
        fprintf(fid, '%s', 'nodim;nodim;secs;um/sec;um/sec;nodim;nodim;nodim;nodim;nodim;nodim');
        fprintf(fid, '%s', 'numberOfBins');
        fprintf(fid, '\n');
        fprintf(fid, '%s', num2str(bins));
        fprintf(fid, '\n');
        fprintf(fid, '%s', 'W');
        fprintf(fid, '\n');
        fprintf(fid, '%s', num2str(1));
        fprintf(fid, '\n');
        fclose(fid);    
        
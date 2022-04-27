function readBioemergencesDetection(outputfile, vtkCenters, t)
    % David Pastor Escuredo 2012
    % Read selection csv to generate seed for tracking
    
    V = readVTKCenters(vtkCenters);
    numOfCent=size(V,1);

   % V(1:5,:)
    count=1;
    fid=fopen(outputfile, 'w');
    for i=1:numOfCent
        fprintf(fid, '%s', [num2str(count) ';' num2str(t) ';' num2str(V(i,1)) ';' num2str(V(i,2)) ';' num2str(V(i,3)) ';0'] );
        fprintf(fid, '\n');
        count=count+1;
    end
    count
    fclose(fid);

end


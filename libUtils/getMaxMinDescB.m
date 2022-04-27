% David Pastor Escuredo 2012
% Mechanics Framework
% (C) All rights reserved

function [maxDesc minDesc]= getMaxMinDesc(values, maxPercentil, minPercentil, lowerLimitVal, upperLimitVal)
 
    if nargin<5
        upperLimitVal=9999999999999999999;
    end
    if nargin<4
        lowerLimitVal=-1;%-9999999999999999999
    end
    if nargin<3
        minPercentil=1;
    end
    if nargin<2
        maxPercentil=99;
    end
    
    % returns a matrix with the max and min for descriptios (order by
    % columns)
    
    %lowerLimitVal and upperLimitVal put the ranges where the values should
    %be considered.
    
    %put exclusion values too
    for i=1:size(values,2)
        value=values(:,i);
        iv=find(value>lowerLimitVal & value<upperLimitVal);
        value=value(iv);
        %TODO: exclusion list!!!
        
        maxDesc(i)=prctile(value,maxPercentil);
        minDesc(i)=prctile(value,minPercentil);

    end    

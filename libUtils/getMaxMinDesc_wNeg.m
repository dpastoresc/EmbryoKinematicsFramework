function [maxDesc minDesc]= getMaxMinDesc_wNeg(values, maxPercentil, minPercentil)
    % Mechanics Framework
    % returns a matrix with the max and min for descriptios (order by
    % columns)
    
    %-1 means not to be processed and returned
    %go through columns
    for i=1:size(values,2)
        value=values(:,i);
        maxDesc(i)=prctile(value,maxPercentil);
        minDesc(i)=prctile(value,minPercentil);

    end    

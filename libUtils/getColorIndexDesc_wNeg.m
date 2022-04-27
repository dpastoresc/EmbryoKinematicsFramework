%DEPRECATED! not use
function cValue= getColorIndexDesc_wNeg(value, maxVal, minVal, mapLevels)
    % David Pastor Escuredo 2012
    % Mechanics Framework
    % returns a matrix of indexes for each value from 1 to mapLevels range
   
    %-1 means not to be processed and returned
    levels=maxVal-minVal;    
    
    %Calculatin the value index regarding the map and the number of
    %levels
    %%100-0 -> 255-0 -> 256-1
    %size(value);
    %this goes between 0 and levels-1
    %for 256 bins: 0-255
    cValue=fix(double((value-minVal)./double(levels)).*(mapLevels-1));%0-based
    %+1 because of matlab
    %cValue=fix(double((value-minVal)./double(levels)).*(mapLevels-1))+1;1-based
    
    %some checking of values, shouldn't be, buuut...
    %check max value
    im=find(cValue>(mapLevels-1));
    cValue(im)=mapLevels-1;
    %check min value
    imm=find(cValue<0);
    cValue(imm)=0;     
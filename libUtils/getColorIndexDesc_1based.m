%DEPRECATED!! not use
function cValue= getColorIndexDesc_1based(value, maxVal, minVal, mapLevels)
    % David Pastor Escuredo 2012
    % Mechanics Framework
    % returns a matrix of indexes for each value from 1 to mapLevels range
   
    %-1 means not to be processed and returned
    levels=maxVal-minVal;
    minVal;
    maxVal;

    %Calculatin the value index regarding the map and the number of
    %levels
    %%100-0 -> 255-0 -> 256-1
    size(value);
    %this goes between 0 and levels-1
    %for 256 bins: 0-255
    %cValue=fix(double((value-minVal)./double(levels)).*(mapLevels-1));%0-based
    %+1 because of matlab
    cValue=fix(double((value-minVal)./double(levels)).*(mapLevels-1))+1;%1-based
    
    %some checking of values, shouldn't be, buuut...
    %check max value
    %check this if i should put it 0 and mapLevels-1
    im=find(cValue>(mapLevels));
    cValue(im)=mapLevels;
    %check min value
    imm=find(cValue<1);
    cValue(imm)=1;
    
    %All values that were -1 should be kept like that
    cValue(value==-1)=-1;
  
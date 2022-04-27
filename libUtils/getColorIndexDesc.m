% David Pastor Escuredo 2012
% Mechanics Framework
% (C) All rights reserved

%FIX!!! ... if it is from 0->N and we have exceptVal=-1, the range is done
    %with -1->N even when we exclude it after... so we are giving it some
    %offset.
    %If we use a lot of bins, this is not important though
    
%ToDo: Add physical values range/threshold.

function cValue= getColorIndexDesc(value, maxVal, minVal, mapLevels, exceptVal, indexBase)

    if nargin<6
        %0 (general) or 1 (matlab)
        %255-0 -> 256-1
        indexBase=0;
    end
    if nargin<5
        %value to exclude
        exceptVal=-1;
    end
    %indexBase
    %exceptVal
    %pause
    
    levels=maxVal-minVal;        
    cValue=fix(double((value-minVal)./double(levels)).*(mapLevels-1))+indexBase;
    
    %security checks to force the range of the colormap
    %check max value
    im=find(cValue>(mapLevels-1+indexBase));
    cValue(im)=mapLevels-1+indexBase;
    %check min value
    imm=find(cValue<indexBase);
    cValue(imm)=indexBase;
    
    %All values that were -1 should be kept like that
    cValue(value==exceptVal)=exceptVal;
  
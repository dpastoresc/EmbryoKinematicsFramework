function [A im]=m_removeOutlayers(A, maxp, minp, exceptValue, minSamples)
    
    if nargin<5 
        minSamples=2;
    end
    if nargin<4 
        exceptValue=-1;
    end
    if nargin<3 
        minp=5;
    end
    if nargin<2 
        maxp=95;
    end
    mx=prctile(A, maxp);
    mn=prctile(A, minp);
    im=find(A<mx & A>mn & A~=exceptValue & ~isnan(A));
    
    %if we dont have enoug samples out of this removal, we dont apply it
    if length(im)>=minSamples
       A=A(im);
    else
       im=find(A~=exceptValue); 
       %length(im)
       if length(im)>0
        A=A(im);
       end    
    end
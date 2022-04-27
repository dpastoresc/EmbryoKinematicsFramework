% /* --------------------------------------------------------------------------------------
%  * File:    orderEigValues_n1.m
%  * Date:    01/06/2016
%  * Author:  David Pastor Escuredo, david.pastor@upm.es
%  * Version: 0.2
%  * License: BSD
%  * --------------------------------------------------------------------------------------
%  Copyright (c) 2015, David Pastor Escuredo - david.pastor@upm.es
%  All rights reserved.%Mechanics Framework


function evalues=orderEigValues_n1(evalin)

aeig=evalin;
aeig(evalin<1)=1/evalin(evalin<1);
aeig2=sort(aeig);
%evalues=zeros(size(evalin));

e1=max(aeig2);
ie1=find(aeig==e1);
evalues(1:length(ie1))=evalin(ie1);
iothers=find(aeig~=e1);

if size(iothers,2)>0
    rest=aeig(iothers);
    e2=max(rest);
    ie2=find(aeig==e2);
    offset=length(ie1)+1;
    evalues(1,offset:offset+length(ie2)-1)=evalin(ie2);    
    offset=size(evalues,2);
    
    if offset<3
       e3=rest(find(rest~=e2));
       offset=offset+1;
       ie3=find(aeig==e3);
       evalues(offset:offset+length(e3)-1)=evalin(ie3);
    end
end

%evalues=zeros(size(evalin));
%evalues(1,1)=evalin(find(aeig)==aeig2(1,1));
%evalues(1,2)=evalin(find(aeig)==aeig2(1,2));
%evalues(1,3)=evalin(find(aeig)==aeig2(1,3));
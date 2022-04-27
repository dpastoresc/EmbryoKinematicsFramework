% /* --------------------------------------------------------------------------------------
%  * File:    removeDiscarded.m
%  * Date:    01/06/2016
%  * Author:  David Pastor Escuredo, david.pastor@upm.es
%  * Version: 0.2
%  * License: BSD
%  * --------------------------------------------------------------------------------------
%  Copyright (c) 2015, David Pastor Escuredo - david.pastor@upm.es
%  All rights reserved.
function [A isel]=removeDiscarded(A, exceptValue)
    
    isel=find(A~=exceptValue & ~isnan(A));
    A=A(isel);
   
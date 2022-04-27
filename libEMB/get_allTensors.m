% /* --------------------------------------------------------------------------------------
%  * File:    get_allTensors.m
%  * Date:    01/11/2018
%  * Author:  David Pastor Escuredo, research@dpastoresc.org
%  * Version: 0.2
%  * License: BSD
%  * --------------------------------------------------------------------------------------
%  Copyright (c) 2015-2019, David Pastor Escuredo

% Code taken from Benoit Lombardot's framwork

function [Eall Rall Fall] = get_allTensors(Fall)


Eall=cell(size(Fall));
Rall=cell(size(Fall));
Eerrall=cell(size(Fall));

I3 = diag(ones(1,3));
F = zeros(3,3);
Ferr = zeros(3,3);

if iscell(Fall)
    if ~isempty(Fall)>0
        for t=1:length(Fall)

            if size(Fall{t},1)==23
                Eall{t} = zeros(9,size(Fall{t},2));
                %Eerrall{t} = zeros(9,size(Fall{t},2));
                Rall{t} = zeros(9,size(Fall{t},2));
                %Fall2{t} = zeros(9,size(Fall{t},2));
                
                for j=1:size(Fall{t},2)
                    F(:) = Fall{t}(1:9,j);
                    if sum(abs(F(:)))>0
                        E =  1/2*(F*F'-I3);
                        Eall{t}(:,j) = E(:);
                        %Ferr(:) = meanMovErr{t}(j)^2*Fall{t}(15:23,j); % covar mat du bruit sur les param de F
                        Fall{t}(15:23,j) = Ferr(:);
                        %Eerr = F*Ferr*F'; % covar mat du bruit sur les params de E
                        %Eerrall{t}(:,j) = Eerr(:);
                        %passage dans la base de diag de C
                        C = 2*E+I3;% construction de C
                        [P,val] = eig(C); % diagonalisation de C
                        rC = sqrt(val); % racine de rC
                        rC = P*rC*P'; %retour dans la base de départ
                        R = F*inv(rC);
                        Rall{t}(:,j) = R(:);
                        %Fall2{t}(:,j) = F(:);
                    end
                end
                Fall{t} = Fall{t}(1:9,:);
            end
        end
    end
end
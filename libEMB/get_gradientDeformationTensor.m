%David Pastor Escuredo and Benoit Lombardot (CNRS, ISC-PIF, France)
%Mechanics Framework 2014
%(C) All rights reserved


function [A res resRel C] = get_gradientDeformationTensor(pos0, dpos0, pos, dpos, w, cas)

% M = [X;Y;Z];
% m = [x;y;z];
% V = [Vx;Vy;Vz];
% x = X+Vx;
% A : matrice recherchée, c'est le tenseur gradient
% mi-m0 = A(Mi-M0)
% xi-x0 = [(d(Vx)/dX)(M0) , (d(Vx)/dY)(M0) , (d(Vx)/dz)(M0) ] * (Mi-M0) ;
% xi-x0 = A(1,:)*(Mi-M0);

% si n points 
% on a un système de n équations pour chaque variable(x,y,z)
% au sens des moindres carrés on a :
% A(1,:) = (xi-x0)*(Mi-M0)'*((Mi-M0)*(Mi-M0)')^(-1);
% plus généralement:
% A = (mi-m0)*(Mi-M0)'*((Mi-M0)*(Mi-M0)')^(-1);

% on peut associé un poid wi à chaque point Mi
% w = [w1,w2,...,wn];
% W = diag(w);
% mi-m0 = A*(Mi-M0);
% A  = (mi-m0)*W*(Mi-M0)'*((Mi-M0)*W*(Mi-M0)')^(-1);
% on remplace (Mi-M0)' par W*(Mi-M0)'

% on peut ajouter un biais à la position initiale
% xi-x0 = A(1,1:3)*(Mi-M0) + A(1,4);
% on remplace Mi-M0 par [Mi-M0;1];
% mi-m0 = A*[Mi-M0;ones(1,n)];
% A = (mi-m0)*(Mi_-M0_)'*((Mi_-M0_)*(Mi_-M0_)')^(-1);

% rq : les point les plus éloigné prenne un poid plus important
%      amplitude proportionnelle à la distance
%      donc amplitude de l'erreur proportionelle à la distance entre M et M0
%      donc dans l'idéale on respecte mieux les hypothèses du théorème de
%      Gauss Markov (les moindres carré donne la solution optimale ...) si
%      on normalise les M_moins_M0

% pour plus de détail article wikipedia "Linear least squares"

precision = 'single';
M_moins_M0 = zeros(size(pos),precision);
m_moins_m0 = zeros(size(pos),precision);
for k=1:3
    M_moins_M0(k,:) = pos(k,:)-pos0(k);
    m_moins_m0(k,:) = M_moins_M0(k,:) + dpos(k,:)-dpos0(k);
end
% % M_moins_M0_Norm = sum(M_moins_M0.^2,1).^0.5;
% % M_moins_M0 = M_moins_M0./(ones(3,1)*max(M_moins_M0_Norm,0.1));
% % m_moins_m0 = m_moins_m0./(ones(3,1)*max(M_moins_M0_Norm,0.1));

if strcmp(cas, 'std' )    
    B = M_moins_M0';
elseif strcmp(cas, 'weight' )
    % on remplace (Mi-M0)' par W*(Mi-M0)'
    % w : n lignes, 1 colonne;
    B = (w*ones(1,3)).*M_moins_M0';
elseif strcmp(cas, 'shift' )
    % on remplace Mi-M0 par [Mi-M0;1];
    M_moins_M0 = [M_moins_M0; ones(1,size(pos,2))];
    B = M_moins_M0';
elseif strcmp(cas, 'shiftweight' )
    M_moins_M0 = [M_moins_M0; ones(1,size(pos,2))];
    B = (w*ones(1,4)).*(M_moins_M0');    
end

C = M_moins_M0 * B;
if det(C)>1e-8
    C = inv( C );
    A = m_moins_m0*B*C;
    if strcmp(cas, 'std' )
        % residual = |dx-dx^| ;
        % resRel = |dx-dX| ;
        res = sum((m_moins_m0 - A * M_moins_M0).^2,1);
        resRel = sum(sum((m_moins_m0 - M_moins_M0).^2,1)).^0.5/length(res);
        res = sum(res).^0.5 / length(res);
        %resRel = sum(resRel).^0.5 / length(res);
    elseif strcmp(cas, 'weight' )
        res = sum((m_moins_m0 - A * M_moins_M0).^2,1);
        resRel = res./sum((m_moins_m0 - M_moins_M0).^2,1);
        res = sum(w.*res).^0.5 / sum(w);
        resRel = sum(w.*resRel).^0.5 / sum(w);
    elseif strcmp(cas, 'shift' )
        res = sum((m_moins_m0 - A * M_moins_M0).^2,1);
        resRel = res./sum((m_moins_m0 - M_moins_M0).^2,1);
        res = sum(res).^0.5 / length(res);
        resRel = sum(resRel).^0.5 / length(res);
    elseif strcmp(cas, 'shiftweight' )
        res = sum((m_moins_m0 - A * M_moins_M0).^2,1);
        resRel = res./sum((m_moins_m0 - M_moins_M0).^2,1);
        res = sum(w.*res).^0.5 / sum(w);
        resRel = sum(w.*resRel).^0.5 / sum(w);
    end
else
    disp(['singular matrix. rq: ' int2str(size(M_moins_M0,2)) ' points'])
    A = zeros(3,3,precision);
    res=0;
    resRel = 0;
    C = zeros(3,3,precision);
end




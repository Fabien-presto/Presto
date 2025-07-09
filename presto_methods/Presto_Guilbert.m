% This file is part of Presto-HF, a matlab toolbox to identify a circuit from
% its response.
%
% SPDX-License-Identifier: AGPL-3.0-or-later
%
% Copyright 2025 by
%   Centre Inria de l'Université Côte d'Azur
%   2004, route des Lucioles
%   06902 Sophia Antipolis Cedex
%
% and by
%   Mines Paris - Armines
%   60, boulevard Saint-Michel
%   75006 Paris
%
% Contributors: Fabien Seyfert, Jean-Paul Marmorat, Martine Olivi
%
% Presto-HF is free software: you can redistribute it and/or modify it under
% the terms of the GNU Affero General Public License as published by the Free
% Software Foundation, either version 3 of the License, or (at your option) any
% later version.
%
% Presto-HF is distributed in the hope that it will be useful, but WITHOUT ANY
% WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
% A PARTICULAR PURPOSE. See the GNU Affero General Public License for more
% details.
%
% You should have received a copy of the GNU Affero General Public License
% along with Presto-HF. If not, see <https://www.gnu.org/licenses/>.
%
function [sys]=Guilbert(num,den)
% Magic numbers 
% Test of symmetrie S12=S21
sym_prec=0.3;
% Test of rank 1 residues
rank_eps=0.01;

Poles=roots(den);
dden=Polyder(den);
B = [];
NP = length(Poles);
% Boucle sur les p�les
for ip = 1:NP
    % Construction des r�sidus
    beta = Poles(ip);
    P = Poles(find( (1:NP) ~= ip ));
    for lt = 1:2
        for ct = 1:2
            q=poly(P);
            R(lt,ct) = polyval(num{lt,ct},beta)/polyval(dden,beta);
        end
    end
    % Test de symm�trie
    nsym = norm(R-transpose(R));
    if nsym > sym_prec 
        mesg = sprintf('\n Non symmetric residue: %d\n',ip);
        disp(R);
        error(mesg);
    end
    % Test du rang
    if rank(R,rank_eps) > 1
        mesg = sprintf('\n Non rank 1 residue: %d\n',ip);
        disp(R);
        disp(ip);
        disp(det(R));
        error(mesg);
    end
    % Factorisation des residus
    a = sqrt(R(1,1));
    % Force residue to be symmetric while remainig of rank 1
    r12corr=sqrt(R(1,2)*R(2,1));
    if max(abs(r12corr-R(1,2)),abs(r12corr-R(2,1)))>max(abs(-r12corr-R(1,2)),abs(-r12corr-R(2,1)))
        r12corr=-r12corr;
    end
    b= r12corr/a;
    %b = sqrt(R(2,2));
    %err1 = abs (a*b-R(1,2));
    %err2 = abs (-a*b-R(1,2));
    %if (err2<=err1)
    %    b=-b;
    %end;
        
    V(1,1)=a;
    V(1,2)=b;
    B = [ B ; V  ];
end

sys.a = diag(Poles);
sys.b = B;
sys.c = transpose(B);
sys.d = zeros(2,2);
for k=1:2
    for j=1:2
        if (length(num{k,j})==length(den))
            sys.d(k,j) = num{k,j}(1)/den(1);
        else
            if (length(num{k,j})>length(den))
                error('Non proper transfert function');
            end
        end
    end
end
% Force D symmetrie 
sys.d=(sys.d+sys.d.')/2;
    

function [pp]=Polyder(p)

pp=polyder(p);
k=length(pp);
n=length(p)-1;
pp=[zeros(1,n-k),pp]; 

    

function [output] = OutOfSampleInversion(ve,H)
% OutOfSampleInversion Computes the ECG inverse solution appliying out of
% sample technique according with Generalization and Regularization for Inverse Cardiac Estimators

%  Copyright (C) 2021  F. M. Melgarejo-Meseguer, E. Everss-Villalba, M. Gutierrez-Fernandez-Calvillo, 
%                      S. Munoz-Romero, F.-J. Gimeno-Blanes, A. Garcia-Alberola, and J.-L. Rojo-Alvarez
%
%    This program is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation version 3.
%
%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.

%    You should have received a copy of the GNU General Public License
%    along with this program.  If not, see <http://www.gnu.org/licenses/>.

% %%%%%%%%% INPUT %%%%%%%%%%%%%
% ve -> (NxS matrix) Voltage measured on torso.
% H ->  (NxM matrix) Direct problem transfer matrix. 
% %%%%%%%%% Output %%%%%%%%%%%%%
% output -> structure 
% output.veEst -> (NxS matrix) Estimated voltage on torso.
% output.vmEst -> (MxS matrix) Estimated voltage on cavitiy.
% output.gamma -> Optimun gamma parameter selected.
% output.time -> processing time in seconds.
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Set FOT parameters and initialization

R = eye(size(H,2)); % FOT regularization matrix.
NG = 30;
mySeed = 19; % for reproducibility
gamma = logspace(-10,0,NG); % set of possible gammas. 
prc = 0.87; % Training signals percentage.

errorTrainVe = zeros(size(gamma));
errorValVe = zeros(size(gamma));
%% Out of sample process
tic
[Ht,~,idxST,idxSV,~,~,Rt,~,Htv] = DescomponH(H,R,prc,mySeed); % H and R matrix decomposition.

HH = H'*H;
RR = R'*R;
HHt = Ht'*Ht;
RRt = Rt'*Rt;
vet = ve(idxST,:); % Training Set.
vev = ve(idxSV,:); % Validation Set.
div = (Ht' * vet);

for g = 1:length(gamma)
    vmEstt =  (HHt + gamma(g) * RRt)\ div; % Regularized least squared
    veEstt = Ht*vmEstt;
    veEstv = Htv*vmEstt;
    errorTrainVe(g) = sqrt(mean((vet(:) - veEstt(:)).^2));
    errorValVe(g) = sqrt(mean((vev(:) - veEstv(:)).^2));
end

%% Computed potentials taking the best gamma given by the method
gamIdxOpt = find(errorValVe==min(errorValVe),1,'first'); % Select the lowest validation error.
div_opt = (H' * ve);
vmEst = (HH + gamma(gamIdxOpt) * RR)\ div_opt;
veEst = H*vmEst;
optGam = gamma(gamIdxOpt);
procTime= toc;

%% Output 
output.veEst = veEst;
output.vmEst = vmEst;
output.gamma = optGam;
output.time = procTime;

end


%% AUX FUNCTIONS %%
function [Ht,Hv,idxST,idxSV,idxNT,idxNV,Lt,Lv,Htv] = DescomponH(H,L,prcTrain,mySeed)
% Generamos los n# de train 
rng(mySeed)
[S,N] = size(H);
sT = floor(prcTrain*S);
nT = floor(prcTrain*N);

% Geranmos los indices de train y valid
idxS = 1:S;
idxN = 1:N;
if length(idxS) == length(idxN)
    idxST = sort(randperm(S,sT)); 
    idxSV= setdiff(idxS,idxST);
    idxNT = idxST;
    idxNV= idxSV;
else
    idxST = sort(randperm(S,sT)); 
    idxSV= setdiff(idxS,idxST);
    idxNT = sort(randperm(N,nT)); 
    idxNV= setdiff(idxN,idxNT);
end


% Creamos las matrices de entrenamiento y validacion
Ht = H(idxST,idxNT); % sensoresTrain x nodosTrain
Hv = H(idxSV,idxNV); % sensoresValid x nodosValid
Lt = L(idxNT,idxNT); % nodosTrain x nodosTrain
Lv = L(idxNV,idxNV); % nodosValid x nodosValid
Htv = H(idxSV,idxNT);%

end

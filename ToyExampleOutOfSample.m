% %%%%% OutOfSampleInversion Toy example %%%%%%

%% Simple case simulation
rng(19); % For reproducibility

% Substrate plane generation
N = 20; % number of samples by side 
xx = linspace(-1,1,N);
yy = xx;
[X,Y] = meshgrid(xx,yy);
nodes = [X(:),Y(:)];
tri = delaunayTriangulation(nodes);
nodes = [tri.Points,zeros(size(X(:)))];

% Sensor plane generation
z = 0.05;
sensors = nodes;
sensors(:,3) = sensors(:,3)+z;

% Gaussian stimulation generation
[nr,nc] = size(X);
vm = zeros(nr*nc,1);
idx = pdist2([0,0,0],nodes)<0.8;
x1 = nodes(idx,1);
y1 = nodes(idx,2);
ex1 = (x1.^2/0.1);
ey1 = (y1.^2/0.1);
vm(idx) = 1*exp(-(ex1+ey1));

% Laplacian computed for H matrix creation
k = 1;
myMesh = tri.ConnectivityList;
for t = 1:length(myMesh) 
    tri1 = [nodes(myMesh(t,1),:);nodes(myMesh(t,1),:);nodes(myMesh(t,2),:)];
    tri2 = [nodes(myMesh(t,2),:);nodes(myMesh(t,3),:);nodes(myMesh(t,3),:)];
    dist =  diag(pdist2(tri1,tri2));
    [~,idxR] = max(dist);
    if idxR == 1
        a1 = [myMesh(t,1),myMesh(t,2)];
        a2 = [myMesh(t,3),myMesh(t,3)];
    elseif idxR == 2
        a1 = [myMesh(t,1),myMesh(t,2)];
        a2 = [myMesh(t,2),myMesh(t,3)];
    else
        a1 = [myMesh(t,1),myMesh(t,1)];
        a2 = [myMesh(t,2),myMesh(t,3)];
    end
        m1(k:k+1) = a1;
        m2(k:k+1) = a2;
        k = k+2;
end
g = graph(m1,m2,'omitselfloops');
A = adjacency(g);
A = A'|A;
G = graph(A);
L = laplacian(G);
L = -L;

% H matrix computation as H = 1/R * Laplacian(Vm)

D = pdist2(sensors,nodes);  
D = 1./D; % Distance matrix.
H = D*L;
ve = D*L*vm(:); % Ve computation, H*vm;

% Out of sample processing
SNR_ve = 20; % in dB
SNR_H = 20; % in dB
ve = awgn(ve,SNR_ve,'measured');
H = awgn(H,SNR_H,'measured');

[output] = OutOfSampleInversion(ve,H);

% Plotting results
figure()
ax1 = subplot(231);
surf(X,Y,reshape(ve,N,N))
title('Real Ve')
xlabel('cm')
ylabel('cm')
zlabel('mV');

ax2 = subplot(232);
surf(X,Y,reshape(output.veEst,N,N))
title('Estimed Ve')
xlabel('cm')
ylabel('cm')
zlabel('mV');

ax3 = subplot(233);
surf(X,Y,abs(reshape(output.veEst-ve,N,N)))
title('Residuals Ve')
xlabel('cm')
ylabel('cm')
zlabel('mV');

ax4 = subplot(234);
surf(X,Y,reshape(vm,N,N))
title('Real Vm')
xlabel('cm')
ylabel('cm')
zlabel('mV');

ax5 = subplot(235);
surf(X,Y,reshape(output.vmEst,N,N))
title('Estimed Vm')
xlabel('cm')
ylabel('cm')
zlabel('mV');

ax6= subplot(236);
surf(X,Y,abs(reshape(output.vmEst-vm,N,N)))
title('Residuals Vm')
xlabel('cm')
ylabel('cm')
zlabel('mV');


linkprop([ax5 ax4 ax6],{'Zlim'});
linkprop([ax1 ax2 ax3 ax4 ax5 ax6],{'CameraPosition'});

%% Real case from EDGAR database
load('ExampleEDGARD.mat')

[output] = OutOfSampleInversion(ve,H);

%Correlation coefficient computation
corVe = zeros(size(ve,1),1);
corVm = zeros(size(vm,1),1);

for i = 1:length(corVe)
    aux = corrcoef(ve(i,:),output.veEst(i,:));
    corVe(i) = aux(1,2);
end

for i = 1:length(corVm)
    aux = corrcoef(vm(i,:),output.vmEst(i,:));
    corVm(i) = aux(1,2);
end


ax1 = subplot(221);
trisurf(TorsoTri,'FacevertexCdata',corVe);
a = colorbar;
a.Label.String = '\rho';
caxis([-1,1])
axis equal

ax2 = subplot(222);
trisurf(TorsoTri,'FacevertexCdata',sqrt(mean((output.veEst-ve).^2,2)));
a = colorbar;
a.Label.String = 'RMSE in mV';
caxis([-2,2])
axis equal

ax3 = subplot(223);
trisurf(CavTri,'FacevertexCdata',corVm);
a = colorbar;
a.Label.String = '\rho';
caxis([-1,1])
axis equal

ax4 = subplot(224);
trisurf(CavTri,'FacevertexCdata',sqrt(mean((output.vmEst-vm).^2,2)));
a = colorbar;
a.Label.String = 'RMSE in mV';
caxis([-2,2])
axis equal

linkprop([ax1 ax2 ax3 ax4],{'CameraPosition'});
%
% Worst signals
figure()

ax1 = subplot(231);
[val,idx] = min(corVe);
plot(1:492,ve(idx,:),1:492,output.veEst(idx,:))
xlabel('Samples')
ylabel('mV')
axis tight
title(['Ve \rho = ' num2str(val)]) 


ax2 = subplot(232);
[val,idx] = min(abs(corVe));
plot(1:492,ve(idx,:),1:492,output.veEst(idx,:))
xlabel('Samples')
ylabel('mV')
axis tight
title(['Ve \rho = ' num2str(val)]) 

ax3 = subplot(233);
[val,idx] = max((corVe));
plot(1:492,ve(idx,:),1:492,output.veEst(idx,:))
xlabel('Samples')
ylabel('mV')
axis tight
title(['Ve \rho = ' num2str(val)]) 


ax4 = subplot(234);
[val,idx] = min(corVm);
plot(1:492,vm(idx,:),1:492,output.vmEst(idx,:))
xlabel('Samples')
ylabel('mV')
axis tight
title(['Vm \rho = ' num2str(val)]) 


ax5 = subplot(235);
[val,idx] = min(abs(corVm));
plot(1:492,vm(idx,:),1:492,output.vmEst(idx,:))
xlabel('Samples')
ylabel('mV')
axis tight
title(['Vm \rho = ' num2str(val)]) 

ax6 = subplot(236);
[val,idx] = max((corVm));
plot(1:492,vm(idx,:),1:492,output.vmEst(idx,:))
xlabel('Samples')
ylabel('mV')
axis tight
title(['Vm \rho = ' num2str(val)]) 
legend({'Real','Estimated'})
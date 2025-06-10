clear;
close all; 
clc;

%% Initializing parameters
var_arr = [0.9 0.3];
x1 = 60; 
x2 = 60; 
nV = x1; 
nSRCS = 7; 
N = 300; 
nIter = 30; %% iterations for dictionary
TR = 1;
Dp = dctbases(N, N); % DCT basis dictionary

%% Constructing spatial and temporal sources
warning('off')
onset_intervals = [60, 35, 40, 30, 40, 10, 60];
durations = [20, 18, 20, 12, 20, 20, 10];
TC = zeros(N, length(onset_intervals)); % Preallocate matrix for efficiency

for i = 1:length(onset_intervals)
    onsets = 0:onset_intervals(i):N-20;
    TC(:, i) = generate_TC(onsets, durations(i), N);
end

% Uncomment the following if you have SPM12 added to the path (this is what we used in the published paper)
% hrf = spm_hrf(TR);

% Create HRF (hemodynamic response function)
hrf_duration = 20; % samples (~40 seconds at TR=2s)
hrf_t = (0:hrf_duration-1)';
% Double gamma HRF model (canonical HRF)
a1 = 6; a2 = 16; b1 = 1; b2 = 1; c = 1/6;
hrf = ((hrf_t/a1).^(a1) .* exp(-(hrf_t-a1)/b1) - ...
    c * (hrf_t/a2).^(a2) .* exp(-(hrf_t-a2)/b2));
hrf = hrf / max(hrf); % normalize
R  = [TC(1,2) zeros(1,length(hrf)-1)];
for i = 1:size(TC, 2)
    TC(:, i) = zscore(toeplitz(TC(:, i), R)*hrf); %
    TC(:, i) = zscore(TC(:, i) + 0.01 * wfbm(1, N)'); 
end

%% CUSTOMISE component 2 (periodic up/down spikes)
task_events = {[30, 80, 140], [60, 110, 180]}; % multiple event types
TC(:,2) = zscore(create_mixed_signal(N, task_events, [1:N], TR));

%%  Custom component-3 (narrow-peak sawtooth)
t       = (0:N-1)';         
Tper    = 48;               
duty    = 0.18;           
rawSaw = sawtooth(2*pi*t/Tper, duty); 
env = 1 + 0.25*sin(2*pi*t/120) + 0.15*sin(2*pi*t/45);
sig3 = env .* rawSaw + 0.03*randn(N,1);
TC(:,3) = zscore(sig3);

%% CUSTOMISE component 4 (periodic up/down spikes)
freq_band = [0.01, 0.1]; 
TC(:,4) = zscore(create_resting_state_signal(N, freq_band, 0.5, TR));

%% CUSTOMISE component 5 (periodic up/down spikes)
period  = 60;                    
peakLag = 10;                        
dipLag  = 30;                         
sig5    = zeros(N,1);
for start = 1:period:N
    p = start + peakLag;  if p <= N, sig5(p) =  6;  end   % positive spike
    d = start + dipLag;   if d <= N, sig5(d) = -6;  end   % negative spike
end
sig5 = smoothdata(sig5,'gaussian',5);
TC(:,5) = zscore(sig5);

%% CUSTOMISE component 3 (periodic up/down spikes)
% % mixed signal
% task_events = {[30, 150], [60, 190]}; % multiple event types
% TC(:,3) = zscore(create_mixed_signal(N, task_events, [1:N], TR));
% 

%% Spatial basis, one-or-many blobs per slice 
oldSize = 21; % hand-drawn grid size
newSizeR = x1; % e.g. 40
newSizeC = x2; % e.g. 40
scaleR = (newSizeR-1)/(oldSize-1);
scaleC = (newSizeC-1)/(oldSize-1);

origBlobs = { ...
 [ 2 9 2 9 9 ; 
 11 19 13 19 16 ];
 [ 2 11 11 20 5 ];
 [ 11 15 11 15 6;
 9 15 2 12 6;
 17 21 8 12 3 ];
 [ 7 15 10 19 4 ];
 [14 19 2 10 12;
 6 12 8 16 3 ];
 [11 19 7 17 6 ];
 [ 2 10 5 15 9 ] };
numSlices = numel(origBlobs);
tmpSM = zeros(x1,x2,numSlices);
SM = zeros(numSlices, x1*x2);
for s = 1:numSlices
    Gslice = zeros(newSizeR,newSizeC);
    
    for row = 1:size(origBlobs{s},1)
        r = origBlobs{s}(row,:);
        r1 = min(newSizeR, round((r(1)-1)*scaleR)+1);
        r2 = min(newSizeR, round((r(2)-1)*scaleR)+1);
        c1 = min(newSizeC, round((r(3)-1)*scaleC)+1);
        c2 = min(newSizeC, round((r(4)-1)*scaleC)+1);
        sig = round(r(5)*mean([scaleR scaleC]));
        
        rows = r2-r1+1; cols = c2-c1+1;
        [Z,Y] = meshgrid(linspace(-(cols-1)/2,(cols-1)/2,cols), ...
                         linspace(-(rows-1)/2,(rows-1)/2,rows));
        R2 = Z.^2 + Y.^2;
        rad = min(rows,cols)/2;
        G1 = exp(-R2/(2*sig^2));                   
        G2 = 0.6 * exp(-R2/(2*(sig*0.7)^2));       
        G3 = 0.3 * exp(-R2/(2*(sig*0.4)^2));      
        G = G1 + G2 + G3;  
        G(R2 > rad^2) = 0; 
        Gslice(r1:r2, c1:c2) = Gslice(r1:r2, c1:c2) + G;
    end
    tmpSM(:,:,s) = Gslice;
    raw_data = reshape(Gslice,1,[]);
    SM(s,:) = zscore(raw_data);
end

%% Generating data
Y = (TC+sqrt(var_arr(1))*randn(N,nSRCS))*(SM+sqrt(var_arr(2))*randn(nSRCS,x1*x2)); 
Y = Y-repmat(mean(Y),size(Y,1),1);


%% PMD
tic;
K1 = 24; K2 = 16;
lambda =0.5; gamma=lambda*ones(1,K1);
tmpX=GPower(Y,gamma,K1,'l1',0);
Y2 = Y*tmpX;
Y2 = Y2*diag(1./sqrt(sum(Y2.*Y2)));
lambda2 = 0.5;
[Wx1,~,~] = pmd_rankK(Y',Y2',K2,lambda2);
X{1} = Wx1';
D{1} = Y*X{1}';
toc

%% ACSDBE
tic
spa = 500; %500
lam = 30;
K = 12;
tStart=tic;
[D{2},X{2},~,~]= my_ACSDBE(Y,Dp(:,1:150),K,spa,lam,nIter,TC,SM); 
toc

%% LSICA
tic
K = nSRCS;
spa = 0.025; %0.015
[D{3},X{3},U,Err]=LSICA(Y,K,spa,nIter,TC,SM); %_mod
toc

%% DPCA
tic;
spa1 = 650; 
spa2 = 20;
spa3 = 80;
K = nSRCS;
[Dq, V, Xq]= svds(Y,K);
[D{4},X{4},~,~,]= my_DPCA_sim_3(Y,Dq,Xq',K,spa1,spa2,spa3,nIter,TC,SM);
toc

%% SDPCA_1
Dp = dctbases(size(Y,1),size(Y,1));
Up = Dp(:,1:150);
PSF_sim3 = spcol(augknt(linspace(0,1,25),4),4,linspace(0,1,x1))';  %20
PSF_sim4 = kron(PSF_sim3, PSF_sim3);
PSF_sim4 = abs(PSF_sim4);
for jjj=1:size(PSF_sim4,1)
    PSF_sim4(jjj,:) =(1/max(PSF_sim4(jjj,:)))*PSF_sim4(jjj,:);
end

Zp = [PSF_sim4]; %; SM  ;
Zp( all(isnan(Zp), 2), : ) = [];

tic
K = nSRCS;
spa = 40; %9 45
[Uq, ~, Zq] = svds(Y, K);
TC2 = TC*diag(1./sqrt(sum(TC.*TC))); Up = [Up TC2(:,[2 3])];
[D{5},X{5},Err,~,U]=SDPCA_1(Y,Uq,Zq',Up,Zp,spa,60,90,0.3,0.3,nIter,TC,SM); %Xp2 250,0.3,
sum(U~=0,2)
toc


%% SDPCA_2
tic
K = nSRCS;
spa = 60; %9 45
[Uq, ~, Zq] = svds(Y, K);
TC2 = TC*diag(1./sqrt(sum(TC.*TC))); 
[D{6},X{6},Err,~,U]=SDPCA_2(Y,Uq,Zq',Up,Zp,spa,60,90,0.3,0.3,nIter,TC,SM); %Xp2 0.3
sum(U~=0,2)
toc


%% Post-Processing
nA=7;
sD{1} = TC;
sX{1} = SM;
for jj =1:nA-1
    %  [sD{jj+1},sX{jj+1},ind]=sort_TSandSM_spatial(TC,SM,D{jj},X{jj},nSRCS);
    [sD{jj+1},sX{jj+1},ind]=sort_TSandSM_temporal(TC,D{jj},X{jj});
    for ii =1:nSRCS
        TCcorr(jj+1,ii,1) =abs(corr(TC(:,ii),D{jj}(:,ind(ii))));
        SMcorr(jj+1,ii,1) =abs(corr(SM(ii,:)',X{jj}(ind(ii),:)'));
    end
end
cTC(1,:) = sum(TCcorr(:,:,1)')
cSM(1,:) = sum(SMcorr(:,:,1)')
mean(cTC)
mean(cSM)
TP = zeros(nA,nSRCS); FP = zeros(nA,nSRCS); FN = zeros(nA,nSRCS); F_score = 0; thr = 0.025;
for i =1:nA-1
    for jjj=1:nSRCS
        SM_gw4(jjj,:) =SM(jjj,:)/norm(SM(jjj,:));
        [~, indd(jjj)] = max(abs(corr(abs(SM_gw4(jjj,:)'),abs(X{i}'))));
        XX{i}(indd(jjj),:) = X{i}(indd(jjj),:)/norm(X{i}(indd(jjj),:));
        TP(i,jjj) = TP(i,jjj) +sum(abs(SM_gw4(jjj,:))>=thr & abs(XX{i}(indd(jjj),:))>=thr);
        FP(i,jjj) = FP(i,jjj) +sum(abs(SM_gw4(jjj,:))<=thr & abs(XX{i}(indd(jjj),:))>=thr);
        FN(i,jjj) = FN(i,jjj) +sum(abs(SM_gw4(jjj,:))>=thr & abs(XX{i}(indd(jjj),:))<=thr);
        F_score(i,jjj) = (2*(TP(i,jjj)))/(2*(TP(i,jjj))+sum(FP(i,jjj))+(FN(i,jjj)));
    end
end
F_score = [zeros(1,nSRCS); F_score];
mean(F_score)

vec = 1:1;
ccTC = mean(TCcorr(:,:,vec),3)
ccSM = mean(SMcorr(:,:,vec),3)
mean(ccTC')
mean(ccSM')
mean(F_score')
clear DS
tau = 0.01;
for k =1:nA-1
    for i =1:nSRCS
        aa = reshape(abs(sX{1}(i,:)),nV,nV); aa = aa/norm(aa); aa(aa<=tau)=0; aa(aa>tau)=1;
        bb = reshape(abs(sX{k+1}(i,:)),nV,nV); bb = full(bb)/norm(full(bb)); bb(bb<=tau)=0; bb(bb>tau)=1;
        DS(i,k) = dice(aa,bb);
    end
end
DS = DS';
DS = [zeros(1,nSRCS); DS];
TT(1,:) = mean(ccTC');
TT(2,:) = mean(ccSM');
TT(3,:) = mean(F_score');
TT(4,:) = mean(DS');
TT

%% Plotting
f = figure; f.Position = [170 120 1600 600]; 
for i = 1:nA  
    row_base = (i-1) * 16;
    subplot_tight(nA, 16, row_base + 1, [0.015 0.01]); 
    if i == 1
        text(0.5, 0.5, '(a)', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'FontSize', 9, 'FontWeight', 'bold', 'Color', 'red');
    else
        letter = sprintf('(%s)', char('a' + i - 1));
        stats = sprintf('\nγ_{T}=%.2f γ_{S}=%.2f\nFS=%.2f DS=%.2f', TT(1,i), TT(2,i), TT(3,i), TT(4,i));
        text(0.5, 0.60, letter, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'FontSize', 9, 'FontWeight', 'bold', 'Color', 'red');
        text(0.5, 0.40, stats, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'FontSize', 9, 'FontWeight', 'bold', 'Color', 'black');
    end
    axis off;
    % Spatial maps (columns 2-8)
    for j = 1:nSRCS
        subplot_tight(nA, 16, row_base + 1 + j, [0.03 0.01]); 
        imagesc(flipdim(reshape(abs(zscore(sX{i}(j,:))), x1, x2), 1));
        colormap('hot');
        set(gca,'XTickLabel','')
        set(gca,'YTickLabel','')
    end
    % Time series (columns 9.5-15.5)
    for j = 1:nSRCS
        subplot_tight(nA, 16, row_base + 8.2 + j, [0.03 0.01]); 
        plot(zscore(sD{i}(:,j)), 'b', 'LineWidth', 0.25);
        axis tight;
        ylim([-3 3]);
        ax = gca;
        ax.XTick = [1 length(sD{i}(:,j))/2 length(sD{i}(:,j))];
        ax.XTickLabel = {'0', num2str(round(length(sD{i}(:,j))/2)), num2str(length(sD{i}(:,j))-1)};
        set(gca,'YTickLabel','')
        grid on;
    end
end

exportgraphics(gcf,'khalid2.png','Resolution',300)

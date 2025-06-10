function [U,Z,R,C,B]= SDPCA_1(Y,Uq,Zq,Up,Zp,lambda,zeta1,zeta2,tau1,tau2,nIter,TC,SM)
%  Uq,Zq,Up,Zp,

K = size(Uq,2);
Psi = eye(K,K);
Phi = zeros(K,K);
U = Uq*Psi;
Z = Phi*Zq;
Lt = size(Up,2);
Ls = size(Zp,1);
sizeU = numel(U);
sizeZ = numel(Z');

Gram_Zp = Zp * Zp';
Gram_Up = Up' * Up;

fprintf('Iteration:     ');
for j= 1:nIter
    Do = U;
    fprintf('\b\b\b\b\b%5i',j);
    F1 = U'*U; G1 = U'*Y;
    iiter = 0;
    Zpp = Z;
    
    while (iiter < nIter)
        iiter = iiter + 1;
        for i =1:K
            xk = 1.0/F1(i,i) * (G1(i,:) - F1(i,:)*Z) + Z(i,:);
            thr = lambda./abs(xk);
            Z(i,:) = sign(xk).*max(0, bsxfun(@minus,abs(xk),thr/2));
        end
        if (norm(Z - Zpp, 'fro')/sizeZ < 1e-6)
            break;
        end
        Zpp = Z;
    end
    
    for k =1:K
        Phi(k,:) = Z(k,:)/Zq;
        Z(k,:) = Phi(k,:)*Zq;
        for ii =1:length(Z(k,:))
            if abs(Z(k,ii))<1
                thr(ii) = lambda./abs(Z(k,ii));
            else
                thr(ii) = lambda./abs(1*(Z(k,ii).^(0.5)));
            end
        end
        Z(k,:) = sign(Z(k,:)).*max(0, bsxfun(@minus,abs(Z(k,:)),thr/2));
    end
    
    all_correlations = abs(Zp * Z');
    Zp_Z_products = Zp * Z';
    B = zeros(K,Ls);
    for k = 1:K
        correlations = all_correlations(:, k);
        maxCorr = max(correlations);
        threshold = tau2 * maxCorr;
        [~, bb] = sort(correlations, 'descend');
        ind = bb(1:zeta2);
        keepInd = ind(correlations(ind) >= threshold);
        B(k, :) = 0;
        if ~isempty(keepInd)
            G_keep = Gram_Zp(keepInd, keepInd);
            b_keep = Zp_Z_products(keepInd, k);
            B(k, keepInd) = (G_keep \ b_keep)';
            Z(k,:) = B(k, :)*Zp;
        end
    end
    
    F2 = Z *Z'; G2 = Y*Z';
    iiter = 0;
    Dpp = U;
    while (iiter < nIter)
        iiter = iiter + 1;
        for jjj = 1: K
            if(F2(jjj,jjj) ~= 0)
                tmp3 = 1.0/F2(jjj,jjj) * (G2(:,jjj) - U*F2(:, jjj)) + U(:,jjj);
                U(:,jjj) = tmp3/(max( norm(tmp3,2),1));
            end
        end
        if (norm(U - Dpp, 'fro')/sizeU < 1e-6)
            break;
        end
        Dpp = U;
    end

    for k =1:K
        Psi(:,k) = Uq\U(:,k);
        U(:,k) = Uq*Psi(:,k);
    end
     
    A = zeros(Lt, K);
    all_U_correlations = abs(Up' * U);
    for k = 1:K
        correlations = all_U_correlations(:, k);
        maxCorr = max(correlations);
        threshold = tau1 * maxCorr;
        [~, bb] = sort(correlations, 'descend');
        ind = bb(1:zeta1);
        keepInd = ind(correlations(ind) >= threshold);
        A(:, k) = 0;
        if ~isempty(keepInd)
            G_keep_U = Gram_Up(keepInd, keepInd);
            b_keep_U = Up(:, keepInd)' * U(:, k);
            A(keepInd, k) = G_keep_U \ b_keep_U;
            A(:,k) = A(:,k)./norm(Up*A(:,k));
            U(:, k) = Up * A(:, k);
        end
    end
    
    R(j) = sqrt(trace((U-Do)'*(U-Do)))/sqrt(trace(Do'*Do));
     
    K2 = size(TC,2);
    [~,~,ind]=sort_TSandSM_spatial(TC,SM,U,Z,K2);
    for ii =1:K2
        TCcorr(ii) =abs(corr(TC(:,ii),U(:,ind(ii))));
        SMcorr(ii) =abs(corr(SM(ii,:)',Z(ind(ii),:)'));
    end
    cTC = sum(TCcorr');
    cSM = sum(SMcorr');
    C(j) =cTC+cSM;

    
end


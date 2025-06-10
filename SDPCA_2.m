function [U,Z,R,C,B]= SDPCA_2(Y,Uq,Zq,Up,Zp,lambda,zeta1,zeta2,tau1,tau2,nIter,TC,SM)
 K = size(Uq,2);
    Psi = eye(K,K);
    Phi = zeros(K,K);
    U = Uq*Psi;
    Z = Phi*Zq;
    Lt = size(Up,2);
    Ls = size(Zp,1);
    A = zeros(Lt,K);
    B = zeros(K,Ls);
    
    Gram_Zp = Zp * Zp';
    Gram_Up = Up' * Up;
    UqTUq_inv = inv(Uq' * Uq); 
    
    fprintf('Iteration: ');
    for iter=1:nIter
        fprintf('\b\b\b\b\b%5i',iter);
        Do = U;
        
        for j =1:K
            Z(j,:) = 0; Psi(:,j) = 0; Phi(j,:) = 0; A(:,j) = 0; B(j,:) = 0;
            E = Y-U*Z;
            xk = U(:,j)'*E;
            thr = lambda./abs(xk);
            zk = sign(xk).*max(0, bsxfun(@minus,abs(xk),thr/2));
            Phi(j,:) = zk/Zq;
            zk = Phi(j,:)*Zq;
            
            for i =1:length(zk)
                if abs(zk(i))<1
                    thr(i) = lambda./abs(zk(i));
                else
                    thr(i) = lambda./abs(1*(zk(i).^(0.5)));
                end
            end
            zk = sign(zk).*max(0, bsxfun(@minus,abs(zk),thr/2));
            
            correlations = abs(Zp * zk');
            maxCorr = max(correlations);
            threshold = tau2 * maxCorr;
            [~, bb] = sort(correlations, 'descend');
            ind = bb(1:zeta2);
            keepInd = ind(correlations(ind) >= threshold);
            
            if ~isempty(keepInd)
                G_keep = Gram_Zp(keepInd, keepInd);
                b_keep = Zp(keepInd, :) * zk';
                B(j, keepInd) = (G_keep \ b_keep)';
                Z(j,:) = B(j, :)*Zp;
            end
            
            rInd = find(Z(j,:));
            if (length(rInd)>1)
                dk = E(:,rInd)*Z(j,rInd)';
                
                Psi(:,j) = UqTUq_inv * (Uq' * dk);
                uk = Uq*Psi(:,j);
                
                correlations = abs(Up' * uk);
                maxCorr = max(correlations);
                threshold = tau1 * maxCorr;
                [~, bb] = sort(correlations, 'descend');
                ind = bb(1:zeta1);
                keepInd = ind(correlations(ind) >= threshold);
                
                if ~isempty(keepInd)
                    G_keep_U = Gram_Up(keepInd, keepInd);
                    b_keep_U = Up(:, keepInd)' * uk;
                    A(keepInd, j) = G_keep_U \ b_keep_U;
                    A(:,j) = A(:,j)./norm(Up*A(:,j));
                    U(:, j) = Up * A(:, j);
                end
            end
        end
        R(iter) = sqrt(trace((U-Do)'*(U-Do)))/sqrt(trace(Do'*Do));
        
        K2 = size(TC,2);
        [~,~,ind]=sort_TSandSM_spatial(TC,SM,U,Z,K2);
        for ii =1:K2
            TCcorr(ii) =abs(corr(TC(:,ii),U(:,ind(ii))));
            SMcorr(ii) =abs(corr(SM(ii,:)',Z(ind(ii),:)'));
        end
        cTC = sum(TCcorr');
        cSM = sum(SMcorr');
        C(iter) =cTC+cSM;
        
    end
end



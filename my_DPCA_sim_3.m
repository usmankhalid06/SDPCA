function [D,X,Err,C]= my_DPCA_sim_3(Y,Dp,Xp,K,spa1,spa2,spa3,nIter,TC,SM)
    K1 = size(Dp,2);
    K2 = size(Xp,1);
    A = eye(K1,K);
    B = zeros(K,K2);
    D = Dp*A;
    X = B*Xp;
    fprintf('Iteration:     ');
    R = Y;
    gridSize = 70;
    minClusterSize = 1000; %- Minimum size of clusters to retain

%         sparseAdj = sparse_knn_adjacency(coords, k);

    for iter=1:nIter
        fprintf('\b\b\b\b\b%5i',iter);
        Do = D;
        for j =1:K
            X(j,:) = 0; A(:,j) = 0; B(j,:) = 0;
            E = R-D*X;
            
            xk = D(:,j)'*E;
            thr1 = spa1./abs(xk); 
            xkk = sign(xk).*max(0, bsxfun(@minus,abs(xk),thr1/2));
            B(j,:) = xkk*Xp'/(Xp*Xp');
            X(j,:) = B(j,:)*Xp;

            X(j,:) = firm_thresholding_nonadaptive(X(j,:), spa2/2, spa3/2);  
            
            rInd = find(X(j,:));
            if (length(rInd)>1)
                tmp3 = E(:,rInd)*X(j,rInd)'; 
                A(:,j)= (Dp'*Dp)\Dp'*tmp3;    
                A(:,j) = A(:,j)./norm(Dp*A(:,j));
                D(:,j) = Dp*A(:,j);
            end          
        end     
        Err(iter) = sqrt(trace((D-Do)'*(D-Do)))/sqrt(trace(Do'*Do));

        K2 = size(TC,2);
        [~,~,ind]=sort_TSandSM_spatial(TC,SM,D,X,K2);
        for ii =1:K2
            TCcorr(ii) =abs(corr(TC(:,ii),D(:,ind(ii))));
            SMcorr(ii) =abs(corr(SM(ii,:)',X(ind(ii),:)'));
        end
        cTC = sum(TCcorr');
        cSM = sum(SMcorr');
        C(iter) =cTC+cSM;
        

    end
end





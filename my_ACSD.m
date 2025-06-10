function [D,X,Err,C]= my_ACSD(Y,D,spa,nIter,TC,SM)
    X = zeros(size(D,2),size(Y,2)); 
%     fprintf('Iteration:     ');
    D = Y(:,1:size(D,2));  D = D*diag(1./sqrt(sum(D.*D))); 
    for iter=1:nIter
%         fprintf('\b\b\b\b\b%5i',iter);
        Dold = D;
        for j =1:size(D,2)
            X(j,:) = 0;
            E = Y-D*X;
            xk = D(:,j)'*E; 
            thr = spa./abs(xk);
            X(j,:) = sign(xk).*max(0, bsxfun(@minus,abs(xk),thr/2));
            rInd = find(X(j,:));
            if (length(rInd)<1)
                D(:,j)= randn(size(D(:,j),1), 1);
                [~,ind]= max(sum(Y-D*X.^2)); 
                D(:,j)= Y(:,ind)/norm(Y(:,ind));
            else
                D(:,j) = E(:,rInd)*X(j,rInd)'./norm(E(:,rInd)*X(j,rInd)');
            end                 
        end      
        Err(iter) = sqrt(trace((D-Dold)'*(D-Dold)))/sqrt(trace(Dold'*Dold));

        
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

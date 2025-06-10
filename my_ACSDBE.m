function [D,X,Err,C]= my_ACSDBE(Y,Dp,K,spa,lam,nIter,TC,SM)
    A = eye(size(Dp,2),K);
%     D = Dp*A;  
    D = Y(:,1:K);   D = D*diag(1./sqrt(sum(D.*D))); 
    X = zeros(size(D,2),size(Y,2)); 
%     fprintf('Iteration:     ');
    for iter=1:nIter
%         fprintf('\b\b\b\b\b%5i',iter);
        Dold = D;
        for j =1:size(D,2)
            X(j,:) = 0; A(:,j) = 0;
            E = Y-D*X;
            xk = D(:,j)'*E; 
            thr = spa./abs(xk);
            X(j,:) = sign(xk).*max(0, bsxfun(@minus,abs(xk),thr/2));
            rInd = find(X(j,:));
            if (length(rInd)<1)
%                 tmp1 = randn(size(Y,1),1);
%                 [~,bb]= sort(abs(Dp'*tmp1),'descend');
%                 ind2 = bb(1:lam);
%                 A(ind2,j)= (Dp(:,ind2)'*Dp(:,ind2))\Dp(:,ind2)'*tmp1;
%                 A(:,j)= A(:,j)./norm(Dp*A(:,j));
%                 D(:,j) = Dp*A(:,j);
            else
                tmp3 = E(:,rInd)*X(j,rInd)'; 
                [~,bb]= sort(abs(Dp'*tmp3),'descend');
                ind = bb(1:lam);
                A(ind,j)= (Dp(:,ind)'*Dp(:,ind))\Dp(:,ind)'*tmp3;
                A(:,j) = A(:,j)./norm(Dp*A(:,j));
                D(:,j) = Dp*A(:,j);
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



            

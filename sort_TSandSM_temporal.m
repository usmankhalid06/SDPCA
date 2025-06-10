function [srtd_Zt,srtd_Zs,ind_Zt]=sort_TSandSM_temporal(TC,Zt,Zs) 
   srcs = size(TC,2);
   for j=1:srcs
       [~, ind_Zt(j)]  = max(abs(corr(TC(:,j),Zt)));
       srtd_Zs(j,:) = Zs(ind_Zt(j),:);
       srtd_Zt(:,j) = sign(corr(TC(:,j),Zt(:,ind_Zt(j))))*Zt(:,ind_Zt(j));
   end  
    
end
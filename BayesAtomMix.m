function [MU,ETA,Nhat,Nhata] = BayesAtomMix( n, d, S, alpha,s,disp_int,full,trunc, nmc  )
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

MU = zeros(nmc,d);
%acc_a = zeros(nmc,d); acc_b = zeros(nmc,d);
acc = zeros(nmc,d);
mu = unifrnd(0,1,[d 1]); eta = 1./d.*ones(d,1);
Nhat = zeros(nmc,1);
ETA = zeros(nmc,d);
if trunc
    alpha0 = alpha;
    nalpha = length(alpha0);
    alpha = 0;    
    Nhata = zeros(nmc,nalpha);
end

for t=1:nmc
    
    ctr = 0;
   for j=1:d
        ctr = ctr + 1;          
        prop_muj = unifrnd(mu(j)-s(j),mu(j)+s(j));
        if prop_muj<0
            prop_muj = -prop_muj;
        end
        if prop_muj>1
            prop_muj = 1-(prop_muj-1);
        end
        prop_mu = mu; prop_mu(j) = prop_muj;
        
        if full
            prop_ll = loglik_full(prop_mu,eta);
            curr_ll = loglik_full(mu,eta);            
        else
            prop_ll = loglik(prop_mu,eta);
            curr_ll = loglik(mu,eta);
        end


        lrr = prop_ll-curr_ll;
        acc(t,j) = rand<exp(lrr);
        if acc(t,j)
            mu(j) = prop_mu(j);            
        end
   end
   
   if full
      lr = zeros(S+1,d);
   else
      lr = zeros(S,d); 
   end
   
   for j=1:d
       wts0 = zeros(d,1);
       wts0(j) = 1;
       if full
           [~,pk] = loglik_full(mu,wts0);
           pr = pk;
       else
           [~,pk] = loglik(mu,wts0);
           pr = pk(2:end)./(1-pk(1));
       end
       
       lr(:,j) = log(pr);       
   end
   lr = bsxfun(@plus,log(eta'),lr);
   lr = bsxfun(@minus,lr,max(lr,[],2));
   lr = exp(lr);
   lcprob = bsxfun(@times,lr,1./sum(lr,2));
   
   if full
       zi = zeros(S+1,d);
       for j=1:S+1
           zi(j,:) = mnrnd(n(j),lcprob(j,:));
       end
   else
       zi = zeros(S,d);
       for j=1:S
           zi(j,:) = mnrnd(n(j),lcprob(j,:));
       end
   end

   
   zcts = sum(zi,1);
   tmp = gamrnd(1./d+zcts,1);
   eta = tmp'./sum(tmp);
   
   
   if full
       [~,pr] = loglik_full(mu,eta);
       rho = 1-pr(1);
       n(1) = nbinrnd(sum(n(2:end)),rho);
       Nhat(t) = sum(n);       
   else
       [~,pr] = loglik(mu,eta);
        Nhat(t) = sum(n)./(1-pr(1));
   end
   
%    if Nhat(t)>3000
%       disp('check'); 
%    end
   
   if trunc
      for l=1:nalpha
        prca = pr_cap(mu,eta,alpha0(l));
        Nhata(t,l) = binornd(sum(n),prca);
      end
   end
   
   if mod(t,disp_int)==0
      %figure(1); plot((1:t)',Nhat(1:t)); drawnow;
      mean(acc(1:t,:),1)
      mean(Nhat(1:t))
      if trunc
         mean(Nhata(1:t,:)) 
      end
   end
   
   ETA(t,:) = eta;
   MU(t,:) = mu;
   
end
    






    function [ll,pk] =loglik(mu,wts)
        pk = zeros(S+1,1);
        for k=0:S
            pk(k+1) = nchoosek(S,k).*wts'*(mu.^k.*(1-mu).^(S-k));
        end
        
        ll = (n*log(pk(2:end))-sum(n).*log(1-pk(1)));
    end

    function [ll,pk] =loglik_full(mu,wts)
        pk = zeros(S+1,1);
        for k=0:S
            pk(k+1) = nchoosek(S,k).*wts'*(mu.^k.*(1-mu).^(S-k));
        end
        
        ll = (n*log(pk));
    end

    function prca = pr_cap(mu,wts,alpha) 
        wts = wts./sum(wts);
%         pr0 = zeros(d,1);
%         for h=1:d
%            pr0(h) = 1-(1-mu(h)).^S; 
%         end
        
        prca = sum(wts(mu>alpha));
        prca(prca<0) = 0;
        prca(prca>1) = 1;
    end

end


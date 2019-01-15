function [A,B,ETA,Nhat,Nhata] = BayesBetaMix( n, d, S, alpha,sa,sb,a0,b0,disp_int,full,trunc, nmc  )
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

A = zeros(nmc,d); B = zeros(nmc,d);
%acc_a = zeros(nmc,d); acc_b = zeros(nmc,d);
acc = zeros(nmc,d);
a = unifrnd(0,1,[d 1]); b = unifrnd(0,1,[d 1]); eta = 1./d.*ones(d,1);
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
        prop_ab = exp(normrnd([log(a(j)) log(b(j))],[sa(j) sb(j)]));
        prop_a = a; prop_b = b;
        prop_a(j) = prop_ab(1); prop_b(j) = prop_ab(2);
        if full
            prop_ll = loglik_full(prop_a,prop_b,eta,alpha);
            curr_ll = loglik_full(a,b,eta,alpha);            
        else
            prop_ll = loglik(prop_a,prop_b,eta,alpha);
            curr_ll = loglik(a,b,eta,alpha);
        end

        log_pr = (a0-1).*log(prop_a(j)./a(j))-b0.*(prop_a(j)-a(j)) + (a0-1).*log(prop_b(j)./b(j))-b0.*(prop_b(j)-b(j));
        log_q = log(prop_a(j)./a(j))+log(prop_b(j)./b(j));

        lrr = prop_ll-curr_ll+log_pr+log_q;
        acc(t,j) = rand<exp(lrr);
        if acc(t,j)
            a(j) = prop_a(j);
            b(j) = prop_b(j);
        end
   end
   
   if d>1
       if full
          lr = zeros(S+1,d);
       else
          lr = zeros(S,d); 
       end

       for j=1:d
           wts0 = zeros(d,1);
           wts0(j) = 1;
           if full
               [~,pk] = loglik_full(a,b,wts0,alpha);
               pr = pk;
           else
               [~,pk] = loglik(a,b,wts0,alpha);
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
               %zi(j,:) = multirnd(n(j),lcprob(j,:));
           end
       else
           zi = zeros(S,d);
           for j=1:S
               zi(j,:) = mnrnd(n(j),lcprob(j,:));
               %zi(j,:) = multirnd(n(j),lcprob(j,:));
           end
       end
       zcts = sum(zi,1);
       tmp = gamrnd(1./d+zcts,1);
       eta = tmp'./sum(tmp);
   else
       eta = 1; 
   end
   
   if full
       [~,pr] = loglik_full(a,b,eta,alpha);
       rho = 1-pr(1);
       n(1) = nbinrnd(sum(n(2:end)),rho);
       Nhat(t) = sum(n);       
   else
       [~,pr] = loglik(a,b,eta,alpha);
        Nhat(t) = sum(n)./(1-pr(1));
   end
   
   if trunc
      for l=1:nalpha
        prca = pr_cap(a,b,eta,alpha0(l));
        %Nhata(t,l) = binornd(sum(n),prca);
        Nhata(t,l) = binrnd(sum(n),prca);
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
   A(t,:) = a;
   B(t,:) = b;
   
end
    
    function x = binrnd(n,p)
       tmp = rand(n,1);
       x = sum(tmp<=p);
    end

    function x = multirnd(n,p)
       p = p./sum(p);
       dp = length(p);
       [ps,idp] = sort(p,'descend');
       x = zeros(dp,1); sx = sum(x);
       for j=1:dp
            if sx < n
                x(j) = binornd(n-sx,ps(j));
                sx = sum(x);
            end
       end
       x = x(idp);
    end


    function prc = capprobs(a,b,S)
        prc = zeros(S+1,1);
        for j=0:S
           prc(j) = nchoosek(S,j)*gamma(a+j)*gamma(b-j+S)./(beta(a,b)*gamma(a+b+S));
        end
    end



    function [ll,pk] =loglik(as,bs,wts,alpha)
        pk = zeros(S+1,1);
        for k=0:S
            pk(k+1) = -nchoosek(S,k).*wts'*((gamma(as).*gamma(bs).*(beta(as+k,bs-k+S).*betainc(alpha,as+k,bs-k+S).*gamma(as+bs+S)-gamma(as+k).*gamma(bs-k+S)))./...
                (beta(as,bs).*(gamma(as).*gamma(bs)-beta(as,bs).*betainc(alpha,as,bs).*gamma(as+bs)).*gamma(as+bs+S)));
        end
        %pk = pk(2:end)./(1-pk(1));
        
        ll = (n*log(pk(2:end))-sum(n).*log(1-pk(1)));
    end

    function [ll,pk] =loglik_full(as,bs,wts,alpha)
        pk = zeros(S+1,1);
        for k=0:S
            pk(k+1) = -nchoosek(S,k).*wts'*((gamma(as).*gamma(bs).*(beta(as+k,bs-k+S).*betainc(alpha,as+k,bs-k+S).*gamma(as+bs+S)-gamma(as+k).*gamma(bs-k+S)))./...
                (beta(as,bs).*(gamma(as).*gamma(bs)-beta(as,bs).*betainc(alpha,as,bs).*gamma(as+bs)).*gamma(as+bs+S)));
        end
        %pk = pk(2:end)./(1-pk(1));
        
        ll = (n*log(pk));
    end

    function prca = pr_cap(as,bs,wts,alpha)        
        prca = wts'*((gamma(as).*gamma(bs) - beta(as,bs).*betainc(alpha, as, bs).*gamma(as + bs))./(gamma(as).*gamma(bs)));                
    end

end


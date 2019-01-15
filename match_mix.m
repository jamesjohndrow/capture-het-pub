clear; 
run_pre = false;
dosims = true;
%[thet_hat,ll] = bino_max([0 0],100);

if run_pre
    S = 7;
    %c = [0.200357, 0.194929, 0.174171, 0.149978, 0.121691, 0.0923489,0.0665257];
    %c = [0.119836, 0.158131, 0.155216, 0.14812, 0.142939, 0.13926, 0.136499];
    c = [0.107202, 0.139955, 0.136878, 0.13039, 0.125666, 0.122309, 0.11979, 0.11781];
    n = round(2000.*c);
    n = n(2:end);
    %S = 5;
    %c = [0.227109, 0.221593, 0.192539, 0.157337, 0.118575, 0.0828472];
    
    alpha = .3;
    c0 = c(1);
    c = c(2:end)./(1-c(1));
    d = 3;
    
    %[thet_hat,loss,ps,wts] = min_mix(zeros(5,1),c',1-(1-alpha)^(1/S),S);
    
    [ thet_hat, loss, ps, wts ] = ml_mix( normrnd(0,1,[2*d-1 1]), 1-(1-alpha)^(1/S), n, S, d );
    
    mmt_check = zeros(S+1,1);
    for j=0:S
        mmt_check(j+1) = wts'*(nchoosek(S,j).*ps.^(j).*(1-ps).^(S-j));
    end
    
    ps
    wts
    [mmt_check(2:end)./(1-mmt_check(1)) c']
    [1-mmt_check(1) 1-c0]
    
    d = 2;
    alpha = .3;
    c = [0.107202, 0.139955, 0.136878, 0.13039, 0.125666, 0.122309, 0.11979, 0.11781];
    %c = [0.119836, 0.158131, 0.155216, 0.14812, 0.142939, 0.13926, 0.136499];
    S = 7;
    n = round(5000.*c);
    n = n(2:end);
    
    [ thet_hat, loss, pk, as, bs, wts ] = ml_mixbeta( normrnd(0,1,[3*d-1 1]), 1-(1-alpha)^(1/S), n, S, d );
    
    pk
    as
    bs
    wts
    
    [A,B,eta,Nhat] = BayesBetaMix( n, d, S, 1-(1-alpha)^(1/S),[.25 .25],[.25 .25],1,1,1000,false, 100000  );
    save('Output/betamix.mat','A','B','eta','Nhat');
    
    c = [0.0769614, 0.104211, 0.127879, 0.149002, 0.161626, 0.159402,0.135965, 0.0849526];
    alpha = .1;
    n = mnrnd(5000,c./sum(c));
    n = n(2:end);
    
    [A,B,eta,Nhat] = BayesBetaMix( n, d, S, 1-(1-alpha)^(1/S),[.25 .25],[.25 .25],1/2,1/2,1000,false, 100000  );
    save('Output/betamix2.mat','A','B','eta','Nhat');
    
    n = mnrnd(5000,c./sum(c));
    [A,B,eta,Nhat] = BayesBetaMix( n, d, S, 1-(1-alpha)^(1/S),[.25 .25],[.25 .25],1/2,1/2,1000,true, 100000  );
    save('Output/betamix2_full.mat','A','B','eta','Nhat');
    
    c = [0.109153, 0.109492, 0.0939661, 0.0986536, 0.124706, 0.157342, 0.171996, 0.134692];
    n = mnrnd(5000,c./sum(c)); d = 2; alpha = .1; S = 7;
    [A,B,eta,Nhat,Nhata] = BayesBetaMix( n, d, S, 1-(1-alpha)^(1/S),[.1 .1],[.1 .1],1/4,1/4,1000,true,true, 100000  );
    save('Output/betamix_full.mat','A','B','eta','Nhat');
    
    d=3;
    [mu,eta,Nhat,Nhata] = BayesAtomMix( n, d, S, 1-(1-alpha)^(1/S),[.05 .05 .05],1000,true,true, 100000  );
end

if dosims
    % actually run for paper
    Cs = [0.219674, 0.231806, 0.195864, 0.177008, 0.120624, 0.0471984, 0.0078256;...
        0.39078, 0.106409, 0.0892353, 0.141162, 0.155763, 0.0933228, 0.0233282;...
        %0.609268, 0.201868, 0.0804118, 0.039413, 0.0272994, 0.0223828, 0.0193565;...
        %0.360956, 0.229352, 0.115296, 0.0791537, 0.0723199, 0.0714915, 0.0714307;...
        %0.142857, 0.142857, 0.142857, 0.142857, 0.142857, 0.142857, 0.142857];
        0.225586, 0.123047, 0.102539, 0.0976563, 0.102539, 0.123047, 0.225586;...
        0.625, 0.25, 0.0892857, 0.0274725, 0.00686813, 0.00124875, 0.000124875;...
        0.142857, 0.142857, 0.142857, 0.142857, 0.142857, 0.142857, 0.142857];
    
    %Cs(:,2:end) = bsxfun(@times,Cs(:,2:end),1./(1-Cs(:,1)));

    nsim = 5; N = 500; S = 6; alpha = [.005 .01 .05 .1];
    nmc = 100000;

    for j=4:5
        disp([num2str(j) ' : beta']);
        n = mnrnd(N,Cs(j,:)./sum(Cs(j,:)));
        d=1;
        [A,B,eta,Nhat,Nhata] = BayesBetaMix( n, d, S, 1-(1-alpha).^(1/S),[.2 .2],[.2 .2],1/4,1/4,10000,true,true, nmc );

        disp([num2str(j) ' : atom']);
        d=3;
        [mu,eta_atom,Nhat_atom,Nhata_atom] = BayesAtomMix( n, d, S, 1-(1-alpha).^(1/S),[.05 .05 .05],10000,true,true, nmc );
        save(strcat('Output/bayes_results_case_',num2str(j),'.mat'),'A','B','n','eta','eta_atom','mu','Nhat','Nhat_atom','Nhata','Nhata_atom');

    end
end

% hare dataq
n = [2 25 22 13 5 1 2];
d=1; nmc = 100000; S = 6; alpha = [.005 .01 .05 .1];
[A,B,eta,Nhat,Nhata] = BayesBetaMix( n, d, S, 1-(1-alpha).^(1/S),[.2 .2],[.2 .2],1/4,1/4,10000,true,true, nmc );
d=3;
[mu,eta_atom,Nhat_atom,Nhata_atom] = BayesAtomMix( n, d, S, 1-(1-alpha).^(1/S),[.05 .05 .05],1000,true,true, nmc );
save('Output/bayes_results_hares.mat','A','B','n','eta','eta_atom','mu','Nhat','Nhat_atom','Nhata','Nhata_atom');







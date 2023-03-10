function out = MVALSE_HN( Y, N, X, Noise, is_prior, model_order, theta_prior, kappa_prior)
% MVALSE algorithm in Heteroscedastic Noise Environment
% INPUTS:
%   Y  - measurement vector of size M*L
%   N  - is the assumed number of complex sinusoids
%   X  - the true signal - used for computing the MSE vs iterations
%   Noise  - pre-assumed noise case 
%   is_prior - indicator determininL if the prior information is used
%           is_prior=0 will not use prior information
%           is_prior=1 will use prior information
%   model_order  - different ways to update the model order
%           '0'  - use the traditional way
%           '1'  - use tht new way
%   theta_prior, kappa_prior  - prior information of frequencies
% OUTPUTS:
%   out - structure
%      .freqs      - vector of frequency estimates
%      .amps       - vector of amplitude estimates
%      .X_estimate - reconstructed siLnal
%      .noise_var  - estimate of the noise variance
%      .iterations - number of iterations until converLence
%      .mse        - evolution of the mse of x_estimate with iterations
%      .K          - evolution of the estimated number of components with iterations
%
% See full paper:
%       "Gridless Variational Direction-of-Arrival Estimation in Heteroscedastic Noise Environment"
%       preprint available at https://arxiv.org/pdf/1912.11738.pdf
%       by Qi Zhang, Jiang Zhu, Yuantao Gu and Zhiwei Xu

%% Initialization of the posterior pdfs of the frequencies
[M,L] = size(Y);
M_vec = (0:1:M-1)';
M_vec_init = (1:1:M-1)';
A     = zeros(M,N);
J     = zeros(N,N,L);
H     = zeros(N,L);
W     = zeros(N,L);
C     = zeros(N,N,L);
T     = 500;       
mse   = zeros(T,1);
Kt    = zeros(T,1);
t     = 1;
if is_prior
    eta_prior = kappa_prior.*exp(1j*theta_prior);
end
    
%% Initialization of the posterior pdfs of the frequencies
res   = Y;
if is_prior
    prior_set = 1:1:length(theta_prior);
end
for n=1:N
    YI = res;
    R  = YI*YI';
    gamma = zeros(M-1,1);
    for i=2:M
        for k=1:i-1
            gamma(i-k) = gamma(i-k) + R(i,k);
        end
    end
    if(n==1)
        eta_toep = toeplitz([trace(Y'*Y);gamma(:,1)])/M;  
        Y_d = sort(real(eig(eta_toep)));
        nu = mean(Y_d(1:floor(N/4)))/L;
        K   = floor(N/2);
        rho = K/N;
        tau = (trace(Y'*Y)/M - nu*L)/(rho*N);
        Nu = nu*ones(M,L);
    end
    eta   = 2*gamma/(M+nu/tau)/nu;
    if (is_prior == 0||isempty(prior_set))
        [H2_theta,H2_k]= Heu(eta, M_vec_init);
    else
        [H2_theta,H2_k,prior_index]= Heu_Prior(eta, M_vec_init, eta_prior(prior_set));
        prior_index_new =  prior_set(prior_index);
        prior_set = setdiff(prior_set,prior_index_new);
    end

    A(:,n) = exp(1j*M_vec*H2_theta).*( besseli(M_vec,H2_k,1)/besseli(0,H2_k,1) );
    % initialize a_hat and w_hat
    for ll = 1:L
        J(1:n-1,n,ll) = A(:,1:n-1)'*diag(1./Nu(:,ll))*A(:,n); J(n,1:n-1,ll) = J(1:n-1,n,ll)'; J(n,n,ll) = sum(1./Nu(:,ll));
        H(n,ll) = A(:,n)'*diag(1./Nu(:,ll))*Y(:,ll);
        C(1:n,1:n,ll) = eye(n)/(J(1:n,1:n,ll) + 1/tau*eye(n));
        W(1:n,ll) = C(1:n,1:n,ll)*H(1:n,ll);
    end
    % save mse and K at initialization
    res = Y - A(:,1:n)*W(1:n,:);
    if n == K 
        X_esti_old    = A(:,1:n)*W(1:n,:);
        mse(t) = norm(X - X_esti_old)^2/norm(X)^2;
        Kt(t)  = K;
    end
end
cont = 1;
while cont
    t = t + 1;
%% Update the support and weights
    s     = false(N,1); % Initialize s
    K     = 0;
    W     = zeros(N,L);
    C     = zeros(N,N,L);
    while(1)
        delta = zeros(N,1);
        for i = 1:N
            if(~s(i))
                for ll = 1:L
                    A4_v = real(1/(sum(1./Nu(:,ll)) + 1/tau - J(s,i,ll)'*C(s,s,ll)*J(s,i,ll)));
                    A4_u = A4_v*(H(i,ll)-J(s,i,ll)'*W(s,ll));
                    delta_ll = log(A4_v/tau) + A4_u*A4_u'/A4_v;
                    delta(i) = delta(i) + delta_ll;
                end
                if model_order == 0
                    delta(i) = real(delta(i)) + log(rho/(1-rho));
                elseif model_order == 1
                    delta(i) = real(delta(i))/L + log(rho/(1-rho));
                end
            else
                for ll = 1:L
                    delta_ll = -log(C(i,i,ll)/tau) - W(i,ll)*W(i,ll)'/C(i,i,ll);
                    delta(i) = delta(i) + delta_ll;
                end
                if model_order == 0
                    delta(i) = real(delta(i)) - log(rho/(1-rho));
                elseif model_order == 1
                    delta(i) = real(delta(i))/L - log(rho/(1-rho));
                end
            end
            if(K == N)
                delta = -1*ones(N,1);
            end
        end
        [~,k] = max(real(delta));
        if(delta(k)<=0)
            break
        else
             if(~s(k))
                for ll = 1:L
                    A4_vmax = 1/((sum(1./Nu(:,ll))) + 1/tau - J(s,k,ll)'*C(s,s,ll)*J(s,k,ll));
                    A4_umax = A4_vmax*(H(k,ll)-J(s,k,ll)'*W(s,ll));
                    W(s,ll) = W(s,ll) - C(s,s,ll)'*J(s,k,ll)*A4_umax;
                    W(k,ll) = A4_umax;
                    C_temp = [C(s,s,ll),zeros(K,1);zeros(1,K+1)] + A4_vmax* [C(s,s,ll)*J(s,k,ll);-1]*[C(s,s,ll)*J(s,k,ll);-1]'; 
                    C(s,s,ll) = C_temp(1:K,1:K);
                    C(s,k,ll) = C_temp(1:K,K+1);
                    C(k,s,ll) = C_temp(K+1,1:K);
                    C(k,k,ll) = C_temp(K+1,K+1);
                end
                s(k) = ~s(k); K = K+1;
            else
                s(k) = ~s(k); K = K-1;
                for ll = 1:L
                    W(s,ll) = W(s,ll) - C(s,k,ll)*W(k,ll)/C(k,k,ll);
                    C(s,s,ll) =  C(s,s,ll) - C(s,k,ll)*C(k,s,ll)/C(k,k,ll);
                end
            end   
            for ll = 1:L
                C(:,:,ll) = (C(:,:,ll)+C(:,:,ll)')/2; % ensure the diagonal is real
            end
        end
    end
%% update v, rho, tau
    if (K>0)
        switch Noise
            case 0
                 nu_temp = 0;
                 for ll = 1:L
                     for i = 1:M
                         nu_temp = nu_temp + norm(Y(i,ll)-A(i,s)*W(s,ll),2)^2 + A(i,s)*C(s,s,ll)*A(i,s)' + sum((abs(W(s,ll)').^2).*(1-abs(A(i,s)).^2));
                     end
                 end
                 nu_temp = real(nu_temp/M/L);
                 Nu = nu_temp*ones(M,L);
            case 1
                 for ll = 1:L
                     nu_temp = 0;
                     for i = 1:M
                         nu_temp = nu_temp + norm(Y(i,ll)-A(i,s)*W(s,ll),2)^2 + A(i,s)*C(s,s,ll)*A(i,s)' + sum((abs(W(s,ll)').^2).*(1-abs(A(i,s)).^2));
                     end
                     nu_temp = nu_temp/M;
                     Nu(:,ll) = real(nu_temp)*ones(M,1);
                 end
            case 2
                for i = 1:M
                    nu_temp = 0;
                    for ll = 1:L
                        nu_temp = nu_temp + norm(Y(i,ll)-A(i,s)*W(s,ll),2)^2 + A(i,s)*C(s,s,ll)*A(i,s)' + sum((abs(W(s,ll)').^2).*(1-abs(A(i,s)).^2));
                    end
                    nu_temp = nu_temp/L;
                    Nu(i,:) = real(nu_temp)*ones(1,L);
                end
            case 3
                for ll = 1:L
                     for i = 1:M
                         Nu(i,ll)  = real(norm(Y(i,ll)-A(i,s)*W(s,ll),2)^2 + A(i,s)*C(s,s,ll)*A(i,s)' + sum((abs(W(s,ll)').^2).*(1-abs(A(i,s)).^2)));
                     end
                 end 
        end
        tau = trace(W(s,:)'*W(s,:));
        for ll = 1:L
            tau = tau + trace(C(s,s,ll));
        end
        tau = real(tau/(K*L));
        if K<N
            rho = K/N;
        else
            rho = (N-1)/N;         % just to avoid the potential issue of log(1-rho) when rho=1
        end
    else
        rho = 1/N;                 % just to avoid the potential issue of log(rho) when rho=0
    end
%% update eta, A, and estimate frequencies
    th_esti = zeros(K,1);
    ka_esti = zeros(K,1);
    if is_prior
        prior_set = 1:1:length(theta_prior);
    end
    s_ind = 1:N; s_ind = s_ind(s); % indices of the non-zero components
    for s_ind_k = 1:K
        i = s_ind(s_ind_k);
        eta = zeros(M,1);
        for ll = 1:L
            eta = eta + 2*(diag(1./Nu(:,ll)))*(Y(:,ll)*W(i,ll)' - A(:,s)*C(s,i,ll) - A(:,s)*W(s,ll)*W(i,ll)'...
                  + A(:,i)*(C(i,i,ll)+trace(W(i,ll)*W(i,ll)')));
        end
        if (is_prior == 0||isempty(prior_set))
            [H2_theta,H2_k] = Heu(eta,M_vec);
        else
            [H2_theta,H2_k,prior_index] = Heu_Prior(eta, M_vec, eta_prior(prior_set));
            prior_index_new =  prior_set(prior_index);
            prior_set = setdiff(prior_set,prior_index_new);
        end
        A(:,i) = exp(1j*M_vec*H2_theta).*( besseli(M_vec,H2_k,1)/besseli(0,H2_k,1) );
        th_esti(s_ind_k) = H2_theta;
        ka_esti(s_ind_k) = H2_k;
    end
    for ll = 1:L
        J(:,:,ll) = A'*(diag(1./Nu(:,ll)))*A(:,:);
        J(:,:,ll) = J(:,:,ll)';
        J(:,:,ll) = J(:,:,ll) - diag(diag(J(:,:,ll))) + sum(1./Nu(:,ll))*eye(N);
        H(:,ll) = A(:,:)'*(diag(1./Nu(:,ll)))*Y(:,ll);
    end
%% stopping criterion:
    % the relative change of the reconstructed signalis below threshold or
    % max number of iterations is reached
    X_esti     = A(:,s)*W(s,:);
    mse(t) = norm(X_esti-X)^2/norm(X)^2;
    Kt(t)  = K;
    if (norm(X_esti-X_esti_old)/norm(X_esti_old)<1e-6) || (norm(X_esti_old)==0&&norm(X_esti-X_esti_old)==0) || (t >= T)
        cont = 0;
        mse(t+1:end) = mse(t);
        Kt(t+1:end)  = Kt(t);
    end
    X_esti_old = X_esti;
end
out = struct('freqs',th_esti,'kappas',ka_esti,'amps',W(s,:),'X_estimate',X_esti,'noise_var',Nu,'iterations',t,'mse',mse,'K',Kt);
end


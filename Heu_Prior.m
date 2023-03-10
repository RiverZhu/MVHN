function [ theta, kappa, prior_index] = Heu_Prior( eta, m, eta_prior)
N     = length(m);
ka    = abs(eta);
A     = besseli(1,ka,1)./besseli(0,ka,1);
kmix  = Ainv( A.^(1./m.^2) );
k     = N;
eta_q = kmix(k) * exp( 1i * ( angle(eta(k)) + 2*pi*(1:m(k)).' )/m(k) );
for k = N-1:-1:1
    if m(k) ~= 0
        phi   = angle(eta(k));
        eta_q = eta_q + kmix(k) * exp( 1i*( phi + 2*pi*round( (m(k)*angle(eta_q) - phi)/2/pi ) )/m(k) );
    end
end

Rec_prior = ones(length(eta_q),1);
eta_record = eta_q;
for i=1:m(end)
    for j=1:length(eta_prior)
        if(abs(eta_prior(j) + eta_record(i))>=abs(eta_q(i)))
            eta_q(i) = eta_prior(j) + eta_record(i);
            Rec_prior(i) = j;
        end
    end
end

[~,in] = max(abs(eta_q));
mu     = angle(eta_q(in));
d1    = real(1j*eta_prior(Rec_prior(in))'*exp(1j*mu) + eta'*(1j*m.*exp(1j*m*mu)));
d2    = -real(eta_prior(Rec_prior(in))'*exp(1j*mu) + eta'*(m.*m.*exp(1j*m*mu)));

if d2<0 
    theta  = mu - d1/d2;
    kappa  = Ainv( exp(0.5/d2) );
else   
    theta  = mu;
    kappa  = abs(eta_q(in));
end
prior_index = Rec_prior(in);
end




function [theta, kappa] = Heu( eta, m )
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
[~,in] = max(abs(eta_q));
mu     = angle(eta_q(in));
d1     = -imag( eta' * ( m    .* exp(1i*m*mu) ) );
d2     = -real( eta' * ( m.^2 .* exp(1i*m*mu) ) );
if d2<0 
    theta  = mu - d1/d2;
    kappa  = Ainv( exp(0.5/d2) );
else   
    theta  = mu;
    kappa  = abs(eta_q(in));
end
end


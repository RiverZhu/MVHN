function [ k ] = Ainv( R )
% Returns the approximate solution of the equation R = A(k),
% where A(k) = I_1(k)/I_0(k) is the ration of modified Bessel functions of
% the first kind of first and zero order
% Uses the approximation from
%       Mardia & Jupp - Directional Statistics, Wiley 2000, pp. 85-86.
%
% When input R is a vector, the output is a vector containing the
% corresponding entries

    k   = R;                  % define A with same dimensions
    in1 = (R<.53);            % indices of the entries < .53
    in3 = (R>=.85);           % indices of the entries >= .85
    in2 = logical(1-in1-in3); % indices of the entries >=.53 and <.85
    R1  = R(in1);             % entries < .53
    R2  = R(in2);             % entries >=.53 and <.85
    R3  = R(in3);             % entries >= .85

    % compute for the entries which are < .53
    if ~isempty(R1)
        t      = R1.*R1;
        k(in1) = R1 .* ( 2 + t + 5/6*t.*t );
    end
    % compute for the entries which are >=.53 and <.85
    if ~isempty(R2)
        k(in2) = -.4 + 1.39*R2 + 0.43./(1-R2);
    end
    % compute for the entries which are >= .85
    if ~isempty(R3)
        k(in3) = 1./( R3.*(R3-1).*(R3-3) );
    end
end


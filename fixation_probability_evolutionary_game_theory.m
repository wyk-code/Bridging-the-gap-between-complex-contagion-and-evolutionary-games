function phi = fixation_probability_evolutionary_game_theory(a,b,c,d,w,i0)
%-------------------------------------------------------------
% Moran process (well-mixed)
% exponential fitness mapping
%
% payoff matrix
%       A      B
% A     a      b
% B     c      d
%
% N   : population size
% w   : selection intensity
% i0  : initial number of mutants (A)
%
% output:
% phi : fixation probability of A
%-------------------------------------------------------------

N = 100;

gamma = zeros(N-1,1);

if i0 == 0
    phi = 0;
    return;
end

if i0 == N
    phi = 1;
    return;
end

%% ---------- Neutral selection ----------

if w == 0
    phi = i0/N;
    return;
end

for j = 1:N-1

    % fraction of A
    x = j/N;

    % expected payoff
    piA = ((j-1)*a + (N-j)*b)/(N-1);
    piB = (j*c + (N-j-1)*d)/(N-1);

    % exponential fitness
    fA = exp(w*piA);
    fB = exp(w*piB);

    % gamma=T-/T+
    gamma(j) = fB/fA;

end

% cumulative products
prod_gamma = cumprod(gamma);

% denominator
den = 1 + sum(prod_gamma);

% numerator
if i0==1
    num = 1;
else
    num = 1 + sum(prod_gamma(1:i0-1));
end

phi = num/den;

end
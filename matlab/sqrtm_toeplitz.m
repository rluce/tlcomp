function [S, flag, res] = sqrtm_toeplitz(T, maxiter, eps_stop, do_scale)
% S = sqrtm_toeplitz(r)
% S = sqrtm_toeplitz(r, maxiter)
% S = sqrtm_toeplitz(r, maxiter, eps_stop)
% S = sqrtm_toeplitz(r, maxiter, eps_stop, do_scale)
%
% [S, flag] = sqrtm_toeplitz(r);
% [S, flag, res] = sqrtm_toeplitz(r);
%
%
% Simplistic sqrtm Newton iteration for an SPD Toeplitz matrix T.
%
% Input:
%   r        -- First row of T
%   maxiter  -- (optional) Iteration limit (default: 10).
%   eps_stop -- (optional) Stopping criterion tolerance.  The Newton iteration
%               is stopped as soon as |S_{k+1}^2 - T|_F <= eps_stop * |T|_F
%               (default: 1e-8).
%   do_scale -- (optional) true/false switch that results in a damped Newtion
%               iteration (default: true).
%               See: "Optimally scaled Newton iterations for the matrix
%               square root", by Bernhard Beckermann, some proceedings.
%
% Passing an empty matrix for any of the optional parameters results in
% default parameter choice.
%
% Output:
%   S        -- An approximation S \approx sqrtm(T).
%   flag     -- Termination flag
%               flag==0 <=> Normal termination
%               flag==1 <=> Iteration limit reached
%               flag==2 <=> Convergence problem (residual norm increaased
%                           from one step to the next.
%   res      -- |S^2 - T|_F of the returned iterate S

if nargin < 4 || isempty(do_scale)
    do_scale = true;
end

if nargin < 3 || isempty(eps_stop)
    eps_stop = 1e-8;
end

if nargin < 2 || isempty(maxiter)
    maxiter = 10;
end

norm_T = norm(T, 'fro');
S = TLMat(T.c, T.r);

fprintf('Input matrix is of size %d, |T|_F=%8.2e\n', size(T,1), norm_T);
if do_scale
    tic;
    mu = get_mus(T, maxiter);
    fprintf('Computing mu took %.2fs\n', toc);
else
    mu = ones(maxiter, 1);
end

do_stop = false;
res = +inf;
flag = -1;
iter = 0;
fprintf('%4s  %5s  %8s  %8s  %8s\n', 'iter', 'drank', '|S^2-T|', 'mu', 'it-time')
while ~do_stop
    t_iter = tic;
    iter = iter + 1;
    S_old = S;
    S = .5 * (mu(iter) * S + (T / S) / mu(iter));
    res_old = res;
    res = gennorm(S^2 - T);
    
    % Determine whether we can stop the iteration
    if res <= eps_stop * norm_T
        % Converged w.r.t given tolerance
        do_stop = true;
        flag = 0;
    elseif res > res_old
        % residual has increased, we force termination
        do_stop = true;
        S = S_old;
        res = res_old;
        flag = 2;
    elseif iter >= maxiter
        % iteration limit reached
        do_stop = true;
        flag = 1;
    end
    fprintf('%4d  %5d  %8.2e  %8.2e  %7.2fs\n', ...
        iter, drank(S), res, mu(iter), toc(t_iter));
end

switch flag
    case 0
        fprintf('Converged in %d iterations to tolerance %8.2e\n', ...
            iter, eps_stop);
    case 1
        fprintf('Iteration limit (%d) reached, quitting.\n', maxiter);
    case 2
        fprintf('Convergence problem in the Newton iteration, quitting.\n');
    otherwise
        % If this happens this is a bug
        assert(false);
end

fprintf('Final relative residual |S^2 - T|_F / |T|_F = %8.2e\n', res/norm_T);

end

function mu = get_mus(T, maxiter)
opts.tol = 1e-5;
opts.isreal = true;
opts.issym = true;
opts.p = 40;
e = eigs(@(x) T*x, size(T,1), 2, 'BE', opts);
e_min = min(e);
e_max = max(e);

fprintf('EV bounds: [%8.2e, %8.2e]\n', e_min, e_max);

mu = zeros(1,maxiter);
kappa = sqrt(e_min/e_max);
mu(1) = (e_min * e_max)^(-1/4);
mu(2) = sqrt( 2 * sqrt(kappa) / (1 + kappa) );
for k = 3:maxiter
    mu(k) = sqrt( 2 * mu(k-1) / (1 + mu(k-1)^2) );
end

end

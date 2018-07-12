function [X, e, time] = select_mse_explicit_energy(A, T, target_mse, lambda, s, C, e0)
[n, m] = size(A);

%% check if CVX is installed
try
    run('cvx_setup');
catch err
    error('CVX problem.');
end

tic;
% weights
W = ones(m, T);
% selection variable
X = zeros(m, T);
% indices of "1" entries
the_zeros = [];
% indices of "0" entries
the_ones = [];

% in each time instance select one of the sensors, in order
for i = 1:T
    the_ones = [the_ones (i-1)*m+i];
end

while(1)
    the_ones_old = the_ones;
    
    %%% solve main optimization problem
    cvx_begin
        cvx_quiet('true');
        cvx_solver('sedumi');
        
        variable X(m, T);
        variable e(m);

        minimize (vec(W)'*vec(X) + lambda*norm(e))
        subject to
                for t = 1:T
                    trace_inv(A*diag(X(:, t))*A') <= target_mse;
                end
                X >= 0;
                X <= 1;
                X(the_zeros) == 0;
                X(the_ones) == 1;
                e >= 0;
                sum(X, 2) >= 1;
                (diag(s) + C)*X*ones(T, 1) <= e0 + e;
    cvx_end

    %%% update weights
    W = 1./(X + 10e-5);

    %%% force entry to 1
    the_ones = find(abs(X)>=0.95);
    if (length(the_ones) <= length(the_ones_old))
        found = 0;
        [~, indices] = sort(vec(X), 'descend');
        for ii = 1:m*T
            if (X(indices(ii)) <= 1-10e-4)
                the_ones = unique([the_ones; indices(ii)]);
                X(indices(ii)) = 1;
                found = 1;
                break;
            end
        end
        if (found == 0)
            ii = randsample(m*T, 1);
            the_ones = unique([the_ones; ii]);
            X(ii) = 1;
        end
    end
    the_zeros = find(abs(X)<=10e-6);

    %%% binary solution reached
    if ((length(the_ones) + length(the_zeros)) == m*T)
        X(the_ones) = 1;
        X(the_zeros) = 0;
        break;
    end
end
time = toc;

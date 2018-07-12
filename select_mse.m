function [x, time] = select_mse(A, target_mse)
[n, m] = size(A);

%% check if CVX is installed
try
    run('cvx_setup');
catch err
    error('CVX problem.');
end

tic;
% weights
w = ones(m, 1);
% selection variable
x = zeros(m, 1);
% indices of "1" entries
the_zeros = [];
% indices of "0" entries
the_ones = [];

while(1)
    the_ones_old = the_ones;
    
    %%% solve main optimization problem
    cvx_begin
        cvx_quiet('true');
        cvx_solver('sedumi');
        
        variable x(m, 1);

        minimize (w'*x)
        subject to
                trace_inv(A*diag(x)*A') <= target_mse;
                x >= 0;
                x <= 1;
                x(the_zeros) == 0;
                x(the_ones) == 1;
    cvx_end

    %%% update weights
    w = 1./(abs(x) + 10e-5);

    %%% force entry to 1
    the_ones = find(abs(x)>=0.95);
    if (length(the_ones) <= length(the_ones_old))
        [~, indices] = sort(x, 'descend');
        for ii = 1:m
            if (x(indices(ii)) <= 1-10e-4)
                the_ones = unique([the_ones; indices(ii)]);
                x(indices(ii)) = 1;
                break;
            end
        end
    end
    the_zeros = find(abs(x)<=10e-7);

    %%% binary solution reached
    if ((length(the_ones) + length(the_zeros)) == m)
        x(the_ones) = 1;
        x(the_zeros) = 0;
        break;
    end
end
time = toc;

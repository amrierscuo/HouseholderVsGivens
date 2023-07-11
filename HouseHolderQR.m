%HouseHolderQR.m
function [Q, R] = HouseHolderQR(A)

    [m, n] = size(A);
    Q = eye(m);
    R = A;

    for i = 1:n
        X = R(i:m, i);
        [Ht, ~] = HouseHolderMatrix(X);
        H = [eye(i-1) zeros(i-1, m-i+1); zeros(m-i+1, i-1) Ht];
        R = H * R;
        Q = Q * H;
    end
end

function [Ht, k] = HouseHolderMatrix(x)
    n = length(x);
    e1 = zeros(n, 1);
    e1(1) = 1;
    k = sign(x(1)) * norm(x);
    v = x + k * e1;
    Ht = eye(n) - 2 * (v * v') / (v' * v);
end

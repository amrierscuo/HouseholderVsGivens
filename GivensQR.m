%GivensQR.m
function [Q, R] = GivensQR(A)
    [m, n] = size(A);
    R = A;
    G = eye(m);

    for k = 1:n
        for i = k+1:m
            [c, s] = GivRot(R(k, k), R(i, k));
            for j = k:n
                t = c * R(k, j) + s * R(i, j);
                R(i, j) = -s * R(k, j) + c * R(i, j);
                R(k, j) = t;
            end
            G = GMatrix(c, s, k, i, m) * G;
        end
    end

    Q = G';
end

function [c, s] = GivRot(a, b)
    r = sqrt(a^2 + b^2);
    c = a / r;
    s = b / r;
end

function G = GMatrix(c, s, k, i, m)
    G = eye(m);
    G(k, k) = c;
    G(i, i) = c;
    G(k, i) = s;
    G(i, k) = -s;
end

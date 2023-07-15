%% File: QRvsSVD\Test2QRvsSVD.m

% Genera dati casuali
m = 20;
n = 10;
times = 20;  % Esegui solo 5 iterazioni

% Inizializza le variabili per salvare i risultati:
% quadrate, sovradimensionate e sottodimensionate
arrayIndex = 1:times;
dimRowSquare = zeros(1, times);
dimColumnSquare = zeros(1, times);
errorSolmQRSquare = zeros(1, times);
errorSolhQRSquare = zeros(1, times);
errorSolgQRSquare = zeros(1, times);
errorSolmSVDSquare = zeros(1, times);
errorSolhSVDSquare = zeros(1, times);
errorSolgSVDSquare = zeros(1, times);

dimRowOverdetermined = zeros(1, times);
dimColumnOverdetermined = zeros(1, times);
errorSolmQROverdetermined = zeros(1, times);
errorSolhQROverdetermined = zeros(1, times);
errorSolgQROverdetermined = zeros(1, times);
errorSolmSVDOverdetermined = zeros(1, times);
errorSolhSVDOverdetermined = zeros(1, times);
errorSolgSVDOverdetermined = zeros(1, times);

dimRowUnderdetermined = zeros(1, times);
dimColumnUnderdetermined = zeros(1, times);
errorSolmQRUnderdetermined = zeros(1, times);
errorSolhQRUnderdetermined = zeros(1, times);
errorSolgQRUnderdetermined = zeros(1, times);
errorSolmSVDUnderdetermined = zeros(1, times);
errorSolhSVDUnderdetermined = zeros(1, times);
errorSolgSVDUnderdetermined = zeros(1, times);

% Inizializza le variabili per salvare le matrici
A_square_all = cell(1, times);
x_square_all = cell(1, times);
b_square_all = cell(1, times);

A_overdetermined_all = cell(1, times);
x_overdetermined_all = cell(1, times);
b_overdetermined_all = cell(1, times);

A_underdetermined_all = cell(1, times);
x_underdetermined_all = cell(1, times);
b_underdetermined_all = cell(1, times);

% Inizializza le variabili per salvare i tempi di esecuzione
tQR_Square = zeros(1, times);
tSVD_Square = zeros(1, times);
tQR_Overdetermined = zeros(1, times);
tSVD_Overdetermined = zeros(1, times);
tQR_Underdetermined = zeros(1, times);
tSVD_Underdetermined = zeros(1, times);

% Loop di iterazione
for i = 1:times
    % Genera dati casuali per matrici quadrate
    A_square = rand(m, m);
    x_square = rand(m, 1);
    b_square = A_square * x_square;

    % Genera dati casuali per sistemi sovradeterminati
    A_overdetermined = rand(m, n);
    x_overdetermined = rand(n, 1);
    b_overdetermined = A_overdetermined * x_overdetermined;

    % Genera dati casuali per sistemi sottodeterminati
    A_underdetermined = rand(m, n);  % Invertire m e n
    x_underdetermined = rand(n, 1);
    b_underdetermined = A_underdetermined * x_underdetermined;

    % Salva le matrici quadrate
    A_square_all{i} = A_square;
    x_square_all{i} = x_square;
    b_square_all{i} = b_square;

    % Salva le matrici sovradeterminate
    A_overdetermined_all{i} = A_overdetermined;
    x_overdetermined_all{i} = x_overdetermined;
    b_overdetermined_all{i} = b_overdetermined;

    % Salva le matrici sottodeterminate
    A_underdetermined_all{i} = A_underdetermined;
    x_underdetermined_all{i} = x_underdetermined;
    b_underdetermined_all{i} = b_underdetermined;

    % Calcola gli errori con il metodo QR per matrici quadrate
    [~, ~, ~, ~, errorSolhQRSquare(i)] = ComputeErrors(@HouseHolderQR, A_square, x_square, b_square);
    [~, ~, ~, ~, errorSolgQRSquare(i)] = ComputeErrors(@GivensQR, A_square, x_square, b_square);
    [errorSolmQRSquare(i), ~, ~, ~, errorSolOriginal] = ComputeErrors(@qr, A_square, x_square, b_square);
    fprintf('Square Matrix (size %d): Error in computed solution, Machine Error (infinity norm): %.2e\n', m, errorSolOriginal);

    % Calcola gli errori con il metodo SVD per matrici quadrate
    [errorSolmSVDSquare(i), ~, ~, ~, ~] = ComputeErrors(@svd, A_square, x_square, b_square);

    % Calcola gli errori con il metodo QR per sistemi sovradeterminati
    [~, ~, ~, ~, errorSolhQROverdetermined(i)] = ComputeErrors(@HouseHolderQR, A_overdetermined, x_overdetermined, b_overdetermined);
    [~, ~, ~, ~, errorSolgQROverdetermined(i)] = ComputeErrors(@GivensQR, A_overdetermined, x_overdetermined, b_overdetermined);
    [errorSolmQROverdetermined(i), ~, ~, ~, errorSolOriginal] = ComputeErrors(@qr, A_overdetermined, x_overdetermined, b_overdetermined);
    fprintf('Overdetermined System (size %dx%d): Error in computed solution, Machine Error (infinity norm): %.2e\n', m, n, errorSolOriginal);

    % Calcola gli errori con il metodo SVD per sistemi sovradeterminati
    [errorSolmSVDOverdetermined(i), ~, ~, ~, ~] = ComputeErrors(@svd, A_overdetermined, x_overdetermined, b_overdetermined);

    % Calcola gli errori con il metodo QR per sistemi sottodeterminati
    [~, ~, ~, ~, errorSolhQRUnderdetermined(i)] = ComputeErrors(@HouseHolderQR, A_underdetermined, x_underdetermined, b_underdetermined);
    [~, ~, ~, ~, errorSolgQRUnderdetermined(i)] = ComputeErrors(@GivensQR, A_underdetermined, x_underdetermined, b_underdetermined);
    [errorSolmQRUnderdetermined(i), ~, ~, ~, errorSolOriginal] = ComputeErrors(@qr, A_underdetermined, x_underdetermined, b_underdetermined);
    fprintf('Underdetermined System (size %dx%d): Error in computed solution, Machine Error (infinity norm): %.2e\n', m, n, errorSolOriginal);

    % Calcola gli errori con il metodo SVD per sistemi sottodeterminati
    [errorSolmSVDUnderdetermined(i), ~, ~, ~, ~] = ComputeErrors(@svd, A_underdetermined, x_underdetermined, b_underdetermined);

    % Calcola i tempi di esecuzione con il metodo QR e SVD per matrici quadrate
    [~, tQR_Square(i)] = ComputeExecutionTime(@qr, A_square, b_square);
    [~, tSVD_Square(i)] = ComputeExecutionTime(@svd, A_square, b_square);

    % Calcola i tempi di esecuzione con il metodo QR e SVD per sistemi sovradeterminati
    [~, tQR_Overdetermined(i)] = ComputeExecutionTime(@qr, A_overdetermined, b_overdetermined);
    [~, tSVD_Overdetermined(i)] = ComputeExecutionTime(@svd, A_overdetermined, b_overdetermined);

    % Calcola i tempi di esecuzione con il metodo QR e SVD per sistemi sottodeterminati
    [~, tQR_Underdetermined(i)] = ComputeExecutionTime(@qr, A_underdetermined, b_underdetermined);
    [~, tSVD_Underdetermined(i)] = ComputeExecutionTime(@svd, A_underdetermined, b_underdetermined);

    % Aggiorna m e n
    m = m + 10;
    n = n + 10;

    dimRowSquare(i) = m;
    dimColumnSquare(i) = m;

    dimRowOverdetermined(i) = m;
    dimColumnOverdetermined(i) = n;

    dimRowUnderdetermined(i) = m;
    dimColumnUnderdetermined(i) = n;
end

%% Salva dati .mat
% Salva i dati nelle tabelle
dimensionByStepQRSquare = table(dimRowSquare, dimColumnSquare, 'VariableNames', {'RowSize', 'ColumnSize'});
dimensionByStepSVDSquare = dimensionByStepQRSquare;

dimensionByStepQROverdetermined = table(dimRowOverdetermined, dimColumnOverdetermined, 'VariableNames', {'RowSize', 'ColumnSize'});
dimensionByStepSVDOverdetermined = dimensionByStepQROverdetermined;

dimensionByStepQRUnderdetermined = table(dimRowUnderdetermined, dimColumnUnderdetermined, 'VariableNames', {'RowSize', 'ColumnSize'});
dimensionByStepSVDUnderdetermined = dimensionByStepQRUnderdetermined;

% Salva le matrici quadrate
save('matrici_quadrate.mat', 'A_square_all', 'x_square_all', 'b_square_all');

% Salva le matrici sovradeterminate
save('matrici_sovradeterminate.mat', 'A_overdetermined_all', 'x_overdetermined_all', 'b_overdetermined_all');

% Salva le matrici sottodeterminate
save('matrici_sottodeterminate.mat', 'A_underdetermined_all', 'x_underdetermined_all', 'b_underdetermined_all');

%Salva i tempi di esecuzione
save('execution_times.mat', 'tQR_Square', 'tSVD_Square', 'tQR_Overdetermined', 'tSVD_Overdetermined', 'tQR_Underdetermined', 'tSVD_Underdetermined');

%% Grafici
% Crea il grafico dei risultati per il metodo QR per matrici quadrate
fprintf('QR Square Average Time: %.4fs\n', mean(tQR_Square));
fprintf('SVD Square Average Time: %.4fs\n', mean(tSVD_Square));
fprintf('QR Overdetermined Average Time: %.4fs\n', mean(tQR_Overdetermined));
fprintf('SVD Overdetermined Average Time: %.4fs\n', mean(tSVD_Overdetermined));
fprintf('QR Underdetermined Average Time: %.4fs\n', mean(tQR_Underdetermined));
fprintf('SVD Underdetermined Average Time: %.4fs\n', mean(tSVD_Underdetermined));

figure;
loglog(dimensionByStepQRSquare.RowSize, errorSolhQRSquare, '-o');
hold on;
loglog(dimensionByStepQRSquare.RowSize, errorSolgQRSquare, '-o');
loglog(dimensionByStepQRSquare.RowSize, errorSolmQRSquare, '-o');
hold off;
legend('Householder QR Method', 'Givens QR Method', 'MATLAB qr Method');
title('Solution Error vs Matrix Size (QR) - Square Matrices');
xlabel('Matrix Size');
ylabel('Solution Error');
grid on;

% Crea il grafico dei risultati per il metodo SVD per matrici quadrate
figure;
loglog(dimensionByStepSVDSquare.RowSize, errorSolmSVDSquare, '-o');
title('Solution Error vs Matrix Size (SVD) - Square Matrices');
xlabel('Matrix Size');
ylabel('Solution Error');
grid on;

% Crea il grafico dei risultati per il metodo QR per sistemi sovradeterminati
figure;
loglog(dimensionByStepQROverdetermined.RowSize, errorSolhQROverdetermined, '-o');
hold on;
loglog(dimensionByStepQROverdetermined.RowSize, errorSolgQROverdetermined, '-o');
loglog(dimensionByStepQROverdetermined.RowSize, errorSolmQROverdetermined, '-o');
hold off;
legend('Householder QR Method', 'Givens QR Method', 'MATLAB qr Method');
title('Solution Error vs Matrix Size (QR) - Overdetermined Systems');
xlabel('Matrix Size');
ylabel('Solution Error');
grid on;

% Crea il grafico dei risultati per il metodo SVD per sistemi sovradeterminati
figure;
loglog(dimensionByStepSVDOverdetermined.RowSize, errorSolmSVDOverdetermined, '-o');
title('Solution Error vs Matrix Size (SVD) - Overdetermined Systems');
xlabel('Matrix Size');
ylabel('Solution Error');
grid on;

% Crea il grafico dei risultati per il metodo QR per sistemi sottodeterminati
figure;
loglog(dimensionByStepQRUnderdetermined.RowSize, errorSolhQRUnderdetermined, '-o');
hold on;
loglog(dimensionByStepQRUnderdetermined.RowSize, errorSolgQRUnderdetermined, '-o');
loglog(dimensionByStepQRUnderdetermined.RowSize, errorSolmQRUnderdetermined, '-o');
hold off;
legend('Householder QR Method', 'Givens QR Method', 'MATLAB qr Method');
title('Solution Error vs Matrix Size (QR) - Underdetermined Systems');
xlabel('Matrix Size');
ylabel('Solution Error');
grid on;

% Crea il grafico dei risultati per il metodo SVD per sistemi sottodeterminati
figure;
loglog(dimensionByStepSVDUnderdetermined.RowSize, errorSolmSVDUnderdetermined, '-o');
title('Solution Error vs Matrix Size (SVD) - Underdetermined Systems');
xlabel('Matrix Size');
ylabel('Solution Error');
grid on;

%% Funzioni

function [x_hat, t] = ComputeExecutionTime(QRmethod, A, b)
    tic;
    if isequal(QRmethod, @svd)
        [U, S, V] = QRmethod(A);
        x_hat = V * (S \ (U' * b));
    else
        [Q, R] = QRmethod(A);
        y = Q' * b;
        x_hat = R \ y;
    end
    t = toc;
end

function [errorSol, t, errorQR, errorQ, errorSolOriginal] = ComputeErrors(QRmethod, A, x, b)
    [m, n] = size(A);
    k = rank(A);
    tic;
    if isequal(QRmethod, @svd)
        [U, S, V] = QRmethod(A);
        if m == n
            % Case of square matrix
            x_hat = V * (S \ (U' * b));
        elseif m > n
            % Case of overdetermined system
            if k == n
                % Full rank
                x_hat = V(:, 1:n) * (S(1:n, 1:n) \ (U(:, 1:n)' * b));
            else
                % Rank deficient
                x_hat = V(:, 1:k) * (S(1:k, 1:k) \ (U(:, 1:k)' * b));
            end
        elseif m < n
            % Case of underdetermined system
            if k == m
                % Full rank
                x_hat = V(:, 1:m) * (S(1:m, 1:m) \ (U(:, 1:m)' * b));
            else
                % Rank deficient
                x_hat = V(:, 1:k) * (S(1:k, 1:k) \ (U(:, 1:k)' * b));
            end
        end
        errorQR = NaN;  % Indicate that errorQR and errorQ are not applicable for SVD
        errorQ = NaN;
    else
        [Q, R] = QRmethod(A);
        y = Q' * b;
        x_hat = R \ y;
        errorQR = norm(A - Q * R, 'fro') / norm(A, 'fro');
        errorQ = norm(Q'*Q - eye(size(Q,2)), 'fro');
    end
    t = toc;

    errorSol = norm(b - A * x_hat) / norm(b);
    errorSolOriginal = norm(x - x_hat, Inf);  % Changed to infinity norm
end

function [Q, R] = GivensQR(A)
    [m, n] = size(A);
    Q = eye(m);
    R = A;

    for j = 1:n
        for i = m:-1:(j + 1)
            G = eye(m);
            [c, s] = GivensRotation(R(i - 1, j), R(i, j));
            G([i - 1, i], [i - 1, i]) = [c, s; -s, c];
            R = G' * R;
            Q = Q * G;
        end
    end
end

function [c, s] = GivensRotation(a, b)
    if b == 0
        c = 1;
        s = 0;
    else
        if abs(b) > abs(a)
            r = a / b;
            s = 1 / sqrt(1 + r^2);
            c = s * r;
        else
            r = b / a;
            c = 1 / sqrt(1 + r^2);
            s = c * r;
        end
    end
end

function [Q, R] = HouseHolderQR(A)
    [m, n] = size(A);
    Q = eye(m);
    R = A;

    for j = 1:n
        [v, beta] = HouseHolderVector(R(j:m, j));
        R(j:m, j:n) = R(j:m, j:n) - beta * v * (v' * R(j:m, j:n));
        Q(j:m, :) = Q(j:m, :) - beta * v * (v' * Q(j:m, :));
    end
end

function [v, beta] = HouseHolderVector(x)
    v = x;
    sigma = v(1) / abs(v(1));
    v = v / sigma;
    v(1) = 1;
    beta = 2 / (v' * v);
end

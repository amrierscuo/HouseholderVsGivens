%% File: QRvsSVD\QRvsSVD.m

% Genera random data
m = 20;
n = 10;
times = 5;  % Esegui solo 5 iterazioni

% Inizializza le variabili per salvare i risultati:
% quadrate, sovradimensionate e sottodimensionate
arrayIndex = 1:times;
dimRowSquare = ones(1, times);
dimColumnSquare = ones(1, times);
errorSolmQRSquare = ones(1, times);
errorSolhQRSquare = ones(1, times);
errorSolgQRSquare = ones(1, times);
errorSolmSVDSquare = ones(1, times);
errorSolhSVDSquare = ones(1, times);
errorSolgSVDSquare = ones(1, times);

dimRowOverdetermined = ones(1, times);
dimColumnOverdetermined = ones(1, times);
errorSolmQROverdetermined = ones(1, times);
errorSolhQROverdetermined = ones(1, times);
errorSolgQROverdetermined = ones(1, times);
errorSolmSVDOverdetermined = ones(1, times);
errorSolhSVDOverdetermined = ones(1, times);
errorSolgSVDOverdetermined = ones(1, times);

dimRowUnderdetermined = ones(1, times);
dimColumnUnderdetermined = ones(1, times);
errorSolmQRUnderdetermined = ones(1, times);
errorSolhQRUnderdetermined = ones(1, times);
errorSolgQRUnderdetermined = ones(1, times);
errorSolmSVDUnderdetermined = ones(1, times);
errorSolhSVDUnderdetermined = ones(1, times);
errorSolgSVDUnderdetermined = ones(1, times);

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
    % Genera random data per matrici quadrate
    A_square = rand(m, m);
    x_square = rand(m, 1);
    b_square = A_square * x_square;

    % Genera random data per sistemi sovradeterminati
    A_overdetermined = rand(m, n);
    x_overdetermined = rand(n, 1);
    b_overdetermined = A_overdetermined * x_overdetermined;

    % Genera random data per sistemi sottodeterminati
    A_underdetermined = rand(m, n);
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
    [errorSolmQRSquare(i), ~, ~, ~, errorSolhQRSquare(i)] = ComputeErrors(@GivensQR, A_square, x_square, b_square);
    [~, ~, ~, ~, errorSolgQRSquare(i)] = ComputeErrors(@HouseHolderQR, A_square, x_square, b_square);

    % Calcola gli errori con il metodo SVD per matrici quadrate
    [errorSolmSVDSquare(i), ~, ~, ~, ~] = ComputeErrors(@svd, A_square, x_square, b_square);

    % Calcola gli errori con il metodo QR per sistemi sovradeterminati
    [errorSolmQROverdetermined(i), ~, ~, ~, errorSolhQROverdetermined(i)] = ComputeErrors(@GivensQR, A_overdetermined, x_overdetermined, b_overdetermined);
    [~, ~, ~, ~, errorSolgQROverdetermined(i)] = ComputeErrors(@HouseHolderQR, A_overdetermined, x_overdetermined, b_overdetermined);

    % Calcola gli errori con il metodo SVD per sistemi sovradeterminati
    [errorSolmSVDOverdetermined(i), ~, ~, ~, ~] = ComputeErrors(@svd, A_overdetermined, x_overdetermined, b_overdetermined);

    % Calcola gli errori con il metodo QR per sistemi sottodeterminati
    [errorSolmQRUnderdetermined(i), ~, ~, ~, errorSolhQRUnderdetermined(i)] = ComputeErrors(@GivensQR, A_underdetermined, x_underdetermined, b_underdetermined);
    [~, ~, ~, ~, errorSolgQRUnderdetermined(i)] = ComputeErrors(@HouseHolderQR, A_underdetermined, x_underdetermined, b_underdetermined);

    % Calcola gli errori con il metodo SVD per sistemi sottodeterminati
    [errorSolmSVDUnderdetermined(i), ~, ~, ~, ~] = ComputeErrors(@svd, A_underdetermined, x_underdetermined, b_underdetermined);
    % Calcola i tempi di esecuzione con il metodo QR e SVD per matrici quadrate
    [~, tQR_Square(i)] = ComputeExecutionTime(@GivensQR, A_square, b_square);
    [~, tSVD_Square(i)] = ComputeExecutionTime(@svd, A_square, b_square);

    % Calcola i tempi di esecuzione con il metodo QR e SVD per sistemi sovradeterminati
    [~, tQR_Overdetermined(i)] = ComputeExecutionTime(@GivensQR, A_overdetermined, b_overdetermined);
    [~, tSVD_Overdetermined(i)] = ComputeExecutionTime(@svd, A_overdetermined, b_overdetermined);

    % Calcola i tempi di esecuzione con il metodo QR e SVD per sistemi sottodeterminati
    [~, tQR_Underdetermined(i)] = ComputeExecutionTime(@GivensQR, A_underdetermined, b_underdetermined);
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
loglog(dimensionByStepQRSquare.RowSize, errorSolmQRSquare, '-o');
hold on;
loglog(dimensionByStepQRSquare.RowSize, errorSolhQRSquare, '-o');
loglog(dimensionByStepQRSquare.RowSize, errorSolgQRSquare, '-o');
hold off;
legend('QR Method', 'Householder QR Method', 'Givens QR Method');
title('Solution Error vs Matrix Size (QR) - Square Matrices');
xlabel('Matrix Size');
ylabel('Solution Error');
grid on;

% Crea il grafico dei risultati per il metodo SVD per matrici quadrate
figure;
loglog(dimensionByStepSVDSquare.RowSize, errorSolmSVDSquare, '-o');
hold on;
loglog(dimensionByStepSVDSquare.RowSize, errorSolhSVDSquare, '-o');
loglog(dimensionByStepSVDSquare.RowSize, errorSolgSVDSquare, '-o');
hold off;
legend('SVD Method');
title('Solution Error vs Matrix Size (SVD) - Square Matrices');
xlabel('Matrix Size');
ylabel('Solution Error');
grid on;

% Crea il grafico dei risultati per il metodo QR per sistemi sovradeterminati
figure;
loglog(dimensionByStepQROverdetermined.RowSize, errorSolmQROverdetermined, '-o');
hold on;
loglog(dimensionByStepQROverdetermined.RowSize, errorSolhQROverdetermined, '-o');
loglog(dimensionByStepQROverdetermined.RowSize, errorSolgQROverdetermined, '-o');
hold off;
legend('QR Method', 'Householder QR Method', 'Givens QR Method');
title('Solution Error vs Matrix Size (QR) - Overdetermined Systems');
xlabel('Matrix Size');
ylabel('Solution Error');
grid on;

% Crea il grafico dei risultati per il metodo SVD per sistemi sovradeterminati
figure;
loglog(dimensionByStepSVDOverdetermined.RowSize, errorSolmSVDOverdetermined, '-o');
hold on;
loglog(dimensionByStepSVDOverdetermined.RowSize, errorSolhSVDOverdetermined, '-o');
loglog(dimensionByStepSVDOverdetermined.RowSize, errorSolgSVDOverdetermined, '-o');
hold off;
legend('SVD Method');
title('Solution Error vs Matrix Size (SVD) - Overdetermined Systems');
xlabel('Matrix Size');
ylabel('Solution Error');
grid on;

% Crea il grafico dei risultati per il metodo QR per sistemi sottodeterminati
figure;
loglog(dimensionByStepQRUnderdetermined.RowSize, errorSolmQRUnderdetermined, '-o');
hold on;
loglog(dimensionByStepQRUnderdetermined.RowSize, errorSolhQRUnderdetermined, '-o');
loglog(dimensionByStepQRUnderdetermined.RowSize, errorSolgQRUnderdetermined, '-o');
hold off;
legend('QR Method', 'Householder QR Method', 'Givens QR Method');
title('Solution Error vs Matrix Size (QR) - Underdetermined Systems');
xlabel('Matrix Size');
ylabel('Solution Error');
grid on;

% Crea il grafico dei risultati per il metodo SVD per sistemi sottodeterminati
figure;
loglog(dimensionByStepSVDUnderdetermined.RowSize, errorSolmSVDUnderdetermined, '-o');
hold on;
loglog(dimensionByStepSVDUnderdetermined.RowSize, errorSolhSVDUnderdetermined, '-o');
loglog(dimensionByStepSVDUnderdetermined.RowSize, errorSolgSVDUnderdetermined, '-o');
hold off;
legend('SVD Method');
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
    tic;
    [Q, R] = QRmethod(A);
    y = Q' * b;
    x_hat = R \ y;
    t = toc;

    errorSol = norm(b - A * x_hat) / norm(b);
    errorQR = norm(A - Q * R, 'fro') / norm(A, 'fro');
    errorQ = norm(Q(:,1:size(A,2))'*Q(:,1:size(A,2)) - eye(size(A,2)), 'fro');
    errorSolOriginal = norm(x - x_hat) / norm(x);
end

function [Q, R] = GivensQR(A)
    [m, n] = size(A);    % Implementazione del metodo di Givens per la fattorizzazione QR.
    Q = eye(m);
    R = A;

    for j = 1:n
        for i = m:-1:(j+1)
            G = eye(m);
            [c, s] = givensrotation(R(i-1,j), R(i,j));
            G([i-1, i],[i-1, i]) = [c -s; s c]; 
            R = G'*R;
            Q = Q*G;
        end
    end
end

function [c, s] = givensrotation(a, b)
    if b == 0
        c = 1;
        s = 0;
    else
        if abs(b) > abs(a)
            r = a / b;
            s = 1 / sqrt(1 + r^2);
            c = s*r;
        else
            r = b / a;
            c = 1 / sqrt(1 + r^2);
            s = c*r;
        end
    end
end

function [Q, R] = HouseHolderQR(A)
    [m, n] = size(A);    % Implementazione del metodo di Householder per la fattorizzazione QR.
    Q = eye(m);
    R = A;

    for k = 1:n
        x = R(k:m,k);
        e = zeros(length(x),1);
        e(1) = norm(x);
        u = sign(x(1))*e + x;
        v = u / norm(u);
        R(k:m,:) = R(k:m,:) - 2 * v * (v' * R(k:m,:));
        Q(:,k:m) = Q(:,k:m) - 2 * (Q(:,k:m) * v) * v';
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

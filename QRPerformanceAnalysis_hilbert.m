% File: QRPerformanceAnalysis.m
% Init
sizes = 10:10:100;  % Dimensioni delle matrici
nSizes = length(sizes);

% Preallocazione
householderTime = zeros(1, nSizes);
givensTime = zeros(1, nSizes);
condNum = zeros(1, nSizes);
normAQR = zeros(2, nSizes);
normQQT = zeros(2, nSizes);
absError = zeros(2, nSizes);

% Ciclo attraverso le dimensioni delle matrici
for i = 1:nSizes
    % Genera una matrice di Hilbert
    A = hilb(sizes(i));
    
    fprintf('\nMatrice A di dimensione %d:\n', sizes(i));
    disp(A);
    
    % Householder QR
    tic;
    [Q, R] = HouseHolderQR(A);
    householderTime(i) = toc;
    
    % Calcoli per Householder
    A_ricostruita = Q * R;
    normAQR(1, i) = norm(A - A_ricostruita, inf);
    normQQT(1, i) = norm(Q*Q' - eye(size(Q,1)), inf);
    absError(1, i) = norm(A - A_ricostruita, 'fro');
    
    % Stampa i risultati per Householder
    fprintf('\nDimensione: %d (Householder)\n', sizes(i));
    fprintf('Tempo di esecuzione: %.16f secondi\n', householderTime(i));
    fprintf('A ricostruita:\n');
    disp(A_ricostruita);
    fprintf('A - A_ricostruita:\n');
    disp(A - A_ricostruita);
    fprintf('Q*Q^T - I:\n');
    disp(Q*Q' - eye(size(Q,1)));
    fprintf('||A-QR||∞: %.16f\n', normAQR(1, i));
    fprintf('||QQ^T-I||∞: %.16f\n', normQQT(1, i));
    fprintf('Errore assoluto di ricostruzione: %.16f\n', absError(1, i));
    
    % Givens QR
    tic;
    [Q, R] = GivensQR(A);
    givensTime(i) = toc;
    
    % Calcoli per Givens
    A_ricostruita = Q * R;
    normAQR(2, i) = norm(A - A_ricostruita, inf);
    normQQT(2, i) = norm(Q*Q' - eye(size(Q,1)), inf);
    absError(2, i) = norm(A - A_ricostruita, 'fro');
    
    % Stampa i risultati per Givens
    fprintf('\nDimensione: %d (Givens)\n', sizes(i));
    fprintf('Tempo di esecuzione: %.16f secondi\n', givensTime(i));
    fprintf('A ricostruita:\n');
    disp(A_ricostruita);
    fprintf('A - A_ricostruita:\n');
    disp(A - A_ricostruita);
    fprintf('Q*Q^T - I:\n');
    disp(Q*Q' - eye(size(Q,1)));
    fprintf('||A-QR||∞: %.16f\n', normAQR(2, i));
    fprintf('||QQ^T-I||∞: %.16f\n', normQQT(2, i));
    fprintf('Errore assoluto di ricostruzione: %.16f\n', absError(2, i));
    
    % Calcolo dell'errore di ricostruzione relativo
    relError = norm(A - A_ricostruita, 'fro') / norm(A, 'fro');
    
    fprintf('Errore relativo di ricostruzione: %.16f\n', relError);
    
    % Calcolo del numero di condizionamento
    condNum(i) = cond(A);
    
    % Analisi delle dimensioni delle matrici
    fprintf('\n--- Analisi delle dimensioni delle matrici ---\n');
    fprintf('Dimensione: %d\n', sizes(i));
    fprintf('Tempo di esecuzione (Householder): %.16f secondi\n', householderTime(i));
    fprintf('Tempo di esecuzione (Givens): %.16f secondi\n', givensTime(i));
    fprintf('Numero di condizionamento: %.16f \n', condNum(i));
    fprintf('Media ||A-QR||∞: %e\n', mean(normAQR(:, i)));
    fprintf('Media ||QQ^T-I||∞: %e\n', mean(normQQT(:, i)));
    fprintf('Media Errore assoluto di ricostruzione: %e\n', mean(absError(:, i)));
    
    % Confronto con altre tecniche
    % Eseguiamo altre tecniche di fattorizzazione QR, decomposizione LU, SVD, ecc.
    % Confrontiamo i risultati ottenuti
    
    % Scalabilità
    % Eseguiamo lo script su matrici di dimensioni sempre maggiori, hilbert?
    % Osserviamo come i tempi di esecuzione e altre misure si comportano

end

% Crea un grafico dei tempi di esecuzione
figure;
plot(sizes, householderTime, 'r', 'LineWidth', 2);
hold on;
plot(sizes, givensTime, 'b', 'LineWidth', 2);
legend('Householder', 'Givens', 'Location', 'NorthWest');
xlabel('Dimensione della matrice');
ylabel('Tempo di esecuzione (secondi)');
title('Tempi di esecuzione Householder e Givens');
grid on;


% Crea un grafico per l'errore di ricostruzione relativo
figure;
plot(sizes, normAQR(1,:), 'r', 'LineWidth', 2);
hold on;
plot(sizes, normAQR(2,:), 'b', 'LineWidth', 2);
legend('Householder', 'Givens', 'Location', 'NorthWest');
xlabel('Dimensione della matrice');
ylabel('Errore di ricostruzione relativo');
title('Errore di ricostruzione relativo Householder e Givens');
grid on;

% Mostra tutti i grafici
hold off;

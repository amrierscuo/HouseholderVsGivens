% Caricamento e analisi dei set di dati
datasets = {'data/conchiglia.mat', 'data/superficie1.mat'};
figure('Name', 'Error Results', 'Position', [100, 100, 1200, 600]);

for i = 1:length(datasets)
    disp(['Loading dataset: ', datasets{i}]);
    load(datasets{i});
    
    % Trasposizione dei dati
    disp('Transposing the data...');
    Ltrasp = L';
    Mtrasp = M';
    
    % Calcolo QR tramite rotazioni di Householder
    disp('Calculating QR decomposition using Householder rotations...');
    tic
    [Qh, Rh] = HouseHolderQR(Ltrasp);
    Nqrtesth = QRSystemResolution(Qh, Rh, Mtrasp);
    timeH = toc;
    disp(['Time taken for Householder method: ', num2str(timeH)]);
    
    % Calcolo QR tramite rotazioni di Givens
    disp('Calculating QR decomposition using Givens rotations...');
    tic
    [Qg, Rg] = GivensQR(Ltrasp);
    Nqrtestg = QRSystemResolution(Qg, Rg, Mtrasp);
    timeG = toc;
    disp(['Time taken for Givens method: ', num2str(timeG)]);
    
    % Calcolo diretto MATLAB
    disp('Calculating using MATLAB direct method...');
    tic
    Nmatlab = Ltrasp \ Mtrasp;
    timeM = toc;
    disp(['Time taken for MATLAB direct method: ', num2str(timeM)]);
    
    % Calcolo tramite SVD
    disp('Calculating using Singular Value Decomposition...');
    tic
    [Usvd, Ssvd, Vsvd] = svd(Ltrasp);
    Nqrtestsvd = QRSystemResolution(Usvd, Ssvd * Vsvd', Mtrasp);
    timeSVD = toc;
    disp(['Time taken for Singular Value Decomposition method: ', num2str(timeSVD)]);
    
    % Calcolo degli errori
    disp('Calculating errors...');
    errorFatt = [norm(Nqrtesth' * L - M, 'fro'); norm(Nqrtestg' * L - M, 'fro'); norm(Nmatlab' * L - M, 'fro'); norm(Nqrtestsvd' * L - M, 'fro')];
    errorN = [norm(Nqrtesth, 'fro'); norm(Nqrtestg, 'fro'); norm(Nmatlab, 'fro'); norm(Nqrtestsvd, 'fro')];
    errorNRel = [norm(Nqrtesth, 'fro') / (norm(Nqrtesth, 'fro') + eps); norm(Nqrtestg, 'fro') / (norm(Nqrtestg, 'fro') + eps); norm(Nmatlab, 'fro') / (norm(Nmatlab, 'fro') + eps); norm(Nqrtestsvd, 'fro') / (norm(Nqrtestsvd, 'fro') + eps)];
    errorNmathRel = [norm(Nmatlab - Nqrtesth, 'fro') / (norm(Nmatlab, 'fro') + eps); norm(Nmatlab - Nqrtestg, 'fro') / (norm(Nmatlab, 'fro') + eps); NaN; norm(Nmatlab - Nqrtestsvd, 'fro') / (norm(Nmatlab, 'fro') + eps)];
    
    timeAll = [timeH; timeG; timeM; timeSVD];

    data = [errorFatt, errorN, errorNRel, errorNmathRel, timeAll];
    
    % Visualizzazione dei risultati
    disp('Displaying results...');
    resultTable = array2table(data, 'VariableNames', {'FactError', 'AbsNormError', 'RelNormError', 'RelNormErrorMatlab', 'ExecutionTime'}, 'RowNames', {'Householder', 'Givens', 'Matlab Direct', 'SVD'});
    disp(resultTable);
    
    % Plotting results
    subplot(2, length(datasets), i);
    loglog(data(:, 1), '-o', 'LineWidth', 2);
    hold on;
    loglog(data(:, 2), '-s', 'LineWidth', 2);
    loglog(data(:, 3), '-d', 'LineWidth', 2);
    loglog(data(:, 4), '-^', 'LineWidth', 2);
    hold off;
    title('Error Results');
    xlabel('Method');
    ylabel('Error');
    xticks(1:length(datasets));
    xticklabels({'Householder', 'Givens', 'Matlab Direct', 'SVD'});
    legend('FactError', 'AbsNormError', 'RelNormError', 'RelNormErrorMatlab');
    grid on;
    set(gca, 'FontSize', 12);
    
    % Displaying all images in a single plot
    subplot(2, length(datasets), length(datasets) + i);
    if strcmp(datasets{i}, 'data/conchiglia.mat')
        montage(reshape(Mtrasp, r, s, []), 'DisplayRange', []);
    else
        montage(reshape(M, r, s, []), 'DisplayRange', []);
    end
    colormap gray;
    axis image;
    title(['Dataset: ', datasets{i}]);
end

% Definizione delle funzioni (riportate le stesse funzioni dal codice originale)
function [Q, R] = HouseHolderQR(A)
    [m, n] = size(A);
    R = A;
    Q = eye(m);
    for k = 1:n
        x = R(k:m,k);
        e = zeros(length(x),1);
        e(1) = 1;
        u = sign(x(1))*norm(x)*e + x;
        v = u / norm(u);
        R(k:m,k:n) = R(k:m,k:n) - 2 * (v * (v' * R(k:m,k:n)));
        Q(k:m,:) = Q(k:m,:) - 2 * (v * (v' * Q(k:m,:)));
    end
    Q = Q';
end

function [Q, R] = GivensQR(A)
    [m,n] = size(A);
    Q = eye(m);
    R = A;
    for j = 1:n
        for i = m:-1:(j+1)
            G = eye(m);
            [c,s] = GivensRotation( R(i-1,j),R(i,j) );
            G([i-1, i],[i-1, i]) = [c -s; s c];
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

function x = QRSystemResolution(Q, R, b)
    y = Q'*b;
    x = R\y;
end

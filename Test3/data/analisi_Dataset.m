% analisi_Dataset.m

% Carica il dataset 'superficie1.mat'
load superficie1.mat

% Estrai le dimensioni delle immagini nel dataset
[n, m] = size(M);
r = 101; % Numero di righe in ogni immagine
s = 101; % Numero di colonne in ogni immagine

disp(['Dimensioni di M (superficie1): ', num2str(n), ' x ', num2str(m)]);
disp(['Dimensioni dell''immagine (superficie1): ', num2str(r), ' x ', num2str(s)]);
disp(['Numero di immagini in superficie1: ', num2str(m)]);

% Visualizza le immagini nel dataset 'superficie1.mat'
figure();
for i = 1:m
    I = mat2gray(reshape(M(:, i), r, s));
    subplot(2, 4, i);
    imshow(I);
end
title('Immagini originali di Superficie1');

% Elabora o modifica i dati nel dataset 'superficie1.mat', se necessario
% Esempio: applica una semplice trasformazione lineare ai dati
M_processed = M * 2;

% Salva le modifiche apportate ai dati in un nuovo file .mat
save('superficie1_processed.mat', 'M_processed');

% Ripeti lo stesso procedimento per il dataset 'conchiglia.mat'
load conchiglia.mat

[n, m] = size(M);
r = 705; % Numero di righe in ogni immagine
s = 885; % Numero di colonne in ogni immagine

disp(['Dimensioni di M (conchiglia): ', num2str(n), ' x ', num2str(m)]);
disp(['Dimensioni dell''immagine (conchiglia): ', num2str(r), ' x ', num2str(s)]);
disp(['Numero di immagini in conchiglia: ', num2str(m)]);

figure();
for i = 1:m
    I = mat2gray(reshape(M(:, i), r, s));
    subplot(5, 4, i);
    imshow(I);
end
title('Immagini originali conchiglia');

% Elabora o modifica i dati nel dataset 'conchiglia.mat', se necessario
% Esempio: applica una semplice trasformazione lineare ai dati
M_processed = M * 2;

% Salva le modifiche apportate ai dati in un nuovo file .mat
save('conchiglia_processed.mat', 'M_processed');

% Elimina i file 'conchiglia_processed.mat' e 'superficie1_processed.mat'
delete('conchiglia_processed.mat');
delete('superficie1_processed.mat');

Utilizzo delle funzioni HouseHolderQR e GivensQR

Queste due funzioni servono per calcolare la fattorizzazione QR di una matrice A utilizzando metodi diversi: Householder e Givens.

Per utilizzare queste funzioni, seguire questi passaggi:
-------------------
1. apri gli script HouseHolderQR e GivensQR in MATLAB nella cartella corrente.
-------------------
2. Crea o importa la matrice A da fattorizzare(r>c). Ad esempio:

   A = [1 2 3; 4 5 6; 7 8 9];
   A =  rand(5);
   A =  rand(6,5);
-------------------
3. Chiama la funzione HouseHolderQR per eseguire la fattorizzazione QR con il metodo di Householder:

   [Q, R] = HouseHolderQR(A);

   Ora Q e R contengono rispettivamente la matrice ortogonale Q e la matrice triangolare superiore R risultanti dalla fattorizzazione QR.
-------------------
4. Oppure chiama la funzione GivensQR per eseguire la fattorizzazione QR con il metodo di Givens:

   [Q, R] = GivensQR(A);

   Anche in questo caso, Q e R contengono rispettivamente la matrice ortogonale Q e la matrice triangolare superiore R risultanti dalla fattorizzazione QR.
-------------------
5. runna questo nella command window per verificare il funzionamento
disp('Matrice Q:')
disp(Q)
disp('Matrice R:')
disp(R)
-------------------
6. verifica
A_ricostruita = Q * R;
differenza = norm(A - A_ricostruita);

disp('Differenza tra la matrice originale A e la matrice ricostruita:');
disp(differenza);
disp(A_ricostruita)

Entrambe le funzioni restituiscono le matrici Q e R, che possiamo usare per i confronti.
Se la differenza è vicina a zero allora la funzione è corretta con un errore numerico accettabile
-------------------
7. stampa tutto
disp('Matrice A:')
disp(A)
disp('Matrice Q:')
disp(Q)
disp('Matrice R:')
disp(R)
disp('Matrice A ricostruita:')
disp(A_ricostruita)
disp('Differenza tra la matrice originale A e la matrice ricostruita:');
disp(differenza);

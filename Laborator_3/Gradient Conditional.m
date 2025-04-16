clc, clear
n = 70;                              % valoarea de redimensionare a img.
A_bar = imread("bunny.jpg");         % citirea img
A_bar = rgb2gray(A_bar);             % convertire in alb-negru
A_bar =im2double(A_bar);             % convertirea in valori double
A_bar = imresize(A_bar,[n n]);       % redimensionarea  pozei  
figure(1)
imshow(A_bar);                           % afisarea pozei initiale 
title("Imaginea originala")

nrintraricunoscute=3000;             % setam un numar de intrari cunoscute 
rPerm = randperm(n*n);                          %generarea random a indicilor pentru intrarile cunoscuti
omega = sort(rPerm(1 : nrintraricunoscute));    %intrarile care se cunosc
A = nan(n); A(omega) = A_bar(omega); 
figure(2)
title("Imagine cu pixeli lipsa")
imshow(A) 

A_k = randn(n, n);
A_k1 = zeros(n, n);
G = zeros(n, n); %gradientul functiei de eroare
cond_oprire = [];

dR = 70;              % coeficient de scalare
eroareMP = 1e-5;          % pentru MP
eroareGC = 1e-5;          % pentru GC
itermax = 10000;       % nr max iteratii
iterr = 1;
eroare_parcurs = 1;

while eroare_parcurs >= eroareGC && iterr <= itermax
    cond_oprire = [cond_oprire, eroare_parcurs];
    A_k = A_k1;
    
    G(omega) = A_k(omega) - A(omega);

    v1 = randn(n, 1); 
    v1 = v1 / norm(v1);
    aproxl = v1; %aproximarea lui v1
    eroare = 1;
    alfa1 = 1;
    TtT = G' * G;

    while eroare >= eroareMP
        v2 = TtT * v1 / norm(TtT * v1);
        alfa2 = norm(G * v1);
        aprox2 = G * v2 / alfa2;
        eroare = abs(aprox2' * G * v2 - alfa1);
        alfa1 = alfa2;
        v1 = v2;
        aproxl = aprox2;
    end

    S = -dR * aproxl * v1';
    a = 2 / (iterr + 2);
    A_k1 = A_k + a * (S - A_k);
    
    eroare_parcurs = norm(A_k1 - A_k, 'fro');
    iterr = iterr + 1;
end

figure(3)
imshow(A_k1)

figure;
semilogy(cond_oprire);
xlabel('Iteratii');
ylabel('norma');
title('GC');
grid on;






clc; clear
n = 70;                              % valoarea de redimensionare a img.

A_bar = imread("bunny.jpg");         % citirea img
A_bar = rgb2gray(A_bar);             % convertire in alb-negru
A_bar = im2double(A_bar);            % convertirea in valori double
A_bar = imresize(A_bar, [n n]);      % redimensionarea  pozei  

figure(1)
imshow(A_bar);                       % afisarea pozei initiale 
title("Imaginea originala")

nrintraricunoscute = 3000;                                  % setam un numar de intrari cunoscute 
rPerm = randperm(n * n);                                    % generarea random a indicilor pentru intrarile cunoscute
omega = sort(rPerm(1 : nrintraricunoscute));                % intrarile care se cunosc

A = nan(n); 
A(omega) = A_bar(omega);                                    

figure(2)
imshow(A) 
title("Imagine cu pixeli lipsa")

epsilon = 1e-3;                      
c = 10;                              % constantă pas α_k = c / k

A_nou = randn(n, n);                 
A_nou(omega) = A(omega);            

iter = 10000;                        
i = 1;
eroare = 1;

conditie_oprire = [];               % vector pentru memorarea normei la fiecare pas

while eroare >= epsilon && i < iter
    conditie_oprire = [conditie_oprire, eroare];  

    A_vechi = A_nou;                            

    [U, ~, V] = svd(A_vechi);                    
    k = i + 1;
    alfa = c / k;                                

    A_nou = A_vechi - alfa * U * V';             
    A_nou(omega) = A(omega);                     % valorile deja cunoscute sunt folosite pentru a calcula norma

    eroare = norm(A_nou - A_vechi, 'fro');       
    i = i + 1;                                  
end

figure(3)
imshow(A_nou)
title("Imagine reconstruita cu Gradient Proiectat")

figure;
semilogy(conditie_oprire);                     
xlabel('Iteratii');
ylabel('Norma');
title('GP');
grid on;

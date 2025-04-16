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

cvx_begin
    variable Y(n, n)
    minimize(norm_nuc(Y)) %norma nucleara
    subject to %supus la
        Y(omega) == A(omega)
cvx_end

figure(3)
imshow(Y)
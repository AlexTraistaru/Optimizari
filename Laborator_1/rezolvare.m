load 'data.mat';
[n, N]=size(R_train); 
R=R_train;
k=100; 
U=rand(n, k); 
M=rand(N, k); 
U_c=U; 
M_c=M;
epsilon=0.01;
alfa=0.003; 
maxim_iteratii=1000; 
iter=0; 
f1=1; 
f2=0; 
functii_obiectiv=zeros(maxim_iteratii, 1); 
eroare_relativa=zeros(maxim_iteratii, 1); 
while abs(f2-f1)>epsilon && iter<maxim_iteratii
    f1=f2;
    f2=0;
    for i=1:n
        elemente=find(R(i, :));
        for j=1:length(elemente)
            grad=R(i, elemente(j))-U(i, :)*(M(elemente(j), :))';
            for p=1:k
                U_c(i, p)=U(i, p)+alfa*grad*M(elemente(j), p);
                M_c(elemente(j), p)=M(elemente(j), p)+alfa*grad*U(i, p);
            end
            f2=f2+(1/2)*(R(i, elemente(j))-U_c(i, :)*(M_c(elemente(j), :))')^2;
        end
    end
    iter=iter+1;
    U=U_c;
    M=M_c;
    %disp(abs(f1-f2));
    functii_obiectiv(iter)=f2;
    f2=f2/2;
    eroare_relativa(iter)=abs(f2-f1);
end



x=1:iter;
figure;
plot(x, functii_obiectiv(1:iter, 1), 'g', 'LineWidth', 2);
grid on;
title('Graficul functiei obiectiv');

figure;
semilogy(x, eroare_relativa(1:iter, 1), 'k', 'LineWidth', 2);
grid on;
title('Graficul erorii relative');


utilizatori1=randi(500);
utilizatori2=randi(500);
produs_scalar=U(utilizatori1, :)*U(utilizatori2, :)';
produs_norm=norm(U(utilizatori1, :))*norm(U(utilizatori2, :));
cosinus=produs_scalar/produs_norm;

fprintf('Similaritatea cosinusoidala intre utilizatorii %d si %d este %.4f\n', utilizatori1, utilizatori2, cosinus);


clc, clear, close all
%% Descrierea sistemului pendulul inversat pe un carucior
M = 0.5;   %Masa cartului
m = 0.2;  %Masa pendulului
b = 0.1;  %coef de frecare pt carucior
I = 0.006;%momentul de inertie
g = 9.8;  %acceleratia gravitationala
l = 0.3;  %lungimea lantului
p = I*(M+m)+M*m*l^2; %denominator for the A and B matrices
A = [0 1 0 0;
     0 -(I+m*l^2)*b/p (m^2*g*l^2)/p 0;
     0 0 0 1;
     0 -(m*l*b)/p m*g*l*(M+m)/p 0];
B = [0;
     (I+m*l^2)/p;
     0;
     m*l/p];
    
Q = [10 0 0 0;
     0 1 0 0;
     0 0 10 0;
     0 0 0 1];
R = 1;
N = 3; % orizontul de predictie
% Constrangerea de tip box pentru intrare u
ub = 3; 
lb = -3;
% Starea initiala a sistemului
z0 = [2; 1 ; 0.3; 0];
% Referinta de urmarit
x_ref = 0;
z_ref = [0;0; 0; 0];
maxIter = 20;
steps = 0;
u =[];
z=z0;
x = zeros(N,1); % warm start pentru metoda bariera

epsilon=0.0001;
sigma=0.6;
alfa=1;

%% Metoda bariera cu Newton
while steps < maxIter  
   [H,q,C,d]= denseMPC(A,B,Q,R,z0,N,ub,lb, z_ref, x_ref); 

   m_c = size(C,1);
   teta=1;
   while m_c*teta>=epsilon 
        gradient=H*x+q; %H e simetrica
        hessiana = H;
        constrangeri=C*x-d;
        for i=1:m_c
                gradient=gradient-teta*(C(i, :))'/constrangeri(i);
                hessiana = hessiana + teta*((C(i, :)' * C(i, :)) / (constrangeri(i))^2);
        end
        % Metoda Newton
        while norm(gradient)>=epsilon
            dir = inv(hessiana) * gradient; %calculez directia de deplasare
            %-inv(hessiana) e directia de coborare
            x = x - alfa * dir;
            constrangeri=C*x-d;
            gradient=H*x+q;
            for i=1:m_c
                gradient=gradient-teta*(C(i, :))'/constrangeri(i);
                hessiana = hessiana + teta*((C(i, :)' * C(i, :)) / (constrangeri(i))^2);
            end
        end
        teta=teta*sigma;
    end 

   u = [u,x(1,1)]; 
   z_new = A* z0 + B *x(1,1); 
   z0 = z_new;
   z = [z z0];
   steps = steps +1;
end

%% Grafic rezultate - Newton
figure(1)
subplot(2,2,1)
plot(z(1,:))
legend('x Pozitia caruciorului');
title('Pozitia caruciorului');

subplot(2,2,2)
plot(z(2,:))
legend('Viteza caruciorului');
title('Viteza caruciorului');

subplot(2,2,3)
plot(z(3,:))
legend('\phi Abaterea pendulului');
title('Abaterea pendulului');

subplot(2,2,4)
plot(z(4,:))
legend('Viteza unghiulara pendul');
title('Viteza unghiulara');

figure(2)
plot(u);
legend('u = Forta aplicata caruciorului');
title('Comanda aplicata');

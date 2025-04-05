
clc
clear
%% Citirea si preprocesare datelor
dataFolder = 'train';
categories = {'cat', 'dog'};

imds = imageDatastore(fullfile(dataFolder, categories), 'LabelSource', 'foldernames', 'IncludeSubfolders', true);

dataFolder = 'test';
imds_t = imageDatastore(fullfile(dataFolder, categories), 'LabelSource', 'foldernames', 'IncludeSubfolders', true);

data=[];                            % date de antrenare
data_t =[];                         % date de test
while hasdata(imds) 
    img = read(imds) ;              % citeste o imagine din datastore
    img = imresize(img, [227 227]);
    %figure, imshow(img); % decomentati pentru a vizualiza pozele din
    %pause                        %baza de dat 
    img = double(rgb2gray(img));
    data =[data, reshape(img, [], 1)];

end
eticheta = double(imds.Labels == 'cat'); % eticheta 1 - pisica, 0 - caine
while hasdata(imds_t) 
    img = read(imds_t) ;              % citeste o imagine din datastore
    img = imresize(img, [227 227]); 
    %figure, imshow(img);    % decomentati pentru a vizualiza pozele din
    % pause                  %baza de date
    img = double(rgb2gray(img));
    data_t =[data_t, reshape(img, [], 1)];
    
end
eticheta_t = double(imds_t.Labels == 'cat'); % eticheta 1 - pisica, 0 - caine
%Reducerea dimensiuni
data_pca = pca(data, 'NumComponents', 45)';
data_pca_t = pca(data_t, 'NumComponents', 45)';

clear categories imds imds_t img dataFolder data data_t

[n, N] = size(data_pca); % n - nr de caracteristici; N- nr de poze

w = randn(n,1);       

vpc = [];
vpi = [];

alpha = 0.9;
epsilon = 1e-8;
maxIter = 10000;

iter = 0;
wpc = w;
wpi = w;
while (norm(grad(wpc,eticheta,data_pca)) > epsilon && iter < maxIter)
    wpc = wpc - alpha * inv(Hess(wpc,eticheta,data_pca)) * grad(wpc,eticheta,data_pca)';
    vpc = [vpc, norm(grad(wpc,eticheta,data_pca))];
    iter = iter + 1;
end

iter = 0;

while (norm(grad(wpi,eticheta,data_pca)) > epsilon && iter < maxIter)
    f_alpha = @(alpha) obj(wpi - alpha * inv(Hess(wpi,eticheta,data_pca)) * (grad(wpi,eticheta,data_pca))', eticheta,data_pca);
    alpha = fminbnd(f_alpha, 0, 1); % caută pasul optim
    wpi = wpi - alpha * inv(Hess(wpi,eticheta,data_pca)) * (grad(wpi,eticheta,data_pca))';

    vpi = [vpi, norm(grad(wpi,eticheta,data_pca))];
    iter = iter + 1;
end

figure; 
plot(1:length(vpc), vpc, 'LineWidth', 2); hold on;
plot(1:length(vpi), vpi, 'LineWidth', 2);


function g = sigmoid(z)
    g = 1.0 ./ (1.0 + exp(-z));
end

function [A] = grad(w,y,x)
    N = size(x,2); 
    h = zeros(N,1);
    for i = 1:N
        h(i) = sigmoid(w'*x(:,i));
    end
    A = ((h-y)' * x') / N;
end

function [h] = Hess(w,y,x)
    m = size(x,2); %nr de exemple - coloane
    n = size(x,1); %nr de caracteristici - linii
    z = zeros(m,1); %in el pun hi(1-hi)
    for i = 1:m
        e = sigmoid(w'*x(:,i));
        z(i) = e*(1-e);
    end
    Q = diag(z);
    h = (x * Q * x') / m;
end

function g = obj(w,y,X)
    N = size(X,2);
    h = zeros(N,1);
    for i = 1:N
        h(i) = sigmoid(w'*X(:,i));
    end
    h = min(max(h, 1e-10), 1 - 1e-10); %din cauza lui log(0) 
    g = (-y' * log(h) - (1 - y)' * log(1 - h)) / N;
end


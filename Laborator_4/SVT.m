function [x] = SVT(A, a)
    [U, S, V] = svd(A, "econ");
    sigma = diag(S); % extrage valorile singulare
    sigma_th = max(sigma - a, 0); %fac maximul cu 0 in cazul in care am numere negative, valorile singulare trebuie sa fie pozitive
    S_th = diag(sigma_th); % reconstruim matricea
    x = U * S_th * V'; % reconstruiesc matricea
end

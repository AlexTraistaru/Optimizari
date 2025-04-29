function [x] = ST(A, a)
    x = sign(A) .* max(abs(A) - a, 0);
end

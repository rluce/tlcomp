function val = is_exact_toeplitz(A)
% val = is_exact_toeplitz(A) -- true/false, indicating wether A is a
% Toeplitz matrix.
diff = A - toeplitz(A(:,1), A(1,:));
val = isempty(nonzeros(diff(:)));

end
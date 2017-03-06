function val = is_exact_toeplitz(A)
% val = is_exact_toeplitz(A) -- true/false, indicating wether A is a
% Toeplitz matrix.

if isempty(A)
    val = true;
    return;
end

diff = A - toeplitz(A(:,1), A(1,:));
val = isempty(nonzeros(diff(:)));

end
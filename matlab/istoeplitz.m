function val = istoeplitz(A)
% val = istoeplitz(A) -- true/false, indicating whether A is an exact, full
% Toeplitz matrix.

if isempty(A)
    val = true;
    return;
end

diff = A - toeplitz(A(:,1), A(1,:));
val = isempty(nonzeros(diff(:)));

end

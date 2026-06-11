%return A.'(:)
function a = MatTranspVec(A)
    AT      = A.';
    a       = AT(:);
end
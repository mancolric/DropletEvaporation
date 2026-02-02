clear
clc

CRStoi          = cell(5,1);

for ii=1:5
    CRStoi{ii}  = (ii-1)*ones(3);
end

a               = CRStoi(cellfun(@(x) any(x, 'all'), CRStoi));

matrizResultado = prod(cat(3, a{:}), 3);
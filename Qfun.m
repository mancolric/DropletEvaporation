function Q=Qfun(model, T, rhoy)

        Q             = cell(model.nDiff,1);

        C             = cellfun(@(r, mw) r ./ (1e6*mw), rhoy, num2cell(model.species_mw(:)), 'UniformOutput', false);   %[mol/cm3]
        MatrixExpo    = model.ROrder((1:model.nFuels) + model.nInerts, min(1:model.nSpecies, model.nInerts+1)).';
        CROrder       = cellfun(@(c, exp) abs(c).^exp, repmat(C, 1, model.nFuels), num2cell(MatrixExpo), 'UniformOutput', false);
     
        idx_Arr = (1:model.nFuels) + model.nInerts;
        A1 = model.Arr(idx_Arr, 1);
        A2 = model.Arr(idx_Arr, 2);
        A3 = model.Arr(idx_Arr, 3);
        calc_omega = @(j) (A2(j) * T.^A3(j) .* exp(-A1(j) ./ (model.R .* T))) .* prod(cat(3, CROrder{:,j}), 3);
        omega = arrayfun(calc_omega, (1:model.nFuels)', 'UniformOutput', false);
       
        Omega3D = cat(3, omega{:});
        
        calc_Q = @(kk) model.species_mw(kk) .* 1e6 .* sum(Omega3D .* reshape(model.DiffStoi((1:model.nFuels)+model.nInerts, kk), 1, 1, []), 3);
        Q(1:model.nSpecies) = arrayfun(calc_Q, (1:model.nSpecies)', 'UniformOutput', false);
        Q{end}                  = zeros(size(Q{end-1}));

end
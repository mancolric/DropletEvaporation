function value = calcula_lambda_liq_mix(gota, Xf_molar, rho_liq_individual, k_liq_individual)

N_comb=size(gota.tipo_combustible,2); 


% Liquid mixture thermal conductivity: Li correlation
% Perry's Chemical Engineers' Handbook, 8th Edition, pp. 2-512
% for id=1:N_comb
% 
%     for jd=1:N_comb
%         sumandos_denominador_phi(jd)=Xf_molar(jd)/rho_liq_individual(jd);
%         sumatorio_sumandos_denominador_phi=sum(sumandos_denominador_phi);
%     end
% 
%     for jd=1:N_comb
%         numerador_phi(jd)=Xf_molar(jd)/rho_liq_individual(jd);
%         phi(jd)=numerador_phi(jd)/sumatorio_sumandos_denominador_phi;
%         numerador_phi(id)=Xf_molar(id)/rho_liq_individual(id);
%         phi(id)=numerador_phi(id)/sumatorio_sumandos_denominador_phi;
%         producto(id,jd)=phi(id)*phi(jd)*2*k_liq_individual(id)*k_liq_individual(jd)/(k_liq_individual(id)+k_liq_individual(jd));
%     end
%     suma_filas(id)=sum(producto(id,:));
% 
% end
% value=sum(suma_filas);
numerador_phi = Xf_molar ./ rho_liq_individual; % 2xN
denominador_phi = sum(numerador_phi, 1);        % 1xN
phi = numerador_phi ./ denominador_phi;         % 2xN

k = zeros(size(k_liq_individual));
t = zeros(N_comb,N_comb,size(k_liq_individual,2));
% Precalcular términos del producto para cada par (i,j) = (1,1), (1,2), (2,1), (2,2)
for id=1:N_comb
    k(id,:)   = k_liq_individual(id,:); % 1xN
    phi(id,:) = phi(id,:);
end

for id=1:N_comb
    for jd=1:N_comb
        t(id,jd,:) = phi(id,:) .* phi(jd,:) .* (2 .* k(id,:) .* k(jd,:) ./ (k(id,:) + k(jd,:)));
    end
end

value = sum(sum(t));
value = reshape(value, 1, []);


return


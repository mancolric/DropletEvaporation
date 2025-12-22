clear all
close all

load 'debug2.mat'

%Assembly full Jacobian:
% N_l     = N_ul+N_vl;
% N_g     = N_ug+N_vg;
% J_iv    = cat(1, J_l.iv, J_lz.iv, ...
%                 J_g.iv+N_l, J_gz.iv+N_l, ...
%                 J_zl.iv+N_l+N_g, J_zg.iv+N_l+N_g, J_z.iv+N_l+N_g);
% J_jv    = cat(1, J_l.jv, J_lz.jv+N_l+N_g, ...
%                 J_g.jv+N_l, J_gz.jv+N_l+N_g, ...
%                 J_zl.jv, J_zg.jv+N_l, J_z.jv+N_l+N_g);
% J_sv    = cat(1, J_l.sv, 1.0*J_lz.sv, ... 
%                 J_g.sv, 1.0*J_gz.sv, ...
%                 1.0*J_zl.sv, 1.0*J_zg.sv, J_z.sv);
% J       = sparse(J_iv(:), J_jv(:), J_sv(:), ...
%                 N_l+N_g+N_L+N_R+1, N_l+N_g+N_L+N_R+1);

 Srinv_factors   = 1.0 ./ [ ...
                                NF_omega_l*NF_l(1:end-1); NF_omega_l*NF_l(1)*NF_tau; ...
                                NF_omega_g*NF_g(1:end-1); NF_omega_g*NF_g(1)*NF_tau; ...
                                NF_g(end)*NF_g(1:end-1)*NF_tau; ...
                                max(DeltaT_0,300)*NF_tau; max(DeltaP_0,1e5)*NF_tau; NF_l(1)*NF_tau; NF_g(1)*NF_tau ];
          
J       = Srinv{3}*A_n{3}*Sy(block_ul(1):end, block_ul(1):end);
J2      = J;

fscaled     = rand(size(J2,1),1);

i0          = block_mesh_g(end);

figure()
spy(J2(block_z([5,6])-i0, :))
% J2(block_z(5)-i0, block_vl-i0)     = 0.0;


display(J2(block_z-i0, block_vl-i0))
display(J2(block_z-i0, block_z(end)-i0))


% J2          = Jscaled;
% for II=1:2
%     j1      = block_ul(1)-1+II*101;
%     j2      = block_z(1)-1+II;
%     J2(:,j1)    = Jscaled(:,j1) + Jscaled(:,j2);
%     J2(:,j2)    = Jscaled(:,j2)*(1e-7/5);
% end
% for II=1:6
%     j1      = block_ug(1)+(II-1)*101;
%     j2      = block_z(1)-1+2+II;
%     J2(:,j1)    = Jscaled(:,j1) + Jscaled(:,j2);
%     J2(:,j2)    = Jscaled(:,j2)*(1e-7/5);
% end

Jscaled_F   = LUFactorization(J2);
x           = LUSolve(Jscaled_F, fscaled);
disp('OK')
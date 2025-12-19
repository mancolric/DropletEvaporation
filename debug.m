clear all
close all

load 'test.mat'

blocks      = { block_ul, block_vl, block_ug, block_vg, block_z };
% blocks      = { block_ul(1:80), block_ul(81:160), block_vl, block_ug, block_vg, block_z(6) };
% blocks      = { block_ul(1:160), block_vl, block_ug, block_vg, block_z(6) };

% for II=length(blocks)
%     for JJ=1:length(blocks)
%         
%         block1      = blocks{II};
%         block2      = blocks{JJ};
% 
%         J1          = Jscaled(block1, block2);
%         J2          = Jest(block1, block2);
% 
%         [ii,jj,ss]  = find(J1);
%         [M,N]       = size(J1);
%         aux         = ii+(jj-1)*M;
% 
%         figure()
%         plot(J1(aux), '-xb')
%         hold on
%         plot(J2(aux), '+g')
%         title([II, JJ])
%         
%     end
% end

fscaled     = rand(size(Jscaled,1),1);

J2          = Jscaled;
for II=1:2
    j1      = block_ul(1)-1+II*101;
    j2      = block_z(1)-1+II;
    J2(:,j1)    = Jscaled(:,j1) + Jscaled(:,j2);
    J2(:,j2)    = Jscaled(:,j2)*(1e-7/5);
end
for II=1:6
    j1      = block_ug(1)+(II-1)*101;
    j2      = block_z(1)-1+2+II;
    J2(:,j1)    = Jscaled(:,j1) + Jscaled(:,j2);
    J2(:,j2)    = Jscaled(:,j2)*(1e-7/5);
end

Jscaled_F   = LUFactorization(J2);
x           = LUSolve(Jscaled_F, fscaled);

clear all
close all

load 'test.mat'

blocks      = { block_ul, block_vl, block_ug, block_vg, block_z };
% blocks      = { block_ul(1:80), block_ul(81:160), block_vl, block_ug, block_vg, block_z(6) };
% blocks      = { block_ul(1:160), block_vl, block_ug, block_vg, block_z(6) };

for II=1:length(blocks)
    for JJ=1:length(blocks)
        
        block1      = blocks{II};
        block2      = blocks{JJ};

        J1          = Jscaled(block1, block2);
        J2          = Jest(block1, block2);

        [ii,jj,ss]  = find(J1);
        [M,N]       = size(J1);
        aux         = ii+(jj-1)*M;

        figure()
        plot(J1(aux), '-xb')
        hold on
        plot(J2(aux), '+g')
        title([II, JJ])
        
    end
end
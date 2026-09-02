function [dout] = mleDecoder(d,constell)
%% 
% input:
% d- input symbols
% constell - constellation map
% output:
% dout - output symbols (projection on constellation)
%%
    M = size(d,1);
    N = size(d,2);
    dout = zeros(size(d));
    for ii = 1:M
        for jj = 1:N
            dis = abs(d(ii,jj) - constell);
            [~,ind] = min(dis);
            dout(ii,jj) = constell(ind);
        end
    end
end
function [Ke, Fe, S_local] = mitc4_element(X0, Xt, Vn0, V10, V20, Vnt, V1t, V2t, thick, C_mat)
% MITC4_ELEMENT Wrapper for classical MITC4 shell element
[Ke, Fe, S_local] = mitc4_plus_element(X0, Xt, Vn0, V10, V20, Vnt, V1t, V2t, thick, C_mat, false);
end

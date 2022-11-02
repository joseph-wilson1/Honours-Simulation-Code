function [Lt,L_index] = L_t(xvec,l_birth)
% Find cells whose length excees l_birth.
N = length(xvec)-1; %Number of cells.
Lt = 0;
L_index = [];
for i=1:N
    if li(xvec,i) >= l_birth
        Lt = Lt + 1;
        L_index = [L_index i];
    end
end
function l = li(xvec,i) 
%Length of cell is right cell border (x_{i+1} - x_{i}).
    l = abs(xvec(i+1) - xvec(i));
end
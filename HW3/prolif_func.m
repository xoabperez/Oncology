function [k] = prolif_function(k_in)
% k_in are 25 values that will be mapped to the correct location using this
% function
kt = reshape(k_in,[5 5]);
kt = imresize(kt,[50 50],'nearest');
k = ones(100,100);
k(25:74,25:74) = 2*kt;

k(k<0) = 0;




end
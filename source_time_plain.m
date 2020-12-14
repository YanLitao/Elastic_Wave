% function f = source_time_plain(z,t)
% 
% %%%% Source time function
% 
% %%%% Quasi-monochromatic wave train
% w = 3;
% f = sin(z + w*t);
% 
function f = source_time_plain(t,half_dur)

%%%% Source time function

%%%% Gaussian function

f = exp(-100*(t-half_dur).^2);


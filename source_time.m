function f = source_time(t,half_dur)

%%%% Source time function

%%%% Gaussian function
f = 2*exp(-200*(t-half_dur).^2);


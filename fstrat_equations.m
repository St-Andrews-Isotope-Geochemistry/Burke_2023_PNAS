function y = fstrat_equations(x,d34M, d33M,fb, d34b,d33b,  d34t, d33t,  l)

% This function allows solving for fstrat (x(1)) and d34strat(x(2)) 
% It is used by the script fstrat_solve.m
% Written by Andrea Burke
% Citation: Burke et al. (2023) "High sensitivity of summer temperatures to stratospheric sulfur
% loading from volcanoes in the Northern Hemisphere." Proceedings of the National Academy of Sciences (PNAS).


y(1) =  fb*d34b+ (1-fb-x(1))*d34t + x(1).*x(2) -d34M;
y(2) = fb*d33b + (1-fb-x(1))*d33t + x(1).*((x(2)/1000+1)^l-1)*1000- d33M; 
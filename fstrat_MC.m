%% This script sets Monte Carlo simulations of the solver for sytem of equations to determine fraction of stratospheric sulfate (fstrat) based on multiple sulfur isotopes
%% Written by Andrea Burke
%% When using cite Burke et al. (2023) "High sensitivity of summer temperatures to stratospheric sulfur
%% loading from volcanoes in the Northern Hemisphere." Proceedings of the National Academy of Sciences (PNAS).

function [save_solutions] = fstrat_MC(d34M, d33M,fb, d34b,d33b,  d34t, d33t,  l, stratmin_input, stratmax_input)

% How many times will the equations be solved?
iterations = length(d34M);

%create variable to store solutions 
save_solutions = [];

% increase number of steps in solver from default
options = optimoptions('fsolve','MaxFunctionEvaluations',1000, 'Display','none');

for i = 1:iterations

%% Call function (fstratsolver) that computes x, the solutions of the equations. In
%% this case x(1) is fstrat and x(2) is d34S_strat

% pass in parameters to the function
F = @(x)fstrat_equations(x,d34M(i), d33M(i),fb(i), d34b(i), d33b(i),  d34t(i), d33t(i),  l(i));

%create variable to store solutions
solutions = nan(2,2);

% set initial conditions
x0 = [0.5 15;0.5 -15];

%solve equations
[x,~,exitflag,~] = fsolve(F,x0(1,:),options);

 if exitflag ==1 % if solution was converged upon
     solutions(1,:) = x;
 end

%solve equations
[x,~,exitflag,~]  = fsolve(F,x0(2,:),options);

 if exitflag ==1 % if solution was converged upon
     solutions(2,:) = x;
 end

 fmax = 1 - fb(i); % this is the maximum fraction fstrat can be, given the background fraction

 if length(stratmin_input) == iterations
     stratmin = stratmin_input(i);
 elseif length(stratmin_input)==1
     stratmin = stratmin_input(1);
 else
     error('unexpected size of stratmin_input')
 end

 if length(stratmax_input) == iterations
     stratmax = stratmax_input(i);
 elseif length(stratmax_input)==1
     stratmax = stratmax_input(1);
 else
     error('unexpected size of stratmax_input')
 end


% Check solutions - were two found?
 if not(isnan(solutions(1,1)) & not(isnan(solutions(2,1)))) && solutions(1,1)~=solutions(2,1)    %if 2 unique solutions were found
    solutionind = find(solutions(:,1)>=0 & solutions(:,1)<=fmax & solutions(:,2) >= stratmin &solutions(:,2) <= stratmax ); %limit possible solutions to those that have fstrat within 0 and fmax
    if length(solutionind)==1 & isreal(solutions(solutionind,:))
        tobesaved = [solutions(solutionind,:) d34t(i)  l(i)];
    save_solutions = [save_solutions; tobesaved]; %save solutions that have fstrat between 0 and 1
    elseif length(solutionind)==2 & isreal(solutions(solutionind,:))
        tobesaved = [solutions(solutionind(1),:) d34t(i)  l(i);solutions(solutionind(2),:) d34t(i)  l(i)] ;
    save_solutions = [save_solutions; tobesaved] ; %save solutions that have fstrat between 0 and 1
    end
 elseif solutions(1,1)~=solutions(2,1) % if only one unique solution was found
    solutionind = find(solutions(:,1)>=0 & solutions(:,1)<=fmax & solutions(:,2) >= stratmin &solutions(:,2) <= stratmax); %limit possible solutions to those that have fstrat within 0 and fmax
    tobesaved = [solutions(solutionind,:) d34t(i)  l(i)];
    if solutionind & isreal(solutions(solutionind,:))
    save_solutions = [save_solutions; tobesaved]; %save solutions that have fstrat between 0 and fmax
    end
 elseif solutions(1,1)==solutions(2,1)  %if same solution was returned twice
    if solutions(1,1)>=0 && solutions(1,1)<=fmax && solutions(:,2) >= stratmin &&solutions(:,2) <= stratmax && isreal(solutions(1,:)) %if it's fstrat is between 0 and fmax

        tobesaved = [solutions(1,:) d34t(i)  l(i)];
    save_solutions = [save_solutions; tobesaved];%save solutions that have fstrat between 0 and fmax
    end 
   
 end

end

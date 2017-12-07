%HW7
%GB comments
1a 70. The equations generate fixed points at 1 and 0. Having 0 of anything won’t affect the system and is therefore a fixed point. Anything greater than 0 is a strong enough perturbation to the system that pushes  the population growth to its maximum capacity, which is a fixed point at 1 in this case. 
1b. 90. It does affect how quick the system reaches its fixed point values. 
1c. 90 I had to correct two mislabeled variables in your function to get it to work correctly. 
1d 70. No axis labels and no explaination of your results
2a 30 equation is not correct. This models activation and not a toggle switch. Use [V/(1+x(2)^4)-x(1); V/(1+x(1)^4)-x(2)];
2b. 100 correctly generated plots, but the equations used are wrong. I will give full credit
2c  100 same as 2b
overall 79

% Problem 1: Modeling population growth
% The simplest model for a growing population assumes that each current
% individual has equal likelihood to divide, which yields a differential
% equation dx/dt = a*x where a is the division rate. It is easy to see that
% this will grow exponentially without bound. A simple modification is to
% assume that the growth rate slows done as the population reaches some
% maximum value N so that dx/dt = a*x*(1-x/N). Defining X = x/N, we have 
% dX/dt = a*X*(1-X).

% Part 1. This equation has two fixed points at 0 and 1. Explain the
% meaning of these two points.

%Based on the wording of the question, a population of 0 or 1 doesn't make sense because there
%need to be 2+ individuals for sexual reproduction.

% Part 2: Evaluate the stability of these fixed points. Does it depend on
% the value of the parameter a?

%The stability of these two fixed points doesn't change. If x is
%substituted with 0, the equation yields an answer of 0 for dx/dt. If x is
%substituted with 1, the equation also yields an answer of 0 as dx/dt = a*
%(1)(1-1) equals 0.

% Part 3: Write a function that takes two inputs - the initial condition x0
% and the a parameter and integrates the equation forward in time. Make
% your code return two variables - the timecourse of X and the time
% required for the population to reach 99% of its maximum value. 

[time_course,time] = population_time(x0,a);

% Part 4: Another possible model is to consider discrete generations
% instead allowing the population to vary continuously. e.g. X(t+1) = a*
% X(t)*(1-X(t)). Consider this model and vary the a parameter in the range 0
% < a <= 4. For each value of a choose 200 random starting points  0 < x0 < 1 
% and iterate the equation forward to steady state. For each final
% value Xf, plot the point in the plane (a,Xf) so that at the end you will
% have produced a bifucation diagram showing all possible final values of
% Xf at each value of a. Explain your results. 

figure;
for a = 0.1:0.5:4
    for n = 1:200
    x0 = rand(); xf = x0;
    xf = a*xf*(1-xf);
    xf = a*xf*(1-xf);
    xf = a*xf*(1-xf);
    xf = a*xf*(1-xf);
    xf = a*xf*(1-xf);
    xf = a*xf*(1-xf);
    plot(a,xf,'.');
    hold on;
    end
end

% Problem 2. Genetic toggle switches. 
% Consider a genetic system of two genes A and B in which each gene
% product represses the expression of the other. Make the following
% assumptions: 
% a. Repression is cooperative:  each promotor region of one gene has 4
% binding sites for the other protein and that all of these need to be
% occupied for the gene to be repressed. 
% b. You can ignore the intermediate mRNA production so that the product of
% the synthesis of one gene can be assumed to directly repress the other
% c. the system is prefectly symmetric so that the degradation
% times, binding strengths etc are the same for both genes. 
% d. You can choose time and concentration scales so that all Michaelis
% binding constants and degradation times are equal to 1. 
%
% Part 1. Write down a two equation model (one for each gene product) for
% this system. Your model should have one free parameter corresponding to the
% maximum rate of expression of the gene, call it V. 
%
%Equation 1: dx/dt = (ku + kb(R/K)^n)/(1+R/K)^n) - x

%Equation 2: dA/dt = V*B^4/((1+B^4))-A, dB/dt = V*A^4/((1+A^4))-B


% Part 2. Write code to integrate your model in time and plot the results for V = 5 for two cases, 
% one in which A0 > B0 and one in which B0 > A0. 

V = 5;
rhs = @(t,x)[((V*x(2)^4)/(1+x(2)^4))-x(1);
    ((V*x(1)^4)/(1+x(1)^4))-x(2)];

%A0 > B0
sol = ode23(rhs, [0 10], [8,4]);
figure; plot(solx, soly(1,:));
hold on;
plot(solx, soly(2,:));

%B0 > A0
sol = ode23(rhs, [0 10],[2,4]);
figure; plot(solx, soly(1,:));
hold on;
plot(solx, soly(2,:));


% Part 3. By any means you want, write code to produce a bifurcation diagram showing all
% fixed points of the system as a function of the V parameter. 

gxfunc2 = @(x,V) ((V*x^4)/(1+x^4))-x;
figure; hold on;
for V = 1:0.1:5
    gxfunc = @(x) gxfunc2(x,V);
    for x0 = 1:0.1:5
        [rt,~,exitflag] = fzero(gxfunc, x0);
        if exitflag == 1
            plot(V, rt, 'k.');
            hold on;
        end
    end
end

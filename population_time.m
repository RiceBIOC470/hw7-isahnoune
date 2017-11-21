function [time_course,time] = population_time(x0,a)
N = 1000*0.99;
rhs = @(t,x) [a*x(1)*(1-x(1))];
x = x0; X = x/N;
sol = ode23(rhs, [x N],X);
time_course = (solx);
time = length(time_course);
   
end

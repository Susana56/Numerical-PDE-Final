% I organized the figures I produced in this main file

%% Figure 1 in the report: plot exact solution to display Boundary layer

ep = 1e-5;
h= 1/100;
x= (0:h:1)';
u = @(x) x- (exp(-(1-x)/ep) - exp(-1/ep))/(1 - exp (-1/ep));

plot(x, u(x))
title("Exact solution of the ODE with \epsilon = "+ ep,'FontSize', 14)
grid on, xlabel x



%% Central Difference scheme: Recreate plot on page 47, oscilations
clc
clear

% -eps*u'' + u' = 0 on (0,1), u(0) = 0, u(1)=1

b = @(x) 1;
c = @(x) 0;
f = @(x) 0;
alpha = 0;
beta = 1;

figure    % Figure 2
h = 0.025;
epsi = 1e-6;
[Y, x_bins, x] = Central_difference(h, epsi, alpha, beta, b, c,  f);

xf = 0:1/1000:1;
yexact = (exp((-(1-xf))/epsi)-exp(-1/epsi))/(1-exp(-1/epsi));

plot(x,abs(Y))
hold on
plot(xf,yexact)
title("2\epsilon << h with \epsilon = " + epsi + " and h = " + h, 'FontSize', 14)
grid on, xlabel x, legend ('finite diff','exact', 'FontSize', 14)

%% C. D. S: 2\epsilon ~ h
% -eps*u'' + u' = 0 on (0,1), u(0) = 0, u(1)=1
clc
clear

b = @(x) 1;
c = @(x) 0;
f = @(x) 0;
alpha = 0;
beta = 1;

figure    % Figure 3
h = 0.02;
epsi = 0.01;
[Y2, x_bins, x] = Central_difference(h, epsi, alpha, beta, b, c,  f);

xf = 0:1/1000:1;
yexact = (exp((-(1-xf))/epsi)-exp(-1/epsi))/(1-exp(-1/epsi));

plot(x,abs(Y2))
hold on
plot(xf,yexact)
title("2\epsilon = h with \epsilon = " + epsi + " and h = " + h, 'FontSize', 14)
grid on, xlabel x, legend ('finite diff','exact', 'FontSize', 14)



%% C. D. S: 2\epsilon > h
clc
clear

b = @(x) 1;
c = @(x) 0;
f = @(x) 0;
alpha = 0;
beta = 1;

figure    % Figure 4
h = 0.0002;
epsi = 0.001;
[Y, x_max, x] = Central_difference(h, epsi, alpha, beta, b, c,  f);

xf = 0:1/1000:1;
yexact = (exp((-(1-xf))/epsi)-exp(-1/epsi))/(1-exp(-1/epsi));

plot(x,abs(Y))
hold on
plot(xf,yexact)
title("2\epsilon > h with \epsilon = " + epsi + " and h = " + h, 'FontSize', 14)
grid on, xlabel x, legend ('finite diff','exact', 'FontSize', 14)


%% Upwind scheme: \epsilon > h

% -eps*u'' + u' = 0 on (0,1), u(0) = 0, u(1)=1

b = @(x) 1;
c = @(x) 0;
f = @(x) 0;
alpha = 0;
beta = 1;

figure    % Figure 5
h = 1/30;
epsi = 1e-1;
[Y, x_bins, x] = Simple_UPW(h, epsi, alpha, beta, b, c,  f);

xf = 0:1/1000:1;
yexact = (exp((-(1-xf))/epsi)-exp(-1/epsi))/(1-exp(-1/epsi));

plot(x,Y)
hold on
plot(xf,yexact)

title("Simple Upwind scheme : \epsilon = " + epsi + " and h = " + h, 'FontSize', 14)
grid on, xlabel x, legend ('finite diff','exact')

%% Upwind scheme: \epsilon < h 

b = @(x) 1;
c = @(x) 0;
f = @(x) 0;
alpha = 0;
beta = 1;

figure    % Figure 6
h = 1/30;
epsi = 1e-5;
[Y, x_bins, x] = Simple_UPW(h, epsi, alpha, beta, b, c,  f);

xf = 0:1/1000:1;
yexact = (exp((-(1-xf))/epsi)-exp(-1/epsi))/(1-exp(-1/epsi));

plot(x,Y)
hold on
plot(xf,yexact)

title("Simple Upwind scheme : \epsilon = " + epsi + " and h = " + h, 'FontSize', 14)
grid on, xlabel x, legend ('finite diff','exact')
%% Upwind scheme: failing example page, page 51

%%-eu''-u'=0, u(0)=0, u(1)=1
% has a boundary layer at x=0

b = @(x) -1;
c = @(x) 0;
f = @(x) 0;
alpha = 0;
beta = 1;

figure    % Figure 6
h = 1/30;
epsi = 1e-3;
[Y, x_bins, x] = Simple_UPW(h, epsi, alpha, beta, b, c,  f);

xf = 0:1/1000:1;
yexact = exp(-xf/epsi)/(exp(-1/epsi)-1) - (1/(exp(-1/eps)-1));

plot(x,Y)
hold on
plot(xf,yexact)

title("Simple Upwind scheme : \epsilon = " + epsi + " and h = " + h, 'FontSize', 14)
grid on, xlabel x, legend ('finite diff','exact')


%% Plot the typical error of the Upwind scheme

% Replicate fig 2. 3page 58
% for -eu''-u'=0 u(0)=u(1)=0

epsi = 1e-7 ;  
N_s = [10,20,40,60,80,100, 200,400,600,800,1000, 2000,4000,6000,8000,10000 , 20000,40000,60000, 80000, 100000, 200000, 400000, 600000, 800000, 1000000, 2000000, 4000000,6000000, 8000000, 10000000]' ;
error_y = zeros(length(N_s),1);


for i = 1:1:length(N_s)
    h= 1/N_s(i);
    r = epsi/ (epsi + h);
    u_1 = (1 - r)/(1-r^N_s(i));
    
    u_exact = (1-exp(-h/epsi)) / 1-exp(-1/epsi);
    
    error_y(i) = abs(u_1 - u_exact);
end

plot(N_s, error_y)
title("Error at the first grid point nearest boundary layer ",'FontSize', 14)
xlabel('N')
ylabel('u(x_1) - u_1')

%% Upwind scheme with Artificial diffusion
%%-eu''+u'=2x , u(0)=u(1)=1
% boundary layer at x=1
clc
clear

b = @(x) 1;
c = @(x) 0;
f = @(x) 2*x;
alpha = 1;
beta = 1;
sigma = @(x) sqrt(1+x^2);  

h = 1/30;
epsi = 1e-2;
[Y, x_bins, x] = UPW_AD(h, epsi, alpha, beta, b, c,  f, sigma);

xf = 0:1/1000:1;
yexact = (1+((2*epsi+1)/(exp(1/epsi)-1))) - ((2*epsi+1)/(exp(1/epsi)-1))*exp(xf/epsi) + (xf).^2 + 2*epsi*xf;

plot(x,Y)
hold on
plot(xf,yexact)

title("A.D. Upwind Scheme : \epsilon = " + epsi + " and h = " + h, 'FontSize', 14)
grid on, xlabel x, legend ('finite diff','exact')

%% Let's check convergence of U.W.S. A.D. at the grid point closest to the boundary layers
%%-eu''+u'=2x , u(0)=u(1)=1
% boundary layer at x=1
b = @(x) 1;
c = @(x) 0;
f = @(x) 2*x;
alpha = 1;
beta = 1;
sigma = @(x) sqrt(1+x^2); 
epsi = 1e-2 ;  

N_s = [10,20,40,60,80,100, 200,400,600,800,1000]' ;
error_y = zeros(length(N_s),1);

for i = 1:1:length(N_s)
    h= 1/N_s(i);
    [Y, x_bins, x] = UPW_AD(h, epsi, alpha, beta, b, c,  f, sigma);
    u_N_1 = Y(x_bins-1);
    u_exact = (1+((2*epsi+1)/(exp(1/epsi)-1))) - ((2*epsi+1)/(exp(1/epsi)-1))*exp((1-h)/epsi) + (1-h)^2 + 2*epsi*(1-h);
    
    
    error_y(i) = abs(u_N_1 - u_exact);
end

plot(N_s, error_y)
title("Error at the first grid point nearest boundary layer ",'FontSize', 14)
xlabel('N')
ylabel('u(x_{N-1}) - u_{N-1}')

%% Midpoint scheme & Gushchin-Shchennnikov scheme : trial with same function as ^ 

%%-eu''+u'=2x , u(0)=u(1)=1
% boundary layer at x=1
clc
clear

b = @(x) 1;
c = @(x) 0;
f = @(x) 2*x;
alpha = 1;
beta = 1;
sigma = @(x) sqrt(1+x^2);  

h = 1/30;
epsi = 1e-2;
[Y, x_bins, x] = MP_UPS(h, epsi, alpha, beta, b, c,  f);
[Y_2, x_2_bins, x_2] = GS_S(h, epsi, alpha, beta, b, c,  f);
xf = 0:1/1000:1;
yexact = (1+((2*epsi+1)/(exp(1/epsi)-1))) - ((2*epsi+1)/(exp(1/epsi)-1))*exp(xf/epsi) + (xf).^2 + 2*epsi*xf;

plot(x,Y)
hold on
plot(x_2, Y_2)
plot(xf,yexact)

title("Extra 2 schemes : \epsilon = " + epsi + " and h = " + h, 'FontSize', 14)
grid on, xlabel x, legend ('Midpoint Upwind','Gushchin-Shchennikov Scheme','exact')

%% Trial: apply all the schemes to one ODE problem 

% -eu''+u' = x for u(0)=u(1)=0

clc
clear

b = @(x) 1;
c = @(x) 0;
f = @(x) x;
alpha = 0;
beta = 0;
sigma = @(x) sqrt(1+x^2);  
h = 0.002;
epsi = 0.001;

% Central difference 
[Y_1, x_1_max, x_1] = Central_difference(h, epsi, alpha, beta, b, c,  f);

%Simple Upwind
[Y_2, x_2_bins, x_2] = Simple_UPW(h, epsi, alpha, beta, b, c,  f);

% Upwind with Artificial Diffusion
[Y_3, x_3_bins, x_3] = UPW_AD(h, epsi, alpha, beta, b, c,  f, sigma);

%Midpoint Upwind scheme
[Y_4, x_4_bins, x_4] = MP_UPS(h, epsi, alpha, beta, b, c,  f);

%Gushchin-Shchennikov scheme
[Y_5, x_5_bins, x_5] = GS_S(h, epsi, alpha, beta, b, c,  f);

xf = 0:1/1000:1;
yexact = ((epsi+(1/2))/(exp(1/epsi)-1)) - ((epsi+(1/2))/(exp(1/epsi)-1))*exp(xf/epsi) + ((xf).^2)/2 + epsi*xf;

plot(x_1,Y_1,'.-')
hold on
plot(x_2,Y_2)
plot(x_3,Y_3,'--')
plot(x_4,Y_4, '+ -')
plot(x_5,Y_5)
plot(xf,yexact,'b')
title("All Scheme : \epsilon = " + epsi + " and h = " + h, 'FontSize', 14)
grid on, xlabel x, legend ('C. D. S','Simple UPW','UPW AD','MP UPW','GS','exact')

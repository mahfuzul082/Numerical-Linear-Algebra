clear all; clc;

%inputs
n = 100;
A = diag(ones(n,1)) + diag(ones(n-1,1),1) + diag(ones(n-1,1),-1);
b = zeros(n, 1);
b(50) = 1;

%Gauss elimination with partial pivoting
x = gauss_pp(A, b);
%built in
x_ = A \ b;

%plot
f = figure(1);
plot(1:n, x, "bo-", LineWidth=1.5, MarkerSize=8, DisplayName="Gaussian elimination w/ partial pivoting"); hold on;
plot(1:n, x_, "rx", LineWidth=1.5, MarkerSize=9, DisplayName="Matlab built-in"); hold off;
set(gca,"TicklabelInterpreter","latex","FontSize", 15);
legend("boxoff");
legend("Location", "northwest", "interpreter", "latex");
xlabel("$j$", "interpreter", "latex");
ylabel("$x_j$", "interpreter", "latex");
pbaspect([3 1 1]);
exportgraphics(f,"lin_solv.pdf");
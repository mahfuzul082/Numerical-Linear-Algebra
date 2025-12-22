clear all; clc;

% inputs
n = 200; % size of A
A = spdiags([-ones(n,1) 2*ones(n,1) -ones(n,1)], -1:1, n, n); % A matrix
b = ones(n,1); % b (rhs)
tol = 1e-3; % tolerance of residual
x_0 = zeros(n,1); % initial guess
max_iter = 1e5; % maximum no. of iterations

% exact solution
x_exact = A \ b;

fprintf("===== Jacobi =====\n");
[x_j, num_iter_j, time_j, res_j] = my_jacobi(A, b, x_0, tol, max_iter);
fprintf("Solve time: %.8f s, No. of iterations: %d\n", time_j, num_iter_j);
err_j = norm(x_j - x_exact);
fprintf("Final error norm = %.8e\n\n", err_j);

fprintf("===== Gauss-Seidel =====\n");
[x_gs, num_iter_gs, time_gs, res_gs] = my_gauss_seidel(A, b, x_0, tol, max_iter);
fprintf("Solve time: %.8f s, No. of iterations: %d\n", time_gs, num_iter_gs);
err_gs = norm(x_gs - x_exact);
fprintf("Final error norm = %.8e\n\n", err_gs);

fprintf("===== GMRES =====\n");
[x_gm, num_iter_gm, time_gm, res_gm] = my_gmres(A, b, tol);
fprintf("Solve time: %.8f s, No. of iterations: %d\n", time_gm, num_iter_gm);
err_gm = norm(x_gm - x_exact);
fprintf("Final error norm = %.8e\n\n", err_gm);

fprintf("===== Conjugate Gradient =====\n");
[x_cg, num_iter_cg, time_cg, res_cg] = my_cg(A, b, tol);
fprintf("Solve time: %.8f s, No. of iterations: %d\n", time_cg, num_iter_cg);
err_cg = norm(x_cg - x_exact);
fprintf("Final error norm = %.8e\n", err_cg);

fig1 = figure(1);
loglog(res_j(:,1), res_j(:,2), "ro-", "LineWidth", 1.5, "DisplayName", "Jacobi"); hold on;
loglog(res_gs(:,1), res_gs(:,2), "bo-", "LineWidth", 1.5, "DisplayName", "Gauss-Seidel"); hold on;
loglog(res_gm(:,1), res_gm(:,2), "o-", "LineWidth", 1.5, "DisplayName", "GMRES"); hold on;
loglog(res_cg(:,1), res_cg(:,2), "o-", "LineWidth", 1.5, "DisplayName", "Conjugate Gradient"); hold off;
set(gca,"TicklabelInterpreter","latex","FontSize", 15);
legend("boxoff");
legend("Location", "northeast", "interpreter", "latex");
ylim([1e-3 1e3]);
xlabel("\textit{No. of iterations}", "interpreter", "latex");
ylabel("$\|b-Ax\|_2$", "interpreter", "latex");
pbaspect([2 1 1]);
exportgraphics(fig1,"lin_solve.pdf");
% MATLAB script to solve a saddle point problem similar to the C++ code

% Parameters
n = 5; % Number of grid points
h = 1.0 / n; % Grid spacing
c_val = 1.0; % Average over volume constraint value

% Assemble the Poisson matrix A
A = spdiags([-1/h*ones(n,1), 2/h*ones(n,1), -1/h*ones(n,1)], -1:1, n, n);

% Assemble the constraint matrix C
C = h * ones(n, 1);

% Assemble the right-hand side vector F
F = h * ones(n, 1);

% Assemble the constraint vector c
c = c_val;

% Assemble the saddle point matrix G
G = [A, C; C', 0];

% Assemble the right-hand side vector b
b = [F; c];

% Solve the linear system G * x = b
x = G \ b;

% Compute the Schur complement S = -C' * inv(A) * C
invA = inv(A);
S = -C' * invA * C;

% Display the results
disp('Poisson Matrix A:');
disp(full(A));

disp('Saddle Point Matrix G:');
disp(full(G));

disp('Combined Vector b:');
disp(full(b));

disp('Solution x:');
disp(full(x));

% Extract the solution components
x_F = x(1:n);
x_c = x(n+1);

disp('Solution component x_F:');
disp(full(x_F));

disp('Solution component x_c:');
disp(full(x_c));

% Display the Schur complement
disp('Schur Complement S:');
disp(full(S));
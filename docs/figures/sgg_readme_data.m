% Dr. Denys Dutykh (Mathematics Department, Khalifa University of Science and Technology, Abu Dhabi, UAE)
%
% Data of the figure of the README: the run of matlab/run_graph.m (three kinks of speed 0.95 on the
% graph of four vertices and six edges), repeated step by step without its graphics. The parameters,
% the initial condition, the vertex rule and the leapfrog scheme are those of run_graph.m.
%
% Run it from this folder with GNU Octave or MATLAB:
%   octave --no-gui sgg_readme_data.m        or, in MATLAB,   run sgg_readme_data
% It writes sgg_frames.mat beside itself (the solution every 0.5 time units, the embedding of the
% graph in the plane and the discrete energy at every step); plot_readme_figure.py then draws
% sgg-kinks.png from it.

here = fileparts(mfilename('fullpath'));

% Parameters of run_graph.m
N  = 500;                     % nodes on each edge
l  = 10.0;                    % length of each edge
x  = linspace(0, l, N);       % coordinate along an edge
dx = x(2) - x(1);
dt = 0.5*l/N;                 % time step
Nt = 3300;                    % number of time steps (final time 33)
t  = 0.0;
dtdx2 = (dt/dx)^2;
e = ones(N, 1);
L = dtdx2*spdiags([e -2*e e], -1:1, N, N);

% Graph: six edges, each oriented from its first vertex to its second
M  = 6;
Nv = 4;
E  = [1 2; 2 3; 2 4; 3 4; 4 1; 1 3];

% Nearest lattice points of the three edges meeting at each vertex: [node, edge]
Vn = zeros(Nv, 3, 2);
Vn(1, :, :) = [2 1; 2 6; N-1 5];
Vn(2, :, :) = [2 2; 2 3; N-1 1];
Vn(3, :, :) = [2 4; N-1 2; N-1 6];
Vn(4, :, :) = [2 5; N-1 3; N-1 4];

% Initial condition: kinks of speed c0 on edges 1, 5 and 6, the constant 2*pi on edges 2, 3 and 4
c0 = 0.95;
u0 = zeros(N, M);
u1 = zeros(N, M);
u2 = zeros(N, M);
kink0 = 4.0*atan(exp((x - x(N/2))/sqrt(1 - c0*c0)));
kink1 = 4.0*atan(exp((x - x(N/2) - c0*dt)/sqrt(1 - c0*c0)));
u0(:, 1) = kink0;  u1(:, 1) = kink1;
u0(:, 2:4) = 2*pi; u1(:, 2:4) = 2*pi;
u0(end:-1:1, 5) = kink0;  u1(end:-1:1, 5) = kink1;
u0(:, 6) = kink0;  u1(:, 6) = kink1;

% Embedding of the graph in the plane, as in run_graph.m (for the drawing only)
th = linspace(0.0, 1.0, N)';
X = zeros(N, 2, M);
X(:, 1, 1) = 0.0;                                X(:, 2, 1) = -l*(1 - th);
X(:, 1, 2) = l*cos(pi/6)*th;                     X(:, 2, 2) = l*sin(pi/6)*th;
X(:, 1, 3) = l*cos(5*pi/6)*th;                   X(:, 2, 3) = l*sin(5*pi/6)*th;
X(:, 1, 4) = l*cos(pi/6*(1-th) + 5*pi/6*th);     X(:, 2, 4) = l*sin(pi/6*(1-th) + 5*pi/6*th);
X(:, 1, 5) = l*cos(5*pi/6*(1-th) + 3*pi/2*th);   X(:, 2, 5) = l*sin(5*pi/6*(1-th) + 3*pi/2*th);
X(:, 1, 6) = l*cos(3*pi/2*(1-th) + 13*pi/6*th);  X(:, 2, 6) = l*sin(3*pi/2*(1-th) + 13*pi/6*th);
V = [0.0, 0.0, l*cos(pi/6), l*cos(5*pi/6); -l, 0.0, l*sin(pi/6), l*sin(5*pi/6)];

% Discrete energy between two consecutive levels a and b (trapezoidal weights on each edge)
w = dx*ones(N, 1);  w([1 N]) = dx/2;
potential = @(u) 0.5*sum(sum(diff(u).^2))/dx + sum(w'*(1 - cos(u)));
energy = @(a, b) 0.5*sum(w'*((b - a)/dt).^2) + 0.5*(potential(a) + potential(b));

keep = 50;                                % store the solution every 50 steps (0.5 time units)
frames = zeros(N, M, Nt/keep + 1);
tkeep = zeros(Nt/keep + 1, 1);
frames(:, :, 1) = u0;
energies = zeros(Nt, 1);                  % energy between the levels t - dt and t

for i = 1:Nt
  t = t + dt;
  for j = 1:M
    % vertex rule of run_graph.m: the value at each end of an edge is the mean of the values at
    % the nearest nodes of the three edges that meet there
    u0(1, j) = mean(u0(sub2ind([N M], Vn(E(j,1), :, 1), Vn(E(j,1), :, 2))));
    u1(1, j) = mean(u1(sub2ind([N M], Vn(E(j,1), :, 1), Vn(E(j,1), :, 2))));
    u0(N, j) = mean(u0(sub2ind([N M], Vn(E(j,2), :, 1), Vn(E(j,2), :, 2))));
    u1(N, j) = mean(u1(sub2ind([N M], Vn(E(j,2), :, 1), Vn(E(j,2), :, 2))));
    u2(:, j) = 2*u1(:, j) - u0(:, j) + L*u1(:, j) - dt*dt*sin(u1(:, j));
  end
  energies(i) = energy(u0, u1);
  u0 = u1; u1 = u2;
  if mod(i, keep) == 0
    frames(:, :, i/keep + 1) = u0;
    tkeep(i/keep + 1) = t;
  end
end

save(fullfile(here, 'sgg_frames.mat'), 'frames', 'tkeep', 'X', 'V', 'E', 'energies', ...
     'dt', 'dx', 'c0', 'u0', 'u1', '-v7');
gamma = 1/sqrt(1 - c0^2);
fprintf('three kinks: 3*8*gamma = %.10f\n', 24*gamma);
fprintf('discrete energy: first %.10f, last %.10f, min %.10f, max %.10f\n', ...
        energies(1), energies(end), min(energies), max(energies));
fprintf('final time %.12f, solution between %.6f and %.6f\n', t, min(u0(:)), max(u0(:)));

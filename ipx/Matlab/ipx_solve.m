%IPX_SOLVE interior point linear programming solver
%
% [result,stats]=ipx_solve(model,params);
%
% solves the linear programming problem
%
%   min c'x subject to Ax{<,>,=}rhs, lb<=x<=ub.
%
% model must be a struct with fields:
%   A       m-by-n sparse matrix
%   rhs     size m double vector, entries must be finite
%   obj     size n double vector, entries must be finite
%   lb      size n double vector, entries can be -inf
%   ub      size n double vector, entries can be +inf
%   sense   size m char vector, entries in {'<','>','='}
%
% params is an optional argument to provide solver parameters. If given, it
% must be a struct; see ipx_defaults() for recognized parameter names.
%
% result is an optional output argument that returns the solution. It is a
% struct that contains at least the fields
%   status, status_ipm, status_crossover, errflag.
% If status_ipm is 1 or 2, then the solution from the interior point method
% is returned as fields
%   x, xl, xu, slack, y, zl, zu.
% If status_crossover is 1 or 2, then the basic solution from the crossover
% method is returned as fields
%  basis_x, basis_slack, basis_y, basis_z, cbasis, vbasis.
% See the reference guide for a definition of these quantities.
%
% stats is an optional output argument that returns solver statistics.

error('ipx_solve mexFunction not found');

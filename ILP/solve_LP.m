function [x_opt] = solve_LP(f, a, b, a_eq, b_eq, lb, ub, linprog_options)
%SOLVE_LP Solves a generic LP
%
%	Version: 1.11
%	Date: 1/06/20
%	Author: Thomas Thayer (tthayer@ucmerced.edu)
%
%	This function solves a generic Linear Program (LP), such that the sum of f * x_opt is minimized and the inequality constraints are upper bounds.
%   Inputs:
%		f: the function being minimized
%		a: the matrix of inequality functions
%		b: the matrix of inequality constraints
%		a_eq: the matrix of equality functions
%		b_eq: the matrix of equality constraints
%		lb: a vector of lower bounds for x_opt
%		ub: a vector of upper bounds for x_opt
%		linprog_options: pass through options for Matlab's linprog(), which if exists, will force usage of Matlab's linear program solver
%	Outputs:
%		x_opt: the vector of optimal x values
%TODO: Not tested with mac or windows, not tested with gurobi

	%% Check for other solvers
	no_other_solver = 0;
	[cplex_flag, cplex_path] = system('command -v cplex');
	[gurobi_flag, gurobi_path] = system('command -v gurobi');
	[scip_flag, scip_path] = system('command -v scip');
	scip_flag = scip_flag | ~exist('SaveMPS') | ~exist('BuildMPS') | ~exist('read_cplexsol');
	if exist('linprog_options') && ~isempty(linprog_options)
		no_other_solver = 1;
	else
		linprog_options = optimoptions(@linprog, 'Preprocess', 'none', 'Display', 'off', 'MaxIter', 10*(size(b,2) + size(b_eq,2) + size(a,2)));
	end
	no_other_solver = no_other_solver || any([~cplex_flag, ~gurobi_flag, ~scip_flag]);
	%% Solve LP with best found solver
	warning('off','all');
	if cplex_flag == 0 && (no_other_solver == 0)
		%addpath('/opt/ibm/ILOG/CPLEX_Studio1271/cplex/matlab/x86-64_linux');
		cplex_path = split(cplex_path, 'bin');
		cplex_path = cplex_path{1};
		addpath([cplex_path(1:end-1), '/matlab/x86-64_linux']);
		x_opt = cplexlp(f, a, b, a_eq, b_eq, lb, ub);
	elseif gurobi_flag == 0 && (no_other_solver == 0)
		gurobi_path = split(gurobi_path, 'bin');
		gurobi_path = gurobi_path{1};
		addpath([gurobi_path(1:end-1), '/examples/matlab']);
		[x_opt, fval, exit_flag, output, lambda] = linprog(f, x_int, a, b, a_eq, b_eq, lb, ub);
		if any(exit_flag == [-2:0])
			error('Linprog was not successful. Check inputs.');
		end
	elseif scip_flag == 0 && (no_other_solver == 0)
		SaveMPS("lp.mps", BuildMPS(a, b, a_eq, b_eq, f, lb, ub));
		system('scip -c "read $PWD/lp.mps optimize write solution $PWD/lp.sol quit" > $PWD/lp.log');
		x_opt = import_SOL_LP('lp.sol');
		x_opt = full(sparse(x_opt(:, 1), ones(size(x_opt, 1), 1), x_opt(:, 2), length(f), 1));
		delete('lp.mps', 'lp.log', 'lp.sol');
	else
		%linprog_options = optimoptions(@linprog, 'Algorithm', interior-point, 'Preprocess', 'none', 'Display', 'off', 'MaxIter', 10*(size(b,2) + size(b_eq,2) + size(a,2)));
		[x_opt, fval, exit_flag, output] = linprog(f, a, b, a_eq, b_eq, lb, ub, linprog_options);
		if any(exit_flag == [-2:0])
			error('Linprog was not successful. Check inputs.');
		end
	end
	warning('on','all');
	%% old
% 	%% Check for other solvers
% 	no_other_solver = 0;
% 	if exist('SaveMPS') && exist('BuildMPS')
% 		[scip_flag, output] = system('command -v scip');
% 		if exist('read_cplexsol')
% 			[cplex_flag, output] = system('command -v cplex');
% 		else
% 			cplex_flag = 1;
% 		end
% 		if ~(scip_flag == 0) && ~(cplex_flag == 0)
% 			no_other_solver = 1;
% 		end
% 	else
% 		no_other_solver = 1;
% 	end
% 	if exist('linprog_options') && ~isempty(linprog_options)
% 		no_other_solver = 1;
% 	else
% 		linprog_options = optimoptions(@linprog, 'Preprocess', 'none', 'Display', 'off', 'MaxIter', 10*(size(b,2) + size(b_eq,2) + size(a,2)));
% 	end
% 	%% Solve LP with best found solver
% 	%no_other_solver = 1; % force matlab solver
% 	%cplex_flag = 1; % dont use cplex solver
% 	if no_other_solver
% 		%linprog_options = optimoptions(@linprog, 'Algorithm', interior-point, 'Preprocess', 'none', 'Display', 'off', 'MaxIter', 10*(size(b,2) + size(b_eq,2) + size(a,2)));
% 		[x_opt, fval, exit_flag, output] = linprog(f, a, b, a_eq, b_eq, lb, ub, linprog_options);
% 		if any(exit_flag == [-2:0])
% 			error('Linprog was not successful. Check inputs.');
% 		end
% 	else
% 		warning('off','all');
% 		if cplex_flag == 0
% 			%system('cplex -c "read $PWD/lp.mps" "optimize" "write $PWD/lp.sol" > lp.log');
% 			%x_opt = read_cplexsol('lp.sol');
% 			addpath('/opt/ibm/ILOG/CPLEX_Studio1271/cplex/matlab/x86-64_linux');
% 			x_opt = cplexlp(f, a, b, a_eq, b_eq, lb, ub);
% 			delete('clone1.log', 'clone2.log', 'cplex.log');
% 		elseif scip_flag == 0
% 			SaveMPS("lp.mps", BuildMPS(a, b, a_eq, b_eq, f, lb, ub));
% 			system('scip -c "read $PWD/lp.mps optimize write solution $PWD/lp.sol quit" > $PWD/lp.log');
% 			x_opt = import_SOL_LP('lp.sol');
% 			x_opt = full(sparse(x_opt(:, 1), ones(size(x_opt, 1), 1), x_opt(:, 2), length(f), 1));
% 		end
% 		delete('lp.mps', 'lp.log', 'lp.sol');
% 		warning('on','all');
% 	end
end


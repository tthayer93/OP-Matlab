function [x_opt] = solve_ILP(f, x_int, a, b, a_eq, b_eq, lb, ub, intlinprog_options)
%SOLVE_ILP Solves a generic ILP
%
%	Version: 2.1
%	Date: 04/2/21
%	Author: Thomas Thayer (tthayer@ucmerced.edu)
%
%	This function solves a generic Integer Linear Program (ILP), such that the sum of f * x_opt is minimized and the inequality constraints are upper bounds.
%	Inputs:
%		f: the function being minimized
%		x_int: the x variables that are contrained as integers
%		a: the matrix of inequality functions
%		b: the vector of inequality constraints
%		a_eq: the matrix of equality functions
%		b_eq: the vector of equality constraints
%		lb: a vector of lower bounds for x_opt
%		ub: a vector of upper bounds for x_opt
%		intlinprog_options: pass through options for Matlab's linprog(), which if exists, will force usage of Matlab's linear program solver
%	Outputs:
%		x_opt: the vector of optimal x values
%
%TODO: Not tested with mac or windows, not tested with gurobi

	%% Check for other solvers
	no_other_solver = 0;
	[cplex_flag, cplex_path] = system('command -v cplex');
	[gurobi_flag, gurobi_path] = system('command -v gurobi');
	[scip_flag, scip_path] = system('command -v scip');
	scip_flag = scip_flag | ~exist('SaveMPS') | ~exist('BuildMPS') | ~exist('read_cplexsol');
	if exist('intlinprog_options') && ~isempty(intlinprog_options)
		no_other_solver = 1;
	else
		intlinprog_options = optimoptions(@intlinprog, 'Display', 'off', 'Heuristics', 'advanced');
	end
	no_other_solver = no_other_solver || any([~cplex_flag, ~gurobi_flag, ~scip_flag]);
	%% Solve LP with best found solver
	warning('off','all');
	if cplex_flag == 0 && (no_other_solver == 0)
		%addpath('/opt/ibm/ILOG/CPLEX_Studio1271/cplex/matlab/x86-64_linux');
		cplex_path = split(cplex_path, 'bin');
		cplex_path = cplex_path{1};
		addpath([cplex_path(1:end-1), '/matlab/x86-64_linux']);
		ctype = zeros(1, length(f));
		ctype(:) = 'C';
		ctype(x_int) = 'I';
		x_opt = cplexmilp(f, a, b, a_eq, b_eq, [], [], [], lb, ub, char(ctype));
	elseif gurobi_flag == 0 && (no_other_solver == 0)
		gurobi_path = split(gurobi_path, 'bin');
		gurobi_path = gurobi_path{1};
		addpath([gurobi_path(1:end-1), '/examples/matlab']);
		[x_opt, fval, exit_flag, output] = intlinprog(f, x_int, a, b, a_eq, b_eq, lb, ub);
		if any(exit_flag == [-2:0])
			error('Linprog was not successful. Check inputs.');
		end
	elseif scip_flag == 0 && (no_other_solver == 0)
		SaveMPS("lp.mps", BuildMPS(a, b, a_eq, b_eq, f, lb, ub, 'Integer', x_int));
		system('scip -c "read $PWD/lp.mps optimize write solution $PWD/lp.sol quit" > $PWD/lp.log');
		x_opt = import_SOL_LP('lp.sol');
		x_opt = full(sparse(x_opt(:, 1), ones(size(x_opt, 1), 1), x_opt(:, 2), length(f), 1));
		delete('lp.mps', 'lp.log', 'lp.sol');
	else
		%intlinprog_options = optimoptions(@intlinprog, 'Display', 'off', 'Heuristics', 'advanced');
		[x_opt, fval, exit_flag, output] = intlinprog(f, x_int, a, b, a_eq, b_eq, lb, ub, [], intlinprog_options);
		if any(exit_flag == [-2:0])
			error('Linprog was not successful. Check inputs.');
		end
	end
	warning('on','all');
	
	
	%% old
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
% 	if exist('linprog_options') && ~isempty(intlinprog_options)
% 		no_other_solver = 1;
% 	else
% 		intlinprog_options = optimoptions(@intlinprog, 'Display', 'off', 'Heuristics', 'advanced');
% 	end
% 	%% Solve LP with best found solver
% 	%no_other_solver = 1; % force matlab solver
% 	%cplex_flag = 1; % dont use cplex solver
% 	if no_other_solver
% 		%intlinprog_options = optimoptions(@intlinprog, 'Display', 'off', 'Heuristics', 'advanced');
% 		[x_opt, fval, exit_flag, output] = intlinprog(f, x_int, a, b, a_eq, b_eq, lb, ub, [], intlinprog_options);
% 		if any(exit_flag == [-2:0])
% 			error('Linprog was not successful. Check inputs.');
% 		end
% 	else
% 		warning('off','all');
% 		if cplex_flag == 0
% 			%system('cplex -c "read $PWD/lp.mps" "optimize" "write $PWD/lp.sol" > lp.log');
% 			%x_opt = read_cplexsol('lp.sol');
% 			addpath('/opt/ibm/ILOG/CPLEX_Studio1271/cplex/matlab/x86-64_linux');
% 			ctype = zeros(1, length(f));
% 			ctype(:) = 'C';
% 			ctype(x_int) = 'I';
% 			x_opt = cplexmilp(f, a, b, a_eq, b_eq, [], [], [], lb, ub, char(ctype));
% 			delete('clone1.log', 'clone2.log', 'cplex.log');
% 		elseif scip_flag == 0
% 			SaveMPS("lp.mps", BuildMPS(a, b, a_eq, b_eq, f, lb, ub, 'Integer', x_int));
% 			system('scip -c "read $PWD/lp.mps optimize write solution $PWD/lp.sol quit" > $PWD/lp.log');
% 			x_opt = import_SOL_LP('lp.sol');
% 			x_opt = full(sparse(x_opt(:, 1), ones(size(x_opt, 1), 1), x_opt(:, 2), length(f), 1));
% 		end
% 		delete('lp.mps', 'lp.log', 'lp.sol');
% 		warning('on','all');
% 	end

end
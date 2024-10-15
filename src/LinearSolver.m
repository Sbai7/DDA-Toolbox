classdef LinearSolver < handle
% LinearSolver  Helper class for sparse iterative linear solvers and 
% preconditioning methods. 
%
% Author: M.A. Sbai, Ph.D.
%
% Copyright (C) 2024 Mohammed Adil Sbai
%
% This program is free software: you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by the
% Free Software Foundation, either version 3 of the License, or (at your
% option) any later version.
%
% This program is distributed in the hope that it will be useful, but 
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
% Public License for more details.
%
% You should have received a copy of the GNU General Public License along
% with this program. If not, see <http://www.gnu.org/licenses/>.


    properties (Access = protected)
        rhs                     % last rhs of this linear solver  
        runs                    % number of calls to solve method 
        solver                  % name of currently used solver
        prec                    % name of currently used preconditioner
    end
    
    properties (Access = public)
        %--- iterative linear solver input parameters 
        tol                     % tolerance for iterative solver convergence
        maxit                   % max. number of linear iterations
        icholtype               % incomplete Cholesky factorization type
        michol                  % modified incomplete Cholesky flag 
        ilutype                 % incomplete LU factorization type 
        droptol                 % drop tolerance for ichol/ilu preconditioner
        gmres_restart           % GMRES's restart parameter 
        
        %--- iterative linear solver output parameters 
        conv_flag               % flag indicating convergence behavior
        rrn                     % relative residual norm 
        iter                    % number of performed iterations 
        res_vec                 % residual vector 
    end
    
    methods 
        %--- class setters ------------------------------------------------                
        function set.ilutype(this,type)
            validateattributes(type, {'char'}, {'nonempty'});
            
            supported_types = {'nofill','crout','ilutp'};
            tf = ismember(type, supported_types);
            if tf 
                this.ilutype = type;
            else 
                error('ilutype field must be one of the following: %s', ... 
                    ['(''' supported_types{1} ''' ''' supported_types{2} ... 
                     ''' ''' supported_types{3} ''')'] );
            end
        end
        
        function set.icholtype(this,type)
            validateattributes(type, {'char'}, {'nonempty'});
            
            supported_types = {'nofill','ict'};
            tf = ismember(type, supported_types);
            if tf 
                this.icholtype = type;
            else 
                error('icholtype field must be one of the following: %s', ... 
                    ['(''' supported_types{1} ''' ''' supported_types{2} ...
                     ''')'] );
            end
        end
        
        function set.michol(this,type) 
            validateattributes(type, {'char'}, {'nonempty'});
            
            supported_types = {'on','off'};
            tf = ismember(type, supported_types);
            if tf 
                this.michol = type;
            else 
                error('michol field must be one of the following: %s', ... 
                    ['(''' supported_types{1} ''' ''' supported_types{2} ...
                     ''')'] );
            end
        end 
    end
    
    methods (Access = public)
        function this = LinearSolver()
            % Class constructor 
            this.runs          = 0;            
            this.icholtype     = 'nofill';
            this.michol        = 'off';
            this.ilutype       = 'nofill';
            this.tol           = 1e-15;
            this.maxit         = 100;
            this.droptol       = 1e-6;
            this.gmres_restart = 3;
        end
        
        %---
        function solution = solve(this,matrix,rhs,varargin)
            % Solves the sparse linear system of equations Ax=b where A is the
            % sparse matrix of this instance of LinearSolver class, and b is
            % its right hand side vector. The solver type (direct or iterative)
            % could be explicitly specified. 
            %
            % varargin is a set of 'key', 'value' pairs.
            %
            % Supported options are:
            % 'Solver'      key to specify the iterative solver. Accepted
            % values for this key are 'mldivide', 'pcg', 'bicgstab',
            % 'bicgstabl', 'tfqmr', and 'gmres'. The direct solver (i.e. 
            % 'mldivide') is selected by default when this option is omitted.
            % 
            % Solution is a vector of the same shape and size of rhs
            %
%            validateattributes(matrix, {'sparse'}, {'square'});
            validateattributes(rhs, {'double'}, {'column'});
            assert(size(matrix,1) == size(rhs,1));
            
            % set default options 
            persistent p
            if isempty(p)
                p = inputParser;
                defaultType = 'mldivide';
                supportedSolvers = {defaultType,'pcg','bicgstab', ...
                    'bicgstabl','tfqmr','gmres'};
                p.addParameter('Solver',defaultType, ... 
                               @(x) any(validatestring(x,supportedSolvers)));
            end
            
            parse(p,varargin{:});
            type = p.Results.Solver;
            
            % Update the rhs member of the solver class
            this.rhs = rhs;
           
            if strcmp(type,'mldivide') 
                % solve linear system using LU factorization with partial 
                % pivoting 
                solution = matrix\rhs;
                
                this.solver = 'mldivide';
                this.runs   = this.runs + 1;
                
            elseif strcmp(type,'pcg')
                % modified incomplete Cholesky preconditioner
                L = ichol(matrix,struct('type',this.icholtype, ... 
                                        'michol',this.michol,  ...
                                        'droptol',this.droptol)); 
                if strcmp(this.michol,'off')
                    this.prec = 'incomplete cholesky';
                else
                    this.prec = 'modified incomplete cholesky';
                end
                
                % solve linear system using preconditionned conjugate
                % gradient method 
                [solution,this.conv_flag,this.rrn,this.iter,this.res_vec] = ... 
                    pcg(matrix,rhs,this.tol,this.maxit,L,L');
                
                this.solver = 'pcg';
                this.runs   = this.runs + 1;
                
            elseif strcmp(type,'bicgstab')
                % Create a preconditioner with ilu, since A is nonsymmetric
                [L,U] = ilu(matrix,struct('type',this.ilutype, ...
                                          'droptol',this.droptol));
                this.prec = 'incomplete factorization';
                
                % solve the system of linear equations by Biconjugate 
                % gradients stabilized method 
                [solution,this.conv_flag,this.rrn,this.iter,this.res_vec] = ...
                    bicgstab(matrix,rhs,this.tol,this.maxit,L,U);
                
                this.solver = 'bicgstab';
                this.runs   = this.runs + 1;
                
            elseif strcmp(type,'bicgstabl')
                % Create a preconditioner with ilu, since A is nonsymmetric
                [L,U] = ilu(matrix,struct('type',this.ilutype, ...
                                          'droptol',this.droptol));
                this.prec = 'incomplete factorization';
                
                % solve the system of linear equations by Biconjugate 
                % gradients stabilized l method 
                [solution,this.conv_flag,this.rrn,this.iter,this.res_vec] = ...
                    bicgstabl(matrix,rhs,this.tol,this.maxit,L,U);
                
                this.solver = 'bicgstabl';
                this.runs   = this.runs + 1;
                
            elseif strcmp(type,'tfqmr')
                % Create a preconditioner with ilu, since A is nonsymmetric
                [L,U] = ilu(matrix,struct('type',this.ilutype, ...
                                          'droptol',this.droptol));
                this.prec = 'incomplete factorization';
                
                % solve the system of linear equations by Transpose free 
                % quasi- minimum residual method 
                [solution,this.conv_flag,this.rrn,this.iter,this.res_vec] = ...
                    tfqmr(matrix,rhs,this.tol,this.maxit,L,U);
                
                this.solver = 'tfqmr';
                this.runs   = this.runs + 1;
                
            elseif strcmp(type,'gmres')
                % Create a preconditioner with ilu, since A is nonsymmetric
                [L,U] = ilu(matrix,struct('type',this.ilutype, ...
                                          'droptol',this.droptol));
                this.prec = 'incomplete factorization';
                
                % solve the system of linear equations by Generalized Minimal
                % Residual method 
                [solution,this.conv_flag,this.rrn,this.iter,this.res_vec] = ...
                    gmres(matrix,rhs,this.gmres_restart,this.tol,this.maxit,L,U);
                
                this.solver = 'gmres';
                this.runs   = this.runs + 1;
                
            else 
                error('Invalid Solver specified.');
            end
        end 

        
        %------------------------------------------------------------------ 
        
        function y = converged(this)
        % Returns true if the linear solver has converged to the desired 
        % tolerance within maximum allowed iterations, and false otherwise.
        % 
        
            assert(isa(this,'LinearSolver')); 
            y = (this.conv_flag==0);
            
        end
        
        %------------------------------------------------------------------
        
        function message = convergenceMessage(this)
        % Returns a message depending on the convergence status as a string  
        % to avoid dealing with specific valuse of this object conv_flag's 
        % property. 
        % 
        
            assert(isa(this,'LinearSolver')); 
            
            if     this.conv_flag==0 
                message = ['iterative solution by ' this.solver ... 
                           ' and ' this.prec ' preconditioning' ... 
                           ' methods converged in ' num2str(this.iter) ... 
                           ' iterations'];
                
            elseif this.conv_flag==1
                message = 'no convergence: increase maximum number of iterations.';
                
            elseif this.conv_flag==2
                message = 'no convergence: ill-conditioned preconditioner.';
                
            elseif this.conv_flag==3
                message = 'no convergence: iteration stagnation.';
                
            elseif this.conv_flag==4
                message = 'no convergence: underflow or overflow occured.';
                
            end
            
        end
        
        function plot(this,varargin)
        % Plots the convergence of iterative solvers supported by this
        % class such as pcg, bicgstab, etc. 
        % 
        % extra arguments in varargin must be the same as those for varagin
        % in plot command. They're used for changing chart line properties
        % such as Line and Markers properties (LineStyle, LineWidth, Color,
        % LineJoin, AlignVertexCenters, Marker, MarkerSize,
        % MarkerEdgeColor, and MarkerFaceColor).
        %
            
            assert(isa(this,'LinearSolver')); 
            
            if this.runs==0 || strcmp(this.solver,'mldivide')
                return, 
            end
            
            if strcmp(this.solver,'bicgstab') || ... 
                   strcmp(this.solver,'tfqmr') 
                % bicgstab and tfqmr algorithms use half iterations to 
                % update the residual vector
                semilogy(0:0.5:this.iter,this.res_vec/norm(this.rhs), ... 
                         varargin{:});
                     
            elseif strcmp(this.solver,'bicgstabl') 
                % bicgstabl algorithm uses quarter iterations to update the
                % residual vector
                semilogy(0:0.25:this.iter,this.res_vec/norm(this.rhs), ... 
                         varargin{:});
                     
            elseif strcmp(this.solver,'gmres') 
                
                % gmres algorithm uses 1/n iterations to update the
                % residual vector where n is the gmres restart parameter
                if this.gmres_restart==0 
                    semilogy([0,1:0.5:this.iter(2)], ... 
                        this.res_vec/norm(this.rhs),varargin{:});
                         
                else
                    last_iter = (length(this.res_vec)+this.gmres_restart-2) ... 
                        /this.gmres_restart;
                    semilogy([0,1:1/(this.gmres_restart):last_iter], ... 
                        this.res_vec/norm(this.rhs),varargin{:});
                    
                end
           
                
            else
                semilogy(0:this.iter,this.res_vec/norm(this.rhs), ... 
                         varargin{:});
                
            end
            xlabel('iteration number');
            ylabel('relative residual');            
        end
        
        %------------------------------------------------------------------ 
      
    end
end

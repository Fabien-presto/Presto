% This file is part of Presto-HF, a matlab toolbox to identify a circuit from
% its response.
%
% SPDX-License-Identifier: AGPL-3.0-or-later
%
% Copyright 2025 by
%   Centre Inria de l'Université Côte d'Azur
%   2004, route des Lucioles
%   06902 Sophia Antipolis Cedex
%
% and by
%   Mines Paris - Armines
%   60, boulevard Saint-Michel
%   75006 Paris
%
% Contributors: Fabien Seyfert, Jean-Paul Marmorat, Martine Olivi
%
% Presto-HF is free software: you can redistribute it and/or modify it under
% the terms of the GNU Affero General Public License as published by the Free
% Software Foundation, either version 3 of the License, or (at your option) any
% later version.
%
% Presto-HF is distributed in the hope that it will be useful, but WITHOUT ANY
% WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
% A PARTICULAR PURPOSE. See the GNU Affero General Public License for more
% details.
%
% You should have received a copy of the GNU Affero General Public License
% along with Presto-HF. If not, see <https://www.gnu.org/licenses/>.
%
function [X,FVAL,EXITFLAG,OUTPUT,LAMBDA,GRAD,HESSIAN] = ...
                            arl2Min(FUN,X,A,B,Aeq,Beq,LB,UB,NONLCON,options,varargin)
defaultopt = optimset('display','final','LargeScale','on', ...
   'TolX',1e-6,'TolFun',1e-6,'TolCon',1e-6,'DerivativeCheck','off',...
   'Diagnostics','off',...
   'GradObj','off','GradConstr','off','Hessian','off','MaxFunEvals','100*numberOfVariables',...
   'DiffMaxChange',1e-1,'DiffMinChange',1e-8,...
   'PrecondBandWidth',0,'TypicalX','ones(numberOfVariables,1)','MaxPCGIter','max(1,floor(numberOfVariables/2))', ...
   'TolPCG',0.1,'MaxIter',400,'HessPattern',[]);
% If just 'defaults' passed in, return the default options in X
if nargin==1 & nargout <= 1 & isequal(FUN,'defaults')
   X = defaultopt;
   return
end

large = 'large-scale';
medium = 'medium-scale';

if nargin < 4, error('arl2Min requires at least four input arguments'); end
if nargin < 10, options=[];
   if nargin < 9, NONLCON=[];
      if nargin < 8, UB = [];
         if nargin < 7, LB = [];
            if nargin < 6, Beq=[];
               if nargin < 5, Aeq =[];
               end, end, end, end, end, end
if isempty(NONLCON) & isempty(A) & isempty(Aeq) & isempty(UB) & isempty(LB)
   error('arl2Min is for constrained problems. Use FMINUNC for unconstrained problems.')
end

if nargout > 4
   computeLambda = 1;
else 
   computeLambda = 0;
end

caller='constr';
lenVarIn = length(varargin);
XOUT=X(:);
numberOfVariables=length(XOUT);

options = optimset(defaultopt,options);
switch optimget(options,'display')
case {'off','none'}
   verbosity = 0;
case 'iter'
   verbosity = 2;
case 'final'
   verbosity = 1;
otherwise
   verbosity = 1;
end

% Set to column vectors
B = B(:);
Beq = Beq(:);



[XOUT,l,u,msg] = checkbounds(XOUT,LB,UB,numberOfVariables);
if ~isempty(msg)
   EXITFLAG = -1;
   [FVAL,OUTPUT,LAMBDA,GRAD,HESSIAN] = deal([]);
   X(:)=XOUT;
   if verbosity > 0
      disp(msg)
   end
   return
end
lFinite = l(~isinf(l));
uFinite = u(~isinf(u));


meritFunctionType = 0;



diagnostics = isequal(optimget(options,'diagnostics','off'),'on');
gradflag =  strcmp(optimget(options,'GradObj'),'on');
hessflag = strcmp(optimget(options,'Hessian'),'on');
if isempty(NONLCON)
   constflag = 0;
else
   constflag = 1;
end
gradconstflag =  strcmp(optimget(options,'GradConstr'),'on');
line_search = strcmp(optimget(options,'largescale','off'),'off'); % 0 means trust-region, 1 means line-search



% Convert to inline function as needed
if ~isempty(FUN)  % will detect empty string, empty matrix, empty cell array
   [funfcn, msg] = fprefcnchk(FUN,'arl2Min',length(varargin),gradflag,hessflag);
else
   errmsg = sprintf('%s\n%s', ...
      'FUN must be a function name, valid string expression, or inline object;', ...
      ' or, FUN may be a cell array that contains these type of objects.');
   error(errmsg)
end

if constflag % NONLCON is non-empty
   [confcn, msg] = fprefcnchk(NONLCON,'arl2Min',length(varargin),gradconstflag,[],1);
else
   confcn{1} = '';
end





[rowAeq,colAeq]=size(Aeq);
% if only l and u then call sfminbx
if ~line_search & isempty(NONLCON) & isempty(A) & isempty(Aeq) & gradflag
   OUTPUT.algorithm = large;
   % if only Aeq beq and Aeq has as many columns as rows, then call sfminle
elseif ~line_search & isempty(NONLCON) & isempty(A) & isempty(lFinite) & isempty(uFinite) & gradflag ...
      & colAeq >= rowAeq
   OUTPUT.algorithm = large;
elseif ~line_search
   warning(['Trust region method does not currently solve this type of problem,',...
         sprintf('\n'), 'switching to line search.'])
   if isequal(funfcn{1},'fungradhess')
      funfcn{1}='fungrad';
      warning('Hessian provided by user will be ignored in line search algorithm')
      
   elseif  isequal(funfcn{1},'fun_then_grad_then_hess')
      funfcn{1}='fun_then_grad';
      warning('Hessian provided by user will be ignored in line search algorithm')
   end    
   hessflag = 0;
   OUTPUT.algorithm = medium;
elseif line_search
   OUTPUT.algorithm = medium;
   if issparse(Aeq) | issparse(A)
      warning('can not do sparse with line_search, converting to full')
   end
   
   % else call nlconst
else
   error('Unrecognized combination of OPTIONS flags and calling sequence.')
end


lenvlb=length(l);
lenvub=length(u);




if isequal(OUTPUT.algorithm,medium)
   CHG = 1e-7*abs(XOUT)+1e-7*ones(numberOfVariables,1);
   i=1:lenvlb;
   lindex = XOUT(i)<l(i);
   if any(lindex),
      XOUT(lindex)=l(lindex)+1e-4; 
   end
   i=1:lenvub;
   uindex = XOUT(i)>u(i);
   if any(uindex)
      XOUT(uindex)=u(uindex);
      CHG(uindex)=-CHG(uindex);
   end
   X(:) = XOUT;
else
   arg = (u >= 1e10); arg2 = (l <= -1e10);
   u(arg) = inf*ones(length(arg(arg>0)),1);
   l(arg2) = -inf*ones(length(arg2(arg2>0)),1);
   if min(min(u-XOUT),min(XOUT-l)) < 0, 
      XOUT = startx(u,l);
      X(:) = XOUT;
   end
end

% Evaluate function
GRAD=zeros(numberOfVariables,1);
HESS = [];



switch funfcn{1}
case 'fun'
   try
      f = feval( funfcn{3}, X,varargin{:});
   catch
      errmsg = sprintf('%s\n%s\n\n%s',...
         'arl2Min cannot continue because user supplied objective function', ...
         ' failed with the following error:', lasterr);
      error(errmsg);
   end
case 'fungrad'
   try
      [f,GRAD(:)] = arl2EvalPsi(  X,varargin{:});
   catch
      errmsg = sprintf('%s\n%s\n\n%s',...
         'arl2Min cannot continue because user supplied objective function', ...
         ' failed with the following error:', lasterr);
      error(errmsg);
   end
case 'fungradhess'
   try
      [f,GRAD(:),HESS] = arl2EvalPsi(  X,varargin{:});
   catch
      errmsg = sprintf('%s\n%s\n\n%s',...
         'arl2Min cannot continue because user supplied objective function', ...
         ' failed with the following error:', lasterr);
      error(errmsg);
   end
otherwise
   error('Undefined calltype in arl2Min');
end



% Evaluate constraints
switch confcn{1}
case 'fun'
   try 
      [ctmp,ceqtmp] = arl2SPModulus(X,varargin{:});
      c = ctmp(:); ceq = ceqtmp(:);
      cGRAD = zeros(numberOfVariables,length(c));
      ceqGRAD = zeros(numberOfVariables,length(ceq));
   catch
      if findstr(xlate('Too many output arguments'),lasterr)
         if isa(confcn{3},'inline')
            errmsg = sprintf('%s%s%s\n%s\n%s\n%s', ...
               'The inline function ',formula(confcn{3}),' representing the constraints',...
               ' must return two outputs: the nonlinear inequality constraints and', ...
               ' the nonlinear equality constraints.  At this time, inline objects may',...
               ' only return one output argument: use an M-file function instead.');
         else
            errmsg = sprintf('%s%s%s\n%s%s', ...
               'The constraint function ',confcn{3},' must return two outputs:',...
               ' the nonlinear inequality constraints and', ...
               ' the nonlinear equality constraints.');
         end
         error(errmsg)
      else
         errmsg = sprintf('%s\n%s\n\n%s',...
            'arl2Min cannot continue because user supplied nonlinear constraint function', ...
            ' failed with the following error:', lasterr);
         error(errmsg);
      end
   end
   
case 'fungrad'
   try
      [ctmp,ceqtmp,cGRAD,ceqGRAD] = arl2SPModulus(X,varargin{:});
      c = ctmp(:); ceq = ceqtmp(:);
   catch
      errmsg = sprintf('%s\n%s\n\n%s',...
         'arl2Min cannot continue because user supplied nonlinear constraint function', ...
         ' failed with the following error:', lasterr);
      error(errmsg);
   end
case ''
   c=[]; ceq =[];
   cGRAD = zeros(numberOfVariables,length(c));
   ceqGRAD = zeros(numberOfVariables,length(ceq));
otherwise
   error('Undefined calltype in arl2Min');
end



non_eq = length(ceq);
non_ineq = length(c);
[lin_eq,Aeqcol] = size(Aeq);
[lin_ineq,Acol] = size(A);
[cgrow, cgcol]= size(cGRAD);
[ceqgrow, ceqgcol]= size(ceqGRAD);

eq = non_eq + lin_eq;
ineq = non_ineq + lin_ineq;

if ~isempty(Aeq) & Aeqcol ~= numberOfVariables
   error('Aeq has the wrong number of columns.')
end
if ~isempty(A) & Acol ~= numberOfVariables
   error('A has the wrong number of columns.')
end
if  cgrow~=numberOfVariables & cgcol~=non_ineq
   error('Gradient of the nonlinear inequality constraints is the wrong size.')
end
if ceqgrow~=numberOfVariables & ceqgcol~=non_eq
   error('Gradient of the nonlinear equality constraints is the wrong size.')
end

if diagnostics > 0
   % Do diagnostics on information so far
   msg = diagnose('arl2Min',OUTPUT,gradflag,hessflag,constflag,gradconstflag,...
      line_search,options,XOUT,non_eq,...
      non_ineq,lin_eq,lin_ineq,l,u,funfcn,confcn,f,GRAD,HESS,c,ceq,cGRAD,ceqGRAD);
end

% call algorithm

if isequal(OUTPUT.algorithm,medium)
   [X,FVAL,lambda,EXITFLAG,OUTPUT,GRAD,HESSIAN]=...
      nlconst(funfcn,X,l,u,full(A),B,full(Aeq),Beq,confcn,options, ...
      verbosity,gradflag,gradconstflag,hessflag,meritFunctionType,...
      CHG,f,GRAD,HESS,c,ceq,cGRAD,ceqGRAD,varargin{:});
   LAMBDA=lambda;
else
   disp(OUTPUT.algorithm);
   error('Bad algorithm');
end



%--------------------------------------------------------------------        fprefcnchk
function [allfcns,msg] = fprefcnchk(funstr,caller,lenVarIn,gradflag,hessflag,constrflag)
%PREFCNCHK Pre- and post-process function expression for FUNCHK.
%   [ALLFCNS,MSG] = PREFUNCHK(FUNSTR,CALLER,lenVarIn,GRADFLAG) takes
%   the (nonempty) expression FUNSTR from CALLER with LenVarIn extra arguments,
%   parses it according to what CALLER is, then returns a string or inline
%   object in ALLFCNS.  If an error occurs, this message is put in MSG.
%
%   ALLFCNS is a cell array: 
%    ALLFCNS{1} contains a flag 
%    that says if the objective and gradients are together in one function 
%    (calltype=='fungrad') or in two functions (calltype='fun_then_grad')
%    or there is no gradient (calltype=='fun'), etc.
%    ALLFCNS{2} contains the string CALLER.
%    ALLFCNS{3}  contains the objective (or constraint) function
%    ALLFCNS{4}  contains the gradient function
%    ALLFCNS{5}  contains the hessian function (not used for constraint function).
%  
%    NOTE: we assume FUNSTR is nonempty.
% Initialize
if nargin < 6
   constrflag = 0;
end
if constrflag
   graderrmsg = 'Constraint gradient function expected (OPTIONS.GradConstr==''on'') but not found.';
   warnstr = ...
      sprintf('%s\n%s\n%s\n','Constraint gradient function provided but OPTIONS.GradConstr==''off'';', ...
      '  ignoring constraint gradient function and using finite-differencing.', ...
      '  Rerun with OPTIONS.GradConstr==''on'' to use constraint gradient function.');
else
   graderrmsg = 'Gradient function expected OPTIONS.GradObj==''on'' but not found.';
   warnstr = ...
      sprintf('%s\n%s\n%s\n','Gradient function provided but OPTIONS.GradObj==''off'';', ...
      '  ignoring gradient function and using finite-differencing.', ...
      '  Rerun with OPTIONS.GradObj==''on'' to use gradient function.');
   
end
msg='';
allfcns = {};
funfcn = [];
gradfcn = [];
hessfcn = [];
if gradflag & hessflag 
   calltype = 'fungradhess';
elseif gradflag
   calltype = 'fungrad';
else
   calltype = 'fun';
end

% {fun}
if isa(funstr, 'cell') & length(funstr)==1
   % take the cellarray apart: we know it is nonempty
   if gradflag
      error(graderrmsg)
   end
   [funfcn, msg] = fcnchk(funstr{1},lenVarIn);
   if ~isempty(msg)
      if constrflag
         msg = ['NONLCON must be a function name.'];
      end
      
      error(msg);
   end
   
   % {fun,[]}      
elseif isa(funstr, 'cell') & length(funstr)==2 & isempty(funstr{2})
   if gradflag
      error(graderrmsg)
   end
   [funfcn, msg] = fcnchk(funstr{1},lenVarIn);
   if ~isempty(msg)
      if constrflag
         msg = ['NONLCON must be a function name.'];
      end
      
      error(msg);
   end  
   
   % {fun, grad}   
elseif isa(funstr, 'cell') & length(funstr)==2 % and ~isempty(funstr{2})
   
   [funfcn, msg] = fcnchk(funstr{1},lenVarIn);
   if ~isempty(msg)
      if constrflag
         msg = ['NONLCON must be a function name.'];
      end
      error(msg);
   end  
   [gradfcn, msg] = fcnchk(funstr{2},lenVarIn);
   if ~isempty(msg)
      if constrflag
         msg = ['NONLCON must be a function name.'];
      end
      
      error(msg);
   end
   calltype = 'fun_then_grad';
   if ~gradflag
      warning(warnstr);
      calltype = 'fun';
   end
   
   
   % {fun, [], []}   
elseif isa(funstr, 'cell') & length(funstr)==3 ...
      & ~isempty(funstr{1}) & isempty(funstr{2}) & isempty(funstr{3})
   if gradflag
      error(graderrmsg)
   end
   if hessflag
      error('Hessian function expected but not found.')
   end
   
   [funfcn, msg] = fcnchk(funstr{1},lenVarIn);
   if ~isempty(msg)
      if constrflag
         msg = ['NONLCON must be a function name.'];
      end
      
      error(msg);
   end  
   
   % {fun, grad, hess}   
elseif isa(funstr, 'cell') & length(funstr)==3 ...
      & ~isempty(funstr{2}) & ~isempty(funstr{3})
   [funfcn, msg] = fcnchk(funstr{1},lenVarIn);
   if ~isempty(msg)
      if constrflag
         msg = ['NONLCON must be a function name.'];
      end
      
      error(msg);
   end  
   [gradfcn, msg] = fcnchk(funstr{2},lenVarIn);
   if ~isempty(msg)
      if constrflag
         msg = ['NONLCON must be a function name.'];
      end
      
      error(msg);
   end
   [hessfcn, msg] = fcnchk(funstr{3},lenVarIn);
   if ~isempty(msg)
      if constrflag
         msg = ['NONLCON must be a function name.'];
      end
      
      error(msg);
   end
   calltype = 'fun_then_grad_then_hess';
   if ~hessflag & ~gradflag
      hwarnstr = sprintf('%s\n%s\n%s\n','Hessian and gradient functions provided ', ...
         '  but OPTIONS.Hessian==''off'' and OPTIONS.GradObj==''off''; ignoring Hessian and gradient functions.', ...
         '  Rerun with OPTIONS.Hessian==''on'' and OPTIONS.GradObj==''on'' to use derivative functions.');
      warning(hwarnstr)
      calltype = 'fun';
   elseif hessflag & ~gradflag
      warnstr = ...
         sprintf('%s\n%s\n%s\n','Hessian and gradient functions provided ', ...
         '  but OPTIONS.GradObj==''off''; ignoring Hessian and gradient functions.', ...
         '  Rerun with OPTIONS.Hessian==''on'' and OPTIONS.GradObj==''on'' to use derivative functions.');
      warning(warnstr)
      calltype = 'fun';
   elseif ~hessflag & gradflag
      hwarnstr = ...
         sprintf('%s\n%s\n%s\n','Hessian function provided but OPTIONS.Hessian==''off'';', ...
         '  ignoring Hessian function,', ...
         '  Rerun with OPTIONS.Hessian==''on'' to use Hessian function.');
      warning(hwarnstr);
      calltype = 'fun_then_grad';
   end
   
   % {fun, grad, []}   
elseif isa(funstr, 'cell') & length(funstr)==3 ...
      & ~isempty(funstr{2}) & isempty(funstr{3})
   if hessflag
      error('Hessian function expected but not found.')
   end
   [funfcn, msg] = fcnchk(funstr{1},lenVarIn);
   if ~isempty(msg)
      if constrflag
         msg = ['NONLCON must be a function name.'];
      end
      
      error(msg);
   end  
   [gradfcn, msg] = fcnchk(funstr{2},lenVarIn);
   if ~isempty(msg)
      if constrflag
         msg = ['NONLCON must be a function name.'];
      end
      
      error(msg);
   end
   calltype = 'fun_then_grad';
   if ~gradflag
      warning(warnstr);
      calltype = 'fun';
   end
   
   % {fun, [], hess}   
elseif isa(funstr, 'cell') & length(funstr)==3 ...
      & isempty(funstr{2}) & ~isempty(funstr{3})
   error('Hessian function given without gradient function.')
   
elseif ~isa(funstr, 'cell')  %Not a cell; is a string expression, function name string or inline object
   [funfcn, msg] = fcnchk(funstr,lenVarIn);
   if ~isempty(msg)
      if constrflag
         msg = ['NONLCON must be a function name.'];
      end
      
      error(msg);
   end   
   if gradflag % gradient and function in one function/M-file
      gradfcn = funfcn; % Do this so graderr will print the correct name
   end  
   if hessflag & ~gradflag
      hwarnstr = ...
         sprintf('%s\n%s\n%s\n','OPTIONS.Hessian==''on'' ', ...
         '  but OPTIONS.GradObj==''off''; ignoring Hessian and gradient functions.', ...
         '  Rerun with OPTIONS.Hessian==''on'' and OPTIONS.GradObj==''on'' to use derivative functions.');
      warning(hwarnstr)
   end
   
else
   errmsg = sprintf('%s\n%s', ...
      'FUN must be a function name or inline object;', ...
      ' or, FUN may be a cell array that contains these type of objects.');
   error(errmsg)
end

allfcns{1} = calltype;
allfcns{2} = caller;
allfcns{3} = funfcn;
allfcns{4} = gradfcn;
allfcns{5} = hessfcn;


%-----------------------------------------------------------------      checkbounds
function [x,lb,ub,msg] = checkbounds(xin,lbin,ubin,nvars)
%CHECKBOUNDS Move the initial point within the (valid) bounds.
%   [X,LB,UB,X,FLAG] = CHECKBOUNDS(X0,LB,UB,nvars) 
%   checks that the upper and lower
%   bounds are valid (LB <= UB) and the same length as X (pad with -inf/inf
%   if necessary); warn if too long.  Also make LB and UB vectors if not 
%   already.
%   Finally, inf in LB or -inf in UB throws an error.

%   Copyright (c) 1990-98 by The MathWorks, Inc.
%   $Revision: 1.6 $  $Date: 2003/12/02 11:49:54 $
%   Mary Ann Branch 5-1-98

msg = [];
% Turn into column vectors
lb = lbin(:); 
ub = ubin(:); 
xin = xin(:);

lenlb = length(lb);
lenub = length(ub);
lenx = length(xin);

% Check maximum length
if lenlb > nvars
   warning('Length of lower bounds is > length(x); ignoring extra bounds');
   lb = lb(1:nvars);   
   lenlb = nvars;
elseif lenlb < nvars
   lb = [lb; -inf*ones(nvars-lenlb,1)];
   lenlb = nvars;
end

if lenub > nvars
   warning('Length of upper bounds is > length(x); ignoring extra bounds');
   ub = ub(1:nvars);
   lenub = nvars;
elseif lenub < nvars
   ub = [ub; inf*ones(nvars-lenub,1)];
   lenub = nvars;
end

% Check feasibility of bounds
len = min(lenlb,lenub);
if any( lb( (1:len)' ) > ub( (1:len)' ) )
   count = full(sum(lb>ub));
   if count == 1
      msg=sprintf(['\nExiting due to infeasibility:  %i lower bound exceeds the' ...
            ' corresponding upper bound.\n'],count);
   else
      msg=sprintf(['\nExiting due to infeasibility:  %i lower bounds exceed the' ...
            ' corresponding upper bounds.\n'],count);
   end 
end
% check if -inf in ub or inf in lb   
if any(eq(ub, -inf)) 
   error('-Inf detected in upper bound: upper bounds must be > -Inf.');
elseif any(eq(lb,inf))
   error('+Inf detected in lower bound: lower bounds must be < Inf.');
end

x = xin;

%--------------------------------------------------------------    compdir
function [SD, dirType] = compdir(Z,H,gf,nvars,f);
% COMPDIR Computes a search direction in a subspace defined by Z. 
%    Helper function for NLCONST.
%    Returns Newton direction if possible.
%    Returns random direction if gradient is small.
%    Otherwise, returns steepest descent direction.  
%    If the steepest descent direction is small it computes a negative
%    curvature direction based on the most negative eigenvalue.
%    For singular matrices, returns steepest descent even if small.

%   Copyright (c) 1990-98 by The MathWorks, Inc.
%   $Revision: 1.6 $  $Date: 2003/12/02 11:49:54 $
%   Mary Ann Branch 10-20-96.

% Define constant strings
Newton = 'Newton';
Random = 'random';
SteepDescent = 'steepest descent';
Eigenvector = 'eigenvector';

%  SD=-Z*((Z'*H*Z)\(Z'*gf));
dirType = [];
% Compute the projected Newton direction if possible
 projH = Z'*H*Z;
 [R, p] = chol(projH);
 if ~p  % positive definite: use Newton direction
   SD = - Z*(R \ ( R'\(Z'*gf)));
   dirType = Newton;
 else % not positive definite
   % If the gradient is small, try a random direction:
      % Sometimes the search direction goes to zero in negative
      % definite problems when the current point rests on
      % the top of the quadratic function. In this case we can move in
      % any direction to get an improvement in the function so 
      % foil search direction by giving a random gradient.
   if norm(gf) < sqrt(eps)
      SD = -Z*Z'*(rand(nvars,1) - 0.5);
      dirType = Random;
   else
     % steepest descent
     stpDesc = - Z*(Z'*gf);
     % check if ||SD|| is close to zero 
     if norm(stpDesc) > sqrt(eps)       
        SD = stpDesc;
        dirType = SteepDescent;
     else
        % Look for a negative curvature direction
        %  Some attempt at efficiency: usually it's
        %  faster to use EIG unless many variables.
        if nvars < 400  
           [VV,DD] = eig(projH);
           [smallRealEig, eigind] = min(diag(DD));
           ev = VV(:,eigind(1));
        else
           options.disp = 0;
           [ev, smallRealEig, flag] = eigs(projH,1,'sr',options);
           if flag  % Call to eigs failed
              [VV,DD] = eig(projH);
              [smallRealEig, eigind] = min(diag(DD));
              ev = VV(:,eigind(1));
           end
        end
        
        if smallRealEig < 0
          % check the sign of SD and the magnitude.
          SDtol = 100*eps*norm(gf); % Note: we know norm(gf) > sqrt(eps)
          Zev = Z*ev;
          if Zev'*gf > SDtol
            SD = -Zev;
            dirType = Eigenvector;
          elseif Zev'*gf < SDtol
            SD = Zev;
            dirType = Eigenvector;
          else % 
            SD = stpDesc;
            dirType = SteepDescent; 
          end
        else % The projected Hessian is singular,i.e., zero direction is ok
          %  -- will propagate thru the algorithm.
          SD = stpDesc;
          dirType = SteepDescent; 
        end % smallRealEig < 0
      end % randSD'*(gf) < -SDtol
   end %  norm(stpDesc) > sqrt(eps)  
 end % ~p

   %--------------------------------------------------------------    diagnose
function msg = diagnose(caller,OUTPUT,gradflag,hessflag,constflag,gradconstflag,line_search,OPTIONS,XOUT,non_eq,...
   non_ineq,lin_eq,lin_ineq,LB,UB,funfcn,confcn,f,GRAD,HESS,c,ceq,cGRAD,ceqGRAD);
% DIAGNOSE prints diagnostic information about the function to be minimized
%    or solved.

%   Copyright (c) 1990-98 by The MathWorks, Inc.
%   $Revision: 1.6 $ $Date: 2003/12/02 11:49:54 $

msg = [];

pstr = sprintf('\n%s\n%s\n',...
   '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%',...
   '   Diagnostic Information ');
disp(pstr)

if ~isempty(funfcn{1})
   if isa(funfcn{3},'inline')
      funformula = formula(funfcn{3});
   else 
      funformula = funfcn{3};
   end
   if isa(funfcn{4},'inline')
      gradformula = formula(funfcn{4});
   else 
      gradformula = funfcn{4};
   end
   
   if isa(funfcn{5},'inline')
      hessformula = formula(funfcn{5});
   else 
      hessformula = funfcn{5};
   end
end

if ~isempty(confcn{1})
   if isa(confcn{3},'inline')
      conformula = formula(confcn{3});
   else 
      conformula = confcn{3};
   end
   if isa(confcn{4},'inline')
      gradcformula = formula(confcn{4});
   else 
      gradcformula = confcn{4};
   end
else
   conformula = '';
   gradcformula = '';
end

disp(['Number of variables: ', int2str(length(XOUT)),sprintf('\n')])
if ~isempty(funfcn{1})
   disp('Functions ')
   switch funfcn{1}
   case 'fun'
      % display 
      disp([' Objective:                            ',funformula]);
      
   case 'fungrad'
      if gradflag
         disp([' Objective and gradient:               ',funformula]);
      else
         disp([' Objective:                            ',funformula]);
         disp( '   (set OPTIONS.GradObj=''on'' to use user provided gradient function)') 
      end
      
   case 'fungradhess'
      if gradflag & hessflag
         disp([' Objective, gradient and Hessian:      ',funformula]);
      elseif gradflag
         disp([' Objective and gradient:               ',funformula]);
         disp( '   (set OPTIONS.Hessian to ''on'' to use user provided Hessian function)') 
      else
         disp([' Objective:                            ',funformula]);
         disp( '   (set OPTIONS.GradObj=''on'' to use user provided gradient function)')
         disp( '   (set OPTIONS.Hessian to ''on'' to use user provided Hessian function)') 
      end
      
      
   case 'fun_then_grad'
      disp([' Objective:                            ',funformula]);
      if gradflag
         disp([' Gradient:                             ',gradformula]);
      end
      if hessflag
         disp('-->Ignoring OPTIONS.Hessian --no user Hessian function provided')
      end
      
   case 'fun_then_grad_then_hess'
      disp([' Objective:                            ',funformula]);
      if gradflag & hessflag
         disp([' Gradient:                             ',gradformula]);
         disp([' Hessian:                              ',hessformula]);
      elseif gradflag
         disp([' Gradient:                             ',gradformula]);
      end   
   otherwise
      
   end
   
   if ~gradflag
      disp(' Gradient:                             finite-differencing')
   end
   % shape of grad
   
   if ~hessflag & (isequal('fmincon',caller) | isequal('fminunc',caller))
      disp(' Hessian:                              finite-differencing (or Quasi-Newton)')
   end
   % shape of hess
   
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~isempty(confcn{1})
   switch confcn{1}
      
   case 'fun'
      disp([' Nonlinear constraints:                ',conformula]);
   case 'fungrad'
      if gradconstflag
         disp([' Nonlinear constraints and gradient:   ',conformula]);
      else
         disp([' Nonlinear constraints:                ',conformula]);
         disp( '   (set OPTIONS.GradConstr to ''on'' to use user provided gradient of constraints function)') 
      end
      
   case 'fun_then_grad'
      disp([' Nonlinear constraints:                ',conformula]);
      if gradconstflag
         disp([' Nonlinear constraints gradient:       ',gradcformula]);
      end
      
   otherwise
      
   end
   
   if ~constflag
      disp(' Nonlinear constraints:                finite-differencing')
   end
   if ~gradconstflag
      
      disp(' Gradient of nonlinear constraints:    finite-differencing')
   end
   disp([sprintf('\n'),'Constraints'])  
   disp([' Number of nonlinear inequality constraints: ',int2str(non_ineq)])
   disp([' Number of nonlinear equality constraints:   ',int2str(non_eq)])
   
elseif isequal(caller,'fmincon') | isequal(caller,'fminimax') | ...
      isequal(caller,'fgoalattain') | isequal(caller,'fseminf')
   disp([sprintf('\n'),'Constraints'])
   disp(' Nonlinear constraints:             do not exist')
   
end

disp(' ')


switch caller
   
case {'fmincon','linprog','quadprog','lsqlin','fminimax','fseminf','fgoalattain'}
   disp([' Number of linear inequality constraints:    ',int2str(lin_ineq)])
   disp([' Number of linear equality constraints:      ',int2str(lin_eq)])
   disp([' Number of lower bound constraints:          ',int2str(nnz(~isinf(LB)))])
   disp([' Number of upper bound constraints:          ',int2str(nnz(~isinf(UB)))])
case {'lsqcurvefit','lsqnonlin'}
   disp([' Number of lower bound constraints:          ',int2str(nnz(~isinf(LB)))])
   disp([' Number of upper bound constraints:          ',int2str(nnz(~isinf(UB)))])
case {'fsolve','fminunc','fsolves'}
otherwise
end

if ~isempty(OUTPUT)
temp = sprintf('\n%s\n   %s\n','Algorithm selected',OUTPUT.algorithm);
disp(temp)
end

pstr = sprintf('\n%s\n%s\n',...
   '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%',...
   ' End diagnostic information ');
disp(pstr)

%--------------------------------------------------------------    formula
function args = formula(fun)
%FORMULA Function formula.
%   FORMULA(FUN) returns the formula for the INLINE object FUN.
%
%   See also INLINE/ARGNAMES, INLINE/CHAR.

%   Copyright (c) 1984-98 by The MathWorks, Inc.
%   $Revision: 1.6 $  $Date: 2003/12/02 11:49:54 $

args = fun.expr;


%--------------------------------------------------------------   graderr
function graderr(finite_diff_deriv, analytic_deriv, evalstr2)
%GRADERR Used to check gradient discrepancy in optimization routines. 

%   Copyright (c) 1990-98 by The MathWorks, Inc.
%   $Revision: 1.6 $  $Date: 2003/12/02 11:49:54 $

finite_diff_deriv = full(finite_diff_deriv); 
analytic_deriv = full(analytic_deriv);
size(finite_diff_deriv);
size(analytic_deriv);
err=max(max(abs(analytic_deriv-finite_diff_deriv)));
disp(sprintf('Maximum discrepancy between derivatives  = %g',err));
if (err > 1e-6*norm(analytic_deriv) + 1e-5) 
    disp('Warning: Derivatives do not match within tolerance')
    disp('Derivative from finite difference calculation:')
    finite_diff_deriv
    disp(['User-supplied derivative, ', evalstr2, ' : '])
    analytic_deriv
    disp('Difference:')
    analytic_deriv - finite_diff_deriv
end

%------------------------------------------------------------------   nlconst
function [x,FVAL,lambda_out,EXITFLAG,OUTPUT,GRADIENT,HESS]= ...
   nlconst(funfcn,x,lb,ub,Ain,Bin,Aeq,Beq,confcn,OPTIONS,...
   verbosity,gradflag,gradconstflag,hessflag,meritFunctionType,...
   CHG,fval,gval,Hval,ncineqval,nceqval,gncval,gnceqval,varargin);
%NLCONST Helper function to find the constrained minimum of a function 
%   of several variables. Called by CONSTR, ATTGOAL. SEMINF and MINIMAX.
%
%   [X,OPTIONS,LAMBDA,HESS]=NLCONST('FUN',X0,OPTIONS,lb,ub,'GRADFUN',...
%   varargin{:}) starts at X0 and finds a constrained minimum to 
%   the function which is described in FUN. FUN is a four element cell array
%   set up by PREFCNCHK.  It contains the call to the objective/constraint
%   function, the gradients of the objective/constraint functions, the
%   calling type (used by OPTEVAL), and the calling function name. 
%
%   Copyright (c) 1990-98 by The MathWorks, Inc.
%   $Revision: 1.6 $  $Date: 2003/12/02 11:49:54 $
%   Andy Grace 7-9-90, Mary Ann Branch 9-30-96.

%   Called by CONSTR, SEMINF, ATTGOAL, MINIMAX.
%   Calls OPTEVAL.
%
%   meritFunctionType==5 for fseminf
%                    ==1 for fminimax & fgoalattain (can use 0, but 1 is default)
%                    ==0 for fmincon

FVAL=[];lambda=[];EXITFLAG =1; OUTPUT=[]; HESS=[];
% Expectations: GRADfcn must be [] if it does not exist.
global OPT_STOP OPT_STEP;
OPT_STEP = 1; 
OPT_STOP = 0; 
% Initialize so if OPT_STOP these have values
lambda = []; HESS = [];
iter = 0;
% Set up parameters.
XOUT=x(:);
% numberOfVariables must be the name of this variable
numberOfVariables = length(XOUT);

Nlconst = 'nlconst';
tolX = optimget(OPTIONS,'tolx');
tolFun = optimget(OPTIONS,'tolfun');
tolCon = optimget(OPTIONS,'tolcon');
DiffMinChange = optimget(OPTIONS,'diffminchange');
DiffMaxChange = optimget(OPTIONS,'diffmaxchange');
DerivativeCheck = strcmp(optimget(OPTIONS,'DerivativeCheck'),'on');
maxFunEvals = optimget(OPTIONS,'maxfunevals');
maxIter = optimget(OPTIONS,'maxIter');
% In case the defaults were gathered from calling: optimset('fminsearch'):
if ischar(maxFunEvals)
  maxFunEvals = 100*numberOfVariables;
end


% Handle bounds as linear constraints
arglb = ~isinf(lb);
lenlb=length(lb); % maybe less than numberOfVariables due to old code
argub = ~isinf(ub);
lenub=length(ub);
boundmatrix = eye(max(lenub,lenlb),numberOfVariables);

if nnz(arglb) > 0     
   lbmatrix = -boundmatrix(arglb,1:numberOfVariables);% select non-Inf bounds 
   lbrhs = -lb(arglb);
else
   lbmatrix = []; lbrhs = [];
end

if nnz(argub) > 0
   ubmatrix = boundmatrix(argub,1:numberOfVariables);
   ubrhs=ub(argub);
else
   ubmatrix = []; ubrhs=[];
end 

bestf = Inf; 
if isempty(confcn{1})
   constflag = 0;
else
   constflag = 1;
end

A = [lbmatrix;ubmatrix;Ain];
B = [lbrhs;ubrhs;Bin];

if isempty(A)
   A = zeros(0,numberOfVariables); B=zeros(0,1);
end
if isempty(Aeq)
   Aeq = zeros(0,numberOfVariables); Beq=zeros(0,1);
end


% Used for semi-infinite optimization:
s = nan; POINT =[]; NEWLAMBDA =[]; LAMBDA = []; NPOINT =[]; FLAG = 2;
OLDLAMBDA = [];

x(:) = XOUT;  % Set x to have user expected size
% Compute the objective function and constraints
if strcmp(funfcn{2},'fseminf')
   f = fval;
   [ncineq,nceq,NPOINT,NEWLAMBDA,OLDLAMBDA,LOLD,s] = ...
      semicon(x,LAMBDA,NEWLAMBDA,OLDLAMBDA,POINT,FLAG,s,varargin{:});
else
   f = fval;
   nceq = nceqval; ncineq = ncineqval;  % nonlinear constraints only
   %c = [Aeq*XOUT-Beq; ceq; A*XOUT-B; c];
end

non_eq = length(nceq);
non_ineq = length(ncineq);
[lin_eq,Aeqcol] = size(Aeq);
[lin_ineq,Acol] = size(A);  % includes upper and lower
eq = non_eq + lin_eq;
ineq = non_ineq + lin_ineq;
nc = [nceq; ncineq];

ncstr = ineq + eq;

if isempty(f)
   error('FUN must return a non-empty objective function.')
end

% Evaluate gradients and check size
if gradflag | gradconstflag %evaluate analytic gradient
   if gradflag
      gf_user = gval;
   end
   
   if gradconstflag
      gnc_user = [  gnceqval, gncval];   % Don't include A and Aeq yet
   else
      gnc_user = [];
   end
   if isempty(gnc_user) & isempty(nc)
      % Make gc compatible
      gnc = nc'; gnc_user = nc';
   end % isempty(gnc_user) & isempty(nc)
end
c = [ Aeq*XOUT-Beq; nceq; A*XOUT-B; ncineq];

OLDX=XOUT;
OLDC=c; OLDNC=nc;
OLDgf=zeros(numberOfVariables,1);
gf=zeros(numberOfVariables,1);
OLDAN=zeros(ncstr,numberOfVariables);
LAMBDA=zeros(ncstr,1);

stepsize=1;

if meritFunctionType==1
   if isequal(funfcn{2},'fgoalattain')
      header = sprintf(['\n                    Attainment                 Directional \n',...
                          ' Iter   F-count       factor      Step-size     derivative    Procedure ']);

   else
   header = sprintf(['\n                       Max                     Directional \n',...
                       ' Iter   F-count  {F,constraints}  Step-size     derivative    Procedure ']);
   end
   formatstr = '%5.0f  %5.0f   %12.4g %12.3g    %12.3g   %s  %s';
else % fmincon is caller
   header = sprintf(['\n                                     max                      Directional \n',...
                       ' Iter   F-count      f(x)         constraint    Step-size      derivative   Procedure ']);
   formatstr = '%5.0f  %5.0f   %12.6g %12.4g %12.3g    %12.3g   %s  %s';
end
if verbosity > 1
   disp(header)
end

HESS=eye(numberOfVariables,numberOfVariables);

numFunEvals=1;
numGradEvals=1;
GNEW=1e8*CHG;
%---------------------------------Main Loop-----------------------------
status = 0; EXITFLAG = 1;
while status ~= 1
   iter = iter + 1;
   
   %----------------GRADIENTS----------------
   
   if ~gradconstflag | ~gradflag | DerivativeCheck
      % Finite Difference gradients (even if just checking analytical)
      POINT = NPOINT; 
      oldf = f;
      oldnc = nc;
      len_nc = length(nc);
      ncstr =  lin_eq + lin_ineq + len_nc;     
      FLAG = 0; % For semi-infinite
      gnc = zeros(numberOfVariables, len_nc);  % For semi-infinite
      % Try to make the finite differences equal to 1e-8.
      CHG = -1e-8./(GNEW+eps);
      CHG = sign(CHG+eps).*min(max(abs(CHG),DiffMinChange),DiffMaxChange);
      OPT_STEP = 1;
      for gcnt=1:numberOfVariables
         if gcnt == numberOfVariables, 
            FLAG = -1; 
         end
         temp = XOUT(gcnt);
         XOUT(gcnt)= temp + CHG(gcnt);
         x(:) =XOUT; 
         if ~gradflag | DerivativeCheck
            if strcmp(funfcn{2},'fseminf')
               f= arl2EvalPsi(x,varargin{3:end});
            else
               
               f = arl2EvalPsi(x,varargin{:});
            end
            
            gf(gcnt,1) = (f-oldf)/CHG(gcnt);
         end
         if ~gradconstflag | DerivativeCheck
            if constflag
               if strcmp(confcn{2},'fseminf')
                  [nctmp,nceqtmp,NPOINT,NEWLAMBDA,OLDLAMBDA,LOLD,s] = ...
                     semicon(x,LAMBDA,NEWLAMBDA,OLDLAMBDA,POINT,FLAG,s,varargin{:});
               else
                  [nctmp,nceqtmp] = arl2SPModulus(x,varargin{:});
               end
               nc = [nceqtmp(:); nctmp(:)];
            end
            % Next line used for problems with varying number of constraints
            if len_nc~=length(nc) & isequal(funfcn{2},'fseminf')
               diff=length(nc); 
               nc=v2sort(oldnc,nc); 
               
            end
            
            if ~isempty(nc)
               gnc(gcnt,:) = (nc - oldnc)'/CHG(gcnt); 
            end
           
         end
         
         OPT_STEP = 0;
         if OPT_STOP
            break;
         end
         XOUT(gcnt) = temp;
         
         if OPT_STOP
            break;
         end
      end % for 
      
      % Gradient check
      if DerivativeCheck == 1 & (gradflag | gradconstflag) % analytic exists
                           
         disp('Function derivative')
         if gradflag
            gfFD = gf;
            gf = gf_user;
            
            if isa(funfcn{4},'inline')
               graderr(gfFD, gf, formula(funfcn{4}));
            else
               graderr(gfFD, gf, funfcn{4});
            end
         end
         
         if gradconstflag
            gncFD = gnc; 
            gnc = gnc_user;
            
            disp('Constraint derivative')
            if isa(confcn{4},'inline')
               graderr(gncFD, gnc, formula(confcn{4}));
            else
               graderr(gncFD, gnc, confcn{4});
            end
         end         
         DerivativeCheck = 0;
      elseif gradflag | gradconstflag
         if gradflag
            gf = gf_user;
         end
         if gradconstflag
            gnc = gnc_user;
         end
      end % DerivativeCheck == 1 &  (gradflag | gradconstflag)
      
      FLAG = 1; % For semi-infinite
      numFunEvals = numFunEvals + numberOfVariables;
      f=oldf;
      nc=oldnc;
   else% gradflag & gradconstflag & no DerivativeCheck 
      gnc = gnc_user;
      gf = gf_user;
   end  
   
   % Now add in Aeq, and A
   if ~isempty(gnc)
      gc = [Aeq', gnc(:,1:non_eq), A', gnc(:,non_eq+1:non_ineq+non_eq)];
   elseif ~isempty(Aeq) | ~isempty(A)
      gc = [Aeq',A'];
   else
      gc = zeros(numberOfVariables,0);
   end
   AN=gc';
   how='';
   OPT_STEP = 2;
   
   %-------------SEARCH DIRECTION---------------
   % For equality constraints make gradient face in 
   % opposite direction to function gradient.
   for i=1:eq 
      schg=AN(i,:)*gf;
      if schg>0
         AN(i,:)=-AN(i,:);
         c(i)=-c(i);
      end
   end
   
   if numGradEvals>1  % Check for first call    
      if meritFunctionType~=5,   
         NEWLAMBDA=LAMBDA; 
      end
      [ma,na] = size(AN);
      try, GNEW=gf+AN'*NEWLAMBDA;
      catch
	disp('The GNEW bug!')
	size(gf)
	size(AN)
	size(NEWLAMBDA)
      end
      GOLD=OLDgf+OLDAN'*LAMBDA;
      YL=GNEW-GOLD;
      sdiff=XOUT-OLDX;
      % Make sure Hessian is positive definite in update.
      if YL'*sdiff<stepsize^2*1e-3
         while YL'*sdiff<-1e-5
            [YMAX,YIND]=min(YL.*sdiff);
            YL(YIND)=YL(YIND)/2;
         end
         if YL'*sdiff < (eps*norm(HESS,'fro'));
            how=' Hessian modified twice';
            FACTOR=AN'*c - OLDAN'*OLDC;
            FACTOR=FACTOR.*(sdiff.*FACTOR>0).*(YL.*sdiff<=eps);
            WT=1e-2;
            if max(abs(FACTOR))==0; FACTOR=1e-5*sign(sdiff); end
            while YL'*sdiff < (eps*norm(HESS,'fro')) & WT < 1/eps
               YL=YL+WT*FACTOR;
               WT=WT*2;
            end
         else
            how=' Hessian modified';
         end
      end
      
      %----------Perform BFGS Update If YL'S Is Positive---------
      if YL'*sdiff>eps
         HESS=HESS ...
            +(YL*YL')/(YL'*sdiff)-((HESS*sdiff)*(sdiff'*HESS'))/(sdiff'*HESS*sdiff);
         % BFGS Update using Cholesky factorization  of Gill, Murray and Wright.
         % In practice this was less robust than above method and slower. 
         %   R=chol(HESS); 
         %   s2=R*S; y=R'\YL; 
         %   W=eye(numberOfVariables,numberOfVariables)-(s2'*s2)\(s2*s2') + (y'*s2)\(y*y');
         %   HESS=R'*W*R;
      else
         how=' Hessian not updated';
      end
      
   else % First call
      OLDLAMBDA=(eps+gf'*gf)*ones(ncstr,1)./(sum(AN'.*AN')'+eps) ;
   end % if numGradEvals>1
   numGradEvals=numGradEvals+1;
   
   LOLD=LAMBDA;
   OLDAN=AN;
   OLDgf=gf;
   OLDC=c;
   OLDF=f;
   OLDX=XOUT;
   XN=zeros(numberOfVariables,1);
   if (meritFunctionType>0 & meritFunctionType<5)
      % Minimax and attgoal problems have special Hessian:
      HESS(numberOfVariables,1:numberOfVariables)=zeros(1,numberOfVariables);
      HESS(1:numberOfVariables,numberOfVariables)=zeros(numberOfVariables,1);
      HESS(numberOfVariables,numberOfVariables)=1e-8*norm(HESS,'inf');
      XN(numberOfVariables)=max(c); % Make a feasible solution for qp
   end
   
   GT =c;
   
   HESS = (HESS + HESS')*0.5;
   [SD,lambda,exitflagqp,outputqp,howqp] ...
      = qpsub(HESS,gf,AN,-GT,[],[],XN,eq,-1, ...
      Nlconst,size(AN,1),numberOfVariables); 
   
   lambda((1:eq)') = abs(lambda( (1:eq)' ));
   ga=[abs(c( (1:eq)' )) ; c( (eq+1:ncstr)' ) ];
   if ~isempty(c)
      mg=max(ga);
   else
      mg = 0;
   end
   
   arl2ExitFlag = arl2CheckConstraints(iter, f, nc); 
   if verbosity>1
      if strncmp(howqp,'ok',2); 
         howqp =''; 
      end
      if ~isempty(how) & ~isempty(howqp) 
         how = [how,'; '];
      end
      if meritFunctionType==1,
         gamma = mg+f;
         CurrOutput = sprintf(formatstr,iter,numFunEvals,gamma,stepsize,gf'*SD,how,howqp); 
         disp(CurrOutput)

      else
        CurrOutput = sprintf(formatstr,iter,numFunEvals,f,mg,stepsize,gf'*SD,how,howqp); 
        disp(CurrOutput)
        % disp([sprintf('%5.0f %12.6g %12.6g ',numFunEvals,f,mg), ...
        %       sprintf('%12.3g  ',stepsize),how, ' ',howqp]);
      end
   end
   LAMBDA=lambda((1:ncstr)');
   OLDLAMBDA=max([LAMBDA';0.5*(LAMBDA+OLDLAMBDA)'])' ;
   
   %---------------LINESEARCH--------------------
   MATX=XOUT;
   MATL = f+sum(OLDLAMBDA.*(ga>0).*ga) + 1e-30;
   infeas = strncmp(howqp,'i',1);
   if meritFunctionType==0 | meritFunctionType == 5
      % This merit function looks for improvement in either the constraint
      % or the objective function unless the sub-problem is infeasible in which
      % case only a reduction in the maximum constraint is tolerated.
      % This less "stringent" merit function has produced faster convergence in
      % a large number of problems.
      if mg > 0
         MATL2 = mg;
      elseif f >=0 
         MATL2 = -1/(f+1);
      else 
         MATL2 = 0;
      end
      if ~infeas & f < 0
         MATL2 = MATL2 + f - 1;
      end
   else
      % Merit function used for MINIMAX or ATTGOAL problems.
      MATL2=mg+f;
   end
   if mg < eps & f < bestf
      bestf = f;
      bestx = XOUT;
      bestHess = HESS;
      bestgrad = gf;
      bestlambda = lambda;
   end
   MERIT = MATL + 1;
   MERIT2 = MATL2 + 1; 
   stepsize=2;
   while  (MERIT2 > MATL2) & (MERIT > MATL) ...
         & numFunEvals < maxFunEvals & ~OPT_STOP
      stepsize=stepsize/2;
      if stepsize < 1e-4,  
         stepsize = -stepsize; 
         
         % Semi-infinite may have changing sampling interval
         % so avoid too stringent check for improvement
         if meritFunctionType == 5, 
            stepsize = -stepsize; 
            MATL2 = MATL2 + 10; 
         end
      end
      XOUT = MATX + stepsize*SD;
      x(:)=XOUT; 
      
      if strcmp(funfcn{2},'fseminf')
         f= arl2EvalPsi(x,varargin{3:end});
         
         [nctmp,nceqtmp,NPOINT,NEWLAMBDA,OLDLAMBDA,LOLD,s] = ...
            semicon(x,LAMBDA,NEWLAMBDA,OLDLAMBDA,POINT,FLAG,s,varargin{:});
         nctmp = nctmp(:); nceqtmp = nceqtmp(:);
         non_ineq = length(nctmp);  % the length of nctmp can change
         ineq = non_ineq + lin_ineq;
         ncstr = ineq + eq;
         
      else
         f = arl2EvalPsi(x,varargin{:});
         if constflag
            [nctmp,nceqtmp] = arl2SPModulus(x,varargin{:});
            nctmp = nctmp(:); nceqtmp = nceqtmp(:);
         else
            nctmp = []; nceqtmp=[];
         end
      end
            
      nc = [nceqtmp(:); nctmp(:)];
      c = [Aeq*XOUT-Beq; nceqtmp(:); A*XOUT-B; nctmp(:)];  
      
      if OPT_STOP
         break;
      end
      
      numFunEvals = numFunEvals + 1;
      ga=[abs(c( (1:eq)' )) ; c( (eq+1:length(c))' )];
      if ~isempty(c)
         mg=max(ga);
      else
         mg = 0;
      end
      
      MERIT = f+sum(OLDLAMBDA.*(ga>0).*ga);
      if meritFunctionType==0 | meritFunctionType == 5
         if mg > 0
            MERIT2 = mg;
         elseif f >=0 
            MERIT2 = -1/(f+1);
         else 
            MERIT2 = 0;
         end
         if ~infeas & f < 0
            MERIT2 = MERIT2 + f - 1;
         end
      else
         MERIT2=mg+f;
      end
   end  % line search loop
   %------------Finished Line Search-------------
   
   if meritFunctionType~=5
      mf=abs(stepsize);
      LAMBDA=mf*LAMBDA+(1-mf)*LOLD;
   end
   % Test stopping conditions (convergence)
   if (max(abs(SD)) < 2*tolX | abs(gf'*SD) < 2*tolFun ) & ...
         (mg < tolCon | (strncmp(howqp,'i',1) & mg > 0 ) )
      if verbosity>0
         if ~strncmp(howqp, 'i', 1) 
            disp('Optimization terminated successfully:')
            if max(abs(SD)) < 2*tolX 
               disp(' Search direction less than 2*options.TolX and')
               disp('  maximum constraint violation is less than options.TolCon')
            else
               disp(' Magnitude of directional derivative in search direction ')
               disp('  less than 2*options.TolFun and maximum constraint violation ')
               disp('  is less than options.TolCon')     
            end
            
            active_const = find(LAMBDA>0);
            if active_const 
               disp('Active Constraints:'), 
               disp(active_const) 
            else % active_const == 0
               disp(' No Active Constraints');
            end 
         end
         
         if (strncmp(howqp, 'i',1) & mg > 0)
            disp('Optimization terminated: No feasible solution found.')
            if max(abs(SD)) < 2*tolX 
               disp(' Search direction less than 2*options.TolX but constraints are not satisfied.')
            else
               disp(' Magnitude of directional derivative in search direction ')
               disp('  less than 2*options.TolFun but constraints are not satisfied.')    
            end
            EXITFLAG = -1;   
         end
      end
      status=1;
   else % continue
      % NEED=[LAMBDA>0] | G>0
      if numFunEvals > maxFunEvals  | OPT_STOP
         XOUT = MATX;
         f = OLDF;
         if ~OPT_STOP
            if verbosity > 0
               disp('Maximum number of function evaluations exceeded;')
               disp('increase OPTIONS.MaxFunEvals')
            end
         end
         EXITFLAG = 0;
         status=1;
      end
      if iter > maxIter
         XOUT = MATX;
         f = OLDF;
         if verbosity > 0
            disp('Maximum number of function evaluations exceeded;')
            disp('increase OPTIONS.MaxIter')
         end
         EXITFLAG = 0;
         status=1;
      end
      if arl2ExitFlag == 0
         XOUT = MATX;
         f = OLDF;
         if verbosity > 0
            disp('Exit on constraints saturation')
         end
         EXITFLAG = 0;
         status=1;
      end
   end 
   
   x(:) = XOUT;
   switch funfcn{1} % evaluate function gradients
   case 'fun'
      ;  % do nothing...will use finite difference.
   case 'fungrad'
      [f,gf_user] = arl2EvalPsi(x,varargin{:});
      gf_user = gf_user(:);
      numGradEvals=numGradEvals+1;
   case 'fun_then_grad'
      gf_user = feval(funfcn{4},x,varargin{:});
      gf_user = gf_user(:);
      numGradEvals=numGradEvals+1;
   otherwise
      error('Undefined calltype in FMINCON');
   end
   numFunEvals=numFunEvals+1;
   
   
   % Evaluate constraint gradients
   switch confcn{1}
   case 'fun'
      gnceq=[]; gncineq=[];
   case 'fungrad'
      [nctmp,nceqtmp,gncineq,gnceq] = arl2SPModulus(x,varargin{:});
      nctmp = nctmp(:); nceqtmp = nceqtmp(:);
      numGradEvals=numGradEvals+1;
   case 'fun_then_grad'
      [gncineq,gnceq] = feval(confcn{4},x,varargin{:});
      numGradEvals=numGradEvals+1;
   case ''
      nctmp=[]; nceqtmp =[];
      gncineq = zeros(numberOfVariables,length(nctmp));
      gnceq = zeros(numberOfVariables,length(nceqtmp));
      
   otherwise
      error('Undefined calltype in FMINCON');
   end
   gnc_user = [gnceq, gncineq];
   gc = [Aeq', gnceq, A', gncineq];
   
   
end % while status ~= 1

% Update 
numConstrEvals = numGradEvals;

% Gradient is in the variable gf
GRADIENT = gf;

% If a better unconstrained solution was found earlier, use it:
if f > bestf 
   XOUT = bestx;
   f = bestf;
   HESS = bestHess;
   GRADIENT = bestgrad;
   lambda = bestlambda;
end

FVAL = f;
x(:) = XOUT;
if (OPT_STOP)
   if verbosity > 0
      disp('Optimization terminated prematurely by user')
   end
end


OUTPUT.iterations = iter;
OUTPUT.funcCount = numFunEvals;
OUTPUT.stepsize = stepsize;
OUTPUT.algorithm = 'medium-scale: SQP, Quasi-Newton, line-search';
OUTPUT.firstorderopt = [];
OUTPUT.cgiterations = [];

[lin_ineq,Acol] = size(Ain);  % excludes upper and lower

lambda_out.lower=zeros(lenlb,1);
lambda_out.upper=zeros(lenub,1);

lambda_out.eqlin = lambda(1:lin_eq);
ii = lin_eq ;
lambda_out.eqnonlin = lambda(ii+1: ii+ non_eq);
ii = ii+non_eq;
lambda_out.lower(arglb) = lambda(ii+1 :ii+nnz(arglb));
ii = ii + nnz(arglb) ;
lambda_out.upper(argub) = lambda(ii+1 :ii+nnz(argub));
ii = ii + nnz(argub);
lambda_out.ineqlin = lambda(ii+1: ii + lin_ineq);
ii = ii + lin_ineq ;
lambda_out.ineqnonlin = lambda(ii+1 : end);

%-------------------------------------------------------------------   qpsub
function [X,lambda,exitflag,output,how]=...
   qpsub(H,f,A,B,lb,ub,X,neqcstr,verbosity,caller,ncstr,numberOfVariables,options)
%QP Quadratic programming subproblem. Handles qp and constrained
%   linear least-squares as well as subproblems generated from NLCONST.
%
%   X=QP(H,f,A,b) solves the quadratic programming problem:
%
%            min 0.5*x'Hx + f'x   subject to:  Ax <= b 
%             x    
%

%   Copyright (c) 1990-98 by The MathWorks, Inc.
%   $Revision: 1.6 $  $Date: 2003/12/02 11:49:54 $
%   Andy Grace 7-9-90. Mary Ann Branch 9-30-96.

% Define constant strings
NewtonStep = 'Newton';
SteepDescent = 'steepest descent';
Conls = 'lsqlin';
Lp = 'linprog';
Qp = 'quadprog';
Qpsub = 'qpsub';
how = 'ok'; 
exitflag = 1;
output = [];
iterations = 0;
if nargin < 13
   options = [];
end

lb=lb(:); ub = ub(:);

msg = nargchk(12,12,nargin);
if isempty(verbosity), verbosity = 1; end
if isempty(neqcstr), neqcstr = 0; end

LLS = 0;
if strcmp(caller, Conls)
   LLS = 1;
   [rowH,colH]=size(H);
   numberOfVariables = colH;
end
if strcmp(caller, Qpsub)
   normalize = -1;
else
   normalize = 1;
end

simplex_iter = 0;
if  norm(H,'inf')==0 | isempty(H), is_qp=0; else, is_qp=1; end



if LLS==1
   is_qp=0;
end

normf = 1;
if normalize > 0
   % Check for lp
   if ~is_qp & ~LLS
      normf = norm(f);
      if normf > 0
         f = f./normf;
      end
   end
end

% Handle bounds as linear constraints
arglb = ~eq(lb,-inf);
lenlb=length(lb); % maybe less than numberOfVariables due to old code
if nnz(arglb) > 0     
   lbmatrix = -eye(lenlb,numberOfVariables);
   
   A=[A; lbmatrix(arglb,1:numberOfVariables)]; % select non-Inf bounds
   B=[B;-lb(arglb)];
end


argub = ~eq(ub,inf);
lenub=length(ub);
if nnz(argub) > 0
   ubmatrix = eye(lenub,numberOfVariables);
   A=[A; ubmatrix(argub,1:numberOfVariables)];
   B=[B; ub(argub)];
end 
ncstr=ncstr + nnz(arglb) + nnz(argub);

% Figure out max iteration count
% For now, don't limit the iterations when qpsub is called from nlconst
maxSQPiters = Inf;
maxiter = optimget(options,'MaxIter',maxSQPiters);

% Used for determining threshold for whether a direction will violate
% a constraint.
normA = ones(ncstr,1);
if normalize > 0 
   for i=1:ncstr
      n = norm(A(i,:));
      if (n ~= 0)
         A(i,:) = A(i,:)/n;
         B(i) = B(i)/n;
         normA(i,1) = n;
      end
   end
else 
   normA = ones(ncstr,1);
end
errnorm = 0.01*sqrt(eps); 

tolDep = 100*numberOfVariables*eps;      
lambda=zeros(ncstr,1);
aix=lambda;
ACTCNT=0;
ACTSET=[];
ACTIND=0;
CIND=1;
eqix = 1:neqcstr; 

%------------EQUALITY CONSTRAINTS---------------------------
Q = zeros(numberOfVariables,numberOfVariables);
R = []; 
indepInd = 1:ncstr; 

if neqcstr>0
   % call equality constraint solver
   [Q,R,A,B,CIND,X,Z,actlambda,how,...
         ACTSET,ACTIND,ACTCNT,aix,eqix,neqcstr,ncstr,remove,exitflag]= ...
      eqnsolv(A,B,eqix,neqcstr,ncstr,numberOfVariables,LLS,H,X,f,normf,normA,verbosity, ...
      aix,how,exitflag);   
   
   if ~isempty(remove)
      indepInd(remove)=[];
      normA = normA(indepInd);
   end
   
   if ACTCNT >= numberOfVariables - 1  
      simplex_iter = 1; 
   end
   [m,n]=size(ACTSET);
   
   if strcmp(how,'infeasible')
      % Equalities are inconsistent, so X and lambda have no valid values
      % Return original X and zeros for lambda.
      output.iterations = iterations;
      return
   end
   
   err = 0;
   if neqcstr > numberOfVariables
      err = max(abs(A(eqix,:)*X-B(eqix)));
      if (err > 1e-8)  % Equalities not met
         how='infeasible';
         % was exitflag = 7; 
         exitflag = -1;
         if verbosity > 0 
            disp('Exiting: The equality constraints are overly stringent;')
            disp('         there is no feasible solution.')
         end
         % Equalities are inconsistent, X and lambda have no valid values
         % Return original X and zeros for lambda.
         output.iterations = iterations;
         return
      else % Check inequalities
         if (max(A*X-B) > 1e-8)
            how = 'infeasible';
            % was exitflag = 8; 
            exitflag = -1;
            if verbosity > 0
               disp('Exiting: The constraints or bounds are overly stringent;')
               disp('         there is no feasible solution.')
               disp('         Equality constraints have been met.')
            end
         end
      end
      if is_qp
         actlambda = -R\(Q'*(H*X+f));
      elseif LLS
         actlambda = -R\(Q'*(H'*(H*X-f)));
      else
         actlambda = -R\(Q'*f);
      end
      lambda(indepInd(eqix)) = normf * (actlambda ./normA(eqix));
      output.iterations = iterations;
      return
   end
   if isempty(Z)
      if is_qp
         actlambda = -R\(Q'*(H*X+f));
      elseif LLS
         actlambda = -R\(Q'*(H'*(H*X-f)));
      else
         actlambda = -R\(Q'*f);
      end
      lambda(indepInd(eqix)) = normf * (actlambda./normA(eqix));
      if (max(A*X-B) > 1e-8)
         how = 'infeasible';
         % was exitflag = 8; 
         exitflag = -1;
         if verbosity > 0
            disp('Exiting: The constraints or bounds are overly stringent;')
            disp('         there is no feasible solution.')
            disp('         Equality constraints have been met.')
         end
      end
      output.iterations = iterations;
      return
   end
   
   
   % Check whether in Phase 1 of feasibility point finding. 
   if (verbosity == -2)
      cstr = A*X-B; 
      mc=max(cstr(neqcstr+1:ncstr));
      if (mc > 0)
         X(numberOfVariables) = mc + 1;
      end
   end
else
   Z=1;
end

% Find Initial Feasible Solution
cstr = A*X-B;
mc=max(cstr(neqcstr+1:ncstr));
if mc>eps
   A2=[[A;zeros(1,numberOfVariables)],[zeros(neqcstr,1);-ones(ncstr+1-neqcstr,1)]];
   quiet = -2;
   [XS,lambdaS,exitflagS,outputS] = qpsub([],[zeros(numberOfVariables,1);1],A2,[B;1e-5], ...
      [],[],[X;mc+1],neqcstr,quiet,Qpsub,size(A2,1),numberOfVariables+1);
   X=XS(1:numberOfVariables);
   cstr=A*X-B;
   if XS(numberOfVariables+1)>eps 
      if XS(numberOfVariables+1)>1e-8 
         how='infeasible';
         % was exitflag = 4; 
         exitflag = -1;
         if verbosity > 0
            disp('Exiting: The constraints are overly stringent;')
            disp('         no feasible starting point found.')
         end
      else
         how = 'overly constrained';
         % was exitflag = 3; 
         exitflag = -1;
         if verbosity > 0
            disp('Exiting: The constraints are overly stringent;')
            disp(' initial feasible point found violates constraints ')
            disp(' by more than eps.');
         end
      end
      lambda(indepInd) = normf * (lambdaS((1:ncstr)')./normA);
      output.iterations = iterations;
      return
   end
end

if (is_qp)
   gf=H*X+f;
   %  SD=-Z*((Z'*H*Z)\(Z'*gf));
   [SD, dirType] = compdir(Z,H,gf,numberOfVariables,f);
   
   % Check for -ve definite problems:
   %  if SD'*gf>0, is_qp = 0; SD=-SD; end
elseif (LLS)
   HXf=H*X-f;
   gf=H'*(HXf);
   HZ= H*Z;
   [mm,nn]=size(HZ);
   if mm >= nn
      %   SD =-Z*((HZ'*HZ)\(Z'*gf));
      [QHZ, RHZ] =  qr(HZ);
      Pd = QHZ'*HXf;
      % Now need to check which is dependent
      if min(size(RHZ))==1 % Make sure RHZ isn't a vector
         depInd = find( abs(RHZ(1,1)) < tolDep);
      else
         depInd = find( abs(diag(RHZ)) < tolDep );
      end  
   end
   if mm >= nn & isempty(depInd) % Newton step
      SD = - Z*(RHZ(1:nn, 1:nn) \ Pd(1:nn,:));
      dirType = NewtonStep;
   else % steepest descent direction
      SD = -Z*(Z'*gf);
      dirType = SteepDescent;
   end
else % lp
   gf = f;
   SD=-Z*Z'*gf;
   dirType = SteepDescent; 
   if norm(SD) < 1e-10 & neqcstr
      % This happens when equality constraint is perpendicular
      % to objective function f.x.
      actlambda = -R\(Q'*(gf));
      lambda(indepInd(eqix)) = normf * (actlambda ./ normA(eqix));
      output.iterations = iterations;
      return;
   end
end

oldind = 0; 

% The maximum number of iterations for a simplex type method is when ncstr >=n:
% maxiters = prod(1:ncstr)/(prod(1:numberOfVariables)*prod(1:max(1,ncstr-numberOfVariables)));

%--------------Main Routine-------------------
while 1 & iterations <= maxiter
   iterations = iterations + 1;
   if isinf(verbosity)
      curr_out = sprintf('Iter: %5.0f, Active: %5.0f, step: %s, proc: %s',iterations,ACTCNT,dirType,how);
      disp(curr_out); 
   end
   
   % Find distance we can move in search direction SD before a 
   % constraint is violated.
   % Gradient with respect to search direction.
   GSD=A*SD;
   
   % Note: we consider only constraints whose gradients are greater
   % than some threshold. If we considered all gradients greater than 
   % zero then it might be possible to add a constraint which would lead to
   % a singular (rank deficient) working set. The gradient (GSD) of such
   % a constraint in the direction of search would be very close to zero.
   indf = find((GSD > errnorm * norm(SD))  &  ~aix);
   
   if isempty(indf) % No constraints to hit
      STEPMIN=1e16;
      dist=[]; ind2=[]; ind=[];
   else % Find distance to the nearest constraint
      dist = abs(cstr(indf)./GSD(indf));
      [STEPMIN,ind2] =  min(dist);
      ind2 = find(dist == STEPMIN);
      % Bland's rule for anti-cycling: if there is more than one 
      % blocking constraint then add the one with the smallest index.
      ind=indf(min(ind2));
      % Non-cycling rule:
      % ind = indf(ind2(1));
   end
   %-----Update X-------------
   
   % Assume we do not delete a constraint
   delete_constr = 0;   
   
   if ~isempty(indf)& isfinite(STEPMIN) % Hit a constraint
      if strcmp(dirType, NewtonStep)
         % Newton step and hit a constraint: LLS or is_qp
         if STEPMIN > 1  % Overstepped minimum; reset STEPMIN
            STEPMIN = 1;
            delete_constr = 1;
         end
         X = X+STEPMIN*SD;
      else
         % Not a Newton step and hit a constraint: is_qp or LLS or maybe lp
         X = X+STEPMIN*SD;          
      end              
   else %  isempty(indf) | ~isfinite(STEPMIN)
      % did not hit a constraint
      if strcmp(dirType, NewtonStep)
         % Newton step and no constraint hit: LLS or maybe is_qp
         STEPMIN = 1;   % Exact distance to the solution. Now delete constr.
         X = X + SD;
         delete_constr = 1;
      else % Not a Newton step: is_qp or lp or LLS
         if is_qp
            % Is it semi-def, neg-def or indef?
            eigoptions.disp = 0;
            ZHZ = Z'*H*Z;
            if numberOfVariables < 400 % only use EIGS on large problems
               [VV,DD] = eig(ZHZ);
               [smallRealEig, eigind] = min(diag(DD));
               ev = VV(:,eigind(1));
            else
               [ev,smallRealEig,flag] = eigs(ZHZ,1,'sr',eigoptions);
               if flag  % Call to eigs failed
                  [VV,DD] = eig(ZHZ);
                  [smallRealEig, eigind] = min(diag(DD));
                  ev = VV(:,eigind(1));
               end
            end
            
         else % define smallRealEig for LLS
            smallRealEig=0;
         end
         
         if (~is_qp & ~LLS) | (smallRealEig < -100*eps) % LP or neg def: not LLS
            % neg def -- unbounded
            if norm(SD) > errnorm
               if normalize < 0
                  STEPMIN=abs((X(numberOfVariables)+1e-5)/(SD(numberOfVariables)+eps));
               else 
                  STEPMIN = 1e16;
               end
               X=X+STEPMIN*SD;
               how='unbounded'; 
               % was exitflag = 5; 
               exitflag = -1;
            else % norm(SD) <= errnorm
               how = 'ill posed';
               % was exitflag = 6; 
               exitflag = -1;

            end
            if verbosity > 0
               if norm(SD) > errnorm
                  disp('Exiting: The solution is unbounded and at infinity;')
                  disp('         the constraints are not restrictive enough.') 
               else
                  disp('Exiting: The search direction is close to zero; ')
                  disp('      the problem is ill-posed.')
                  disp('      The gradient of the objective function may be zero')
                  disp('         or the problem may be badly conditioned.')
               end
            end % if verbosity > 0
            output.iterations = iterations;
            return
         else % singular: solve compatible system for a solution: is_qp or LLS
            if is_qp
               projH = Z'*H*Z; 
               Zgf = Z'*gf;
               projSD = pinv(projH)*(-Zgf);
            else % LLS
               projH = HZ'*HZ; 
               Zgf = Z'*gf;
               projSD = pinv(projH)*(-Zgf);
            end
            
            % Check if compatible
            if norm(projH*projSD+Zgf) > 10*eps*(norm(projH) + norm(Zgf))
               % system is incompatible --> it's a "chute": use SD from compdir
               % unbounded in SD direction
               if norm(SD) > errnorm
                  if normalize < 0
                     STEPMIN=abs((X(numberOfVariables)+1e-5)/(SD(numberOfVariables)+eps));
                  else 
                     STEPMIN = 1e16;
                  end
                  X=X+STEPMIN*SD;
                  how='unbounded'; 
                  % was exitflag = 5;
                  exitflag = -1;
               else % norm(SD) <= errnorm
                  how = 'ill posed';
                  %was exitflag = 6;
                  exitflag = -1;
               end
               if verbosity > 0
                  if norm(SD) > errnorm
                     disp('Exiting: The solution is unbounded and at infinity;')
                     disp('         the constraints are not restrictive enough.') 
                  else
                     disp('Exiting: The search direction is close to zero; ')
                     disp('      the problem is ill-posed.')
                     disp('      The gradient of the objective function may be zero')
                     disp('         or the problem may be badly conditioned.')
                  end
               end % if verbosity > 0
               output.iterations = iterations;
               return
            else % Convex -- move to the minimum (compatible system)
               SD = Z*projSD;
               dirType = 'singular';
               % First check if constraint is violated.
               GSD=A*SD;
               indf = find((GSD > errnorm * norm(SD))  &  ~aix);
               if isempty(indf) % No constraints to hit
                  STEPMIN=1;
                  delete_constr = 1;
                  dist=[]; ind2=[]; ind=[];
               else % Find distance to the nearest constraint
                  dist = abs(cstr(indf)./GSD(indf));
                  [STEPMIN,ind2] =  min(dist);
                  ind2 = find(dist == STEPMIN);
                  % Bland's rule for anti-cycling: if there is more than one 
                  % blocking constraint then add the one with the smallest index.
                  ind=indf(min(ind2));
               end
               if STEPMIN > 1  % Overstepped minimum; reset STEPMIN
                  STEPMIN = 1;
                  delete_constr = 1;
               end
               X = X + STEPMIN*SD; 
            end
         end % if ~is_qp | smallRealEig < -eps
      end % if strcmp(dirType, NewtonStep)
   end % if ~isempty(indf)& isfinite(STEPMIN) % Hit a constraint
   
   if delete_constr
      % Note: only reach here if a minimum in the current subspace found
      if ACTCNT>0
         if ACTCNT>=numberOfVariables-1, 
            % Avoid case when CIND is greater than ACTCNT
            if CIND <= ACTCNT
               ACTSET(CIND,:)=[];
               ACTIND(CIND)=[]; 
            end
         end
         if is_qp
            rlambda = -R\(Q'*(H*X+f));
         elseif LLS
            rlambda = -R\(Q'*(H'*(H*X-f)));
            % else: lp does not reach this point
         end
         actlambda = rlambda;
         actlambda(eqix) = abs(rlambda(eqix));
         indlam = find(actlambda < 0);
         if (~length(indlam)) 
            lambda(indepInd(ACTIND)) = normf * (rlambda./normA(ACTIND));
            output.iterations = iterations;
            return
         end
         % Remove constraint
         lind = find(ACTIND == min(ACTIND(indlam)));
         lind=lind(1);
         ACTSET(lind,:) = [];
         aix(ACTIND(lind)) = 0;
         [Q,R]=qrdelete(Q,R,lind);
         ACTIND(lind) = [];
         ACTCNT = ACTCNT - 2;
         simplex_iter = 0;
         ind = 0;
      else % ACTCNT == 0
         output.iterations = iterations;
         return
      end
      delete_constr = 0;
   end
   
   % Calculate gradient w.r.t objective at this point
   if is_qp
      gf=H*X+f;
   elseif LLS % LLS
      gf=H'*(H*X-f);
      % else gf=f still true.
   end
   
   
   % Update X and calculate constraints
   cstr = A*X-B;
   cstr(eqix) = abs(cstr(eqix));
   % Check no constraint is violated
   if normalize < 0 
      if X(numberOfVariables,1) < eps
         output.iterations = iterations;
         return;
      end
   end
   
   if max(cstr) > 1e5 * errnorm
      if max(cstr) > norm(X) * errnorm 
         if ( verbosity > 0 ) & ( exitflag == 1 )
            disp('Note: The problem is badly conditioned;')
            disp('         the solution may not be reliable') 
            % verbosity = 0;
         end
         how='unreliable'; 
         % exitflag = 2;
         exitflag = -1;
         if 0
            X=X-STEPMIN*SD;
            output.iterations = iterations;
            return
         end
      end
   end
   
   if ind % Hit a constraint
      aix(ind)=1;
      ACTSET(CIND,:)=A(ind,:);
      ACTIND(CIND)=ind;
      [m,n]=size(ACTSET);
      [Q,R] = qrinsert(Q,R,CIND,A(ind,:)');
   end
   if oldind 
      aix(oldind) = 0; 
   end
   if ~simplex_iter
      % Z = null(ACTSET);
      [m,n]=size(ACTSET);
      Z = Q(:,m+1:n);
      ACTCNT=ACTCNT+1;
      if ACTCNT == numberOfVariables - 1, simplex_iter = 1; end
      CIND=ACTCNT+1;
      oldind = 0; 
   else
      rlambda = -R\(Q'*gf);
      
      if isinf(rlambda(1)) & rlambda(1) < 0 
         fprintf('         Working set is singular; results may still be reliable.\n');
         [m,n] = size(ACTSET);
         rlambda = -(ACTSET + sqrt(eps)*randn(m,n))'\gf;
      end
      actlambda = rlambda;
      actlambda(eqix)=abs(actlambda(eqix));
      indlam = find(actlambda<0);
      if length(indlam)
         if STEPMIN > errnorm
            % If there is no chance of cycling then pick the constraint 
            % which causes the biggest reduction in the cost function. 
            % i.e the constraint with the most negative Lagrangian 
            % multiplier. Since the constraints are normalized this may 
            % result in less iterations.
            [minl,CIND] = min(actlambda);
         else
            % Bland's rule for anti-cycling: if there is more than one 
            % negative Lagrangian multiplier then delete the constraint
            % with the smallest index in the active set.
            CIND = find(ACTIND == min(ACTIND(indlam)));
         end
         
         [Q,R]=qrdelete(Q,R,CIND);
         Z = Q(:,numberOfVariables);
         oldind = ACTIND(CIND);
      else
         lambda(indepInd(ACTIND))= normf * (rlambda./normA(ACTIND));
         output.iterations = iterations;
         return
      end
   end %if ACTCNT<numberOfVariables
      
   if (is_qp)
      Zgf = Z'*gf; 
      if ~isempty(Zgf) & (norm(Zgf) < 1e-15) 
         SD = zeros(numberOfVariables,1); 
      else
         [SD, dirType] = compdir(Z,H,gf,numberOfVariables,f);
      end
   elseif (LLS)
      Zgf = Z'*gf;
      HZ = H*Z;
      if (norm(Zgf) < 1e-15)
         SD = zeros(numberOfVariables,1);
      else
         HXf=H*X-f;
         gf=H'*(HXf);
         [mm,nn]=size(HZ);
         if mm >= nn
            [QHZ, RHZ] =  qr(HZ);
            Pd = QHZ'*HXf;
            % SD = - Z*(RHZ(1:nn, 1:nn) \ Pd(1:nn,:));
            % Now need to check which is dependent
            if min(size(RHZ))==1 % Make sure RHZ isn't a vector
               depInd = find( abs(RHZ(1,1)) < tolDep);
            else
               depInd = find( abs(diag(RHZ)) < tolDep );
            end  
         end
         if mm >= nn & isempty(depInd) % Newton step
            SD = - Z*(RHZ(1:nn, 1:nn) \ Pd(1:nn,:));
            dirType = NewtonStep;
         else % steepest descent direction
            SD = -Z*(Z'*gf);
            dirType = SteepDescent;
         end
      end
   else % LP
      if ~simplex_iter
         SD = -Z*(Z'*gf);
         gradsd = norm(SD);
      else
         gradsd = Z'*gf;
         if  gradsd > 0
            SD = -Z;
         else
            SD = Z;
         end
      end
      if abs(gradsd) < 1e-10  % Search direction null
         % Check whether any constraints can be deleted from active set.
         % rlambda = -ACTSET'\gf;
         if ~oldind
            rlambda = -R\(Q'*gf);
         end
         actlambda = rlambda;
         actlambda(1:neqcstr) = abs(actlambda(1:neqcstr));
         indlam = find(actlambda < errnorm);
         lambda(indepInd(ACTIND)) = normf * (rlambda./normA(ACTIND));
         if ~length(indlam)
            output.iterations = iterations;
            return
         end
         cindmax = length(indlam);
         cindcnt = 0;
         newactcnt = 0;
         while (abs(gradsd) < 1e-10) & (cindcnt < cindmax)
            cindcnt = cindcnt + 1;
            if oldind
               % Put back constraint which we deleted
               [Q,R] = qrinsert(Q,R,CIND,A(oldind,:)');
            else
               simplex_iter = 0;
               if ~newactcnt
                  newactcnt = ACTCNT - 1;
               end
            end
            CIND = indlam(cindcnt);
            oldind = ACTIND(CIND);
            
            [Q,R]=qrdelete(Q,R,CIND);
            [m,n]=size(ACTSET);
            Z = Q(:,m:n);
            
            if m ~= numberOfVariables
               SD = -Z*Z'*gf;
               gradsd = norm(SD);
            else
               gradsd = Z'*gf;
               if  gradsd > 0
                  SD = -Z;
               else
                  SD = Z;
               end
            end
         end
         if abs(gradsd) < 1e-10  % Search direction still null
            output.iterations = iterations;
            return;
         end
         lambda = zeros(ncstr,1);
         if newactcnt 
            ACTCNT = newactcnt;
         end
      end
   end
   
   if simplex_iter & oldind
      % Avoid case when CIND is greater than ACTCNT
      if CIND <= ACTCNT
         ACTIND(CIND)=[];
         ACTSET(CIND,:)=[];
         CIND = numberOfVariables;
      end
   end 
end % while 1
if iterations > maxiter
   exitflag = 0;
   how = 'ill-conditioned';   
end

output.iterations = iterations;


%-------------------------------------------------------------------   eqnsolv
function[Q,R,A,B,CIND,X,Z,actlambda,how,...
      ACTSET,ACTIND,ACTCNT,aix,eqix,neqcstr,ncstr,remove,exitflag]= ...
   eqnsolv(A,B,eqix,neqcstr,ncstr,numberOfVariables,LLS,H,X,f,normf,normA,verbosity, ...
   aix,how,exitflag)
% EQNSOLV Helper function for QPSUB.
%    Finds a feasible point with respect to the equality constraints.
%    If the equalities are dependent but not consistent, warning
%    messages are given. If the equalities are dependent but consistent, 
%    the redundant constraints are removed and the corresponding variables 
%    adjusted.

% set tolerances
tolDep = 100*numberOfVariables*eps;      
tolCons = 1e-10;

actlambda = [];
aix(eqix)=ones(neqcstr,1);
ACTSET=A(eqix,:);
ACTIND=eqix;
ACTCNT=neqcstr;
CIND=neqcstr+1;
Z=[]; Anew=[]; Bnew=[]; remove =[];

% See if the equalities form a consistent system:
%   QR factorization of A
[Qa,Ra,Ea]=qr(A(eqix,:));
% Now need to check which is dependent
if min(size(Ra))==1 % Make sure Ra isn't a vector
   depInd = find( abs(Ra(1,1)) < tolDep);
else
   depInd = find( abs(diag(Ra)) < tolDep );
end
if neqcstr > numberOfVariables
   depInd = [depInd; ((numberOfVariables+1):neqcstr)'];
end      

if ~isempty(depInd)
   if verbosity > 0
      disp('The equality constraints are dependent.')
   end
   how='dependent';
   exitflag = 1;
   bdepInd =  abs(Qa(:,depInd)'*B(eqix)) >= tolDep ;
   
   if any( bdepInd ) % Not consistent
      how='infeasible';   
      exitflag = 9;exitflag = -1;
      if verbosity > 0
         disp('The system of equality constraints is not consistent.');
         if ncstr > neqcstr
            disp('The inequality constraints may or may not be satisfied.');
         end
         disp('  There is no feasible solution.')
      end
   else % the equality constraints are consistent
      numDepend = nnz(depInd);
      % delete the redundant constraints:
      % By QR factoring the transpose, we see which columns of A'
      %   (rows of A) move to the end
      [Qat,Rat,Eat]=qr(ACTSET');        
      [i,j] = find(Eat); % Eat permutes the columns of A' (rows of A)
      remove = i(depInd);
      if verbosity > 0
         disp('The system of equality constraints is consistent. Removing');
         disp('the following dependent constraints before continuing:');
         disp(remove)
      end
      A(eqix(remove),:)=[];
      B(eqix(remove))=[];
      neqcstr = neqcstr - nnz(remove);
      ncstr = ncstr - nnz(remove);
      eqix = 1:neqcstr;
      aix=[ones(neqcstr,1); zeros(ncstr-neqcstr,1)];
      ACTIND = eqix;
      ACTSET=A(eqix,:);
      
      CIND = neqcstr+1;
      ACTCNT = neqcstr;
   end % consistency check
end % dependency check

if ~strcmp(how,'infeasible')
   % Find a feasible point
   if max(abs(A(eqix,:)*X-B(eqix))) > tolCons
      X = A(eqix,:)\B(eqix);  
   end
end

[Q,R]=qr(ACTSET');
Z = Q(:,neqcstr+1:numberOfVariables);

% End of eqnsolv.m



%-------------------------------------------------------------------   semicon
function  [nctmp,nceqtmp,NPOINT,NEWLAMBDA,OLDLAMBDA,LOLD,s] = ...
            semicon(x,LAMBDA,NEWLAMBDA,OLDLAMBDA,POINT,FLAG,s,varargin)
error('Semicon should not be called');
%-------------------------------------------------------------------   startx
function xstart = startx(u,l);
%STARTX	Box-centered point
%
% xstart = STARTX(u,l) returns centered point.

%   Copyright (c) 1990-98 by The MathWorks, Inc.
%   $Revision: 1.6 $  $Date: 2003/12/02 11:49:54 $

n = length(u);
onen = ones(n,1);
arg = (u > 1e12);
u(arg) = inf*onen(arg);
xstart = zeros(n,1);
arg1 = (u<inf)&(l==-inf); arg2 = (u== inf)&(l > -inf);
arg3 = (u<inf)&(l>-inf);  arg4 = (u==inf)&(l==-inf);
%
w = max(abs(u),ones(n,1));
xstart(arg1) = u(arg1) - .5*w(arg1);
%
ww = max(abs(l),ones(n,1));
xstart(arg2) = l(arg2) + .5*ww(arg2);
%
xstart(arg3)=(u(arg3)+l(arg3))/2;
xstart(arg4)=ones(length(arg4(arg4>0)),1);




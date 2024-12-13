function [val,w,g,Hes,error_to_be_min]=EvaluateCausalLagrangian(lmult,L,t_flag)

% Evaluates the langrangian at a dual point lmult. L is a structure
% containing all neccessary parameters
% This is a special version of EvaluateLagrangian with a few optimized features
% and independent of it so it can be modified if another type of contrained probleme
% is needed

switch (L.cflag)
case 1
    % Min. causality under modulus contrains
    l2=lmult;
    DD=(L.D)'*diag(l2)*L.D;
    M=L.BB+DD;
    %eig(M)
    %pause;
    w=M\(-1*transpose(conj(L.B))*L.f);
    error_to_be_min=norm(L.B*w+L.f);
    err2=L.D*w;
    g=abs(err2).^2-(L.mc).^2;
    val=error_to_be_min^2+transpose(l2)*g;
    errt=transpose(conj(L.D))*diag(err2);
    Hes=-2*real(transpose(conj(errt))*(M\errt));    
case 2
    % Min. l2 error under moudlus and causality constraints
    l1=lmult(1);
    l2=lmult(2:length(lmult));
    
    DD=transpose(conj(L.D))*diag(l2)*L.D;
    M=L.AA+l1*L.BB+DD;
    Minv=inv(M);
    w=Minv*(-l1*transpose(conj(L.B))*L.f+transpose(conj(L.A))*L.y);
    err1=L.B*w+L.f;
    err2=L.D*w;
    g=norm(err1)^2-(L.c)^2;
    g=[g(1);abs(err2).^2-(L.mc).^2];
    error_to_be_min=norm(L.A*w-L.y);
    val=(error_to_be_min)^2+l1*g(1)+transpose(l2)*g(2:length(g));
    errb=transpose(conj(L.B))*err1;
    errd=transpose(conj(L.D))*diag(err2);
    errt=[errb,errd];
    Hes=-2*real(transpose(conj(errt))*Minv*errt);
end

if (t_flag)
	kk=length(lmult);
        H_t=[];
	eps_t=0.000000001;
	for j=1:kk
        lmult_t=lmult;
	[val_t,v_t,g_t0]=EvaluateCausalLagrangian(lmult_t,L,0);
        lmult_t(j)=lmult_t(j)+eps_t;
	[val_t,v_t,g_t1]=EvaluateCausalLagrangian(lmult_t,L,0);
        H_t=[H_t,(g_t1-g_t0)/eps_t];
	end
	fprintf('Error in hessian %f %% Error in sym: %f %%',norm(Hes-H_t)/norm(Hes)*100,norm(Hes-transpose(Hes))/norm(Hes)*100);
end	

% $Id: EvaluateCausalLagrangian.m,v 1.4 2006/05/16 16:01:08 fseyfert Exp $
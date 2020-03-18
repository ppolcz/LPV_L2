function [gam,Pvars]=lpvL2gain(lpvsys,pargrd,Pbase,Pbase_der,gammin)
% Computes an upper bound for the induced L2 gain of an LPV system. 
% Syntax: [gam,Pvars,chkres]=LPVL2GAIN(lpvsys,pargrd,Pbase,Pbase_der,gammin) 
% Arguments: 
%     lpvsys:     function evaluating the system matrices at p. 
%                 format: function [A,B,C,D]=lpvsys(p)  
%     pargrd:     grid for each parameter and their derivatives,
%                 e.g. if the LPV system depends on two parameters p1 and
%                 p2, where p2 is the time derivative of p1 and we choose
%                 parameter-dependent lyapunov function P(p1,p2) then the time
%                 derivative of P will depend on p1, p2 and p3=d(p2)/dt. So
%                 pargrd should be defined like this:
%                 pargrd={linspace(-1,1,10),-2:0.2:2, [-1 1]}, where the 
%                 linspace(-1,1,10) is for p1, -2:0.2:2 is for p2 and [-1 1] 
%                 is for p3. 
%     Pbase:      basis functions determining the parameter dependence of the Lyapunov function
%                 e.g. Pbase=@(p) [1 p(1) p(2)]; this means P(p1,p2)=P0+p1*P1+p2*P2
%     Pbase_der:  time derivative of Pbase, i.e.
%                 Pbase_der=@(p)= [0 p(2) p(3)]; 
%     gammin:     for numerical reasons it is advisable to give some lower
%                 bound for the gamma. The default value is 1e-5.
%
% Algorithm: the classic method of evaluating the Lyapunov function and the
%            bounded real lemma condition at each grid point and solving
%            the huge set of LMIs simultaneously.  
%
% Outputs: 
%     gam:        the gamma value
%     Pvars:      the components of the lyapunov matrix s.t.
%                 P(p1,p2)=Pvars{1}+p1*Pvars{2}+p2*Pvars{3}



%% Demo
if nargin==0
    clear all; close all;
    lpvsys=@(p) lpvsys_harald(p(1));
    
    pargrd={linspace(2,7,100),[-1,1]}; %Contains both the parameters and their derivatives
    %Pbase    =@(p) [1 p(1) p(1)^2 p(1)^3 p(1)^4 p(1)^5 p(1)^6 1/p(1) 1/p(1)^2 1/p(1)^3];
    %Pbase_der=@(p) [0 1 2*p(1) 3*p(1)^2 4*p(1)^3 5*p(1)^4 6*p(1)^5 -1/p(1)^2 -2*p(1)/p(1)^4 -3*p(1)^2/p(1)^6]*p(2);
    
    Pbase    =@(p) [1];
    Pbase_der=@(p) [0];
    
    gammin=1;
    
    [gam,~]=lpvL2gain(lpvsys,pargrd,Pbase,Pbase_der,gammin);
    return
end;
    
%% Initialization

if ~exist('gammin','var'), gammin=1e-5; end;

% generating the parameter grid
pgrd=multigrid(pargrd{:})
sizeofgrd=size(pgrd,1); 

pgrd(1,:)
Pbase
Pbase(pgrd(1,:))

% initialization of some variables
nbf=length(Pbase(pgrd(1,:)));
[~,B,~,~]=lpvsys(pgrd(1,:));
[nx,nw]=size(B); 

%% Main program

posdef=@(LMI) (LMI>=1e-8*eye(size(LMI)));

% defining the optimization variables
yalmip('clear');
Pvars=cell(1,nbf); 
for i=1:nbf
    Pvars{i}=sdpvar(nx);
end;
gam2=sdpvar(1);
Cons=posdef(gam2-gammin*gammin);

% constructing the LMIs for Bounded Real Lemma 
         
fprintf('Collecting the LMI constraints.........');    
for k=1:sizeofgrd
    fprintf('\b\b\b\b\b\b\b%6.2f%%',100*k/sizeofgrd); 
    [A,B,C,D]=lpvsys(pgrd(k,:));
    bf=Pbase(pgrd(k,:));
    bfdot=Pbase_der(pgrd(k,:));
    P=0; Pdot=0;
    for i=1:nbf
        P=P+bf(i)*Pvars{i};
        Pdot=Pdot+bfdot(i)*Pvars{i};
    end;
    BRL=-[Pdot+P*A+A'*P+C'*C P*B+C'*D; B'*P+D'*C D'*D-gam2*eye(nw)];
    Cons=[Cons, posdef(P), posdef(BRL)];
end;

% calling the solver
fprintf('\nCalling the solver..............');
options=sdpsettings('verbose',0,'solver','mosek');
tic;
sol=solvesdp(Cons,gam2,options);
telapsed=toc;
if sol.problem ~= 0
    fprintf('\n    Something went wrong, but we check the results though!');
else
    fprintf('OK (T=%6.2f)',telapsed);
end;
gam=sqrt(double(gam2)); 
for i=1:nbf
    Pvars{i}=double(Pvars{i});
end;


fprintf('\nChecking the result.............');
[pres,dres]=check(Cons);
msg='OK';
if any(pres<0) ||  any(dres<0), 
    for k=1:sizeofgrd
        [A,B,C,D]=lpvsys(pgrd(k,:));
        bf=Pbase(pgrd(k,:));
        bfdot=Pbase_der(pgrd(k,:));
        P=0; Pdot=0;
        for i=1:nbf
            P=P+bf(i)*Pvars{i};
            Pdot=Pdot+bfdot(i)*Pvars{i};
        end;
        BRL=-[Pdot+P*A+A'*P+C'*C P*B+C'*D; B'*P+D'*C D'*D-gam^2*eye(nw)];
        if any(eig(P)<0) || any(eig(BRL)<0), 
            msg='FAILED'; 
            break;
        end;
    end;
end;
fprintf('%s\n',msg);
fprintf('Gamma value:  %6.4f\n',gam);


    
 

function [A,B,C,D]=lpvsys_harald(p)
% Function to compute the LPV state-space matrices for Harald's system 
% (LPV plant+gain scheduling PI controller)
% p is 1-dimensional with range [2,7] and rate [-1,1]


% Parameter-dependent gain and time-constant for plant
tau=sqrt(133.6-16.8*p);
K=sqrt(4.8*p-8.6);

% PI gains for controller
Kp=(2*0.7*0.25*tau-1)/K;
Ki=0.25^2*tau/K;

% State-space matrices for closed-loop sensitivity function
A=[-1/tau*(1+Kp*K) 1/tau; -Ki*K 0];
B=[1/tau*Kp;Ki];
C=[-K 0]; 
D=1;

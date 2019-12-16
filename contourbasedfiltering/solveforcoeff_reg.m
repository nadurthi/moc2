function [lam,err]=solveforcoeff_reg(A,b,Aineq,bineq,Nlam2ndorder,maxrelerr_percent,method)
% Solve Ac=b with Aineq*c<= bineq
% g is gradient
% J is hessian
% onecol: col-number in A that has all 1s
% 1 to Nlam2ndorder is the number of lam coeff for second order polynomials
% keyboard
% if rcond(A)<10*eps
%     keyboard
%
%     Am = zeros(size(A));
%     bm = zeros(size(b));
%
%
%     for i=1:size(A)
%         mx = max(abs(A(i,:)));
%         mn = min(abs(A(i,:)));
%         mm = mean([mx,mn]);
%
%         Am(i,:) = A(i,:)/mm;
%         bm(i) = b(i)/mm;
%     end
%
%     Aineqm = zeros(size(Aineq));
%     bineqm = zeros(size(bineq));
%     for i=1:size(Aineq)
%         mx = max(abs(Aineq(i,:)));
%         mn = min(abs(Aineq(i,:)));
%         mm = mean([mx,mn]);
%
%         Aineqm(i,:) = Aineq(i,:)/mm;
%         bineqm(i) = bineq(i)/mm;
%     end
%
%
%     A = Am;
%     b = bm;
%     Aineq = Aineqm;
%     bineq = bineqm;
% end

% options = optimoptions('fminunc','Algorithm','trust-region','SpecifyObjectiveGradient',true,'HessianFcn','objective');
% options = optimoptions('fminunc','Algorithm','trust-region',...
%     'SpecifyObjectiveGradient',true,'HessianFcn','objective');

% options = optimoptions('fminunc','Algorithm','trust-region','Display','iter');


[Neq,lamdim]=size(A);
[Nineq,~]=size(Aineq);

Cvec = [0.01,0.1,25,100,250,300,500,800,1200,1900,2500,3500,5000,7500,10000];
%%


% x=[lam,teq];


%% inf norm on error
if strcmp(method,'teqinfnorm_lam1norm')
options = optimoptions('fmincon','MaxIterations',1e6,'Display','iter');
c0 = zeros(lamdim,1);
C=100;
nonlcon=[];
lb= [];
ub= [];
nonlcon = [];
lam = fmincon(@(lam)C*norm(lam,1)+norm(A*lam-b,inf),c0,Aineq,bineq,A,b,lb,ub,nonlcon,options);
end

%% 2 norm on error
if strcmp(method,'teq2norm_lam1norm')
c0 = zeros(lamdim,1);
C=100;
lam = fmincon(@(lam)C*norm(lam,1)+norm(A*lam-b,2),c0,Aineq,bineq,A,b,lb,ub,nonlcon,options);
end
%% LP formulation of 1-norm on lam and inf-norm on teq
if strcmp(method,'LPguru')
    options = optimoptions('linprog','Algorithm','interior-point','Display','iter');
    % x = [lams,tlams,tmax]';
    lb = [];
    ub = [];
    C=100;
    
    
    
    SII=sparse(eye(lamdim));
    Zlamdim = sparse(zeros(lamdim,1));
    Zineq1 = sparse(zeros(Nineq,1));
    Zineqlamdim = sparse(zeros(Nineq,lamdim));
    Zeqlamdim =  sparse(zeros(Neq,lamdim));
    Zlamdimlamdim = sparse(zeros(lamdim,lamdim));
    
    AAineq=[SII,-SII,Zlamdim;
        -SII,-SII,Zlamdim;
        A,Zeqlamdim,-ones(Neq,1);
        -A,Zeqlamdim,-ones(Neq,1);
        Aineq,Zineqlamdim,Zineq1;
        Zlamdimlamdim,-SII,Zlamdim];
    BBineq = [Zlamdim;
        Zlamdim;
        b;
        -b;
        bineq;
        Zlamdim];
    
    % II=eye(lamdim);
    % AAineq=[II,-II,zeros(lamdim,1);
    %         -II,-II,zeros(lamdim,1);
    %         A,zeros(Neq,lamdim),-ones(Neq,1);
    %         -A,zeros(Neq,lamdim),-ones(Neq,1);
    %         Aineq,zeros(Nineq,lamdim),zeros(Nineq,1);
    %         zeros(lamdim,lamdim),-II,zeros(lamdim,1)];
    % BBineq = [zeros(lamdim,1);
    %             zeros(lamdim,1);
    %             b;
    %             -b;
    %             bineq;
    %             zeros(lamdim,1)];
    % x = linprog(f,AAineq,BBineq,[],[],lb,ub,options);
    prevlam=0;
    cnt=0;
    for C=Cvec
        f=[zeros(lamdim,1);C/lamdim*ones(lamdim,1);10];
        [x,fval,exitflag,output,lambda] = linprog_guru(f,AAineq,BBineq,[],[],lb,ub,options);
        lam = x(1:lamdim);
        if exitflag<0
            continue
        end
        
        if cnt==0
            prevlam = lam;
        end
        
        err=A*lam-b;
        [C, 100*abs( err(1)/(b(1)) ),exitflag]
        if 100*abs( err(1)/(b(1)) )>maxrelerr_percent && exitflag>=0
            break
        end
        prevlam = lam;
        cnt=cnt+1;
    end
    lam = prevlam;
end
%% Finalize


err=A*lam-b;

end
%
% function [f, g, H] = rosenboth(x)
% % Calculate objective f
% f = 100*(x(2) - x(1)^2)^2 + (1-x(1))^2;
%
% if nargout > 1 % gradient required
%     g = [-400*(x(2)-x(1)^2)*x(1)-2*(1-x(1));
%         200*(x(2)-x(1)^2)];
%
%     if nargout > 2 % Hessian required
%         H = [1200*x(1)^2-400*x(2)+2, -400*x(1);
%             -400*x(1), 200];
%     end
%
% end
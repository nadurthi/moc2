function pdfnorm = get_interp_pdf_0I_2D_gobackget_L1(model,priorfullnormpdf,Xk1,probsk1,Nm,Tstepk1,Tk1,dtkk1,Xmctest,Xtruth)
% known prior pdf priorpdfnormk at time k
% compute pdfnorm or pdfnormk1 at time k+1
% Tk1 is time at k+1,
% dtkk1 is deltat from k to k+1

%%
[N,dim] =size(Xk1);

[mXk1,PXk1]=MeanCov(Xk1,probsk1/sum(probsk1));


dsX = DataSet(Xk1,probsk1,'TrueState');

[a,b] = max(probsk1);
mXk1=Xk1(b,:)';
dsX.AddMeanCov_to_OI_Trasform(mXk1,5^2*PXk1);
% dsX.AddHyperCubeTrasform(-0.9*ones(dim,1),0.9*ones(dim,1));

LB=-1.0*ones(dim,1);
UB=1.0*ones(dim,1);

% plot3(dsX.X(:,1),dsX.X(:,2),dsX.p,'ro')

% keyboard

% indd = dsX.p<1e-100;
% dsX.p(indd) = 1e-100;

dsX.SortByProb('descend');

logpn = dsX.getlogProb();

if isempty(Xmctest)==0
    [Xnmctest,~] = dsX.ApplyAffineTransform_Original2Final(Xmctest, zeros(size(Xmctest,1),1) );
end
if isempty(Xtruth)==0
    [Xntruth,~] =  dsX.ApplyAffineTransform_Original2Final(Xtruth, zeros(size(Xtruth,1),1) );
end

[mt,Pt] = MeanCov(dsX.X,dsX.p/sum(dsX.p));


transForms_atk = priorfullnormpdf.transForms;
normpdf_atk = priorfullnormpdf.func;
normpdfexp_atk = priorfullnormpdf.expfunc;
transForms_atk1 = dsX.GetTrasnformers();

%% plottinmg
for i=0:0
    figure(80+i)
    dsX.PlotPointsProbs3D([i+1,i+2],'ro');
    hold on
    if isempty(Xmctest)==0
        plot(Xnmctest(:,i+1),Xnmctest(:,i+2),'bs')
    end
    if isempty(Xtruth)==0
        plot(Xntruth(:,i+1),Xntruth(:,i+2),'k*')
    end
    title(['true points and MC: time step = ',num2str(Tstepk1)])
    hold off
end


%% Checking points

% first get the test points
Xcheckk1 = mvurnd(1.3*LB,1.3*UB,dim*1000);
probscheck1 = zeros(size(Xcheckk1,1),1);
logprobscheck1 = zeros(size(Xcheckk1,1),1);

for i=1:size(Xcheckk1,1)
    xnormk1 = Xcheckk1(i,:);
    %     [probscheck1(i),logprobscheck1(i)] = getbacknormprobs(model,Tk1,dtkk1,normpdf_atk,normpdfexp_atk,xnormk1,transForms_atk1,transForms_atk);
    [probscheck1(i),logprobscheck1(i),~,~] = getbacknormprobs(xnormk1,Tk1,dtkk1,priorfullnormpdf,transForms_atk1,model);
end

% logprobscheck1=log(probscheck1+1e-12);
figure(71)
plot3(Xcheckk1(:,1),Xcheckk1(:,2),probscheck1,'bo')
xlabel('x')
xlabel('y')
title('probscheck vs Xcheck')
axis([-1.2,1.2,-1.2,1.2])


%% sanity points

% first get the test points
Xsanity = mvurnd(1.5*LB,1.5*UB,dim*10000);

%% outside points
% between 1 and 1.3 are considered outside points
[Xout,~]=smolyak_sparse_grid_modf(zeros(dim,1),0.8^2*eye(dim),dim,8,'gh');
indxout = isinbox(Xout,1.5*LB,1.5*UB);
Xout = Xout(~indxout,:);

%% go back and get interpolation points


[Xall,~]=smolyak_sparse_grid_modf(zeros(dim,1),eye(dim),dim,8,'gh');
Xall=0.85*Xall/max(max(Xall));

probsk1qi = zeros(size(Xall,1),1);
logprobsk1qi = zeros(size(Xall,1),1);
xkvec = zeros(size(Xall,1),dim);
xknormvec = zeros(size(Xall,1),dim);

for i=1:size(Xall,1)
    %     [probsk1qi(i),logprobsk1qi(i),xkvec(i,:),xknormvec(i,:)] = getbacknormprobs(model,Tk1,dtkk1,normpdf_atk,normpdfexp_atk,Xall(i,:)',transForms_atk1,transForms_atk);
    [probsk1qi(i),logprobsk1qi(i),xkvec(i,:),xknormvec(i,:)] = getbacknormprobs(Xall(i,:)',Tk1,dtkk1,priorfullnormpdf,transForms_atk1,model);
end

Fall = probsk1qi;
logFall = logprobsk1qi;

[yy,ind] = sort(Fall,'descend');
Xall = Xall(ind,:);
Fall = Fall(ind,:);
logFall = logFall(ind,:);

% number of points inside at back time  k
insidept = isinbox(xknormvec,priorfullnormpdf.LB,priorfullnormpdf.UB);
disp('xnorm_atk : inside vs out')
[sum(insidept),length(insidept)]


figure(86)
plot3(Xall(:,1),Xall(:,2),Fall,'bo')
title('Xall and Fall')

FuncTablelog = [Xall,logFall];

cutofflogprob = -9;

indeq = logFall>=cutofflogprob;
Xtrain = Xall(indeq,:);
Xineq = Xall(~indeq,:);

Pf=Basis_polyND(dim,Nm);
lamdim = length(Pf);

Beqtrain = logFall(indeq);
Bineq = cutofflogprob*ones(size(Xineq,1),1);

Nout = size(Xout,1);

A=zeros(size(Xall,1),lamdim);
Aout=zeros(size(Xout,1),lamdim);
%     Aineq=zeros(size(Xextra,1),lamdim);
for ib=1:length(Pf)
    A(:,ib)=evaluate_PolyforLargetSetX(Pf{ib},Xall);
    Aout(:,ib)=evaluate_PolyforLargetSetX(Pf{ib},Xout);
end
Aeqtrain = A(indeq,:);
Aineq = A(~indeq,:);

% Aineq = [Aineq;Aout];
% Bineq = [Bineq;2*cutofflogprob*ones(size(Xout,1),1)];
Bout = normpdfexp(Xout,zeros(1,dim),1^2*eye(dim));
% Bout =  -4*ones(size(Xout,1),1);


Neq = size(Aeqtrain,1);
prevmxentpoly = 0;
% keyboard
for C = [0.01,0.1,25,100,250,300,500,800,1200,1900,2500,3500,5000,7500,10000] % 4,8,
    cvx_begin
    variables teq(Neq)  lam(lamdim)
    minimize( C*norm(lam,1)/lamdim+10*norm_largest(teq,floor(0.75*lamdim))  )   %*norm(teq,2)/Neq  berhu huber norm_largest sum(huber(teq,1e-5))
    subject to
    Aeqtrain*lam==Beqtrain+teq;
    %         Aeqtrain(1:10,:)*lam==Beqtrain(1:10);
    Aineq*lam <= Bineq;
    Aout*lam <= Bout;
    cvx_end
    
    [max(abs(lam)),mean(abs(teq)),100*abs(teq(1))/Beqtrain(1)]
    
    lamsol = lam;
    
    mxentpoly=zeros(1,dim+1);
    for i=1:lamdim
        mxentpoly=add_sub_polyND(mxentpoly, scalar_multiply_polyND(lamsol(i),Pf{i}),'add');
    end
    mxentpoly=simplify_polyND(mxentpoly);
    
    c0 = get_coeff_NDpoly(mxentpoly,zeros(1,dim));
    mxentpoly = update_or_insert_coeff_NDpoly(mxentpoly,zeros(1,dim),c0);
    
    Fsanity = exp(evaluate_polyND_highprec(mxentpoly,(Xsanity)));
    if any(isinf(Fsanity))
        continue
    end
    [C,100*abs(teq(1))/Beqtrain(1) , contains(cvx_status,'Solved')]
    keyboard
    if 100*abs(teq(1))/Beqtrain(1)>=0.01 && contains(cvx_status,'Solved')
        break
    end
    
    prevmxentpoly = mxentpoly;
    
end
mxentpoly = prevmxentpoly;

%     C = 8; % 4,8,
%     cvx_begin
%         variables teq(Neq)  lam(lamdim)
%         minimize( C*norm(lam,1)/lamdim+10*norm_largest(teq,floor(1/2*lamdim))  )   %*norm(teq,2)/Neq  berhu huber norm_largest sum(huber(teq,1e-5))
%         subject to
%         Aeqtrain*lam==Beqtrain+teq;
% %         Aeqtrain(1:10,:)*lam==Beqtrain(1:10);
%         Aineq*lam <= Bineq;
%         Aout*lam <= Bout;
%     cvx_end
%
%     lamsol = lam;
%     max(abs(lam))
%
%     mxentpoly=zeros(1,dim+1);
%     for i=1:lamdim
%         mxentpoly=add_sub_polyND(mxentpoly, scalar_multiply_polyND(lamsol(i),Pf{i}),'add');
%     end
%     mxentpoly=simplify_polyND(mxentpoly);
%     Fsanity = exp(evaluate_polyND_highprec(mxentpoly,(Xsanity)));

Festlog=evaluate_polyND_highprec(mxentpoly,(Xall));
Flogcheck = evaluate_polyND(mxentpoly,Xcheckk1);
Fcheck = exp(Flogcheck);


figure(55)
plot3(Xall(:,1),Xall(:,2),Fall,'bo',Xall(:,1),Xall(:,2),exp(Festlog),'r+')
xlabel('x')
ylabel('y')
title('All points used for interpolation')
axis([-1.5,1.5,-1.5,1.5,0,10])
legend('true','est')

figure(59)
plot3(Xcheckk1(:,1),Xcheckk1(:,2),probscheck1,'bo',Xcheckk1(:,1),Xcheckk1(:,2),Fcheck,'r+')
xlabel('x')
ylabel('y')
title('Checking points')
axis([-1.5,1.5,-1.5,1.5,0,10])
legend('true','est')

figure(60)
plot3(Xsanity(:,1),Xsanity(:,2),Fsanity,'bo')
xlabel('x')
ylabel('y')
title('Sanity points')
axis([-1.5,1.5,-1.5,1.5,0,10])


mxerr = max(abs(probscheck1-Fcheck));
inn = probscheck1>1e-10;
avgrelerr = 100*mean(abs(probscheck1(inn)-Fcheck(inn))./(probscheck1(inn)+eps));
[mxerr,avgrelerr]






mxentpoly=zeros(1,dim+1);
for i=1:lamdim
    mxentpoly=add_sub_polyND(mxentpoly, scalar_multiply_polyND(lamsol(i),Pf{i}),'add');
end
mxentpoly=simplify_polyND(mxentpoly);





%%
keyboard
%% truth test
probest_truth = exppdf_withbounds(dsX.X,mxentpoly,LB,UB,exp(-10));
% [probest_truth,dsX.p]
figure(89)
plot3(dsX.X(:,1),dsX.X(:,2),probest_truth,'bo',dsX.X(:,1),dsX.X(:,2),dsX.p,'r+')
xlabel('x')
ylabel('y')
title('true probs and interpolated probs for prop dsX.X')
legend('estprob','trueprob')
axis([-1.2,1.2,-1.2,1.2])


%% normalize
xl = LB;
xu =  UB;
vv = prod(xu-xl);
for Ninteg=11
    [xint,wint] = GLeg_pts(Ninteg*ones(1,dim), xl, xu);
    pp=exppdf_withbounds(xint,mxentpoly,LB,UB,exp(-10));
    cc=vv*wint(:)'*pp(:);
end
disp(['norm constant is = ',num2str(cc)])

c0 = get_coeff_NDpoly(mxentpoly,zeros(1,dim));
mxentpoly = update_or_insert_coeff_NDpoly(mxentpoly,zeros(1,dim),c0-log(cc));
%%
pdfnorm.LB = LB;
pdfnorm.UB = UB;

pdfnorm.dim =dim;
pdfnorm.mxentpoly_norm = mxentpoly;

pdfnorm.func=@(x)exppdf_withbounds(x,mxentpoly,LB,UB,exp(-10));
pdfnorm.expfunc=@(x)expofaexppdf_withbounds_2D(x,mxentpoly,LB,UB,-10);

pdfnorm.transForms = dsX.GetTrasnformers();

pdfnorm.info = 'norm-pdf';
pdfnorm.pdftype = 'ExpPdf';

% % [Xqi,pts1Dqi,interpPoly1Dqi]=getsparsePts(dim+3,dim);
% [Xqi,~] = GH_points(zeros(dim,1),(1)^2*eye(dim),4);
% % mm=max(max(Xqi));
% % Xqi = 1.1*Xqi/mm;
% normprobs = pdfnorm.func(Xqi);
% Xnew=pdfnorm.transForms.normX2trueX(Xqi);
% probsnew=pdfnorm.transForms.normprob2trueprob(normprobs);
%
%
% figure(100)
% plot3(Xnew(:,1),Xnew(:,2),probsnew,'bo',Xk1(:,1),Xk1(:,2),probsk1,'r+')
%

% pdfnorm.RigTree = RigTree;

% normconst = integratorFuncTrueX_usingpdfnorm(pdfnorm,@(x)constantfunc(x,1),'RegTreeBoxIntegrator');
%
% pdfnorm.func=@(x)(1/normconst)*exp(evaluate_polyND(mxentpoly_norm,x));

%

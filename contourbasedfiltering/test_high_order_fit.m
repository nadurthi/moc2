




clc
close all
clear

format longg

digits(50)

%% time


time.t0=0;
time.tf=10;
time.dt=1;

time.Tvec=time.t0:time.dt:time.tf;
time.Ntsteps=length(time.Tvec);

%% models

model.f=@(dt,tk,xk)duff_prop_model(dt,xk);
model.fback=@(dt,tk,xk)duff_prop_model_back(dt,xk);

model.fn=2;
model.Q = diag([(1e-6)^2,(1e-6)^2]);
model.fstates = {'x1','x2','v1','v2'};

model.h=@(x)duff_meas_model(x);
model.hn=1;
model.hvec=@(x)sensmodel2vec(x,model.h,model.hn);


model.hstates = {'th'};
model.R=diag([(5*pi/180)^2]);


x0=[5,5]';
P0=0.5^2*eye(2);

% xf0=mvnrnd(x0(:)',P0);
xf0 = x0;
Pf0 = P0;


[X1,~] = GH_points(xf0(:),Pf0,4);
probsk1 = mvnpdf(X1,xf0(:)',Pf0);
dsX = DataSet(X1,probsk1,'TrueState');
dsX.AddMeanCov_to_OI_Trasform(xf0(:),Pf0);
% dsX.AddHyperCubeTrasform(-1*ones(dim,1),1*ones(dim,1));
transForms = dsX.GetTrasnformers();

pdfnorm.func=@(x)mvnpdf(x,zeros(1,model.fn),eye(model.fn));
pdfnorm.expfunc=@(x)normpdfexp(x,zeros(1,model.fn),eye(model.fn));
pdfnorm.transForms = dsX.GetTrasnformers();

pdfnorm.info = 'true-0I';
pdfnorm.pdftype = 'ExpPdf';

% pdfnorm.GMMHull = GMMHull;
pdfnorm.LB = -3*ones(1,model.fn);
pdfnorm.UB = 3*ones(1,model.fn);



[probscheck1(i),logprobscheck1(i)] = getbacknormprobs(model,Tk1,dtkk1,normpdf_atk,normpdfexp_atk,xnorm,transForms_atk1,transForms_atk);




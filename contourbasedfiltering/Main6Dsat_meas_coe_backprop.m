%%
% Satellite problem

clc
close all
clear

format longg

digits(50)
%% constants
constants.radii=[6378.137,6378.137,6378.137];

constants.mu      = 3.986004418e5;     % Gravitational Const
constants.Re      = constants.radii(1);          % Earth radius (km)

constants.g0   = 9.8065;            % Sea-level acceleration, (m/s^2)

% Canonical Units
constants.muCan   = 1;
constants.RU      = constants.Re;
constants.TU      = sqrt(constants.RU^3 / constants.mu);
constants.VU      = constants.RU/constants.TU;

constants.trueA2normA=(constants.TU^2/constants.RU);
constants.normA2trueA=(constants.RU/constants.TU^2);

constants.trueV2normV=(constants.TU/constants.RU);
constants.normV2trueV=(constants.RU/constants.TU);

constants.trueX2normX=(1/constants.RU);
constants.normX2trueX=(constants.RU);

constants.trueT2normT=(1/constants.TU);
constants.normT2trueT=(constants.TU);

%% time


time.t0=0 *constants.trueT2normT;
time.tf=48*60*60 *constants.trueT2normT;
time.dt=6*60*60 *constants.trueT2normT;
time.dtplot=0.1*60*60 *constants.trueT2normT;

time.Tvec=time.t0:time.dt:time.tf;
time.Ntsteps=length(time.Tvec);

time.TvecPlots=time.t0: time.dtplot :time.tf;

%% models

model.f=@(dt,tk,xk)processmodel_2bodycoe(dt,1,tk,xk);
model.fback=@(dt,tk,xk)processmodel_2body_backwardcoe(dt,1,tk,xk);

% [a,e,i,om,Om,M]
model.fn=6;
model.Q = diag([1e-18^2,1e-18^2,1e-18^2,1e-18^2,1e-18^2,1e-18^2]);
model.Q(1:3,1:3)=model.Q(1:3,1:3)*constants.trueX2normX^2;
model.Q(4:6,4:6)=model.Q(4:6,4:6)*constants.trueV2normV^2;
model.sQ = sqrtm(model.Q);
model.fstates = {'a','e','i','\Omega','\omega','M'};
% model.fstates = {'x','y','z','vx','vy','vz'};

model.h=@(x)radmodel_coe(x,1);
model.hn=2;
model.hvec=@(x)sensmodel2vec(x,model.h,model.hn);
model.R=diag([(1*pi/180)^2,(1*pi/180)^2]);
model.hstates = {'azi','el'};
% model.R=diag([(0.1/constants.Re)^2]);
model.z_pdf =  @(z,x)mvnpdf(z,model.h(x),model.R);



%% generate truth

x0=[10000,10,3000,-0.1,7.4,2]';
OE = cart2orbelem(x0,constants.mu)
% [ r, v, Ehat ] = FnG(0, time.dt, x0(1:3), x0(4:6), 1);
% [ r1, v1, Ehat ] = FnG(0, time.dt, r, -v, 1);

P0=diag([1^2,1^2,1^2,0.001^2,0.001^2,0.001^2]);

x0(1:3)=x0(1:3)*constants.trueX2normX;
x0(4:6)=x0(4:6)*constants.trueV2normV;
P0(1:3,1:3)=P0(1:3,1:3)*constants.trueX2normX^2;
P0(4:6,4:6)=P0(4:6,4:6)*constants.trueV2normV^2;

X=mvnrnd(x0(:)',P0,1000);
nn=size(X,1);
Xcoe=zeros(size(X));
for i=1:size(X,1)
    Xcoe(i,:) = rv2elm_a(X(i,:),1);
end
[x0coe,P0coe]=MeanCov(Xcoe,1/nn*ones(nn,1));

x0=x0coe;
P0=P0coe;

Xtruth = zeros(time.Ntsteps,model.fn);
Xtruth(1,:)=x0;
for k=2:time.Ntsteps
    Xtruth(k,:)=model.f(time.dt,time.Tvec(k-1),Xtruth(k-1,:));
end
Xplot = zeros(time.Ntsteps,model.fn);
Xplot(1,:)=x0;
for k=2:length(time.TvecPlots)
    Xplot(k,:)=model.f(time.dtplot,time.TvecPlots(k-1),Xplot(k-1,:));
end

figure
plot(Xtruth(:,1),Xtruth(:,2),'ro')
hold on
plot(Xplot(:,1),Xplot(:,2),'b')
axis equal
axis square

% plotting the propagatin of MC
Nmc=5000;
XMC=zeros(Nmc,model.fn,time.Ntsteps);
% XMCcoe=zeros(Nmc,model.fn,time.Ntsteps);
XMC(:,:,1)=mvnrnd(x0',P0,Nmc);
% for i=1:Nmc
% %     [p,a,ecc,incl,omega,argp,nu,m,arglat,truelon,lonper ] = rv2coe(XMC(i,1:3,1),XMC(i,4:6,1), 1);
% %     XMCcoe(i,:,1)=[p,a,ecc,incl,omega,argp];
% end
for i=1:Nmc
    i
    for k=2:time.Ntsteps
        XMC(i,:,k)=model.f(time.dt,time.Tvec(k-1),XMC(i,:,k-1));
%         [p,a,ecc,incl,omega,argp,nu,m,arglat,truelon,lonper ] = rv2coe(XMC(i,1:3,k),XMC(i,4:6,k), 1);
%         XMCcoe(i,:,k)=[p,a,ecc,incl,omega,argp];
    end
end
pMC = mvnpdf(XMC(:,:,1),x0',P0);

%%
% close all
% kk=5;
% states=[5,6];
% figure
% plot(XMCcoe(:,states(1),kk),XMCcoe(:,states(2),kk),'bo')
% figure
% plot(XMC(:,states(1),kk),XMC(:,states(2),kk),'bo')
%%
close all
if 0
    for k=1:time.Ntsteps
        figure(1)
        clf
        subplot(1,2,1)
        plot3(XMC(:,1,k),XMC(:,2,k),log(pMC),'ro')
        hold on
        plot3(Xplot(:,1),Xplot(:,2),Xplot(:,3),'b-')
        title(['k = ',num2str(k)])
        axis equal
        axis square
        
        subplot(1,2,2)
        XXX = zeros(Nmc,model.fn);
        XXXcoe = zeros(Nmc,model.fn);
        XXX = XMC(:,:,k);
        XXXcoe = XMCcoe(:,:,k);
        [m,P] = MeanCov(XXX,ones(Nmc,1)/Nmc);
        A=sqrtm(inv(P));
        for i=1:Nmc
            XXX(i,:) = A*(XXX(i,:)'-m(:));
        end
        plot3(XXX(:,1),XXX(:,2),log(pMC),'ro')
        title(['k = ',num2str(k)])
        axis equal
        axis square
        view([-66,90])
        keyboard
    end
end
%% comparing with UKF and particle filter
xf0=mvnrnd(x0(:)',P0);
% xf0 = x0;
Pf0 = P0;

xfquad = xf0;
Pfquad = P0;

Npf = 15000; %paricle filter points


GMMinitial=getInitialGMM(xf0(:),0.9^2*Pf0,4);
model.gmmmethod='ut';

% generate points on contours for characterisitc solutions
model.pointscalefactor = 0.9;
model.pointGenerator = @(mx,Px)hybridpoints(mx,Px,3,0.2,5,0.6);
% model.pointGenerator = @(mx,Px)mvnrnd(mx(:)',Px,20000);


% [Xut,wut] = UT_sigmapoints(zeros(model.fn,1),eye(model.fn),2);
% 
% [X,w] = model.pointGenerator(zeros(model.fn,1),eye(model.fn));


pf.no_particles =10000;
pf.no_bins = 100;               % when computing the histogram

pf.resample = true;             % resampling after measurement update
pf.neff = 0.8;         
pf.regularizepf = 0;


[X,w] = GH_points(zeros(model.fn,1),0.9^2*eye(model.fn),4);
% [X2,w] = GH_points(zeros(model.fn,1),2^2*eye(model.fn),4);
% X2 = 6*sphere6Dm(6); %3906
% % [X,w] = mvnrnd(zeros(4,1),0.5*eye(4),500);
% X=[X1;X2];

A=sqrtm(Pf0);
for i=1:size(X,1)
    X(i,:) = A*X(i,:)'+xf0(:);
end
probs = mvnpdf(X,xf0(:)',Pf0);


figure
plot(X(:,1),X(:,2),'r+')


Xinitial = X;
probsinitial = probs;

model.quadfunc=@(x,P)UT_sigmapoints(x,P,2);

[Xquad_initial,wquad_initial]=model.quadfunc(xf0(:),Pf0);
probs_quad = mvnpdf(Xquad_initial,xf0(:)',Pf0);


Npt=size(Xinitial,1);
Xpttraj=zeros(Npt,model.fn,time.Ntsteps);
Xpttraj(:,:,1)=Xinitial;
for i=1:Npt
    i
    for k=2:time.Ntsteps
        Xpttraj(i,:,k)=model.f(time.dt,time.Tvec(k-1),Xpttraj(i,:,k-1));
    end
end
% p=4;
% Xr = zeros(Npt,model.fn);
% for ii=1:Npt
%     Xr(ii,:) = Xpttraj(ii,:,p);
% end
% 
% figure
% plot3(Xr(:,1),Xr(:,2),Xr(:,3),'bo')





%% run filter
close all
clc

X_pf=mvnrnd(xf0(:)',P0,pf.no_particles)';
w_pf = 1/pf.no_particles * ones(pf.no_particles,1);

X = Xinitial;
probs = probsinitial;

Xquad=Xquad_initial;
wquad=wquad_initial;

xfquad = xf0;
Pfquad = P0;


pdftraj = cell(time.Ntsteps,2); % 1 is prior and 2 is post
maxBnd= ones(1,model.fn);
% mean_transform0I2true, cov_transform0I2true, maxbnd_transform0I211cube,normpdf, normpdf_expoenet

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

pdftraj{1,1} = pdfnorm;

meas_freq_steps = 1;


plotfolder='simulations/SAT6Dsim1_meas_bringback';
mkdir(plotfolder)

savePriorProps.plotfolder=[plotfolder,'/prior'];
mkdir(savePriorProps.plotfolder)
savePriorProps.saveit=1;

savePostProps.plotfolder=[plotfolder,'/post'];
mkdir(savePostProps.plotfolder)
savePostProps.saveit=1;

EstMOCfilter_mu =   zeros(time.Ntsteps,model.fn);
EstMOCfilter_P =    zeros(time.Ntsteps,model.fn^2);

EstQuadfilter_mu =   zeros(time.Ntsteps,model.fn);
EstQuadfilter_P =    zeros(time.Ntsteps,model.fn^2);

[mX,PX]=MeanCov(X,probs/sum(probs));
EstMOCfilter_mu(1,:) = mX;
EstMOCfilter_P(1,:) = reshape(PX,1,model.fn^2);

EstQuadfilter_mu(1,:) = xfquad;
EstQuadfilter_P(1,:) = reshape(Pfquad,1,model.fn^2);

priorfullnormpdf=pdfnorm;

histXprior{1,1}=X;
histXprior{1,2}=probs;
histXprior{1,3}=priorfullnormpdf;
histXprior{1,4}=priorfullnormpdf;

histXprior{1,5}=xfquad;
histXprior{1,6}=Pfquad;

histXprior{1,7}=xfquad;
histXprior{1,8}=Pfquad;

histXpost{k,9}=X_pf;
histXpost{k,10}=w_pf;

histXpost{k,11}=X_pf;
histXpost{k,12}=w_pf;

for k=2:time.Ntsteps
%     close all
    
    k
    Xmctest = zeros(size(XMC,1),model.fn);
    for ii=1:size(XMC,1)
        Xmctest(ii,:) = XMC(ii,:,k);
    end
    Xmctest=[];
    
    disp([' k = ',num2str(k)])
    
    [X,probs]=propagate_character(X,probs,time.dt,time.Tvec(k),model);
    [Xquad,wquad]=propagate_character(Xquad,wquad,time.dt,time.Tvec(k),model);
    [xfquad,Pfquad]=QuadProp(xfquad,Pfquad,time.dt,time.Tvec(k),model,'ut');

%     GMM=prior_prop_GMM(GMM,time.dt,time.Tvec(k),model);
    
     X_pf = pftimeupdate_custommoc(X_pf, model,time.Tvec(k),k,time.dt);
     w_pf=w_pf;
 
 
    [mX,PX]=MeanCov(X,probs/sum(probs));
    

    disp(['cond = ',num2str(cond(PX))])

%     [priorfullnormpdf,X,probs] = get_interp_pdf_0I_6D_gobackget_coe(model,priorfullnormpdf,X,probs,4,k,time.Tvec(k),time.dt,Xmctest,Xtruth(k,:));
%     pdftraj{k,1} = priorfullnormpdf;


%     fullnormpdf=get_interp_pdf_0I_6Dsat(X,probs,mX,PX,4,3,k,Xmctest,Xtruth(k,:)); %Xtruth(k,:)

%     priorpdfnorm=fullnormpdf;
    
%     plotpdfs_prior_6D([1,2],k,priorfullnormpdf,X,probs,xfquad,Pfquad,Xmctest,Xtruth(k,:),model,savePriorProps)
%     plotpdfs_prior_6D([1,6],k,priorfullnormpdf,X,probs,xfquad,Pfquad,Xmctest,Xtruth(k,:),model,savePriorProps)
%     plotpdfs_prior_6D([2,6],k,priorfullnormpdf,X,probs,xfquad,Pfquad,Xmctest,Xtruth(k,:),model,savePriorProps)
%     plotpdfs_prior_6D([3,6],k,priorfullnormpdf,X,probs,xfquad,Pfquad,Xmctest,Xtruth(k,:),model,savePriorProps)
%     plotpdfs_prior_6D([4,6],k,priorfullnormpdf,X,probs,xfquad,Pfquad,Xmctest,Xtruth(k,:),model,savePriorProps)
%     plotpdfs_prior_6D([5,6],k,priorfullnormpdf,X,probs,xfquad,Pfquad,Xmctest,Xtruth(k,:),model,savePriorProps)

%     keyboard
    
%     histXprior{k,1}=X;
%     histXprior{k,2}=probs;
%     histXprior{k,3}=priorfullnormpdf;
%     
%     histXprior{k,4}=xfquad;
%     histXprior{k,5}=Pfquad;
    
%    keyboard
    pause(1)
    
    
    
    
    % do measurement update
    if k>=2
        if rem(k,meas_freq_steps)==0
            disp("doing meas update")
            zk = model.h(Xtruth(k,:)')+sqrtm(model.R)*randn(model.hn,1);
            zk
            
            % %%%%%%%%%%%%%%%%% MEAS UPDATE %%%%%%%%%%%%%%%%%%%%%%
%             [X,probs,fullnormpdf] =       MeasUpdt_character_modf(X,probs,3,k,zk,Xtruth(k,:),model,Xmctest,11);
%             [X,probs,fullnormpdf] = MeasUpdt_character_modf_6Dsat(X,probs,priorpdfnorm,4,k,zk,Xtruth(k,:),model,Xmctest,5);
             [priorfullnormpdf,postfullnormpdf,X,probs] = get_interp_pdf_0I_6D_gobackget_measupt_coe(zk,model,priorfullnormpdf,X,probs,4,k,time.Tvec(k),time.dt,Xmctest,Xtruth(k,:));
             
             
%             [X,probs,fullnormpdf]=MeasUpdt_character_modf(fullnormpdf,X,probs,4,k,zk,Xtruth,model,Xmctest);
            xfquadprior=xfquad;
            Pfquadprior=Pfquad;
            [xfquad,Pfquad]=QuadMeasUpdt(xfquad,Pfquad,zk,time.dt,time.Tvec(k),model,'ut');
            
%             GMM=post_GMM(GMM,zk,time.dt,time.Tvec(k),model);
            
            X_pfprior=X_pf;
            w_pfprior=w_pf;
            [X_pf, w_pf] = pfmeasupdate_custommoc(X_pf, w_pf, model,time.Tvec(k),k,time.dt, pf, zk);
            
%             keyboard
%             plotpdfs_post_6D([1,2],k,priorfullnormpdf,postfullnormpdf,model,X,probs,xfquad,Pfquad,xfquadprior,Pfquadprior,X_pf, w_pf,X_pfprior,w_pfprior,Xmctest,Xtruth(k,:),zk,savePostProps,0)
%             pause(1)
            plotpdfs_post_6D([1,6],k,priorfullnormpdf,postfullnormpdf,model,X,probs,xfquad,Pfquad,xfquadprior,Pfquadprior,X_pf, w_pf,X_pfprior,w_pfprior,Xmctest,Xtruth(k,:),zk,savePostProps,0)
            pause(1)
%             plotpdfs_post_6D([2,6],k,priorfullnormpdf,postfullnormpdf,model,X,probs,xfquad,Pfquad,xfquadprior,Pfquadprior,X_pf, w_pf,X_pfprior,w_pfprior,Xmctest,Xtruth(k,:),zk,savePostProps,0)
%             plotpdfs_post_6D([3,6],k,priorfullnormpdf,postfullnormpdf,model,X,probs,xfquad,Pfquad,xfquadprior,Pfquadprior,X_pf, w_pf,X_pfprior,w_pfprior,Xmctest,Xtruth(k,:),zk,savePostProps,0)
%             plotpdfs_post_6D([4,6],k,priorfullnormpdf,postfullnormpdf,model,X,probs,xfquad,Pfquad,xfquadprior,Pfquadprior,X_pf, w_pf,X_pfprior,w_pfprior,Xmctest,Xtruth(k,:),zk,savePostProps,0)
%             plotpdfs_post_6D([5,6],k,priorfullnormpdf,postfullnormpdf,model,X,probs,xfquad,Pfquad,xfquadprior,Pfquadprior,X_pf, w_pf,X_pfprior,w_pfprior,Xmctest,Xtruth(k,:),zk,savePostProps,0)
            

            % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            histXpost{k,1}=X;
            histXpost{k,2}=probs;
            histXpost{k,3}=priorfullnormpdf;
            histXpost{k,4}=postfullnormpdf;
            
            histXpost{k,5}=xfquad;
            histXpost{k,6}=Pfquad;
            
            histXpost{k,7}=xfquadprior;
            histXpost{k,8}=Pfquadprior;
            
            histXpost{k,9}=X_pf;
            histXpost{k,10}=w_pf;
            
            histXpost{k,11}=X_pfprior;
            histXpost{k,12}=w_pfprior;
            
            histXpost{k,13}=zk;
            
            priorfullnormpdf = postfullnormpdf;
        end
    end

%     keyboard

%     pause(1)
    
    
    [mestX,PestX]=MeanCov(X,probs/sum(probs));

    
    EstMOCfilter_mu(k,:) = mestX;
    EstMOCfilter_P(k,:) = reshape(PestX,1,model.fn^2);

    EstQuadfilter_mu(k,:) = xfquad;
    EstQuadfilter_P(k,:) = reshape(Pfquad,1,model.fn^2);

end

save([plotfolder,'/sim1.mat'])
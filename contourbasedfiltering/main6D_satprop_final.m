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
time.dt=0.5*60*60 *constants.trueT2normT;
time.dtplot=0.1*60*60 *constants.trueT2normT;

time.Tvec=time.t0:time.dt:time.tf;
time.Ntsteps=length(time.Tvec);

time.TvecPlots=time.t0: time.dtplot :time.tf;

%% models

model.f=@(dt,tk,xk)processmodel_2body(dt,1,tk,xk);
model.fback=@(dt,tk,xk)processmodel_2body_backward(dt,1,tk,xk);

% [a,e,i,om,Om,M]
model.fn=6;
model.Q = diag([0.00001^2,0.00001^2,0.00001^2,0.0000001^2,0.0000001^2,0.0000001^2]);
model.Q(1:3,1:3)=model.Q(1:3,1:3)*constants.trueX2normX^2;
model.Q(4:6,4:6)=model.Q(4:6,4:6)*constants.trueV2normV^2;
model.fstates = {'a','e','i','Om','om','M'};
% model.fstates = {'x','y','z','vx','vy','vz'};

model.h=@(x)radmodel(xk,1);
model.hn=2;
model.hvec=@(x)sensmodel2vec(x,model.h,model.hn);
model.R=diag([(0.1*pi/180)^2,(0.1*pi/180)^2]);
model.hstates = {'azi','el'};
% model.R=diag([(0.1/constants.Re)^2]);
model.z_pdf =  @(z,x)mvnpdf(z,model.h(x),model.R);



%% generate truth

x0=[10000,10,3000,-1,7.5,3.5]';
OE = cart2orbelem(x0,constants.mu)
% [ r, v, Ehat ] = FnG(0, time.dt, x0(1:3), x0(4:6), 1);
% [ r1, v1, Ehat ] = FnG(0, time.dt, r, -v, 1);

P0=diag([0.01^2,0.01^2,0.01^2,0.00001^2,0.00001^2,0.00001^2]);

x0(1:3)=x0(1:3)*constants.trueX2normX;
x0(4:6)=x0(4:6)*constants.trueV2normV;
P0(1:3,1:3)=P0(1:3,1:3)*constants.trueX2normX^2;
P0(4:6,4:6)=P0(4:6,4:6)*constants.trueV2normV^2;

% x0=[-6.01508104006094,-4.59324087112725,-4.12840886105591,0.333796709270969, ...
% 0.00757260259438741,0.103806185473391]';

X=mvnrnd(x0(:)',P0,1000);
nn=size(X,1);


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



% plotting the propagatin of MC
Nmc=10000;
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
    end
end
pMC = mvnpdf(XMC(:,:,1),x0',P0);



%% 

% Xmctest=zeros(size(XMC,1),6);
% for k=1:time.Ntsteps
%     for ii=1:size(XMC,1)
%         Xmctest(ii,:) = XMC(ii,:,k);
%     end
%     [m,p]=MeanCov(Xmctest,ones(Nmc,1)/Nmc);
%     A=sqrtm(inv(p));
% 
%     for ii=1:size(XMC,1)
%         Xmctest(ii,:) = A*(Xmctest(ii,:)-m(:)')';
%     end
%     
%     figure(1)
%     plot(Xmctest(:,1)*constants.normX2trueX,Xmctest(:,2)*constants.normX2trueX,'r.')
%     title([num2str(k),' of ',num2str(time.Ntsteps)])
%     pause(0.5)
%     k
% end
%%
close all
if 0
    for k=1:time.Ntsteps
        figure(2)
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
        XXX = XMC(:,:,k);
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
% xf0=mvnrnd(x0(:)',P0);
xf0 = x0;
Pf0 = P0;

xfquad = xf0;
Pfquad = P0;

Npf = 5000; %paricle filter points

GMMinitial=getInitialGMM(xf0(:),0.7^2*Pf0,7);
model.gmmmethod='ut';

% generate points on contours for characterisitc solutions

model.pointGenerator = @(mx,Px)GH_points(mx,Px,11);


% model.pointGenerator = @(mx,Px)mvnrnd(mx,Px,200);
model.pointscalefactor = 0.4;
[X,w] = model.pointGenerator(zeros(model.fn,1),model.pointscalefactor^2*eye(model.fn));
% [X,w] = mvnrnd(zeros(model.fn,1),1^2*eye(model.fn),100);

% [X,w] = mvnrnd(zeros(4,1),0.5*eye(4),500);

A=sqrt(Pf0);
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





%% run filter
close all

X = Xinitial;
probs = probsinitial;

Xquad=Xquad_initial;
wquad=wquad_initial;

xfquad = xf0;
Pfquad = P0;

GMM = GMMinitial;

pdftraj = cell(time.Ntsteps,2); % 1 is prior and 2 is post
maxBnd= ones(1,model.fn);
% mean_transform0I2true, cov_transform0I2true, maxbnd_transform0I211cube,normpdf, normpdf_expoenet

[X1,~] = GH_points(xf0(:),Pf0,4);
probsk1 = mvnpdf(X1,xf0(:)',Pf0);
dsX = DataSet(X1,probsk1,'TrueState');
dsX.AddMeanCov_to_OI_Trasform(xf0(:),Pf0);
transForms = dsX.GetTrasnformers();

pdfnorm.func=@(x)mvnpdf(x,zeros(1,model.fn),eye(model.fn));
pdfnorm.expfunc=@(x)normpdfexp(x,zeros(1,model.fn),eye(model.fn));
pdfnorm.transForms = dsX.GetTrasnformers();

pdfnorm.info = 'true-0I';
pdfnorm.pdftype = 'ExpPdf';

% pdfnorm.GMMHull = GMMHull;
pdfnorm.LB = -5*ones(1,model.fn);
pdfnorm.UB = 5*ones(1,model.fn);

pdftraj{1,1} = pdfnorm;

meas_freq_steps = 1000000;

histXprior=cell(time.Ntsteps,5);
histXpost=cell(time.Ntsteps,5);

histXprior{1,1} = X;
histXprior{1,2} = probs;

histXpost{1,1} = X;
histXpost{1,2} = probs;
teststeps = [];

plotfolder='simulations/sat6D_prop_backprop_L1';
mkdir(plotfolder)

savePriorProps.plotfolder=[plotfolder,'/prior'];%             [X,probs,fullnormpdf]=MeasUpdt_character_modf_duff(X,probs,priorpdfnorm,4,k,zk,Xtruth(k,:),model,Xmctest,11);
% %             [X,probs,fullnormpdf]=MeasUpdt_character_modf(fullnormpdf,X,probs,4,k,zk,Xtruth,model,Xmctest);
%             [xfquad,Pfquad]=QuadMeasUpdt(xfquad,Pfquad,zk,time.dt,time.Tvec(k),model,'ut');
%             
%             GMM=post_GMM(GMM,zk,time.dt,time.Tvec(k),model);
%             
%             plotpdfs_post_2D(k,priorpdfnorm,fullnormpdf,model,Xtruplot,X,probs,xfquad,Pfquad,GMM,Xmctest,Xtruth(k,:),zk,savePostProps)

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

for k=2:time.Ntsteps
%     close all
    
    k
    Xmctest = zeros(size(XMC,1),model.fn);
    for ii=1:size(XMC,1)
        Xmctest(ii,:) = XMC(ii,:,k);
    end
%     Xmctest=[];
    disp([' k = ',num2str(k)])
    
    [X,probs]=propagate_character(X,probs,time.dt,time.Tvec(k),model);
%     [Xquad,wquad]=propagate_character(Xquad,wquad,time.dt,time.Tvec(k),model);
    [xfquad,Pfquad]=QuadProp(xfquad,Pfquad,time.dt,time.Tvec(k),model,'ut');

%     GMM=prior_prop_GMM(GMM,time.dt,time.Tvec(k),model);
    

    
    

    disp(['cond = ',num2str(cond(PX))])

    [mX,PX]=MeanCov(X,probs/sum(probs));
%     fullnormpdf=get_interp_pdf_0I_duff(X,probs,mX,PX,4,3,k,Xmctest,Xtruth(k,:)); %Xtruth(k,:)
    priorfullnormpdf=get_interp_pdf_0I_6D_gobackget_L1(model,priorfullnormpdf,X,probs,15,k,time.Tvec(k),time.dt,Xmctest,Xtruth(k,:));
%     priorfullnormpdf=get_interp_pdf_0I_2D_gobackget(model,priorfullnormpdf,X,probs,4,k,time.Tvec(k),time.dt,Xmctest,Xtruth(k,:));
    %     fullnormpdf=get_interp_pdf_0I_2D(X,probs,mX,PX,4,k,[],Xtruth(k,:),plotsconf); %Xtruth(k,:)

    
    plotpdfs_prior_2D(k,priorfullnormpdf,X,probs,xfquad,Pfquad,GMM,Xmctest,Xtruth(k,:),savePriorProps)
    
        keyboard
    
%     figure(12)
%     hold on
%     [mG,PG]=GMMcell2array(GMM);
%     plot(mG(:,1),mG(:,2),'b*','MarkerSize',10)
%     hold off
    
%     X=histXprior{k,1};
%     probs=histXprior{k,2};
%     fullnormpdf=histXprior{k,3};
    
    
    histXprior{k,1}=X;
    histXprior{k,2}=probs;
    histXprior{k,3}=priorfullnormpdf;
    
    histXprior{k,4}=xfquad;
    histXprior{k,5}=Pfquad;
    
   
    pause(1)
    
    
%     X=histXprior{k,1};
%     probs=histXprior{k,2};
    
    % do measurement update
    if k>=2
        if rem(k,meas_freq_steps)==0
            disp("doing meas update")
            zk = model.h(Xtruth(k,:)')+sqrtm(model.R)*randn(model.hn,1);
            zk
            
            % %%%%%%%%%%%%%%%%% MEAS UPDATE %%%%%%%%%%%%%%%%%%%%%%
%             [X,probs,fullnormpdf]=MeasUpdt_character_modf_duff(X,probs,priorpdfnorm,4,k,zk,Xtruth(k,:),model,Xmctest,11);
% %             [X,probs,fullnormpdf]=MeasUpdt_character_modf(fullnormpdf,X,probs,4,k,zk,Xtruth,model,Xmctest);
%             [xfquad,Pfquad]=QuadMeasUpdt(xfquad,Pfquad,zk,time.dt,time.Tvec(k),model,'ut');
%             
%             GMM=post_GMM(GMM,zk,time.dt,time.Tvec(k),model);
%             
%             plotpdfs_post_2D(k,priorpdfnorm,fullnormpdf,model,Xtruplot,X,probs,xfquad,Pfquad,GMM,Xmctest,Xtruth(k,:),zk,savePostProps)
%             plotpdfs_post_2D(Tk,pdfnormprior,pdfnorm,model,XtruthAll,X,probs,mquadf,Pquadf,GMM,Xmc,Xtruth,zk,saveprops)
            
            % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            histXpost{k,1}=X;
            histXpost{k,2}=probs;
            histXpost{k,3}=fullnormpdf;

            histXpost{k,4}=xfquad;
            histXpost{k,5}=Pfquad;
            histXpost{k,6}=zk;
            
            keyboard
            
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

%%
save([plotfolder,'/sim1.mat'])


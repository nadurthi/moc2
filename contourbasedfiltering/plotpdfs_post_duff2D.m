% function plotpdfs_post_duff2D(Tk,pdfnormprior,pdfnorm,model,XtruthplotAll,X,probs,mquadf,Pquadf,mquadfprior,Pquadfprior,X_pf, w_pf,X_pfprior,w_pfprior,GMM,Xmc,Xtruth,zk,saveprops)
function plotpdfs_post_duff2D(Tk,pdfnormprior,pdfnorm,model,XtruthplotAll,X,probs,mquadfprior,Pquadfprior,mquadf,Pquadf,GMM,Xmc,Xtruth,zk,saveprops)
nametag='post';
states=[1,2];
states2keep=states;
Xn = pdfnorm.transForms.trueX2normX(X);
pn = pdfnorm.transForms.trueprob2normprob(probs);

dim = model.fn;

if isempty(Xmc)==0
    Xmcnorm = pdfnorm.transForms.trueX2normX(Xmc) ;
end
if isempty(Xtruth)==0
    Xtruthnorm = pdfnorm.transForms.trueX2normX(Xtruth) ;
end
% if isempty(mquadf)==0
%     [x,w] = UT_sigmapoints(mquadf,Pquadf,2);
%     x = pdfnorm.transForms.trueX2normX(x) ;
%     [mquadfnorm,Pquadfnorm]=MeanCov(x,w);
% end
% %% resample pf to make equal weightd
% Npf=length(w_pf);
% I = myresampling(w_pf);
% I = round(I);
% X_old=X_pf;
% for j = 1 : Npf
%     X_pf(:,j) = X_old(:,I(j));
% end
% w_pf = ones(1,Npf)/Npf;
% 
% I = myresampling(w_pfprior);
% I = round(I);
% X_old=X_pfprior;
% for j = 1 : Npf
%     X_pfprior(:,j) = X_old(:,I(j));
% end
% w_pfprior = ones(1,Npf)/Npf;
% 
% X_pfnorm=zeros(size(X_pf));
% X_pfpriornorm=zeros(size(X_pf));
% for i=1:Npf
%  X_pfnorm(:,i)=pdfnorm.transForms.trueX2normX(X_pf(:,i)');
%  X_pfpriornorm(:,i)=pdfnorm.transForms.trueX2normX(X_pfprior(:,i)');
% end
% bwpost=zeros(1,dim);
% bwprior=zeros(1,dim);
% 
% Npf = length(w_pf);
% for i=1:dim
%     s=std(X_pfnorm(:,i));
%     bwpost(i) = s*(4/((dim+2)*Npf))^(1/(dim+4));
%     
%     s=std(X_pfpriornorm(:,i));
%     bwprior(i) = s*(4/((dim+2)*Npf))^(1/(dim+4));
%     
% end



%%
[Xxnormbase,Xynormbase]=meshgrid(linspace(-1.2,1.2,50),linspace(-1.2,1.2,50) );
Xx_priortrue = zeros(size(Xxnormbase));
Xy_priortrue = zeros(size(Xxnormbase));
Xx_posttrue = zeros(size(Xxnormbase));
Xy_posttrue = zeros(size(Xxnormbase));

Xx_prior2post_norm = zeros(size(Xxnormbase));
Xy_prior2post_norm = zeros(size(Xxnormbase));

pdfprior2postprobs_norm = zeros(size(Xxnormbase));
QuadFilprior2postprobs_norm = zeros(size(Xxnormbase));


pdfpriorprobs_true = zeros(size(Xxnormbase));
QuadFilpriorprobs_true = zeros(size(Xxnormbase));

pdfpostprobs_true = zeros(size(Xxnormbase));
QuadFilpostprobs_true = zeros(size(Xxnormbase));
QuadFilpostprobs_norm = zeros(size(Xxnormbase));

pdfpostprobs_norm = zeros(size(Xxnormbase));
pdfpriorprobs_norm = zeros(size(Xxnormbase));

for i=1:size(Xxnormbase,1)
    for j=1:size(Xxnormbase,2)
        norm_pt = [Xxnormbase(i,j),Xynormbase(i,j)];
        posttrue_pt = pdfnorm.transForms.normX2trueX(norm_pt);
        priortrue_pt = pdfnormprior.transForms.normX2trueX(norm_pt);
        prior2post_norm_pt = pdfnorm.transForms.trueX2normX(priortrue_pt);
        
        Xx_posttrue(i,j) = posttrue_pt(1);
        Xy_posttrue(i,j) = posttrue_pt(2);
        
        Xx_priortrue(i,j) = priortrue_pt(1);
        Xy_priortrue(i,j) = priortrue_pt(2);
        
        Xx_prior2post_norm(i,j) = prior2post_norm_pt(1);
        Xy_prior2post_norm(i,j) = prior2post_norm_pt(2);
        
        pdfpostprobs_norm(i,j) = pdfnorm.func(norm_pt);
        pdfpriorprobs_norm(i,j) = pdfnormprior.func(norm_pt);
        
        pdfpostprobs_true(i,j) = pdfnorm.transForms.normprob2trueprob( pdfpostprobs_norm(i,j) );
        pdfpriorprobs_true(i,j) = pdfnormprior.transForms.normprob2trueprob( pdfpriorprobs_norm(i,j) );
        
        pdfprior2postprobs_norm(i,j) = pdfnorm.transForms.trueprob2normprob( pdfpriorprobs_true(i,j) );
        
        QuadFilpostprobs_true(i,j) = mvnpdf(posttrue_pt,mquadf(:)',Pquadf);
        QuadFilpostprobs_norm(i,j) = pdfnorm.transForms.trueprob2normprob( QuadFilpostprobs_true(i,j) );
        
        QuadFilpriorprobs_true(i,j) = mvnpdf(priortrue_pt,mquadfprior(:)',Pquadfprior);
        QuadFilprior2postprobs_norm(i,j) = pdfnorm.transForms.trueprob2normprob( QuadFilpriorprobs_true(i,j) );
        

    end
end



%% Get the z-probs p(z)
% [x,w] = UT_sigmapoints(mquadf,Pquadf,2);
dopz = 0;
if dopz == 1
    Z=zeros(size(X,1),model.hn);
    for i=1:length(probs)
        Z(i,:) = model.h(X(i,:));
    end
    [mz,Pz] = MeanCov(Z,probs/sum(probs));
    Pz=Pz+model.R;
    
    sqPz = sqrtm(Pz);
    % keyboard
    if model.hn==1
        Zx = linspace(mz-3*sqPz,mz+3*sqPz,100);
        %     Pzprobs_true = zeros(size(Zx));
        %     for i=1:1:length(Zx)
        %         I = integratorFuncTrueX_usingpdfnorm(pdfnormprior,@(x)mvnpdf(repmat(Zx(i),size(x,1),1),model.h(x),model.R),'RegTreeBoxIntegrator');
        %         Pzprobs_true(i) = I;
        %     end
        Pzprobs_true = zeros(size(Zx));
        %     J=cell(size(Zx,1),1);
        LB=-1.1*ones(1,dim);
        UB=1.1*ones(1,dim);
        voll = prod(UB-LB);
        [XX,WW]=GLgn_pts(LB,UB,9);
        %         keyboard
        for i=1:1:length(Zx)
            i
            
            zz=Zx(i);
            pnorm = pdfnormprior.func(XX);
            Xtr=pdfnormprior.transForms.normX2trueX(XX);
            pz=mvnpdf(repmat(zz,size(Xtr,1),1),model.hvec(Xtr),model.R);
            Pzprobs_true(i)=sum(WW.*pnorm.*pz);
            
            
        end
        
        Pzk = interp1(Zx,Pzprobs_true,zk);
        figure(93)
        plot(Zx,Pzprobs_true,'b',zk,0,'k^',[zk,zk],[0,Pzk],'k--','linewidth',2,'MarkerSize',6)
        xlabel('z')
        ylabel('p(z)')
        hold off
        if saveprops.saveit==1
            saveas(gcf,[saveprops.plotfolder,'/PZtrueSurf_',nametag,'_',num2str(Tk)],'png')
            saveas(gcf,[saveprops.plotfolder,'/PZtrueSurf_',nametag,'_',num2str(Tk)],'fig')
        end
        
    end
    if model.hn==2
        [Zx,Zy] = meshgrid(linspace(mz(1)-3*sqPz(1,1),mz(1)+3*sqPz(1,1),25),linspace(mz(2)-3*sqPz(2,2),mz(2)+3*sqPz(2,2),25));
        Pzprobs_true = zeros(size(Zx));
        J=cell(size(Zx,1),1);
        for i=1:1:size(Zx,1)
            J{i}=zeros(1,size(Zx,2));
            for j=1:1:size(Zx,2)
                zz=[Zx(i,j),Zy(i,j)];
                I = integratorFuncTrueX_usingpdfnorm(pdfnormprior,@(x)mvnpdf(repmat(zz,size(x,1),1),model.h(x),model.R),'RegTreeBoxIntegrator');
                %             Pzprobs_true(i,j) = I;
                J{i}(j) = I;
            end
        end
        for i=1:1:size(Zx,1)
            for j=1:1:size(Zx,2)
                Pzprobs_true(i,j) = J{i}(j);
            end
        end
        figure(93)
        contour(Zx,Zy,Pzprobs_true,15)
        hold on
        plot(zk(1),zk(2),'k*','linewidth',2,'MarkerSize',6)
        title('Pz true contour')
        xlabel('z_1')
        ylabel('z_2')
        hold off
        if saveprops.saveit==1
            saveas(gcf,[saveprops.plotfolder,'/PZtrueContour_',nametag,'_',num2str(Tk)],'png')
            saveas(gcf,[saveprops.plotfolder,'/PZtrueContour_',nametag,'_',num2str(Tk)],'fig')
        end
        
        
        figure(94)
        surf(Zx,Zy,Pzprobs_true,'FaceColor','green','EdgeColor','none','FaceAlpha',0.7);
        camlight right; lighting phong
        alpha 0.4
        hold on
        plot(zk(1),zk(2),'k*','linewidth',2,'MarkerSize',6)
        title('Pz true surf')
        xlabel('z_1')
        ylabel('z_2')
        hold off
        if saveprops.saveit==1
            saveas(gcf,[saveprops.plotfolder,'/PZtrueSurf_',nametag,'_',num2str(Tk)],'png')
            saveas(gcf,[saveprops.plotfolder,'/PZtrueSurf_',nametag,'_',num2str(Tk)],'fig')
        end
        
    end
end

%% contour levels
mm = max(max(pdfpostprobs_norm));
clevelsnorm = exp( linspace(-5,log(mm),15) );

mm = max(max(pdfpostprobs_true));
clevelstrue = exp( linspace(-5,log(mm),15) );

axislims_norm = [-1.2,1.2,-1.2,1.2];
%% postpdf norm space reg contours
% keyboard
figure(41)
contour(Xxnormbase,Xynormbase,pdfprior2postprobs_norm,clevelsnorm,'LineColor','b')
grid on
box off
hold on
contour(Xxnormbase,Xynormbase,pdfpostprobs_norm,clevelsnorm,'LineColor','r')

if isempty(Xmc)==0
    plot(Xmcnorm(:,1),Xmcnorm(:,2),'r.')
end
if isempty(Xtruth)==0
    plot(Xtruthnorm(:,1),Xtruthnorm(:,2),'k*','linewidth',2)
end

% title(['time step = ',num2str(Tk),' cond = ',num2str(cond(Pquad))])
xlabel('x_1')
ylabel('x_2')
title('Prior and Post Reg Contour norm')
axis equal
axis square
hold off
axis(axislims_norm)

if saveprops.saveit==1
    saveas(gcf,[saveprops.plotfolder,'/NormPriorPostContour_',nametag,'_',num2str(Tk)],'png')
    saveas(gcf,[saveprops.plotfolder,'/NormPriorPostContour_',nametag,'_',num2str(Tk)],'fig')
end
%%
figure(42)
surf(Xxnormbase,Xynormbase,pdfprior2postprobs_norm,'FaceColor','blue','EdgeColor','none','FaceAlpha',0.6);
camlight right; lighting phong
alpha 0.4
hold on
surf(Xxnormbase,Xynormbase,pdfpostprobs_norm,'FaceColor','red','EdgeColor','none','FaceAlpha',0.6);
camlight right; lighting phong
alpha 0.4


% plot3(Xn(:,1),Xn(:,2),pn,'bo')

if isempty(Xmc)==0
    plot(Xmcnorm(:,1),Xmcnorm(:,2),'r.')
end
if isempty(Xtruthnorm)==0
    plot(Xtruthnorm(:,1),Xtruthnorm(:,2),'k*','linewidth',2)
end

% title(['time step = ',num2str(Tk),' cond = ',num2str(cond(Pquad))])
xlabel('x_1')
ylabel('x_2')
title('Prior and Post Reg Surf norm')
axis equal
axis square
hold off
axis(axislims_norm)
if saveprops.saveit==1
    saveas(gcf,[saveprops.plotfolder,'/NormPriorPostSurf_',nametag,'_',num2str(Tk)],'png')
    saveas(gcf,[saveprops.plotfolder,'/NormPriorPostSurf_',nametag,'_',num2str(Tk)],'fig')
end
%%
%% postpdf norm space quad-gauss contours
% keyboard
figure(43)
contour(Xxnormbase,Xynormbase,QuadFilprior2postprobs_norm,clevelsnorm,'LineColor','b')
grid on
box off
hold on
contour(Xxnormbase,Xynormbase,QuadFilpostprobs_norm,clevelsnorm,'LineColor','r')

if isempty(Xmc)==0
    plot(Xmcnorm(:,1),Xmcnorm(:,2),'r.')
end
if isempty(Xtruth)==0
    plot(Xtruthnorm(:,1),Xtruthnorm(:,2),'k*','linewidth',2)
end

% title(['time step = ',num2str(Tk),' cond = ',num2str(cond(Pquad))])
xlabel('x_1')
ylabel('x_2')
title('Prior and Post Quad Contour norm')
axis equal
axis square
hold off
axis(axislims_norm)

if saveprops.saveit==1
    saveas(gcf,[saveprops.plotfolder,'/QuadNormPriorPostContour_',nametag,'_',num2str(Tk)],'png')
    saveas(gcf,[saveprops.plotfolder,'/QuadNormPriorPostContour_',nametag,'_',num2str(Tk)],'fig')
end
%%
figure(44)
surf(Xxnormbase,Xynormbase,QuadFilprior2postprobs_norm,'FaceColor','blue','EdgeColor','none','FaceAlpha',0.6);
camlight right; lighting phong
alpha 0.4
hold on
surf(Xxnormbase,Xynormbase,QuadFilpostprobs_norm,'FaceColor','red','EdgeColor','none','FaceAlpha',0.6);
camlight right; lighting phong
alpha 0.4


% plot3(Xn(:,1),Xn(:,2),pn,'bo')

if isempty(Xmc)==0
    plot(Xmcnorm(:,1),Xmcnorm(:,2),'r.')
end
if isempty(Xtruthnorm)==0
    plot(Xtruthnorm(:,1),Xtruthnorm(:,2),'k*','linewidth',2)
end

% title(['time step = ',num2str(Tk),' cond = ',num2str(cond(Pquad))])
xlabel('x_1')
ylabel('x_2')
title('Prior and Post Quad Surf norm')
axis equal
axis square
hold off
axis(axislims_norm)
if saveprops.saveit==1
    saveas(gcf,[saveprops.plotfolder,'/QuadNormPriorPostSurf_',nametag,'_',num2str(Tk)],'png')
    saveas(gcf,[saveprops.plotfolder,'/QuadNormPriorPostSurf_',nametag,'_',num2str(Tk)],'fig')
end
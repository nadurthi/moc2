function plotpdfs_prior_2D(Tk,pdfnorm,X,probs,mquadf,Pquadf,GMM,Xmc,Xtruth,saveprops)
nametag='prior';

Xn = pdfnorm.transForms.trueX2normX(X);
pn = pdfnorm.transForms.trueprob2normprob(probs);


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

[Xx,Xy]=meshgrid(linspace(-1.2,1.2,50),linspace(-1.2,1.2,50) );
Xx_true = zeros(size(Xx));
Xy_true = zeros(size(Xx));

pdfprobs_norm = zeros(size(Xx));
QuadFilprobs_norm = zeros(size(Xx));
GMMprobs_norm = zeros(size(Xx));

pdfprobs_true = zeros(size(Xx));
QuadFilprobs_true = zeros(size(Xx));
GMMprobs_true = zeros(size(Xx));

for i=1:size(Xx,1)
    for j=1:size(Xx,2)
        norm_pt = [Xx(i,j),Xy(i,j)];
        pdfprobs_norm(i,j) = pdfnorm.func(norm_pt);
        
        true_pt = pdfnorm.transForms.normX2trueX(norm_pt);
        Xx_true(i,j) = true_pt(1);
        Xy_true(i,j) = true_pt(2);
        
        pdfprobs_true(i,j) = pdfnorm.transForms.normprob2trueprob( pdfprobs_norm(i,j) );
        
        
        ptrue = mvnpdf(true_pt,mquadf(:)',Pquadf);
        QuadFilprobs_true(i,j) = ptrue;
        QuadFilprobs_norm(i,j) = pdfnorm.transForms.trueprob2normprob(ptrue);
        
        ptrue = evalGMM(GMM,true_pt);
        GMMprobs_true(i,j) = ptrue;
        GMMprobs_norm(i,j) = pdfnorm.transForms.trueprob2normprob(ptrue);

    end
end
%% contour levels
mm = max(max(pdfprobs_norm));
clevelsnorm = exp( linspace(-5,log(mm),15) );

mm = max(max(pdfprobs_true));
clevelstrue = exp( linspace(-5,log(mm),15) );

axislims_norm = [-1.2,1.2,-1.2,1.2];
%% norm contour, Xmc and Xtruth
figure(1)
contour(Xx,Xy,pdfprobs_norm,clevelsnorm)
grid on
box off
hold on
if isempty(Xmc)==0
    plot(Xmcnorm(:,1),Xmcnorm(:,2),'r.')
end
if isempty(Xtruth)==0
    plot(Xtruthnorm(:,1),Xtruthnorm(:,2),'k*','linewidth',2)
end

% title(['time step = ',num2str(Tk),' cond = ',num2str(cond(Pquad))])
xlabel('x_1')
ylabel('x_2')
title('Prior: Reg norm contour')
axis equal
axis square
hold off
axis(axislims_norm)
if saveprops.saveit==1
    saveas(gcf,[saveprops.plotfolder,'/NormContour_',nametag,'_',num2str(Tk)],'png')
    saveas(gcf,[saveprops.plotfolder,'/NormContour_',nametag,'_',num2str(Tk)],'fig')
end
%% norm surf, Xmc and Xtruth
figure(2)
surf(Xx,Xy,pdfprobs_norm,'FaceColor','green','EdgeColor','none','FaceAlpha',0.7);
camlight right; lighting phong
alpha 0.4
hold on
plot3(Xn(:,1),Xn(:,2),pn,'bo')

if isempty(Xmc)==0
    plot(Xmcnorm(:,1),Xmcnorm(:,2),'r.')
end
if isempty(Xtruthnorm)==0
    plot(Xtruthnorm(:,1),Xtruthnorm(:,2),'k*','linewidth',2)
end

% title(['time step = ',num2str(Tk),' cond = ',num2str(cond(Pquad))])
xlabel('x_1')
ylabel('x_2')
title('Prior: Reg norm surf')
axis equal
axis square
hold off
axis(axislims_norm)
if saveprops.saveit==1
    saveas(gcf,[saveprops.plotfolder,'/NormSurf_',nametag,'_',num2str(Tk)],'png')
    saveas(gcf,[saveprops.plotfolder,'/NormSurf_',nametag,'_',num2str(Tk)],'fig')
end
%% true contour, Xmc and Xtruth
figure(3)
[~,h]=contour(Xx_true,Xy_true,pdfprobs_true ,clevelstrue);
grid on
box off
% h.ContourZLevel = 1;
% view([26,43])
hold on
if isempty(Xmc)==0
    plot(Xmc(:,1),Xmc(:,2),'r.')
end
if isempty(Xtruth)==0
    plot(Xtruth(:,1),Xtruth(:,2),'k*','linewidth',2)
end

% title(['time step = ',num2str(Tk),' cond = ',num2str(cond(Pquad))])
xlabel('x_1')
ylabel('x_2')
title('Prior: Reg true contour')
axis equal
axis square
hold off
if saveprops.saveit==1
    saveas(gcf,[saveprops.plotfolder,'/TrueContour_',nametag,'_',num2str(Tk)],'png')
    saveas(gcf,[saveprops.plotfolder,'/TrueContour_',nametag,'_',num2str(Tk)],'fig')
end
%% true surf, Xmc and Xtruth
figure(4)
surf(Xx_true,Xy_true,pdfprobs_true,'FaceColor','r','EdgeColor','none','FaceAlpha',0.7);
camlight right; lighting phong
alpha 0.7
hold on
plot3(X(:,1),X(:,2),probs,'bo')

if isempty(Xmc)==0
    plot(Xmc(:,1),Xmc(:,2),'r.')
end
if isempty(Xtruth)==0
    plot(Xtruth(:,1),Xtruth(:,2),'k*','linewidth',2)
end

% title(['time step = ',num2str(Tk),' cond = ',num2str(cond(Pquad))])
xlabel('x_1')
ylabel('x_2')
title('Prior: Reg true surf')
axis equal
axis square
hold off
if saveprops.saveit==1
    saveas(gcf,[saveprops.plotfolder,'/TrueSurf_',nametag,'_',num2str(Tk)],'png')
    saveas(gcf,[saveprops.plotfolder,'/TrueSurf_',nametag,'_',num2str(Tk)],'fig')
end

%%
% -------------------- Quad Filter Plots------------------------------------------------------------------------
%%
%%  Quad Fil norm contour, Xmc and Xtruth
figure(5)
contour(Xx,Xy,QuadFilprobs_norm,clevelsnorm)
grid on
box off
hold on
if isempty(Xmc)==0
    plot(Xmcnorm(:,1),Xmcnorm(:,2),'r.')
end
if isempty(Xtruth)==0
    plot(Xtruthnorm(:,1),Xtruthnorm(:,2),'k*','linewidth',2)
end

% title(['time step = ',num2str(Tk),' cond = ',num2str(cond(Pquad))])
xlabel('x_1')
ylabel('x_2')
title('Prior: Quad norm contour')
axis equal
axis square
hold off
axis(axislims_norm)
if saveprops.saveit==1
    saveas(gcf,[saveprops.plotfolder,'/QuadFilNormContour_',nametag,'_',num2str(Tk)],'png')
    saveas(gcf,[saveprops.plotfolder,'/QuadFilNormContour_',nametag,'_',num2str(Tk)],'fig')
end
%%  Quad Fil norm surf, Xmc and Xtruth
figure(6)
surf(Xx,Xy,QuadFilprobs_norm,'FaceColor','green','EdgeColor','none','FaceAlpha',0.7);
camlight right; lighting phong
alpha 0.4
hold on
% plot3(Xn(:,1),Xn(:,2),pn,'bo')

if isempty(Xmc)==0
    plot(Xmcnorm(:,1),Xmcnorm(:,2),'r.')
end
if isempty(Xtruth)==0
    plot(Xtruthnorm(:,1),Xtruthnorm(:,2),'k*','linewidth',2)
end

% title(['time step = ',num2str(Tk),' cond = ',num2str(cond(Pquad))])
xlabel('x_1')
ylabel('x_2')
title('Prior: Quad norm surf')
axis equal
axis square
hold off
axis(axislims_norm)
if saveprops.saveit==1
    saveas(gcf,[saveprops.plotfolder,'/QuadFilNormSurf_',nametag,'_',num2str(Tk)],'png')
    saveas(gcf,[saveprops.plotfolder,'/QuadFilNormSurf_',nametag,'_',num2str(Tk)],'fig')
end
%%  Quad Fil true contour, Xmc and Xtruth
figure(7)
[~,h]=contour(Xx_true,Xy_true,QuadFilprobs_true ,clevelstrue);
grid on
box off
% h.ContourZLevel = 1;
% view([26,43])
hold on
if isempty(Xmc)==0
    plot(Xmc(:,1),Xmc(:,2),'r.')
end
if isempty(Xtruth)==0
    plot(Xtruth(:,1),Xtruth(:,2),'k*','linewidth',2)
end

% title(['time step = ',num2str(Tk),' cond = ',num2str(cond(Pquad))])
xlabel('x_1')
ylabel('x_2')
title('Prior: Quad true contour')
axis equal
axis square
hold off
if saveprops.saveit==1
    saveas(gcf,[saveprops.plotfolder,'/QuadFilTrueContour_',nametag,'_',num2str(Tk)],'png')
    saveas(gcf,[saveprops.plotfolder,'/QuadFilTrueContour_',nametag,'_',num2str(Tk)],'fig')
end
%%  Quad Fil true surf, Xmc and Xtruth
figure(8)
surf(Xx_true,Xy_true,QuadFilprobs_true,'FaceColor','r','EdgeColor','none','FaceAlpha',0.7);
camlight right; lighting phong
alpha 0.7
hold on
% plot3(X(:,1),X(:,2),probs,'bo')

if isempty(Xmc)==0
    plot(Xmc(:,1),Xmc(:,2),'r.')
end
if isempty(Xtruth)==0
    plot(Xtruth(:,1),Xtruth(:,2),'k*','linewidth',2)
end

% title(['time step = ',num2str(Tk),' cond = ',num2str(cond(Pquad))])
xlabel('x_1')
ylabel('x_2')
title('Prior: Quad true surf')
axis equal
axis square
hold off
if saveprops.saveit==1
    saveas(gcf,[saveprops.plotfolder,'/QuadFilTrueSurf_',nametag,'_',num2str(Tk)],'png')
    saveas(gcf,[saveprops.plotfolder,'/QuadFilTrueSurf_',nametag,'_',num2str(Tk)],'fig')
end
%%




% -------------------- GMM Filter Plots------------------------------------------------------------------------
%%
%%  GMM Fil norm contour, Xmc and Xtruth
figure(9)
contour(Xx,Xy,GMMprobs_norm,clevelsnorm)
grid on
box off
hold on
if isempty(Xmc)==0
    plot(Xmcnorm(:,1),Xmcnorm(:,2),'r.')
end
if isempty(Xtruth)==0
    plot(Xtruthnorm(:,1),Xtruthnorm(:,2),'k*','linewidth',2)
end

% title(['time step = ',num2str(Tk),' cond = ',num2str(cond(Pquad))])
xlabel('x_1')
ylabel('x_2')
title('Prior: GMM norm contour')
axis equal
axis square
hold off
axis(axislims_norm)
if saveprops.saveit==1
    saveas(gcf,[saveprops.plotfolder,'/GMMNormContour_',nametag,'_',num2str(Tk)],'png')
    saveas(gcf,[saveprops.plotfolder,'/GMMNormContour_',nametag,'_',num2str(Tk)],'fig')
end
%%  GMM Fil norm surf, Xmc and Xtruth
figure(10)
surf(Xx,Xy,GMMprobs_norm,'FaceColor','green','EdgeColor','none','FaceAlpha',0.7);
camlight right; lighting phong
alpha 0.4
hold on
% plot3(Xn(:,1),Xn(:,2),pn,'bo')

if isempty(Xmc)==0
    plot(Xmcnorm(:,1),Xmcnorm(:,2),'r.')
end
if isempty(Xtruth)==0
    plot(Xtruthnorm(:,1),Xtruthnorm(:,2),'k*','linewidth',2)
end

% title(['time step = ',num2str(Tk),' cond = ',num2str(cond(Pquad))])
xlabel('x_1')
ylabel('x_2')
title('Prior: GMM norm surf')
axis equal
axis square
hold off
axis(axislims_norm)
if saveprops.saveit==1
    saveas(gcf,[saveprops.plotfolder,'/GMMNormSurf_',nametag,'_',num2str(Tk)],'png')
    saveas(gcf,[saveprops.plotfolder,'/GMMNormSurf_',nametag,'_',num2str(Tk)],'fig')
end
%%  GMM Fil true contour, Xmc and Xtruth
figure(11)
[~,h]=contour(Xx_true,Xy_true,GMMprobs_true ,clevelstrue);
grid on
box off
% h.ContourZLevel = 1;
% view([26,43])
hold on
if isempty(Xmc)==0
    plot(Xmc(:,1),Xmc(:,2),'r.')
end
if isempty(Xtruth)==0
    plot(Xtruth(:,1),Xtruth(:,2),'k*','linewidth',2)
end

% title(['time step = ',num2str(Tk),' cond = ',num2str(cond(Pquad))])
xlabel('x_1')
ylabel('x_2')
title('Prior: GMM true contour')
axis equal
axis square
hold off
if saveprops.saveit==1
    saveas(gcf,[saveprops.plotfolder,'/GMMTrueContour_',nametag,'_',num2str(Tk)],'png')
    saveas(gcf,[saveprops.plotfolder,'/GMMTrueContour_',nametag,'_',num2str(Tk)],'fig')
end
%%  GMM Fil true surf, Xmc and Xtruth
figure(12)
surf(Xx_true,Xy_true,GMMprobs_true,'FaceColor','r','EdgeColor','none','FaceAlpha',0.7);
camlight right; lighting phong
alpha 0.7
hold on
% plot3(X(:,1),X(:,2),probs,'bo')

if isempty(Xmc)==0
    plot(Xmc(:,1),Xmc(:,2),'r.')
end
if isempty(Xtruth)==0
    plot(Xtruth(:,1),Xtruth(:,2),'k*','linewidth',2)
end

% title(['time step = ',num2str(Tk),' cond = ',num2str(cond(Pquad))])
xlabel('x_1')
ylabel('x_2')
title('GMM true surf')
axis equal
axis square
hold off
if saveprops.saveit==1
    saveas(gcf,[saveprops.plotfolder,'/GMMTrueSurf_',nametag,'_',num2str(Tk)],'png')
    saveas(gcf,[saveprops.plotfolder,'/GMMTrueSurf_',nametag,'_',num2str(Tk)],'fig')
end

end
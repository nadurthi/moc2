function c=treebased_integration(F,LB,UB,Xinitial,abserr)
[N,dim]=size(Xinitial);

%% train the tree
tic
Finitial = F(Xinitial);
states=cell(dim,1);
for i=1:dim
    states{i} = num2str(i);
end

tree = fitrtree(Xinitial,Finitial,'MinParentSize',2*dim,'MaxNumSplits',dim*1000,'MinLeafSize',2*dim,...
    'PredictorNames',states,'ResponseName','F');

boxes=getTree2Boxes(tree,LB,UB);

boxesbnds=zeros(size(boxes,1),2*dim);
for j=1:1:size(boxes,1)
    lb = boxes{j,1};
    ub = boxes{j,2};
    boxesbnds(j,:)=[lb(:)',ub(:)'];
end
disp(['time for tree-fit is ',num2str(toc)])

disp(['#boxes in tree-fit is ',num2str(size(boxes,1))])
figure(99)
clf
hold off
plot2Dboxes_modf(boxes)
% keyboard
%% loop over the boxed
Np=8;
P=cell(Np,2);
for i=1:Np
[xint,wint] = GLeg_pts((1+i)*ones(1,dim), -ones(1,dim), ones(1,dim));
P{i,1} = xint;
P{i,2} = wint;
end

        
Ib=zeros(size(boxes,1),1);
tic
parfor j=1:1:size(boxes,1)
%     lb = boxes{j,1};
%     ub = boxes{j,2};
%     boxdiag = norm(ub-lb);
%     volbox = prod(ub-lb );
    lb = boxesbnds(j,1:dim);
    ub = boxesbnds(j,dim+1:end);
    volbox = prod(ub-lb );
    prevII=100000;
    
    for ip=1:1:Np
        xint = P{ip,1};
        wint = P{ip,2};
        [xint,~]=transform_domain(xint,-ones(1,dim), ones(1,dim),lb,ub);

        II = volbox*sum(F(xint).*wint(:));
        EE = abs(II-prevII);
        if EE<abserr
            Ib(j)=II;
            break
        end
        prevII=II;
    end
    [j,ip,EE]
end
disp(['time for integration is ',num2str(toc)])

c = sum(Ib);

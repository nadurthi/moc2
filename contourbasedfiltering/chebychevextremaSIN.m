function X=chebychevextremaSIN(m)
X=zeros(m,1);
if m==1
    X=0;
    return
end

% for j=1:m
%    X(j) = sin(-pi/2*(j-1)/(m-1));   
% end
X = sin(linspace(-pi/2,pi/2,m));
X=round(X*1e14)/1e14;
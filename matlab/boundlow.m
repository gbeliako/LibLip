function res = boundlow(x,dim,p)
% lower bound
res = max(0,1- p *((1-x(1))+(1-x(2))));
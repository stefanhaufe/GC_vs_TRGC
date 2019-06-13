function [data, Arsig, x, lambdamax] = gen_ar_biv(N, P)
% Guido Nolte, 2006-2014
% Stefan Haufe, 2011-2014

% generates data according to bivariate AR model with unidirectional flow
% from first to second time series

% N: number of data-points
% P: order of AR-model

M = 2; %number of channels;
sigma = 0.2; %scale of random AR-parameters

N0=10000; %length of ignored start 

lambdamax=10;
while lambdamax > 1 || lambdamax < 0.7
  Arsig=[];
  for k=1:P
    aloc = zeros(M);
    aloc([1 2 4]) = randn(3, 1)*sigma;
    Arsig=[Arsig,aloc];
  end
  E=eye(M*P);AA=[Arsig;E(1:end-M,:)];lambda=eig(AA);lambdamax=max(abs(lambda));
end
% lambdamax

x=randn(M,N+N0);
y=x;
for i=P+1:N+N0;
    yloc=reshape(fliplr(y(:,i-P:i-1)),[],1);
    y(:,i)=Arsig*yloc+x(:,i);
end
data=y(:,N0+1:end);

x = x(:, N0+1:end);

Arsig = reshape(Arsig, M, M, P);

% figure; psd(data(1, :)')
% figure; psd(data(2, :)')






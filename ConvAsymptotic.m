% This program is to verify the asymptotic behaviour using MC simulations

clear all
close all

N = 100000; %sample size

% First distribution (Pareto)

% Parameters

%alpha = 0.09; %shape
alpha = 0.7;
%x0 = 0.0003; %scale
x0 = 1;

%RV1 = gamrnd(a1,b1,1,N);

RV1 = x0*(1-unifrnd(0,1,1,N)).^(-1/alpha);

% Second distribution (ML)

% Parameters

%th = 0.09; %shape
th = 0.7;
%nu = 6.1831; %scale
nu = 1/x0;

%RV2 = gamrnd(a2,b2,1,N);

U1 = unifrnd(0,1,1,N);
V1 = unifrnd(0,1,1,N);

RV2 = - ((1/nu)^(1/th)) * log(U1).*(sin(th*pi)./tan(th*pi*V1)-cos(th*pi)).^(1/th);

U2 = unifrnd(0,1,1,N);
V2 = unifrnd(0,1,1,N);

RV3 = - ((1/nu)^(1/th)) * log(U2).*(sin(th*pi)./tan(th*pi*V2)-cos(th*pi)).^(1/th);

% Comparison between RV1 and RV3

[ff1 xx1] = ecdf(RV1);
[ff2 xx2] = ecdf(RV3);

figure(1)

loglog(xx1,1-ff1,'or')
hold on
loglog(xx2,1-ff2,'ob')

ffth2 = 1./xx2.^th;
loglog(xx2,ffth2,'g')

% Sum of random variables

SUMRV = RV1 + RV2;

%SUMRV = RV2 + RV3;

%cutoff for display
index = find(SUMRV<1000);
SUMRV = SUMRV(index);

% Normalized histogram

[ff xx] = histnorm(SUMRV,300);
%edges = [0.003:0.001:max(SUMRV)];
%h = histogram(SUMRV,edges,'Normalization','pdf');
%xx = h.BinEdges;
%ff = h.Values;

figure(2)

loglog(xx,ff,'ob')

hold on

% Exponent -1

ffth1 = 1./xx;

loglog(xx,ffth1)

% integral

ffth2 = (xx.^(2*th-1)).*ml(-xx.^th,th,2*th,2);

loglog(xx,ffth2)

% here we check the asymptotic behaviour of the function

%clear all

%close all

%th = 0.7;

xxx = [10 50 100 500 1000 5000 10000 50000 100000 500000 1000000 5000000 10000000 50000000 100000000];

ffth3 = (xxx.^(2*th-1)).*ml(-xxx.^th,th,2*th,2);

figure (3)

loglog(xxx,ffth3)

hold on

loglog(xxx,1./(xxx.^(th+1)))

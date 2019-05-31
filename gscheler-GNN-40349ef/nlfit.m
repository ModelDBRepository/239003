
%
% fit histogram against lognormal
% using: histogram of log's fitted against Gauss
%
% uses EZFIT
%
function [mu_s,s_s,E] =nlfit(x)

mu=mean(log(x));
s=std(log(x));

f=figure;

hist(log(x),50);
%  fit=ezfit(sprintf('y=b*exp(-(log(x)-mu)^2/(2*s2))/(x*sqrt(2*pi*s2));b=%d;mu=%d',b0,mu0),'dispfitlegend','off','dispeqboxmode','off');

ft=ezfit('gauss');

E = ft.r;

close(f);

mu_s = exp(mu);
s_s = exp(s);


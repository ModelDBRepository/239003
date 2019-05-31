

function data = sim5_winit(data, pars)

switch(pars.W.init),

case 'learnt',
	'LL'
	data.W=zeros(pars.N,pars.N);
	
	for i=1:pars.N,
		for j=1:pars.N,
			data.W(i,j) = pars.lambda*data.I(i)*data.J(j);
		end;
	end;

case 'koulakov',
	'KK'
	W0=lognrnd(0.1,0.2,pars.N,pars.N);
	data.W=zeros(pars.N,pars.N);
	
	for i=1:pars.N,
		for j=1:pars.N,
			data.W(i,j) = pars.lambda*W0(i,j)*data.J(j);
		end;
	end;

	% This makes PRE_SYN-correlated: W=W';
%	data.W=data.W';

case 'lognormal_G',
	W0=lognrnd(pars.W.mean,pars.W.std,pars.N,pars.N);
	data.W=zeros(pars.N,pars.N);
	
	for i=1:pars.N,
		for j=1:pars.N,
			data.W(i,j) = pars.lambda*W0(i,j)*data.G(j);
		end;
	end;

case 'lognormal',
	'NL'
	data.W=lognrnd(pars.W.mean,pars.W.std,pars.N,pars.N);

case 'gauss',
	% random Gauss setting
	data.W= pars.W.mean + pars.W.std*randn(pars.N,pars.N);

case 'const',
	'CC'
	data.W= pars.W.mean + zeros(pars.N,pars.N);
otherwise, disp('wrong selector: W');
end;

switch(pars.W_RL.init),

case 'none',
	data.W_RL=zeros(pars.N,pars.N);
case 'learnt',
	'LL'
	data.W_RL=zeros(pars.N,pars.N);
	
	for i=1:pars.N,
		for j=1:pars.N,
			data.W_RL(i,j) = pars.lambda_RL*data.J(i)*data.J(j);
		end;
	end;

case 'lognormal',
	data.W_RL=lognrnd(pars.W_RL.mean,pars.W_RL.std,pars.N,pars.N);


case 'koulakov',
	'KK'
	W0=lognrnd(0.1,0.2,pars.N,pars.N);
	data.W_RL=zeros(pars.N,pars.N);
	
	for i=1:pars.N,
		for j=1:pars.N,
			data.W_RL(i,j) = pars.lambda*W0(i,j)*data.J(j);
		end;
	end;
	% This makes PRE_SYN-correlated: W=W';
%	data.W_RL=data.W_RL';


case 'gauss',
	% random Gauss setting
	data.W_RL= pars.W_RL.mean + pars.W_RL.std*randn(pars.N,pars.N);

case 'const',
	'CC'
	data.W_RL= pars.W_RL.mean + zeros(pars.N,pars.N);
otherwise, disp('wrong selector');
end;


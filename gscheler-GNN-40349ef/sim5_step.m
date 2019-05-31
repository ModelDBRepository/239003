

function data = sim4_step(data, pars)


EPSP=zeros(pars.N,1);

	%
	%   I---> J
	%
for j=1:pars.N,
   for k=data.I_act,
		%
		% only add connected neurons
		%
	if (data.C(k,j)),
		EPSP(j) = EPSP(j) + data.I(k)*data.W(k,j);
	end;
   end;
end;

	%
	%   J---> J
	%
for j=1:pars.N,
   for k=data.RL_act,
		%
		% only add connected neurons
		%
	if (data.C_RL(k,j)),
		EPSP(j) = EPSP(j) + data.out_J(k)*data.W_RL(k,j);
	end;
   end;
end;

data.out_EPSP = EPSP;

	% prepare IPSC as input-driven and Gauss background

if (pars.H.strength > 0),
	if (strcmp(pars.H.type,'normal')),
		data.IPSP= pars.H.mu + pars.H.sigma*data.rnd.NR;
	%	data.IPSP= pars.H.mu + pars.H.sigma*randn(pars.N,1);
	else
		data.IPSP=lognrnd(pars.H.mu,sqrt(pars.H.s2),pars.N,1);
	end;
else
	data.IPSP=zeros(pars.N,1);
end;
	
IPSP_step=data.IPSP;

	%
	% EPSP - IPSP no negative rates
	%
for j=1:pars.N,
	if (EPSP(j) > IPSP_step(j)),
		EPSP(j) = EPSP(j) -IPSP_step(j);
	else
		EPSP(j) = 1e-30;
	end;
end;


	% convert Hz --> mV
EPSP = pars.f_mV*EPSP;

	% using: 0.3mV = 40pA
EPSC = EPSP *0.04/0.3;    % into nA

data.out_input_pA = EPSC;

	%
	% go through neurons
	%
outJ=zeros(pars.N,1);
for j=1:pars.N,
	outJ(j) = data.G(j)*EPSC(j);
end;

data.out_J = outJ;


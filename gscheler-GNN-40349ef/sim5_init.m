%

function data = sim4_init(data,pars)

		% Hrom
		% these only used for weight calculation
	 data.I=lognrnd(pars.I.mu,sqrt(pars.I.s2),pars.N,1);
	 data.J=lognrnd(pars.J.mu,sqrt(pars.J.s2),pars.N,1);
	 data.out_J=data.J;


	data.G=lognrnd(pars.G.mu,sqrt(pars.G.s2),pars.N,1);

	data.rnd.R=randperm(pars.N);
	data.rnd.RL=randperm(pars.N);
	data.rnd.NR =randn(pars.N,1);



	data.W=zeros(pars.N,pars.N);
	data.W_RL=zeros(pars.N,pars.N);
	
	% active neurons and connections

	%--------------------------------------------------------------------
	% calculate connectivity
	%
		% excit I->J
	data.C = zeros(pars.N,pars.N);
	for j=1:pars.N,
		data.rnd.R=randperm(pars.N);
		data.C(data.rnd.R(1:pars.N*pars.N_perc_conn/100),j) = 1;
	end;
	
		% local recurrence
	data.C_RL = zeros(pars.N,pars.N);
	for j=1:pars.N,
		data.rnd.RL=randperm(pars.N);
		data.C_RL(data.rnd.RL(1:pars.N*pars.N_perc_RL_conn/100),j) = 1;
	end;
	
	%
	% set active neurons
	%
	data.I_act = 1:pars.N;
	data.RL_act = 1:pars.N;


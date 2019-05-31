%
%
%{

	Adjustment of gains


* various initial conditions

* G learn: hebb
* G learn: homeo

%}
%

global data;
global sim;

rand('state',2345);
randn('state',12345);

%-----------------------------------------------------------------
% EXPERIMENT
pars.version = 'psim20';
data.name = 'experiment';

	% number of simulation runs
pars.N_runs = 50;
	% for homeo
pars.N_runs = 80;

%-----------------------------------------------------------------
% ETC
	% fixed params
pars.f_mV = 1e-3;
pars.G.mu = 3;
pars.G.s2 = 0.31;

	% bins for the histogram
pars.bins = 20:2:120;

	% show histogram every other iteration
pars.disp_hist = 2;

%-----------------------------------------------------------------
% NEURONS AND CONNECTIVITY

	% number neurons
pars.N=1000;

	% percentage connectivity
pars.N_perc_conn = 10;
	% percentage recurrent connectivity
pars.N_perc_RL_conn = 0;
	% number of patterns
pars.N_pat = 1;

	% Inihibition
	% 0=no inhib; <>0: inhibition
pars.H.strength = 0;
pars.H.type = 'normal';
pars.H.sigma = 2;

	% that is adjusted wrt recurrence
pars.H.mu = 10;

	% adjustment of lambda parameters
pars.lambda = 0.01;
pars.lambda_RL = 0;


%-----------------------------------------------------------------
% LEARNING

	% adjust gains
pars.learn.gain 	= 1;
	% adjust weights
pars.learn.weight 	= 0;

	% selection of the learning rule (mixing)
	% 100 = only Hebb; 0 = only homeostatic
pars.mix.G = 100;
	% homeostatic
pars.mix.G = 0;

pars.mix.W = 0;


%-----------------------------------------------------------------
% INITIAL I, J populations
%
	%lognormal (RL)
pars.I.mu = 1.6;	% Hrom data
pars.I.s2 = 0.47;
pars.J.mu = 1.6;	% Hrom data
pars.J.s2 = 0.47;

	% Gaussian (RG)
%pars.I.mu = 1.6;	
%pars.I.s2 = 0.27;
%pars.J.mu = 1.6;	
%pars.J.s2 = 0.27;

%pars.I.mu = 1.6;	
%pars.I.s2 = 0.07;
%pars.J.mu = 1.6;	
%pars.J.s2 = 0.07;


%-----------------------------------------------------------------
% INITIAL Weights
%

	% not used here
%pars.W.init = 'koulakov';

	% lognormal (WL)
pars.W.init = 'lognormal';
pars.W.mean = log(1);
pars.W.std = 1.5;

	% Gaussian (WG)
pars.W.init = 'gauss';
pars.W.mean = 1;
pars.W.std = 0.2;  %0.5 produces neg values

	% Uniform (uses Gaussian with very small sigma)
%pars.W.std = 1e-6; 


	% no recurrent connections
pars.W_RL.init = 'none'; 
pars.W_RL.mean = 0;
pars.W_RL.std = 0;


%=================================================
% fixed initialization
%=================================================
data = sim5_init(data, pars);
data = sim5_winit(data, pars);

%-----------------------------------------------------------------
% INITIAL Gains (must be done after sim5 initialization
%
	% Gaussian (GN)
data.G = 27 + 2*randn(pars.N,1);

	% uniform (GU)
%data.G = 50 + 1e-6*randn(pars.N,1);

	% log-normal (GL)
%data.G=lognrnd(3.3,sqrt(0.1),pars.N,1);


%=================================================
% run the simulator
%


	%
	% initial figure (histogram of gains)
	%
figure;
hist(data.G,pars.bins);
set(gca,'Fontsize',16);
xlabel('gain','fontsize',16);
ylabel('number','fontsize',16);
drawnow;


%-----------------------------------------------------------------
	% simulation loop over updates
for run=1:pars.N_runs

	data.pattern_out = zeros(pars.N,pars.N_pat);
	for p=1:pars.N_pat,

		%--------------------------------------------
		% execute one step
		%--------------------------------------------
		data = sim5_step(data,pars);
		rand;
		data.pattern_out(:,p) = data.out_J;
	end;

		%--------------------------------------------
		% MIX at each iteration
		%--------------------------------------------
	data.mix.G = zeros(pars.N,1);
	dd=randperm(pars.N);
	data.mix.G(dd(1:floor(pars.mix.G/100*pars.N))) = 1;

	data.mix.W = zeros(pars.N*pars.N,1);
	dd=randperm(pars.N*pars.N);
	data.mix.W(dd(1:floor(pars.mix.W/100*pars.N*pars.N))) = 1;

		%---------------------------------
		% weight update W
		%---------------------------------
	if (pars.learn.weight && (run > 1)),
	  M=mean(data.out_J)^2;
	  for i=1:pars.N,
		for j=1:pars.N
			if (i~=j),
			  if (data.mix.W((i-1)*pars.N+j))
					% ++ learning
				data.W(i,j) = data.W(i,j)* ...
			(1 + 0.005*(data.out_J_old(i)*data.out_J(j)-5)); 

			  else 
					% -- learning
				data.W(i,j) = data.W(i,j)* ...
			(1 - 0.005*(data.out_J_old(i)*data.out_J(j)-M)); 
			  end;

			end;
		end;
	  end;
	end;

	%---------------------------------
	% gain update 
	%---------------------------------
	if (pars.learn.gain==1),
	    for j=1:pars.N,
		if (data.mix.G(j)),
			% "+" adapt
		    data.G(j)= data.G(j)* (1+ 0.01*(data.out_J(j) -4.35));
		else
			% "-" adapt
		    data.G(j)= data.G(j)*(1 - 0.01*(data.out_J(j) - 2.6));
		end;
	end;
	end;


		%
		% analyze the results
		%
	fprintf('run=%g\n',run);

	[mu_s_G,s_s_G,E_G] = nlfit(data.G);
	all_s_s_G(run) = s_s_G;
	all_mu_s_G(run) = mu_s_G;
	fprintf('MU*_G=%g SIGMA*_G=%g E=%g\n',mu_s_G,s_s_G,E_G);

	[mu_s,s_s,E] = nlfit(reshape(data.W,pars.N*pars.N,1));
	all_s_s(run) = s_s;
	all_mu_s(run) = mu_s;
	fprintf('MU*_W=%g SIGMA*_W=%g E=%g\n',mu_s,s_s,E);

	[mu_s_RJ,s_s_RJ,E_RJ] = nlfit(data.out_J);
	all_s_s_RJ(run) = s_s_RJ;
	all_mu_s_RJ(run) = mu_s_RJ;
	fprintf('MU*_RJ=%g SIGMA*_RJ=%g E=%g\n',mu_s_RJ,s_s_RJ,E_RJ);

	data.out_J_old = data.out_J;

		%
		% plot new figure every other update
		%
	if (mod(run,pars.disp_hist)==0),
		figure;
		xx= data.G(find(data.G < 120));
		hist(xx,pars.bins);
		set(gca,'Fontsize',16);
		xlabel('gain','fontsize',16);
		ylabel('number','fontsize',16);
		title(sprintf('iteration %g',run));
		drawnow;
	end;

end;

%saveas(gcf,'p34_G_GN_WN_RIJN_1','fig')
% based on psim34_G_ad_upd.m

% Toy model of unemployment
% gauthier@vermandel.fr

%----------------------------------------------------------------
% 0. Housekeeping (close all graphic windows)
%----------------------------------------------------------------

close all;

%----------------------------------------------------------------
% 1. Defining variables
%----------------------------------------------------------------

var rr c n u w y k i lb mc pi r q x v_H v_P e mu g tau gy_obs gc_obs gi_obs pi_obs r_obs u_obs varrho;
var e_a e_g e_c e_m e_i e_r e_t;


varexo eta_a eta_g eta_c eta_m eta_i eta_r eta_t;
 
parameters beta delta alpha sigmaC sigmaL delta_N chi phi gy b  Gam eta gamma epsilon kappa rho phi_y phi_pi xi
			tau0 y0 sig theta1 theta2 varphi A xi piss  rho_a rho_g rho_c rho_m rho_i rho_r rho_t;
            
            
%----------------------------------------------------------------
% 2. Calibration
%----------------------------------------------------------------
delta_N = .1;		% separation rate
eta		= .5;		% negotiation share
phi		= 0.05;		% shape hiring cost function
beta 	= 0.993; 	% Discount factor firms
delta 	= 0.025;	% Depreciation rate
alpha 	= 0.30;		% Capital share
gy 		= 0.2;   	% Public spending in GDP
sigmaC 	= 1;		% Consumption risk aversion
sigmaL 	= 2; 		% Elasticity of labor
epsilon = 10;		% Elasticity between goods
rho 	= .8;		% Monetary policy smoothing
phi_y	= 0.1;		% Monetary policy reaction to output
phi_pi	= 1.5;		% Monetary policy reaction to inflation
xi 		= 80;		% Adjustment costs on prices
kappa	= 4;		% adjustment costs on investment
gamma	= .85;		% unemployment insurance as % of real wage
varphi	= 0.2;		% elasticity of emission to GDP
piss	= 1.005;	% steady state inflation

% value of main variables:
tau0 	= 50 /1000;	% value of carbon tax ($/ton)
sig		= 0.2; 		% Carbon intensity USA 0.2 Gt / Trillions USD
y0	 	= 25;		% trillions usd PPA https://data.worldbank.org/indicator/NY.GDP.MKTP.CD
theta1  = 0.05;		% level of abatement costs
theta2  = 2.6;		% curvature abatement cost
Hss		= 1/3;		% labor supply in ss

% autoregressive roots parameters
rho_a	= 0.95;
rho_g	= 0.95;
rho_c	= 0.95;
rho_m  	= 0.95;
rho_i	= 0.95;
rho_r	= 0.40;
rho_t	= 0.8;

%----------------------------------------------------------------
% 3. Model
%----------------------------------------------------------------

model; 
	%% Household
	[name='marginal utility of consumption']
	lb = e_c*c^-sigmaC;
	[name='Euler']
	beta*lb(+1)/lb*rr = 1;
	[name='Labor Supply']
	v_H = w - b - chi/lb + beta*(1-delta_N)*v_H(+1);
	[name='Employment law of motion']
	n = (1-delta_N)*n(-1) + x*u;

	%% Production
	[name='technology']
	y = A*e_a*(k(-1)^alpha)*(n^(1-alpha));
	[name='Capital law of motion']
	i*e_i*(1-(kappa/2)*(i/i(-1)-1)^2) = k-(1-delta)*k(-1);
	[name='FOC k']
	((1-delta)*q(+1)+alpha*varrho(+1)*y(+1)/k) = q*rr;
	[name='FOC i']
 	e_i*q = 1 + e_i*q*(kappa/2)*( 1 + ( 3*i/i(-1)-4 )*i/i(-1) )
			+ beta*lb(+1)/lb*e_i(+1)*q(+1)*kappa*(1-i(+1)/i)*(i(+1)/i)^2;
	[name='FOC n']
	v_P = (1-alpha)*varrho*y/n - w + beta*lb(+1)/lb*(1-delta_N)*v_P(+1);
	[name='FOC x']
	Gam*x^phi = v_P;
	[name='NKPC']
	(1-epsilon) + epsilon*mc - xi*pi*(pi-steady_state(pi)) + xi*beta*lb(+1)/lb*y(+1)/y*pi(+1)*(pi(+1)-steady_state(pi));
	[name='wage']
	(1-e_m*eta)*v_H = e_m*eta*v_P;
	[name='FOC y']
	varrho = mc - theta1*mu^theta2 - tau*(1-varphi)*sig*(1-mu)*y^-varphi;
	[name='FOC mu']
	((tau*sig*y^(1-varphi))/(theta2*theta1))^(1/(theta2-1)) = mu;
	
	
	%% AGGREGATION
	[name='Resources Constraint']
	y = c + i + g  + xi/2*(pi-steady_state(pi))^2*y + Gam/(1+phi)*x^(1+phi)*u;
	[name='Unemployment rate']
	u = 1-n;	
	[name='Fisherian equation']
	rr = r/pi(+1);
	[name='Total emissions']
	e = sig*(1-mu)*y^(1-varphi);
	
	%%% Policy instruments
	[name='Monetary Policy rule']
	r = r(-1)^rho * (steady_state(r)*(pi/steady_state(pi))^phi_pi*(y/steady_state(y))^phi_y)^(1-rho) * e_r;
	[name='Public spending']
	g = gy*steady_state(y)*e_g;
	[name='Carbon tax']
	tau = tau0*e_t;
	
	
	%% Observable variables 
	[name='measurement GDP']
	gy_obs = log(y/y(-1));
	[name='measurement consumption']
	gc_obs = log(c/c(-1));
	[name='measurement investment']
	gi_obs = log(i/i(-1));
	[name='measurement inflation']
	pi_obs = pi - steady_state(pi);
	[name='measurement interest rate']
	r_obs  = r  - steady_state(r);
	[name='measurement unemployment']
	u_obs  = u  - steady_state(u);
	
	[name='shocks']
	log(e_a) = rho_a*log(e_a(-1))+eta_a;
	log(e_g) = rho_g*log(e_g(-1))+eta_g;
	log(e_c) = rho_c*log(e_c(-1))+eta_c;
	log(e_i) = rho_i*log(e_i(-1))+eta_i;
    log(e_m) = rho_m*log(e_m(-1))+eta_m;
    log(e_r) = rho_r*log(e_r(-1))+eta_r;  
    log(e_t) = rho_t*log(e_t(-1))+eta_t;  
end;

    
%----------------------------------------------------------------
% 4. Computation
%----------------------------------------------------------------
steady_state_model;
	pi		= piss;
	rr		= 1/beta;
	r		= rr*pi;
	tau 	= tau0;
	y		= y0;
	mu		= (tau*sig*y^(1-varphi)/(theta2*theta1))^(1/(theta2-1));
	mc		= (epsilon-1)/epsilon;
	varrho 	= mc - theta1*mu^theta2 - tau*(1-varphi)*sig*(1-mu)*y^(-varphi);
	q		= 1;
	u 		= .075;
	n 		= 1-u;
%	k		= n*(A*(rr-(1-delta))/(alpha*varrho))^(1/(alpha-1));
	k		=	alpha*varrho*y/(rr-(1-delta));

	A		= y/(k^alpha*n^(1-alpha));
	i		= delta*k;
	x 		= delta_N*n/u;	
	Gam 	= 0.01*y / (x^(1+phi)*u/(1+phi));
	v_P 	= Gam*x^phi;
	v_H 	= eta/(1-eta)*v_P;
	w 		= (1-alpha)*varrho*y/n - v_P  + beta*((1-delta_N-0*x)*v_P+0*phi/2*x^2);
	c  		= (1-gy)*y-i-Gam/(1+phi)*x^(1+phi)*u;
	lb 		= c^-sigmaC;
	chi 	= ((1-gamma)*w - v_H*(1-beta*(1-delta_N-0*x)))*lb;
	g 		= gy*y;
	b		= gamma*w;
	e 		= sig*(1-mu)*y^(1-varphi);
	e_a 	= 1;
	e_g 	= 1;
	e_c 	= 1;
	e_m 	= 1;
	e_i 	= 1;
	e_r 	= 1;
	e_t 	= 1;
	gy_obs = 0; gc_obs = 0; gi_obs = 0; pi_obs = 0; r_obs = 0; u_obs = 0; 
end;


%%% SIMULATIONS
shocks;
	var eta_a;	stderr 0.01;
	var eta_g;	stderr 0.01;
	var eta_c;	stderr 0.01;
	var eta_m;	stderr 0.01;
	var eta_i;	stderr 0.01;
	var eta_r;	stderr 0.01;
end;
	
resid(1);
check;

stoch_simul(irf=30,order=1) y c i pi r u x ;

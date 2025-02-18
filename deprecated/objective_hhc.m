% Industrial Network, see: `Harmonic Optimal Power Flow with Transformer
% Excitation` by F. Geth and T. Van Acker, pg. 8, § IV.B.

function mpc = industrial_network
mpc.version = '2';
mpc.baseMVA =  100.00;
mpc.principle = 'maximin'

%% bus data
%	bus_id	type    Pd      Qd	    Gs	    Bs	    area	Vm	    Va	    baseKV  zone	Vmax	Vmin
mpc.bus = [
    0       3       0.000   0.000   0.00    0.00    1       1.00    0.00    155.0   1       1.00    1.00;
    1       1       1.000   1.000   0.00    0.00    1       1.00    0.00    155.0   1       1.06    0.94;
    2       1       1.00    1.000   0.00    0.00    1       1.00	0.00    155.0    1       1.06    0.94;
];

%% bus harmonic data 
%column_names%  ref_angle   standard 
mpc.bus_harmonics = [
                0.0         'Clean Bus';                % 0
                0.0         'IEC61000-2-4:2002, Cl. 2'; % 1
                0.0         'IEC61000-2-4:2002, Cl. 2'; % 2
]

%% branch data
%   f_bus	t_bus	r	        x	        b	        rateA	rateB	rateC	ratio	angle	status	angmin	angmax
mpc.branch = [ 
    0       1       0.000000000 0.050000000 0.000000000 9999999 1000    1000    1       0       1       -60     60; % Ssc=2.1e9, XRr=20   
    0       2       0.000000000 0.0250000000 0.000000000 9999999 1000    1000    1       0       1       -60     60; % Ssc=2.1e9, XRr=20    
];

%% generator data
%   bus     Pg      Qg      Qmax    Qmin    Vg      mBase   status  Pmax    Pmin
mpc.gen = [
    0       0.0     0.0     100.00  -100.00 1.05    100.0   1       100.00  0;
];
%column_names%  inf     gsc     bsc
mpc.gen_imp = [
                1       0.0     0.0;
]
mpc.gencost = [
	2	    0.0	    0.0	    3	    0.0     1.0     0.0;
];
% Two Branch Network, designed to test the validity of the rotated second-order
% cone formulation of the load current magnitude constraint.

function mpc = two_branch_network
mpc.version = '2';
mpc.baseMVA =  100.00;

%% bus data
%	bus_id	type    Pd      Qd	    Gs	    Bs	    area	Vm	    Va	    baseKV  zone	Vmax	Vmin
mpc.bus = [
    0       3       0.000   0.000   0.00    0.00    1       1.00    0.00    150.0   1       1.00    1.00;
    1       1       0.000   0.000   0.00    0.00    1       1.00    0.00    150.0   1       1.06    0.94;
    2       1       20.00   5.000   0.00    0.00    1       1.00	0.00    150.0   1       1.06    0.94;
    3       1       20.00   5.000   0.00    0.00    1       1.00    0.00    150.0   1       1.06    0.94;
];

%% bus harmonic data 
%column_names%  ref_angle   angle_range standard 
mpc.bus_harmonics = [
                0.0         0.174532925 'Clean Bus';                % 0
                0.0         0.174532925 'IEC61000-3-6:2008';        % 1
                0.785398164 0.174532925 'IEC61000-2-4:2002, Cl. 2'; % 2
                0.785398164 0.174532925 'IEC61000-2-4:2002, Cl. 2'; % 3
]

%% generator data
%   bus     Pg      Qg      Qmax    Qmin    Vg      mBase   status  Pmax    Pmin
mpc.gen = [
    0       0.0     0.0     100.00  -100.00 1.05    100.0   1       100.00  0;
];
mpc.gencost = [
	2	    0.0	    0.0	    3	    0.0     1.0     0.0;
];
%column_names%  isfilter
mpc.gen_extra = [
	0;
];

%% branch data
%   f_bus	t_bus	r	        x	        b	        rateA	rateB	rateC	ratio	angle	status	angmin	angmax
mpc.branch = [ 
    0       1       0.0024      0.048       0.00000     1000    1000    1000    1       0       1       -60     60; % Ssc=2.1e9, XRr=20     
    1       2       0.060       0.53        0.00000     125     1000    1000    1       0       1       -60     60; 
    1       3       0.060       0.53        0.00000     125     1000    1000    1       0       1       -60     60; 
];

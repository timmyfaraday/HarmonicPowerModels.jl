% Case to test space based matlab matrix
% And other hard to parse cases
% also test data without a generator cost model

function mpc = case_xfmr_Yy0
mpc.version = '2';
mpc.baseMVA =  100.00;

%% bus data
%	bus_id	type    Pd      Qd	    Gs	    Bs	    area	Vm	    Va	    baseKV  zone	Vmax	Vmin
mpc.bus = [
    1       3       0.000   0.000   0.00    0.00    1       1.00    0.00    12.47   1       1.10    0.90;
    2       1       1.500   1.00   0.00    0.00    1       1.00	0.00    4.16    1       1.10    0.90;
];

%% bus harmonic data 
%column_names%  nh_1    nh_3      thdmax
mpc.bus_harmonics = [
                1.000   0.0000   0.5
                1.000   0.4000   0.5
]

% we need to add a filter on the load side to avoid infeasibility - the transformer blocks 3rd harmonics
%% generator data
%   bus     Pg      Qg      Qmax    Qmin    Vg      mBase   status  Pmax    Pmin
mpc.gen = [
    1       0.0     0.0     100.00  -100.00 1.05    100.0   1       100.00  -0.0;
    2       0.0     0.0     100.00  -100.00 1.05    100.0   1       100.00  -0.0;
];

mpc.gencost = [
	2	 0.0	 0.0	 3	   0.00000	   1.000000	   0.000000;
	2	 0.0	 0.0	 3	   0.00000	   2.000000	   0.000000;
];

%column_names%  isfilter
mpc.gen_extra = [
	0;
	1;
];


%% branch data
%   f_bus	t_bus	r	    x	    b	    rateA	rateB	rateC	ratio	angle	status	angmin	angmax
mpc.branch = [ 
];

%% transformer data
%column_names%  f_bus   t_bus   xsc     r1      r2      vg      gnd1    gnd2    re1     xe1     re2     xe2   
mpc.xfmr = [
                1       2       5     1     1     'Yd11'  0       0       0.0     0.0     0.0     0.0
]

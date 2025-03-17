% two bus network with single branch w/o shunt 

function mpc = two_bus_example
mpc.version = '2';
mpc.baseMVA =  1.00;

%% bus data
%	bus_id	type    Pd      Qd	    Gs	    Bs	    area	Vm	    Va	    baseKV  zone    Vmax	Vmin
mpc.bus = [
    1       4       0.000   0.000   0.00    0.00    1       1.00    0.00    10      1       1.10    0.90;
    2       1       1.000   0.500   0.00    0.00    1       1.00	0.00    10      1       1.10    0.90;
];

%% bus harmonic data 
%column_names%  nh_1    nh_2    nh_3    std
mpc.bus_harmonics = [
                1.000   0.000   0.000   'IEC61000-2-4:2002, Cl. 2'
                1.000   0.050   0.050   'IEC61000-2-4:2002, Cl. 2'
]

%% generator data
%   bus     Pg      Qg      Qmax    Qmin    Vg      mBase   status  Pmax    Pmin
mpc.gen = [
    1       0.0     0.0     9.00    -9.00   1.05    100.0   1       10.00   -10.0;
];

%% branch data
%   f_bus	t_bus	r	    x	    b	    rateA	rateB	rateC	ratio	angle	status	angmin	angmax
mpc.branch = [ 
    1       2       0.005   0.01    0       10      10      10      1       0       1       -60     60;
];


% Minimal L-C example of the robust harmonic hosting capacity: one inductive branch from an
% ideal source (clean bus 0) to bus 1, which holds a single harmonic unit and the branch
% capacitance (the to-side half of the line charging; the from-side half sits at the clean
% bus and is shorted at the harmonics). With x = 0.1 pu and b_to = 0.2778 pu the parallel
% resonance seen from bus 1 is at h = 1/sqrt(x*b_to) = 6. A single unit makes the injection
% angle irrelevant, so the hosting capacity at harmonic h is sqrt(limit_h)/|Z_11(h)|.

function mpc = lc_minimal
mpc.version = '2';
mpc.baseMVA =  100.00;
mpc.principle = 'absolute equality'

%% bus data (type 4: clean bus, i.e., ideal source at the fundamental and ideal ground at harmonics)
%	bus_id	type    Pd      Qd	    Gs	    Bs	    area	Vm	    Va	    baseKV  zone	Vmax	Vmin
mpc.bus = [
    0       4       0.000   0.000   0.00    0.00    1       1.00    0.00    33.0    1       1.00    1.00;
    1       1       50.00   0.000   0.00    0.00    1       1.00    0.00    33.0    1       1.10    0.90;
];

%% bus harmonic data
%column_names%  ref_angle   standard
mpc.bus_harmonics = [
                0.0         'Clean Bus';
                0.0         'IEC61000-2-4:2002, Cl. 2';
]

%% branch data
%   f_bus	t_bus	r	        x	        b	        rateA	rateB	rateC	ratio	angle	status	angmin	angmax
mpc.branch = [
    0       1       0.010000    0.100000    0.555556    1000    1000    1000    1       0       1       -60     60;
];

%% generator data
%   bus     Pg      Qg      Qmax    Qmin    Vg      mBase   status  Pmax    Pmin
mpc.gen = [
    0       0.0     0.0     100.00  -100.00 1.00    100.0   1       100.00  0;
];
%column_names%  inf     gsc     bsc
mpc.gen_imp = [
                1       0.0     0.0;
]
mpc.gencost = [
	2	    0.0	    0.0	    3	    0.0     1.0     0.0;
];

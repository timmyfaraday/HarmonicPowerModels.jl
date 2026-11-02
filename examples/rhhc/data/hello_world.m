% Hello-world network for LineCableModels + HarmonicPowerModels: one 33 kV voltage level, no
% transformers. A 150 MVA source feeds a busbar; two XLPE cable sections (15 km + 10 km) and one
% overhead line (10 km) supply three loads. The ~2.5 Mvar of cable charging against the source
% inductance gives a parallel resonance near the 7th-8th harmonic.
%
% The branch r, x, b below are 50 Hz placeholders. examples/rhhc/lcm_network_model.jl replaces the
% cable and line branches at every harmonic order with LineCableModels results (exact pi model).

function mpc = network
mpc.version = '2';
mpc.baseMVA =  100.00;
mpc.principle = 'maximin'

%% bus data (type 4: clean bus, i.e., ideal source at the fundamental and ideal ground at harmonics)
%	bus_id	type    Pd      Qd	    Gs	    Bs	    area	Vm	    Va	    baseKV  zone	Vmax	Vmin
mpc.bus = [
    0       4       0.000   0.000   0.00    0.00    1       1.00    0.00    33.0    1       1.00    1.00;
    1       1       0.000   0.000   0.00    0.00    1       1.00    0.00    33.0    1       1.10    0.90;
    2       1       8.000   2.000   0.00    0.00    1       1.00    0.00    33.0    1       1.10    0.90;
    3       1       5.000   1.500   0.00    0.00    1       1.00    0.00    33.0    1       1.10    0.90;
    4       1       3.000   1.000   0.00    0.00    1       1.00    0.00    33.0    1       1.10    0.90;
];

%% bus harmonic data
%column_names%  ref_angle   standard
mpc.bus_harmonics = [
                0.0         'Clean Bus';                % 0, ideal source behind the source impedance
                0.0         'IEC61000-2-4:2002, Cl. 2'; % 1
                0.0         'IEC61000-2-4:2002, Cl. 2'; % 2
                0.0         'IEC61000-2-4:2002, Cl. 2'; % 3
                0.0         'IEC61000-2-4:2002, Cl. 2'; % 4
]

%% branch data
%   f_bus	t_bus	r	        x	        b	        rateA	rateB	rateC	ratio	angle	status	angmin	angmax
mpc.branch = [
    0       1       0.066340    0.663400    0.000000    1000    1000    1000    1       0       1       -60     60; % source, Ssc = 150 MVA, X/R = 10
    1       2       0.073200    0.134400    0.015310    30      30      30      1       0       1       -60     60; % XLPE 18/30 kV 1x630 Al trefoil, 15 km
    2       3       0.048800    0.089600    0.010210    30      30      30      1       0       1       -60     60; % XLPE 18/30 kV 1x630 Al trefoil, 10 km
    1       4       0.083500    0.316800    0.000366    20      20      20      1       0       1       -60     60; % overhead line, Al ~300 mm2, flat 1.5 m at 10 m, 10 km
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

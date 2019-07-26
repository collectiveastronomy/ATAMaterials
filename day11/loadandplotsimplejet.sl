()=evalfile("simplejet.sl");
fit_fun("simplejet");
list_par;
eval_counts;

()=evalfile("plotscript.sl");

%now cache model on grid
variable loval = _A(10^[-9:3:0.001]);
variable hival = make_hi_grid (loval);
cache_fun("simplejet",loval,hival);
fit_fun("simplejet_cache");
list_par
eval_counts;


evalfile("plotscript.sl");

%%load in some good starting parameters
load_par("sj_cache_start.par");
eval_counts;
evalfile("plotscript.sl);

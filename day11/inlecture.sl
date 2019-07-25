.load load_data
list_data;
xlog;ylog;
set_plot_widths(;m_width=5, r_width=5, re_width=5, d_width=5);
fancy_plot_unit("hz","mjy");
plot_unfold({rid,pid,hid};dsym={4,4,4},dcol={6,4,8}, decol={6,4,8},xrange={3e8,1e20},res=0);
exclude(rid);
fit_fun("constant(Isis_Active_Dataset)*phabs(1)*highecut(1)*(bknpower(1)+diskbb(1)+ gaussian(1))");
list_par;
set_par("phabs(1).nH",0.6,1);
set_par("constant(1).factor",1,1);
list_par;
set_par(14,1,0.9,1.1);
set_par(14,1,0,0.9,1.1);
list_par;
eval_counts;
fit_counts;
fancy_plot_unit("kev","ergs");
plot_unfold({pid,hid};dsym={4,4},dcol={4,8}, decol={4,8},res=1);;
plot_unfold({pid,hid};dsym={4,4},dcol={4,8},decol={4,8},rcol={4,8},recol={4,
8},xrange={2,200},res=1);
list_par;
fit_counts;
()=evalfile("simplejet.sl");
fit_fun("simplejet");
list_par;
fit_fun("constant(Isis_Active_Dataset)*phabs(1)*reflect(1,
fit_fun("constant(Isis_Active_Dataset)*phabs(1)*reflect(1,(nthcomp(1)+gaussian(1)))");
list_par;
freeze(10);
fit_counts;
plot_unfold({pid,hid};dsym={4,4},dcol={4,8},decol={4,8},rcol={4,8}
,recol={4,8},xrange={2,200},res=1);
list_par;
include(rid);
usr_grid([rid,pid],-9,3,0.001,1);
usr_grid([hid],0,3,0.001);
evalfile("simplejet.sl");
fit_fun("simplejet");
list_par;
eval_counts;
list_data;
loval = _A(10^[-9:3:0.001]);
hival = make_hi_grid (loval);
cache_fun("simplejet",loval,hival);
fit_fun("simplejet_cache");
list_par;
eval_counts;
fancy_plot_unit("hz","mjy");
.load plotscript
load_par("sj_cache_start.par")
;
list_par;
eval_counts;
.load plotscript
fit_fun("constant(Isis_Active_Dataset)*phabs(1)*highecut(1)*(bknpowe
fit_fun("constant(Isis_Active_Dataset)*phabs(1)*reflect(1,(simplejet_cache(1)+gaussian(1)+diskbb(1)))");
list_par;
set_par(23,1e-3);
freeze(23);
eval_counts;
.load plotscript
save_par("combinedjet_start.par");
set_par(22,100);
eval_counts;
set_par(22,1e5);
eval_counts;
set_par(22,1e7);
eval_counts;
list_par;
set_par(23,1.e3,0,1.e-4, 1000);
list_par;
set_par(23,1.e-3,0,1.e-4, 1000);
list_par;
eval_counts;
.load plotscript
set_par(23,1.e-3,1,1.e-4, 1000);
set_par(22,1.e8);
eval_counts;
save_par("combinedjet_start2.par");
save_input("inlecture.sl");

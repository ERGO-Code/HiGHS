% Makes bar plot of the fraction of the IPX runtime spend in different parts
% of the algorithm. Uses arrays time_total, time_starting_basis,
% time_kkt_factorize, time_kkt_solve and time_crossover from Matlab session.

% only use models that were solved, but took more than 10s
p = find(isfinite(time_total) & (time_total>10));
time_total = time_total(p);
time_starting_basis = time_starting_basis(p);
time_kkt_factorize = time_kkt_factorize(p);
time_kkt_solve = time_kkt_solve(p);
time_crossover = time_crossover(p);

% sort by total runtime
[time_total,p] = sort(time_total);
time_starting_basis = time_starting_basis(p);
time_kkt_factorize = time_kkt_factorize(p);
time_kkt_solve = time_kkt_solve(p);
time_crossover = time_crossover(p);

% make bar plot
figwidth = 398.3;
figheight = 254.0;
figure('units','points','position',[0,0,figwidth,figheight]);
y = [time_starting_basis./time_total ...
     time_kkt_factorize./time_total ...
     time_kkt_solve./time_total ...
     time_crossover./time_total];
b = bar(y, 'stacked');
b(1).FaceColor = 'k';
b(4).FaceColor = 'k';
for i=1:4
  b(i).EdgeColor = 'w';
end
ax = gca();
ax.TickDir = 'out';

% set x labels to total time
x1 = find(time_total>1e1,1);
x2 = find(time_total>1e2,1);
x3 = find(time_total>1e3,1);
x4 = find(time_total>1e4,1);
xticks([x1 x2 x3 x4]);
xticklabels({'10s', '100s', '1000s', '10000s'});

set(gca,'fontsize',8);
crop_margins();

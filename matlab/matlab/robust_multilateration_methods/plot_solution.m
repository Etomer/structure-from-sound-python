function evalres = plot_solution(asol,gt);

evalres = evaluate_solution(asol,gt);

T = 7791842/96000;
tn = size(asol.s,2);
tt = T*(0:(tn-1))/(tn-1);

clf; 
subplot(3,1,1);
hold off;
plot(tt(find(evalres.sok)),evalres.s_new(:,find(evalres.sok)),'*');
hold on;
plot(tt,gt.sgt_resamp);
subplot(3,1,2);
hold off;
plot(tt,evalres.s_new,'*');
hold on;
plot(tt,gt.sgt_resamp);
axis([0 90 -5 5]);
subplot(3,1,3);
plot(evalres.hh,'*');



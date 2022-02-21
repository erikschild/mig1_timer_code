% variance analysis of bar-1 data
% Brown-Forsythe ANOVA test
% start time based on each strain's rightmost position datapoint

clear all

dd = './';
M1 = load([dd 'fig5_data_ctrl.txt']);
M2 = load([dd 'fig5_data_gof.txt']);
M3 = load([dd 'fig5_data_lof.txt']);

% Control
for i = 1:5
    ind{i} = find(M1(:,3) == i);
end
p1b = M1(ind{2},2);
p1c = M1(ind{3},2);
p1 = [p1b; p1c];
p1max = max(p1);
t1b = p1max-p1b;
t1c = p1max-p1c;
t1 = p1max-p1;
m1b = M1(ind{2},1);
m1c = M1(ind{3},1);
m1 = [m1b; m1c];

% Gain-of-function
for i = 1:5
    ind{i} = find(M2(:,3) == i);
end
p2b = M2(ind{2},2);
p2c = M2(ind{3},2);
p2 = [p2b; p2c];
p2max = max(p2);
t2b = p2max-p2b;
t2c = p2max-p2c;
t2 = p2max-p2;
m2b = M2(ind{2},1);
m2c = M2(ind{3},1);
m2 = [m2b; m2c];

% Loss-of-function
for i = 1:5
    ind{i} = find(M3(:,3) == i);
end
p3b = M3(ind{2},2);
p3c = M3(ind{3},2);
p3 = [p3b; p3c];
p3max = max(p3);
t3b = p3max-p3b;
t3c = p3max-p3c;
t3 = p3max-p3;
m3b = M3(ind{2},1);
m3c = M3(ind{3},1);
m3 = [m3b; m3c];

% Variance
mmin = 10;
mmax = 25;
mthresh = linspace(mmin,mmax,16);

for k = 1:length(mthresh)
    i1 = find(m1 >= mthresh(k));
    i2 = find(m2 >= mthresh(k));
    i3 = find(m3 >= mthresh(k));
    p1_ = p1(i1);
    p2_ = p2(i2);
    p3_ = p3(i3);
    t1_ = t1(i1);
    t2_ = t2(i2);
    t3_ = t3(i3);
    c1 = ones(length(i1),1);
    c2 = 2*ones(length(i2),1);
    c3 = 3*ones(length(i3),1);
    pCGp(k) = vartestn([p1_;p2_],[c1;c2],...
        'testtype','brownforsythe','display','off');
    pCLp(k) = vartestn([p1_;p3_],[c1;c3],...
        'testtype','brownforsythe','display','off');
    pCGt(k) = vartestn([t1_;t2_],[c1;c2],...
        'testtype','brownforsythe','display','off');
    pCLt(k) = vartestn([t1_;t3_],[c1;c3],...
        'testtype','brownforsythe','display','off');
    p3p(k) = vartestn([p1_;p2_;p3_],[c1;c2;c3],...
        'testtype','brownforsythe','display','off');
    p3t(k) = vartestn([t1_;t2_;t3_],[c1;c2;c3],...
        'testtype','brownforsythe','display','off');
end

figure(1); clf
h = plot([mmin mmax],.05*[1 1],'k:',...
    mthresh,pCGp,'b.-',mthresh,pCGt,'bo',mthresh,2*pCGp,'b--',...
    mthresh,pCLp,'r.-',mthresh,pCLt,'ro',mthresh,2*pCLp,'r--',...
    mthresh,p3p,'g.-',mthresh,p3t,'go',...
    'linewidth',2,'markersize',10);
xlabel('minimum mig-1 number')
ylabel('Brown-Forsythe p value')
legend(h(2:end),{'Control-GOF, position','Control-GOF, time',...
    'Bonferoni-corrected (x2)','Control-LOF, position',...
    'Control-LOF, time','Bonferoni-corrected (x2)',...
    'All 3, position','All 3, time'},...
    'location','nw')
set(gca,'fontsize',16)

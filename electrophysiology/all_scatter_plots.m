%% scatter plots Fig. 3A, 3E, 4G and supplementary figure A,F,K
% run this code after rinning the 1st section of the noise correlation code
pp=[];
figure
for kp=1:3
    csi=nanmean(all_m2{1,1}(:,:,kp),2);% all_m1 for main figures, all_m2 for supplementary
    cci=nanmean(all_m2{1,2}(:,:,kp),2);% all_m1 for main figures, all_m2 for supplementary

    ic=[csi cci];
    fn=find(isnan(sum(ic,2))==1);
    ic(fn,:)=[];
    sn=sign(ic);
    subplot(1,3,kp)
    plot([-1.05 1.05],[0 0],'k');
    hold on
    plot([0 0],[-1.05 1.05],'k');
    hold on
    fns=intersect(find(sn(:,1)==1),find(sn(:,2)==1));
    scatter(ic(fns,1),ic(fns,2),10,'r','filled')
    hold on
    fns=intersect(find(sn(:,1)==-1),find(sn(:,2)==1));
    scatter(ic(fns,1),ic(fns,2),10,'g','filled')
    hold on
    fns=intersect(find(sn(:,1)==-1),find(sn(:,2)==-1));
    scatter(ic(fns,1),ic(fns,2),10,'b','filled')
    hold on
    fns=intersect(find(sn(:,1)==1),find(sn(:,2)==-1));
    scatter(ic(fns,1),ic(fns,2),10,'m','filled')
    [b,~,~,~,stats]=regress(ic(:,2),[ones(length(ic(:,2)),1) ic(:,1)]);
    hold on
    plot([-1:0.1:1],[ones(length([-1:0.1:1]),1) [-1:0.1:1]']*b,'k')
    axis image
    [rho,p]=corrcoef(ic(:,1),ic(:,2))
    pp=[pp p(1,2)];
end


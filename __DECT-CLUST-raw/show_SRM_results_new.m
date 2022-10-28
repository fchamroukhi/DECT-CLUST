function show_SRM_results_new(Curves, mixmodel, mixstats, model)

% plot the results of a EM for spatial Regression Mixtures (polynomial,
% Spline, B-Spline, RM with/out mixed random effects)
% FC

t = Curves.abscissas;
T = 1:length(Curves.abscissas);
Y = Curves.ordinates;
%[n, m] = size(Y);
K = size(mixmodel.Betak, 2);
% original data

% figure, plot(t,Y','linewidth',0.001);  %'r-',
% ylabel('y')
% xlabel('x')
% xlim([min(t) max(t)]);
% %set(gca,'ytick',[0.4:0.2:1.4])
% box on;
% % title([' SRM-EM clustering : iteration ', int2str(length(stored_loglik)), '; K = ', int2str(K)]);
% title(['Original data']);
%
colors = {'r','b','g','m','c','k','y','r--','b--','g--','m--','c--','k--','y--',...
    'r-x','b-x','g-x','m-x','c-x','k-x','y-x', 'r.', 'b.','g.','m.','c.','k.','y.'};
% for k=1:K
%     figure
%     
%     sigmak2 = sqrt(mixmodel.Sigmak2(k));
%     Ic_k = [mixstats.Muk(:,k)-2*sigmak2 mixstats.Muk(:,k)+2*sigmak2];
%     
%     
%     hold on,
%     plot(t,Y(mixstats.klas==k,:)','r-','linewidth',0.001);
%     hold on
%     plot(t,mixstats.Muk(:,k),'k-','linewidth',5);
%     
%     hold on
%     plot(t,Ic_k,'k--','linewidth',1);
%     
%     
%     
%     hold on
%     
%     ylabel('y')
%     xlabel('x')
%     xlim([min(t) max(t)]);
%     title(['SRM-EM cluster', int2str(k)]);%, 'iteration ', int2str(length(stored_loglik)), '; K = ', int2str(K)]);
%     box on
% end
% % ylabel('y')
% % xlabel('x')
% % xlim([min(x) max(x)]);
% %set(gca,'xtick',[0:0.2:1])
% %set(gca,'ytick',[0.4:0.2:1.4])
% box on;
% % title(['SRM-EM clustering : iteration ', int2str(length(stored_loglik)), '; K = ', int2str(K)]);
% %


%%%%%%%
% figure,
% for k=1:K%min(K,7)
%     sigmak2 = sqrt(mixmodel.Sigmak2(k));
%     Ic_k = [mixstats.Muk(:,k)-2*sigmak2 mixstats.Muk(:,k)+2*sigmak2];
%     hold on,
%     %plot(x,X(mixstats.klas==k,:)','r','linewidth',0.001);
%     %hold on
%     plot(t,mixstats.Muk(:,k),colors{k},'linewidth',5);
%     %hold on
%     %plot(x,Ic_k,[colors{k},'--'],'linewidth',1);
% end
% ylabel('y')
% xlabel('x')
% xlim([min(t) max(t)]);
% box on;
% title(['SRM-EM clustering (cluster centers) : iter ', int2str(length(mixstats.stored_loglik)), '; K = ', int2str(K)]);

%%%%%%%
figure('units','normalized','outerposition',[0 0 0.5 0.7]);

if(K==2)
    ha = tight_subplot(1,2, [0.02 0.06]);%,[.05 .05],[.05 .05],[.05 .05]);
elseif (K==3)
    ha = tight_subplot(1,3, [0.02 0.06]);%,[.05 .05],[.05 .05],[.05 .05]);
elseif (K==4)
    ha = tight_subplot(2,2, [0.02 0.06]);%,[.05 .05],[.05 .05],[.05 .05]);
elseif (K==5||K==6)
    ha = tight_subplot(3,2, [0.02 0.06]);%,[.05 .05],[.05 .05],[.05 .05]);
elseif (K==7||K==8)
    ha = tight_subplot(4,2, [0.02 0.06]);%,[.05 .05],[.05 .05],[.05 .05]);
elseif (K==9)
    ha = tight_subplot(3,3, [0.02 0.06]);%,[.05 .05],[.05 .05],[.05 .05]);
elseif (K>9 && K<=12)
    ha = tight_subplot(4,3, [0.02 0.06]);%,[.05 .05],[.05 .05],[.05 .05]);
elseif (K>12 && K<=16)
    ha = tight_subplot(4,4, [0.02 0.06]);%,[.05 .05],[.05 .05],[.05 .05]);
elseif (K>6 && K<=20)
    ha = tight_subplot(4,5, [0.02 0.06]);%,[.05 .05],[.05 .05],[.05 .05]);
elseif (K>20 && K<=25)
    ha = tight_subplot(5,5, [0.02 0.06]);
elseif (K>25 && K<=30)
    ha = tight_subplot(6,5, [0.02 0.06]);
elseif (K>30 && K<=36)
    ha = tight_subplot(6,6, [0.02 0.06]);
elseif (K>36 && K<=42)
    ha = tight_subplot(7,6, [0.02 0.06]);
elseif (K>35 && K<=48)
    ha = tight_subplot(8,6, [0.02 0.06]);
elseif (K>48 && K<=54)
    ha = tight_subplot(9,6, [0.02 0.06]);
elseif (K>54 && K<=63)
    ha = tight_subplot(9,7, [0.02 0.06]);
else
    warning('Please add the plot case in show_SRM_results for this many clusters');
end
%title(['SRM-EM clustering : iteration ', int2str(length(mixstats.stored_loglik)), '; K = ', int2str(K)]);
%ylabel('y')

for k=1:K %min(K,7)
    sigmak2 = sqrt(mixmodel.Sigmak2(k));
    Ic_k = [mixstats.Muk(:,k)-2*sigmak2 mixstats.Muk(:,k)+2*sigmak2];
    axes(ha(k));
    plot(T,Y(mixstats.klas==k,:)','r','linewidth',0.001);
    hold on
    plot(T,mixstats.Muk(:,k),'k','linewidth',5);
    hold on
    plot(T,Ic_k,'k--','linewidth',1);
    
    xlim([min(T) max(T)]);
    ylabel(['Cluster ',num2str(k)],'FontWeight','bold');
%     ylim([min(Y,[],'all') max(Y,[],'all')]);
end
box on;
suptitle(['Curve details per cluster - ', char(model)])

for k = (K+1):length(ha)
    axes(ha(k));
    axis off
end

%%%%%%%
figure
plot(mixstats.stored_loglik,'b-');
xlabel('SRM-EM iteration number');
ylabel('observed data log-likelihood');box on;
suptitle(['Convergence for ',char(model),' method'])

%


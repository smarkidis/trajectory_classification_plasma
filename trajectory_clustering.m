% close plots
close all;

%% Read CSV files
disp('Reading CSV files ...');
M0 = csvread("species0_2800-3100cycles.csv", 1, 0);
M1 = csvread("species0_2600-2900cycles.csv", 1, 0);
M3 = csvread("species0_2400-2700cycles.csv", 1, 0);
M4 = csvread("species0_2200-3000cycles.csv", 1, 0);
%Mo = csvread("10000particles_12by6_2200to3000cycle.csv", 1, 0);
disp('Done with reading');
M = [M0  M1 M3 M4(1:301,:)]; %Mo

%% Simulation Parameters
Lx = 40.0;
Ly = 20.0;
np = 40000 % Number of particles to plot 10000
n = np - 1;
% first electrons exiting after approx 280 cycles
% ions are slow, we can use all the cycles
cycles = 300;


% extract x and y for plotting: need to check if any particle goes out of
% box
x = M(1:cycles,1:3:(1+n*3));
y = M(1:cycles,2:3:(2+n*3));

% disp('Plotting dataset ...');
% figure()
% plot(x,y,'.')
% xlabel('x');
% ylabel('y');
% title('40,000 Electron Trajectories During Magnetic Reconnection');
% grid on
% set(gca, 'FontName', 'Times New Roman')
% set(gca, 'FontSize', 16)

%% Preprocessing data: apply FFT
% X - FFT
x = x - M(1,1:3:(1+n*3));
x = x./max(abs(x));
x = abs(fft(x)/cycles);
%x = abs(fft(x).^2); % not working as good as FFT^2 - compress the
%differences
x = x./max(abs(x));

% Y - FFT
y = y - M(1,2:3:(2+n*3));
y = y./max(abs(y));
y = abs(fft(y)/cycles);
%y = abs(fft(y).^2); % not working as good as FFT^2 - compress the
%differences
y = y./max(y);

% energy: we take only mean as energy information is also encoded in
% the trajectory
energy = M(1:cycles,3:3:(3+n*3));
energy = mean(energy);



%% PCA for dimensionality reduction
num_spectral_modes = cycles/2 % can try just with a part of spectrum: cycles/4, cycles/8

traj = [x(1:num_spectral_modes,:) ; y(1:num_spectral_modes,:)];
[coefs,score,~,~,explained] = pca(traj');


% plot PCA
% figure()
% subplot(2,2,1)
% plot(score(:,1),score(:,2),'.');
% %hold on
% %plot(score(1,1),score(1,2),'r*');
% %hold off
% xlabel('PC I')
% ylabel('PC II')
% grid on
% set(gca, 'FontName', 'Times New Roman')
% set(gca, 'FontSize', 18)



% choose the reduced number of PCA components to use in the analysis
num_pca_components = 20
pca_expl = sum(explained(1:num_pca_components))

%% Clustering
disp('Clustering ...')

% for reproducibility
rng(12.1)
% k-means
num_clusters = 12; % need to be square of a number, e.g., 16, 25, 36 or 64
[idx,C] = kmeans(score(:,1:num_pca_components),num_clusters,'Distance','cosine','Replicates',50,'MaxIter',1000,'Display','final');
% cityblock, cosine, sqeuclidean  


% subplot(2,2,2)
% scatter(score(:,1),score(:,2),num_clusters,idx)
% title('K-means - Cosine Distance')
% xlabel('PC I')
% ylabel('PC II')
% grid on
% set(gca, 'FontName', 'Times New Roman')
% set(gca, 'FontSize', 16)



% subplot(2,2,3)
% pareto(explained)
% set(gca, 'FontName', 'Times New Roman')
% set(gca, 'FontSize', 16)
% xlabel('PC')
% ylabel('Variance Explained (%)')
% grid on
% set(gca, 'FontName', 'Times New Roman')
% set(gca, 'FontSize', 16)

% count the number of clusters per idx
figure()
counts = hist(idx,num_clusters)
xlabel('Cluster ID')
ylabel('# Samples')
grid on
set(gca, 'FontName', 'Times New Roman')
set(gca, 'FontSize', 16)
% sort it out
[out,ii] = sort(counts,'descend')


% calculate how good is clustering
% we use silh2 to eliminate samples that do not belong to cluster
disp('Calculating Silhouette ...')
[silh2,h] = silhouette(score(:,1:num_pca_components),idx,'cosine');
bad_samples_in_clusters = sum(silh2<0)




%% Plot trajectories for different categories
figure('Renderer', 'painters', 'Position', [10 10 1200 700])
%figure()
% get particles position
xp = M(1:cycles,1:3:(1+n*3));
yp = M(1:cycles,2:3:(2+n*3));
% plot
plotted_particles = 25






for i=1:num_clusters
    %subplot(sqrt(num_clusters),sqrt(num_clusters),i)
    subplot(4,3,i)

    
    %%% SELECTING 20 PARTICLES %%%%
    xpi = xp(1:cycles,idx==ii(i) & silh2 > 0.99*max(silh2(idx==ii(i))));
    ypi = yp(1:cycles,idx==ii(i) & silh2 > 0.99*max(silh2(idx==ii(i))));
    a = size(xpi);
    if (a(2) < plotted_particles)
        xpi = xp(1:cycles,idx==ii(i) & silh2 > 0.98*max(silh2(idx==ii(i))));
        ypi = yp(1:cycles,idx==ii(i) & silh2 > 0.98*max(silh2(idx==ii(i))));
    end
    a = size(xpi);
    if (a(2) < plotted_particles)
        xpi = xp(1:cycles,idx==ii(i) & silh2 > 0.97*max(silh2(idx==ii(i))));
        ypi = yp(1:cycles,idx==ii(i) & silh2 > 0.97*max(silh2(idx==ii(i))));
    end
    a = size(xpi);
    if (a(2) < plotted_particles)
        xpi = xp(1:cycles,idx==ii(i) & silh2 > 0.96*max(silh2(idx==ii(i))));
        ypi = yp(1:cycles,idx==ii(i) & silh2 > 0.96*max(silh2(idx==ii(i))));
    end
    a = size(xpi);
    if (a(2) < plotted_particles)
        xpi = xp(1:cycles,idx==ii(i) & silh2 > 0.95*max(silh2(idx==ii(i))));
        ypi = yp(1:cycles,idx==ii(i) & silh2 > 0.95*max(silh2(idx==ii(i))));
    end
    a = size(xpi);
    if (a(2) < plotted_particles)
        xpi = xp(1:cycles,idx==ii(i) & silh2 > 0.94*max(silh2(idx==ii(i))));
        ypi = yp(1:cycles,idx==ii(i) & silh2 > 0.94*max(silh2(idx==ii(i))));
    end
    a = size(xpi);
    if (a(2) < plotted_particles)
        xpi = xp(1:cycles,idx==ii(i) & silh2 > 0.93*max(silh2(idx==ii(i))));
        ypi = yp(1:cycles,idx==ii(i) & silh2 > 0.93*max(silh2(idx==ii(i))));
    end
    a = size(xpi);
    if (a(2) < plotted_particles)
        xpi = xp(1:cycles,idx==ii(i) & silh2 > 0.92*max(silh2(idx==ii(i))));
        ypi = yp(1:cycles,idx==ii(i) & silh2 > 0.92*max(silh2(idx==ii(i))));
    end
    a = size(xpi);
    if (a(2) < plotted_particles)
        xpi = xp(1:cycles,idx==ii(i) & silh2 > 0.91*max(silh2(idx==ii(i))));
        ypi = yp(1:cycles,idx==ii(i) & silh2 > 0.91*max(silh2(idx==ii(i))));
    end
    a = size(xpi);
    if (a(2) < plotted_particles)
        xpi = xp(1:cycles,idx==ii(i) & silh2 > 0.9*max(silh2(idx==ii(i))));
        ypi = yp(1:cycles,idx==ii(i) & silh2 > 0.9*max(silh2(idx==ii(i))));
    end
    a = size(xpi);
    if (a(2) < plotted_particles)
        xpi = xp(1:cycles,idx==ii(i) & silh2 > 0.85*max(silh2(idx==ii(i))));
        ypi = yp(1:cycles,idx==ii(i) & silh2 > 0.85*max(silh2(idx==ii(i))));
    end
    a = size(xpi);
    if (a(2) < plotted_particles)
        xpi = xp(1:cycles,idx==ii(i) & silh2 > 0.8*max(silh2(idx==ii(i))));
        ypi = yp(1:cycles,idx==ii(i) & silh2 > 0.8*max(silh2(idx==ii(i))));
    end
    a = size(xpi);
    if (a(2) < plotted_particles)
        xpi = xp(1:cycles,idx==ii(i) & silh2 > 0.75*max(silh2(idx==ii(i))));
        ypi = yp(1:cycles,idx==ii(i) & silh2 > 0.75*max(silh2(idx==ii(i))));
    end
    a = size(xpi);
    if (a(2) < plotted_particles)
        xpi = xp(1:cycles,idx==ii(i) & silh2 > 0.7*max(silh2(idx==ii(i))));
        ypi = yp(1:cycles,idx==ii(i) & silh2 > 0.7*max(silh2(idx==ii(i))));
    end
    a = size(xpi);
    if (a(2) < plotted_particles)
        xpi = xp(1:cycles,idx==ii(i) & silh2 > 0.65*max(silh2(idx==ii(i))));
        ypi = yp(1:cycles,idx==ii(i) & silh2 > 0.65*max(silh2(idx==ii(i))));
    end
    a = size(xpi);
    if (a(2) < plotted_particles)
        xpi = xp(1:cycles,idx==ii(i) & silh2 > 0.6*max(silh2(idx==ii(i))));
        ypi = yp(1:cycles,idx==ii(i) & silh2 > 0.6*max(silh2(idx==ii(i))));
    end
    a = size(xpi);
    if (a(2) < plotted_particles)
        xpi = xp(1:cycles,idx==ii(i) & silh2 > 0.55*max(silh2(idx==ii(i))));
        ypi = yp(1:cycles,idx==ii(i) & silh2 > 0.55*max(silh2(idx==ii(i))));
    end
    %%% FINISHED SELECTING 20 PARTICLES %%%%
    % Plot sub-panel
    perc = counts(ii)/40000*100;
 
    plot(xpi(:,1:plotted_particles),ypi(:,1:plotted_particles),'.')
    %str_title = sprintf('Class %d',i);
    str_title = sprintf('%d - %3.1f%%',i,perc(i));
    title(str_title)
    xlabel('x (d_i)')
    ylabel('y (d_i)')
    grid on
    axis([12 28 8 12]);
    set(gca, 'FontName', 'Times New Roman')
    set(gca, 'FontSize', 16)
    
    if (i == 5)
        % record a specific class
        xpi4 = xpi(:,1:plotted_particles);
        ypi4 = ypi(:,1:plotted_particles);
    end
    
end


%% Plot some classes in large figures
figure()

xp4 = xp(1:cycles,idx==4 & silh2 > 0.13);
yp4 = yp(1:cycles,idx==4 & silh2 > 0.13);
plot(xp4(:,1:20),yp4(:,1:20),'.')
xlabel('x')
ylabel('y')
grid on
axis([12 28 9 11]);
set(gca, 'FontName', 'Times New Roman')
set(gca, 'FontSize', 18)

figure()

xp3 = xp(1:cycles,idx==3 & silh2 > 0.13);
yp3 = yp(1:cycles,idx==3 & silh2 > 0.13);
plot(xp3(:,1:plotted_particles),yp3(:,1:plotted_particles),'.')
xlabel('x')
ylabel('y')
grid on
axis([12 28 9 11]);
set(gca, 'FontName', 'Times New Roman')
set(gca, 'FontSize', 18)


% figure()
% xp2 = xp(1:cycles,idx==2 & silh2 > 0.13);
% yp2 = yp(1:cycles,idx==2 & silh2 > 0.13);
% plot(xp2(:,1:plotted_particles),yp2(:,1:plotted_particles),'.')
% xlabel('x')
% ylabel('y')
% grid on
% axis([12 28 9 11]);
% set(gca, 'FontName', 'Times New Roman')
% set(gca, 'FontSize', 18)
% 
% 
% figure()
% xp1 = xp(1:cycles,idx==1 & silh2 > 0.13);
% yp1 = yp(1:cycles,idx==1 & silh2 > 0.13);
% plot(xp1(:,1:plotted_particles),yp1(:,1:plotted_particles),'.')
% xlabel('x')
% ylabel('y')
% grid on
% %axis([12 28 9 11]);
% axis([10 30 8 12]);
% set(gca, 'FontName', 'Times New Roman')
% set(gca, 'FontSize', 18)

% %% Calculate Clustering with Gaussian Mixture
% disp('Calculating Clustering with Gaussian Mixture ...')
% options = statset('Display','final');
% gm = fitgmdist(score(:,1:num_pca_components),12,'Options',options)
% idx = cluster(gm,score(:,1:num_pca_components));
% % calculate silhouette
% [silh2,h] = silhouette(score(:,1:num_pca_components),idx,'cosine');
% bad_samples_in_clusters = sum(silh2<0)
% 
% 
% % plot
% figure('Renderer', 'painters', 'Position', [10 10 1200 700])
% plotted_particles = 25
% for i=1:12
%     subplot(4,3,i)
%     xpi = xp(1:cycles,idx==i);
%     ypi = yp(1:cycles,idx==i);
%     plot(xpi(:,1:plotted_particles),ypi(:,1:plotted_particles),'.')
%     xlabel('x (d_i)')
%     ylabel('y (d_i)')
%     grid on
%     axis([12 28 9 11]);
% end
% set(gca, 'FontName', 'Times New Roman')
% set(gca, 'FontSize', 18)


% %% Calculate Clustering with Affinity Propagation: This is very slow
% [index,nclusters] = ClusterPointsAffinity(score(:,1:num_pca_components),'cosine',50);
% plotted_particles = 25
% nclusters
% for i=1:nclusters
%     subplot(sqrt(nclusters),sqrt(nclusters),i)
%     xpi = xp(1:cycles,idx==i);
%     ypi = yp(1:cycles,idx==i);
%     plot(xpi(:,1:plotted_particles),ypi(:,1:plotted_particles),'.')
%     xlabel('x (d_i)')
%     ylabel('y (d_i)')
%     grid on
%     axis([12 28 9 11]);
% end



%% T-SNE: very slow
% Y = tsne(traj','Standardize',true,'Perplexity',50);
% figure()
% gscatter(Y(:,1),Y(:,2))
% grid on
% title('t-SNE')
% set(gca,'fontsize', 18)

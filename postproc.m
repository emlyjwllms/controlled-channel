%clear all
clc

% A0 - uncontrolled
% A1 - controlled

filenames0 = string({dir('data/uncontrolled/*').name});
filenames1 = string({dir('data/controlled/*').name});

filenames0 = filenames0(3:end);
filenames1 = filenames1(3:end);
j=1;
b=1;
zz=1;
% for zz = 1:33

for f = 1:length(filenames0)
    
    [x,y,z,xm,ym,zm,U0,V0,W0,P0,nu_t] = read_field("data/uncontrolled/" + filenames0(f));
    
    x = x(1:2:end);
    y = y(1:2:end);
    z = z(1:2:end);
    U0 = U0(1:2:end,2:2:end,2:2:end);
    V0 = V0(2:2:end,1:2:end,2:2:end);
    W0 = W0(2:2:end,2:2:end,1:2:end);
    
    
    U0mean(:,f) = mean(U0(:,:,zz),1);
    V0mean(:,f) = mean(V0(:,:,zz),1);
    W0mean(:,f) = mean(W0(:,:,zz),1);
    
        
    % stack streamwise slices (x)
    
    % y = 0.2240
    %i = 13;
    for i = 1:length(x)
        U_0(:,:,j) = [ squeeze(U0(i,:,:))' squeeze(V0(i,:,:))' squeeze(W0(i,:,:))' ];       
        j = j+1;
    end
        
    
        
end

Umean = mean(U0mean,2);
Vmean = mean(V0mean,2);
Wmean = mean(W0mean,2);


U_0_hat = fft(reshape(squeeze(U_0(zz,:,:)),[length(y),3,j-1]));
U_0_save = reshape(squeeze(U_0(zz,:,:)),[length(y),3,j-1]);

flops
[Mu0(:,:,:,zz),Cu0,Eu0] = compute_POD(U_0_hat,2);
flops

%save('Mu0.mat','Mu0','z','y','x','L')
% end

%%
clear all
clc

load Mu0
load Mu1_xz

Mu0_1 = (squeeze(real(ifft(Mu0(:,1,1,:)))));
Mu0_2 = (squeeze(real(ifft(Mu0(:,1,2,:)))));

Mu1_1 = (squeeze(real(ifft(Mu1(:,1,1,:)))));
Mu1_2 = (squeeze(real(ifft(Mu1(:,1,2,:)))));

%%
set(0, 'defaultAxesTickLabelInterpreter','latex');

figure(1)
plot(Mu0_1(:,16)./max(Mu0_1(:,16)),y,'LineWidth',3)
xlabel('$u/u_c$','interpreter','latex')
ylabel('$y/\delta$','interpreter','latex')
xlim([0,1.1])
ylim([0,1])
yticks([0.2,0.4,0.6,0.8])
set(gcf, 'PaperPositionMode', 'manual')
set(gcf, 'Color', [1 1 1])
set(gca, 'Color', [1 1 1])
set(gcf, 'Position', [0 0 600 350]);
set(gca,'FontSize',20)
set(gca,'FontName','Times New Roman')
exportgraphics(gcf,'A0_POD1_z_pi_4.png','Resolution',300)

figure(2)
plot(Mu0_2(:,16),y,'LineWidth',3)
xlabel("$u'$",'interpreter','latex')
ylabel('$y/\delta$','interpreter','latex')
ylim([0,1])
yticks([0.2,0.4,0.6,0.8])
xticks([-0.01,0,0.01,0.02,0.03,0.04,0.05])

set(gcf, 'PaperPositionMode', 'manual')
set(gcf, 'Color', [1 1 1])
set(gca, 'Color', [1 1 1])
set(gcf, 'Position', [0 0 600 350]);
set(gca,'FontSize',20)
set(gca,'FontName','Times New Roman')
exportgraphics(gcf,'A0_POD2_z_pi_4.png','Resolution',300)
%%
figure(5)
% first mode should be mean velocity profile
contourf(z,y,Mu0_1./max(Mu0_1)), hold on
%xlim([-0.2,1.2])
hcb = colorbar('northoutside'); % u/U_0
colorTitleHandle = get(hcb,'Title');
% titleString = '$u/u_c$';
% set(colorTitleHandle ,'String',titleString,'interpreter','latex');
caxis([0,1])
yticks([0.5,1,1.5,2])
xlabel('$z/\delta$','interpreter','latex')
ylabel('$y/\delta$','interpreter','latex')
set(hcb,'TickLabelInterpreter','latex','FontSize',16)
%title({'Uncontrolled','POD(1)'},'interpreter','latex')
% title('POD(1)','interpreter','latex')
set(gcf, 'PaperPositionMode', 'manual')
set(gcf, 'Color', [1 1 1])
set(gca, 'Color', [1 1 1])
set(gcf, 'Position', [0 0 600 350]);
set(gca,'FontSize',20)
set(gca,'FontName','Times New Roman')
exportgraphics(gcf,'POD1_A0.png','Resolution',300)

figure(6)
% second mode should be fluctuations
contourf(z,y,Mu0_2)
hcb = colorbar('northoutside'); % u'
colorTitleHandle = get(hcb,'Title');
% titleString = "$u'$";
% set(colorTitleHandle ,'String',titleString,'interpreter','latex');
caxis([-0.06,0.06])
yticks([0.5,1,1.5,2])

xlabel("$z/\delta$",'interpreter','latex')
ylabel('$y/\delta$','interpreter','latex')
set(hcb,'TickLabelInterpreter','latex','FontSize',16)
% title('POD(2)','interpreter','latex')
set(gcf, 'PaperPositionMode', 'manual')
set(gcf, 'Color', [1 1 1])
set(gca, 'Color', [1 1 1])
set(gcf, 'Position', [0 0 600 350]);
set(gca,'FontSize',20)
set(gca,'FontName','Times New Roman')
exportgraphics(gcf,'POD2_A0.png','Resolution',300)

%%
figure(7)
contourf(x,z,permute(Mu0_2,[1 2]))
%xlim([-0.2,1.2])
hcb = colorbar('northoutside');
colorTitleHandle = get(hcb,'Title');
% titleString = "$u'$";
% set(colorTitleHandle ,'String',titleString,'interpreter','latex');
set(hcb,'TickLabelInterpreter','latex','FontSize',16)
% caxis([-0.06,0.06])
xlabel('$x/\delta$','interpreter','latex')
ylabel('$z/\delta$','interpreter','latex')
yticks([0.2,0.4,0.6,0.8,1,1.2,1.4])
set(gcf, 'PaperPositionMode', 'manual')
set(gcf, 'Color', [1 1 1])
set(gca, 'Color', [1 1 1])
set(gcf, 'Position', [0 0 600 350]);
set(gca,'FontSize',20)
set(gca,'FontName','Times New Roman')
%exportgraphics(gcf,'A0_POD2_xz.png','Resolution',300)


figure(8)
contourf(x,z,permute(Mu1_2,[1 2]))
hcb = colorbar('northoutside'); % u'
colorTitleHandle = get(hcb,'Title');
%titleString = "$u'$";
%set(colorTitleHandle,'String',titleString,'interpreter','latex','FontName','Times New Roman');
set(hcb,'TickLabelInterpreter','latex','FontSize',16)
% caxis([-0.06,0.06])
xlabel("$x/\delta$",'interpreter','latex')
ylabel('$z/\delta$','interpreter','latex')
yticks([0.2,0.4,0.6,0.8,1,1.2,1.4])
%title('POD(2)','interpreter','latex')
set(gcf, 'PaperPositionMode', 'manual')
set(gcf, 'Color', [1 1 1])
set(gca, 'Color', [1 1 1])
set(gcf, 'Position', [0 0 600 350]);
set(gca,'FontSize',20)
set(gca,'FontName','Times New Roman')
%exportgraphics(gcf,'A1_POD2_xz.png')




%%

% DMD

clear all
close all
clc

% A0 - uncontrolled
% A1 - controlled

filenames0 = string({dir('data/uncontrolled/*').name});
filenames1 = string({dir('data/controlled/*').name});

filenames0 = filenames0(3:end);
filenames1 = filenames1(3:end);
j=1;
b=1;
% for zz = 1:33

for f = 1:length(filenames1)
    
    [x,y,z,xm,ym,zm,U1,V0,W0,P0,nu_t] = read_field("data/controlled/" + filenames1(f));
    
    x = x(1:2:end);
    y = y(1:2:end);
    z = z(1:2:end);
    U1 = U1(1:2:end,2:2:end,2:2:end);
    
%     
%     U0mean(:,f) = mean(U0(:,:,zz),1);
%     V0mean(:,f) = mean(V0(:,:,zz),1);
%     W0mean(:,f) = mean(W0(:,:,zz),1);
    
        
    % stack streamwise slices (x)
    
    %y=0.2
    i = 13;
    %for i = 1:length(x)
        %U_0(:,:,f) = [ squeeze(U0(i,:,:))' squeeze(V0(i,:,:))' squeeze(W0(i,:,:))' ];     
        U_1(:,:,f) = [ squeeze(U1(:,i,:)) ];
    %end
        

end

for f = 1:length(filenames0)
    
    [x,y,z,xm,ym,zm,U0,V0,W0,P0,nu_t] = read_field("data/uncontrolled/" + filenames0(f));
    
    x = x(1:2:end);
    y = y(1:2:end);
    z = z(1:2:end);
    U0 = U0(1:2:end,2:2:end,2:2:end);
    
%     
%     U0mean(:,f) = mean(U0(:,:,zz),1);
%     V0mean(:,f) = mean(V0(:,:,zz),1);
%     W0mean(:,f) = mean(W0(:,:,zz),1);
    
        
    % stack streamwise slices (x)
    
    %y=0.2
    i = 13;
    %for i = 1:length(x)
        %U_0(:,:,f) = [ squeeze(U0(i,:,:))' squeeze(V0(i,:,:))' squeeze(W0(i,:,:))' ];     
        U_0(:,:,f) = [ squeeze(U0(:,i,:)) ];
    %end
        

end

U_0 = reshape(U_0,[33*33,f]);
U_1 = reshape(U_1,[33*33,f]);


[~,Phi0] = compute_DMD(U_0(:,1:end-1),U_0(:,2:end),f-1);
phi0 = reshape(Phi0,[33,33,f-1]);

[~,Phi1] = compute_DMD(U_1(:,1:end-1),U_1(:,2:end),f-1);
phi1 = reshape(Phi1,[33,33,f-1]);




figure(8)
subplot(1,2,1)

contourf(x,z,real(phi0(:,:,1))./max(real(phi0(:,:,1))))
%xlim([-0.2,1.2])
hcb = colorbar('northoutside');
colorTitleHandle = get(hcb,'Title');
% titleString = "$u'$";
% set(colorTitleHandle ,'String',titleString,'interpreter','latex');
set(hcb,'TickLabelInterpreter','latex','FontSize',16)
% caxis([-0.06,0.06])
xlabel('$x/\delta$','interpreter','latex')
ylabel('$z/\delta$','interpreter','latex')
yticks([0.2,0.4,0.6,0.8,1,1.2,1.4])
set(gcf, 'PaperPositionMode', 'manual')
set(gcf, 'Color', [1 1 1])
set(gca, 'Color', [1 1 1])
set(gcf, 'Position', [0 0 1200 350]);
set(gca,'FontSize',20)
set(gca,'FontName','Times New Roman')
%exportgraphics(gcf,'A0_POD2_xz.png','Resolution',300)


subplot(1,2,2)
contourf(x,z,abs(abs(real(phi1(:,:,1)))./max(abs(real(phi1(:,:,1)))))) % controlled
hcb = colorbar('northoutside'); % u'
colorTitleHandle = get(hcb,'Title');
%titleString = "$u'$";
%set(colorTitleHandle,'String',titleString,'interpreter','latex','FontName','Times New Roman');
set(hcb,'TickLabelInterpreter','latex','FontSize',16)
% caxis([-0.06,0.06])
xlabel("$x/\delta$",'interpreter','latex')
ylabel('$z/\delta$','interpreter','latex')
yticks([0.2,0.4,0.6,0.8,1,1.2,1.4])
%title('POD(2)','interpreter','latex')
set(gcf, 'PaperPositionMode', 'manual')
set(gcf, 'Color', [1 1 1])
set(gca, 'Color', [1 1 1])
set(gcf, 'Position', [0 0 1200 350]);
set(gca,'FontSize',20)
set(gca,'FontName','Times New Roman')
exportgraphics(gcf,'DMD1.png')

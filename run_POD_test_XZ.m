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
zz=33;
% for zz = 1:33 -- you can put this in a loop, but it takes a long time for
% some reason

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
    i = 13;
    %for i = 1:length(x)
        U_0(:,:,j) = [ squeeze(U0(:,i,:))' squeeze(V0(:,i,:))' squeeze(W0(:,i,:))' ];       
        j = j+1;
    %end
        
    
        
end

Umean = mean(U0mean,2);
Vmean = mean(V0mean,2);
Wmean = mean(W0mean,2);


U_0_hat = fft(reshape(squeeze(U_0(zz,:,:)),[length(x),3,j-1]));
U_0_save = reshape(squeeze(U_0(zz,:,:)),[length(x),3,j-1]);
[Mu0(:,:,:,zz),Cu0,Eu0] = compute_POD(U_0_hat,2);

%%

Mu0_1 = (squeeze(real(ifft(Mu0(:,1,1,:)))));
Mu0_2 = (squeeze(real(ifft(Mu0(:,1,2,:)))));


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
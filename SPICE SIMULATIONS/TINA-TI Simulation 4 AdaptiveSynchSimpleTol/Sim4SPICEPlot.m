function Sim4SPICEPlot(raw_data_txt)

filename = raw_data_txt(1:end-4);

try
    load(strcat(filename,'.mat'))
catch err
    data_struct = importdata(raw_data_txt);
    data = data_struct.data;
    save(filename,'data');
end  %end of try catch


p = path;
path(p,strcat(pwd,'\MATLAB res\SubAxis'));
c = onCleanup(@()path(p));

%*Time 	vhq	vq	Vc	vh1	vh2	v1	v2

t = data(:,1);
iLh = data(:,2)/1000;  % VQh =  -500ih_L
iL = data(:,3)/1000;   % VQ =  -500i_L
Vc = data(:,4);
v1h = data(:,5);
v2h = data(:,6);
v1 = data(:,7);
v2 = data(:,8);
e_v1 = v1h - v1;
e_v2 = v2h - v2;
e_iL = iLh-iL;

clear data

twoNorm = sqrt(sum(abs([e_v1 e_v2 e_iL]').^2,1)); 

figure(1)
subaxis(1,1,1,'ML',0.10,'MR',0.02,'MT',0.02,'MB',0.11,'SpacingVert',0.06)
plot(iL(18000:end),v2(18000:end))
xlabel('$i_L$','interpreter','latex','fontsize',20)
ylabel('$v_2$','interpreter','latex','fontsize',20)
grid;
axis ([-.004 .004 -1 1])

figure(2)
subaxis(2,1,1,'ML',0.14,'MR',0.05,'MT',0.02,'MB',0.09,'SpacingVert',0.06)
plot(t,Vc)
grid
%xlabel('t, (seconds)')
ylabel('$V_c$ (V)','interpreter','latex','fontsize',20)
xlim([.045 0.085])
yL = get(gca,'YLim');
line([0.06 0.06],yL,'Color','r','LineStyle','-.');


subaxis(2,1,2)
plot(t, twoNorm)
xlabel('t (seconds)')
ylabel('$ e_{\mathrm{norm}}$','interpreter','latex','fontsize',20)
grid
xlim([.045 0.085])
yL = get(gca,'YLim');
line([0.06 0.06],yL,'Color','r','LineStyle','-.');

window = 50;

try
    load(strcat(filename,'LinearFit'));
catch err


% Line Fitting  

iLvec = iL;   
v2vec = v2;
iLhvec = iLh;  
v2hvec = v2h;

[m,n] = size(iLvec);

Lvec = zeros(1,m-window);
Lhvec = zeros(1,m-window);
Rvec = zeros(1,m-window);
Rhvec = zeros(1,m-window);
tvec = t(1:end)';

    for i = 2:(m-window)
        iLdotH1 = diff(iLhvec(i-1:(i-1+window)))./diff(tvec(i-1:(i-1+window))');
        iLdotH2 = diff(iLhvec(i:(i+window)))./diff(tvec(i:(i+window))');
        iLdotH = (iLdotH1 + iLdotH2)/2;
        pH = pinv([iLhvec(i:(i+window-1)) v2hvec(i:(i+window-1))])*iLdotH;
        Lhvec(i) = -1/pH(2);
        Rhvec(i) = pH(1)/pH(2);
    end

    for i = 2:(m-window)
        iLdot1 = diff(iLvec(i-1:(i-1+window)))./diff(tvec(i-1:(i-1+window))');
        iLdot2 = diff(iLvec(i:(i+window)))./diff(tvec(i:(i+window))');
        iLdot = (iLdot1 + iLdot2)/2;
        pH = pinv([iLvec(i:(i+window-1)) v2vec(i:(i+window-1))])*iLdot;
        Lvec(i) = -1/pH(2);
        Rvec(i) = pH(1)/pH(2);
    end
  
 save(strcat(filename,'LinearFit'),'Lvec', 'Rvec', 'Lhvec', 'Rhvec', 'window')   

%***************End For line Fitting B
end %end try catch


%%Plot Linear Fit
figure(3)
subaxis(2,1,1,'ML',0.14,'MR',0.05,'MT',0.02,'MB',0.09,'SpacingVert',0.07)
plot(t(1:end-window),Rhvec,'--',t(1:end-window),Rvec);
ylabel('Resistance ($\Omega$)','interpreter','latex','fontsize',20)
title(['window = ', num2str(window), ' , experimental linear fit w/ 3 CFOA' ])
l=legend('$\widetilde R_{0\mathrm{est}}$','$R_{0\mathrm{est}}$');
set(l,'interpreter','latex') 
grid
xlim([.045 0.085])
ylim([50 130])
subaxis(2,1,2)
plot(t(1:end-window),Lhvec,'--',t(1:end-window),Lvec);
xlabel('t (seconds)')
ylabel('Inductance (H)','interpreter','latex','fontsize',20)
l=legend('$\widetilde L_{\mathrm{est}}$','$L_{\mathrm{est}}$');
set(l,'interpreter','latex') 
grid
xlim([.045 0.085])

%%=====Zoomed****images

I1 = find(t> 0.075,1,'first');
I2 = find(t> 0.085,1,'first');
figure(4)
subaxis(2,1,1,'ML',0.14,'MR',0.05,'MT',0.02,'MB',0.09,'SpacingVert',0.1)
%subplot(211)
plot(t(I1:I2),Vc(I1:I2))
xlim([-inf inf])
grid
%xlabel('t, (seconds)')
ylabel('$V_c$ (V)','interpreter','latex','fontsize',20)

subaxis(2,1,2)
%subplot(212)
plot(t(I1:I2),twoNorm(I1:I2))
xlim([-inf inf])
xlabel('t (seconds)')
ylabel('$ e_{\mathrm{norm}}$','interpreter','latex','fontsize',20)
grid


%%Plot Linear Fit
figure(5)
subaxis(2,1,1,'ML',0.14,'MR',0.05,'MT',0.02,'MB',0.09,'SpacingVert',0.07)
plot(t(I1:I2),Rhvec(I1:I2),'--',t(I1:I2),Rvec(I1:I2));
ylabel('Resistance ($\Omega$)','interpreter','latex','fontsize',20)
xlim([t(I1) t(I2)])
l=legend('$\widetilde R_{0\mathrm{est}}$','$R_{0\mathrm{est}}$');
set(l,'interpreter','latex') 
grid
subaxis(2,1,2)
plot(t(I1:I2),Lhvec(I1:I2),'--',t(I1:I2),Lvec(I1:I2));
xlim([t(I1) t(I2)])
xlabel('t (seconds)')
ylabel('Inductance (H)','interpreter','latex','fontsize',20)
l=legend('$\widetilde L_{\mathrm{est}}$','$L_{\mathrm{est}}$');
set(l,'interpreter','latex') 
grid



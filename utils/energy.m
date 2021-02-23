clear all;
% Simulation parameters ----------------
wpiwci = 6000;
betai = 1;
betae= 0.5;
gamma = 5/3;
dtwci = 0.025;
B0 = 1/wpiwci;
VA = 1/wpiwci;
Ti0 = 0.5*betai/wpiwci^2;
Te0 = 0.5*betae/wpiwci^2;
nct = 256^2*768;
% --------------------------------------

% Pozitive cross-helicity -------------
nvars = 9;

f=fopen("../data/energy.dat","r");
en = fread(f,"double");
m = length(en)/nvars;

en = reshape(en,[nvars m]);
it = [1:m];

R3.update = 20;
R3.it = it;
% total energy
R3.et  = sum(en(1:8,:),1)/nct;
% delta E = (E-E_0)/E_0
R3.det = (R3.et - R3.et(1))/R3.et(1);
% magnetic energy per unit volume normalized to B0^2
R3.b2  = (en(1,:)+en(2,:)+en(3,:))/B0^2/nct - 0.5;
% flow energy per unit volume normalized to VA^2
R3.v2  = (en(4,:)+en(5,:)+en(6,:))/VA^2/nct;
% electron thermal energy normalized to Te0 per unit volume
R3.pe  = en(7,:)/Te0/nct;
% electron thermal energy normalized to Te0 per unit volume
R3.pi  = en(8,:)/nct;
R3.dpi = (R3.pi - R3.pi(1))/R3.et(1);
% time
R3.t = R3.update*R3.it*dtwci;


% Plot data ---------------------------

lw = 5;

% Figure 1: ion thermal energy
figure(1)
clf; hold on;
plot(R3.t,R3.pi,'-r','linewidth',lw)
plot(R3.t,R3.pe,'-b','linewidth',lw)
xlabel('t{\Omega_{ci}}')
legend('Tr(P_i)/(2 n_0 T_{i0})','p_e/(\gamma-1)/(n_0 T_{e0})')
print -dpdf -F:14 "p_energy.pdf"

				% Figure 2: magnetic field
figure(2)
clf; hold on;
%plot(R1.t,R1.b2,'-r','linewidth',lw)
%plot(R4.t,R4.b2,'-r','linewidth',lw)
%plot(R2.t,R2.b2,'--b','linewidth',lw)
plot(R3.t,R3.b2,'-r','linewidth',lw)
%legend('H_c=0','H_c<0','H_c>0')
xlabel('t{\Omega_{ci}}')
ylabel('B^2/2')
print -dpdf -F:14 "b2_energy.pdf"

% Figure 3: total energy 
figure(3)
clf; hold on;
%plot(R4.t,R4.det,'-r','linewidth',lw)
%plot(R2.t,R2.det,'--b','linewidth',lw)
plot(R3.t,R3.det,'-r','linewidth',lw)
plot(R3.t,R3.dpi,'-b','linewidth',lw)
%plot(R0.t,R0.det,':k','linewidth',lw)
%legend('H_c=0','H_c<0','H_c>0','control')
xlabel('t{\Omega_{ci}}')
ylabel('\delta E_{tot}/E(0)')
print -dpdf -F:14 "total_energy.pdf"

% Figure 4: v^2+b^2 
figure(4)
clf; hold on;
%plot(R4.t,R4.b2+R4.v2,'-r','linewidth',lw)
%plot(R2.t,R2.b2+R2.v2,'--b','linewidth',lw)
plot(R3.t,R3.b2+R3.v2,'-r','linewidth',lw)
%legend('H_c=0','H_c<0','H_c>0')
xlabel('t{\Omega_{ci}}')
ylabel('(n V^2+B^2)/2')
print -dpdf -F:14 "decaying_energy.pdf"

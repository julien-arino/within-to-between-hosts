using DifferentialEquations, Plots

function Immune!(du,u,h,p,t)
    ##
    hist1 = h(p, t-tau)[1] #V
    hist2 = h(p, t-tau)[2] #S
    hist3 = h(p, t-tau)[3] #I
    hist4 = h(p, t-tau)[4] #R
    ###
    du[1] = phat*u[3]-dV*u[1];#V
    du[2] = lamS*(1-(u[2]+u[3]+u[4]+u[5])/Smax)*u[2]-beta*u[2]*u[1];#S
    du[3] = (beta/(1+u[7]/eps_F_I))*hist1*hist2*u[8]-dI*u[3];#I
    du[4] = lamS*(1-(u[2]+u[3]+u[4]+u[5])/Smax)*u[4]+(beta*hist1*hist2)*u[9]/(1+eps_F_I/u[7]);#R
    du[5] = dI*u[3]-dD*u[5];#D
    du[6] = psi_F_prod+p_F_I*u[3]/(u[3]+eta_F_I)-k_lin_F*u[6]-k_B_F*((T_star+u[3])*A_F-u[7])*u[6]+k_U_F*u[7];#FU
    du[7] = -k_int_F*u[7]+k_B_F*((T_star+u[3])*A_F-u[7])*u[6]-k_U_F*u[7];#FB
    du[8] = u[8]*delta*(hist3-u[3])#A1
    du[9] = u[9]*delta*(hist4-u[4])#A2
end

tspan = (0.0,200);
#####################################
beta = 0.3;
eta_F_I = 0.022328;
tau = 0.17;
Smax = 0.16;
lamS = 0.74; 
dI = 0.1; 
dD = 8;
phat = 394;
dV = 8.4;
psi_F_prod = 0.25;
k_lin_F = 16.635;   
k_B_F = 0.0107;
T_star = 1.104*1e-4;
k_U_F = 6.072;
k_int_F = 16.968;
R_F_T = 1000;
R_F_I = 1300;
MM_F = 19000;
avo = 6.02214E23;
A_F = (MM_F/avo)*(R_F_I+R_F_T)*(1/5000)*(10^9*1e12);
delta = 0.1;
lags = [tau];
u0 = [1; 0.16;0.0 ;0.0 ; 0.0;0.015 ;1.1e-8 ;1.0 ;1.0];
#####################################################################
#1
p_F_I = 2.8235*1e4*1e-4;
eps_F_I = 2*1e-4;
###
p = [tau,phat,dV,beta,Smax,lamS,dI,eps_F_I,dD,psi_F_prod, p_F_I,eta_F_I,k_lin_F,k_B_F,T_star,k_U_F,k_int_F,A_F,k_B_F,delta];
h(p,t) = [1; 0.16; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0]
z = DDEProblem(Immune!,u0,h,tspan,p;constant_lags=lags);
alg = MethodOfSteps(Tsit5())
sol = solve(z,saveat=0:0.01:3,alg,abstol=1e-8,reltol=1e-8);
###########################################
#V,S,I,R,D,Fu,Fb
V0 = [(t, u[1]) for (u,t) in tuples(sol)]
V01 = [(u[7], u[1]) for (u,t) in tuples(sol)]
V02 = [(u[7], 39581.1*u[7]+0.424384) for (u,t) in tuples(sol)]
S0 = [(t, u[2]) for (u,t) in tuples(sol)]
I0 = [(t, u[3]) for (u,t) in tuples(sol)]
R0 = [(t, u[4]) for (u,t) in tuples(sol)]
D0 = [(t, u[5]) for (u,t) in tuples(sol)]
Fu0 = [(t, u[6]) for (u,t) in tuples(sol)]
Fb0 = [(t, u[7]) for (u,t) in tuples(sol)]
betaave = [((u[7], ((u[1])^3.623)/((u[1])^3.623+7.68^3.623))) for (u,t) in tuples(sol)]
# betaave = [(t, (1.2*(u[1]^2))/(u[1]^2+15.21)) for (u,t) in tuples(sol)]
# betaave = [(u[7], (1.2*(u[1]^2))/(u[1]^2+15.21)) for (u,t) in tuples(sol)]
betaave2 = [(u[1], 1.2*(((39581.1*u[7]+0.424384)^2))/((39581.1*u[7]+0.424384)^2+15.21)) for (u,t) in tuples(sol)]
#betas2 = [(u[7],u[1], 1.2*(((55195.90786608686*u[7]+0.07082905855419551)^2))/((55195.90786608686*u[7]+0.07082905855419551)^2+15.21)) for (u,t) in tuples(sol)]

###
using DelimitedFiles
writedlm("Vave.txt", V0)

%%  S. Cervisiae Data
P.pesocel=1e-12*42.65; %g/cel
P.denscel=1.1126; % g/ml
P.vcel= P.pesocel/P.denscel; %ml
P.Vcel=P.vcel/1000; %L
P.ec=0.01; %mM %% Assumed concentration for internal enzimes
P.tazaacetato_mut=0.031878465;

%% Datos volumétricos de compartimentos 
P.Vtot=0.05; %L
P.biomass0=2.2; %g/L
P.Vcels=P.biomass0*P.Vtot/(P.denscel*1000);
P.vp=1/10*P.Vcels; %L
P.vm=P.vp; %L

%% Time data
Tend=60*40;
P.time=[0:1:Tend]; %min

%% Kinetic parameters

%%% ADH : Ethanol+ NAD <--> Acetaldehído(aca) +NADH 
%%% Etoh , NAD , aca , NADH

P.ADH_Vm=100; %mmol/min
P.ADH_Km_aca=1.5;
P.ADH_Km_NADH=0.122;
P.ADH_Km_NAD=0.059;
P.ADH_Keq=0.0001;
P.ADH_Km_etoh=7.75;

%%% ALDH : NAD+Aca -> Acetate +NADH
%%% NAD , aca , Acetate , NADH

P.ALDH_Vm=300; %mmol/min
P.ALDH_Kia=1;
P.ALDH_Kb=8;
P.ALDH_Ka=0.005;
P.ALDH_Ki_Ald=3.5;
P.ALDH_Ki_Ald_NADH=1;
P.ALDH_Ki_NADH=0.1;

%%% Diffusion between Citoplasm and perox
%%%% Acetate , Acetate(p)

P.Dp=1;

%%%% Acetyl CoA Syntases  Acetate + ATP + CoA = AMP + acCoA + Bifosfato
%%% Peroxisomal

Z.V1=1;
Z.V2=1;
Z.Kia=15;
Z.Kib=1;
Z.Kic=1;
Z.Kip=15;
Z.Kiq=2.7;
Z.Ka=0.35;
Z.Kb=0.65;
Z.Kc=0.16;
Z.Kr=1;
Z.Kir=1;
Z.Keq=1;
Z.Kq=1.8;
Z.Kp=1;
%%% Citoplasmic
W.V1=100;
W.V2=10;
W.Kia=15;
W.Kib=1;
W.Kic=1;
W.Kip=15;
W.Kiq=2.7;
W.Ka=0.35;
W.Kb=0.65;
W.Kc=0.16;
W.Kr=1;
W.Kir=1;
W.Keq=1;
W.Kq=10;
W.Kp=1;

%%%Difusión exterior citoplasm
P.Difex=0.02;
P.kd=0.165;
P.n=200;
P.citcrit=0.174;
P.Acet_coeff=0.05;
P.Acet_coeff_2=2;
P.Acet_coeff_3=4;
P.Acet_coeff_4=0.2;
P.Acet_coeff_5=0.001; 

%V_Acet_diffx=-P.Difex*(Acetatex-P.Acet_coeff*Acetate)*P.Acet_coeff_2*Acetate/(P.Acet_coeff_3+P.Acet_coeff_4*Acetatex)*(1-P.f)+(P.f)*P.kd*Acetate/(P.Acet_coeff_5+Acetate); 

%%%  ATP exchange 

P.ATP_V=0.1; %mmol/min
P.ATP_Km=5; %mM

%%%  AMP exchange

P.AMP_V=0.1; %mmol/min
P.AMP_Km=5.5; %mM

%%% Citrate Syntase perox

P.CS_V=300; %mM/min
P.CS_Km_oaa=0.001; %mM
P.CS_Km_acCoA=0.0045; %mM

%%% Export Citrate p->Ctirate intra

P.Cexp_Vm=50; %mmol/min
P.Cexp_Km=0.011; %mM

%%% Aconitase

P.AC_V=20; %mM/min
P.AC_Km_cit=2.7; %mM
P.AC_Km_iscit=2.5; %mM
P.AC_Keq=0.1; %mM

%%% Isocitrate lyase

P.iscly_Vf=10; %mM/min
P.iscly_Keq=1;
P.iscly_Kms_iscit=0.05;  %mM
P.iscly_Kip_succ=1.8; %mM
P.iscly_Vr=1; %mM/min
P.iscly_Kmq_glio=0.13; %mM
P.iscly_Kmp_succ=0.59; %mM

%%% Malate Syntase

P.MS_V1=100; %mM/min
P.MS_V2=1; %mM/min
P.MS_Ki_AcCoA=1; %mM
P.MS_Kb_glio=0.1; %mM
P.MS_Ka_AcCoA=0.083; %mM
P.MS_Keq=1;
P.MS_Kq_mal=1;
P.MS_Kp_CoA=1;
P.MS_Kiiq_mal=1;
P.MS_Kiip_CoA=15.15;

%%% Malate deshydrogenase

P.MSH_V=100; %mM/min
P.MSH_Km_mal=0.2; %mM
P.MSH_Km_NAD=0.2; %mM
P.MSH_Km_oaa=0.045; %mM
P.MSH_Km_NADH=0.215; %mM
P.MSH_Keq=1;

%%% Transport: oaa to peroxisome

P.oaat_V=100; %mmol/Min
P.oaat_K1=0.1; %mM
P.oaat_K2=1; %mM

%%% Transport: succinate to mitochondria
P.sucfum_V=100;
P.sucfum_Kia_fumm=1;
P.sucfum_Kb_suc=1;
P.sucfum_Ka_fumm=1;

   %% Cineticas Kummer
   
%%% Transport oxalacetate to mitochondria
% oaa->oaam
P.OAC_Vm=1; %mmol/min
P.OAC_Km_oaa=1; %mM
P.OAC_Km_oaam=0.1; %mM

%%% CAK_Conversion:  ATP, AMP and ADP inside citoplasm
% ATP + AMP <--> 2ADP
P.CAK_Vm=1; %mM/min
P.CAK_Km_ATP=1; %mM
P.CAK_Km_AMP=0.5; %mM
P.CAK_Km_ADP=0.55; %mM
P.CAK_Keq=2.86;
%%% CAK_Conversion ATP, AMP and ADP inside  mitochondria
% ATPm + AMPm <--> 2ADPm
P.CAKm_Vm=1; %mM/min
P.CAKm_Km_ATP=0.5; %mM
P.CAKm_Km_AMP=0.5; %mM
P.CAKm_Km_ADP=0.5; %mM
P.CAKm_Keq=1;
%%% MPC: Piruvate Transport to  mitochondria
% pir<->pirm
P.MPC_Vm=1; %mM/min
P.MPC_Km_pir=0.11; %mM
P.MPC_Km_pirm=1.125; %mM
%%%  Transport:  ATP to mitchondria
%ATPm<-->ATP  
P.ATPtm_Vm=1; %mmol/min
P.ATPtm_Km=5; %mM

%%% Transport: NADH and NAD trough mitochondrial membrane 
P.RES_Vm=1; %mmol/min
P.RES_Km_NADH=1; %mM
P.RES_Km_NADHm=8.5  ;%mM
P.RES_Km_NAD=1   ;%mM
P.RES_Km_NADm=8.5  ;%mM
P.RES_Keq=1; 


%%%  Transport: Atp trough mitochondrial membrane
% ADPm<-->ADP
P.ADPtm_Vm=1; %mmol/min
P.ADPtm_Km=5.5; %mM

%%% PDH: m-pyr + CoA + mNAD → m-acCoA + mNADH
P.PDH_Vm=1; %mM/min
P.PDH_Km_pirm=0.075; %mM
P.PDH_Km_CoAm=0.07; %mM
P.PDH_Km_NADm=0.16; %mM
P.PDH_Ki_AcCoAm=0.025;% %mM
P.PDH_Ki_NADHm=0.03; %mM

%%%MAE : mmal+NADm->pirm +NADHm

P.MAE_Vm=1; %mM/min
P.MAE_Km_malm=25; %mM
P.MAE_Km_NADm=1.05; %mM

%%%MDH(m) : malm + NADm -> oaam + NADHm

P.MDHm_Vm=1; %mM/min
P.MDHm_Km_malm=0.35;%mM
P.MDHm_Km_NADm=0.325; %mM
P.MDHm_Keq=0.000015; %mM
P.MDHm_Km_oaam=0.005; %mM
P.MDHm_Km_NADHm=0.5; %mM

%%%FUM(m) : fumm <-> malm
P.FUMm_Vm=1; %mM/min 
P.FUMm_Km_fumm=0.7; %mM
P.FUMm_Km_malm=0.75; %mM
P.FUMm_Ki_ATPm=2.5; %mM
P.FUMm_Keq=5; %mM

%%%SDH succm <-> fum

P.SDH_Vm=1; %mM/min 
P.SDH_Km_succ=4.5; %mM 
P.SDH_Km_fum=3.2; %mM
P.SDH_Keq=1.65;

%%% SCS succCoAm +ADPm <->succm+CoAm+ATPm
P.SCS_Km_SucCoA=0.32; %mM
P.SCS_Km_ADPm=4.35; %mM
P.SCS_Km_succm=27.5; %mM
P.SCS_Km_CoAm=0.15; %mM
P.SCS_Km_ATPm=2.15; %mM
P.SCS_Keq=55; %mM
P.SCS_Vm=1; %mM/min

%%% AKGDH akgm + CoAm + NADm ->sucCoA+NADHm
P.AK_Vm=1; %mM/min
P.AK_Km_akg=2.1; %mM
P.AK_Km_CoAm=0.62; %mM
P.AK_Km_NADm=1; %mM
P.AK_Ki_NADHm=5; %mM
P.AK_Ka_ADPm=5; %mM

%%% IDH icitm +NADm -> akgm +NADHm

P.IDH_Vm=1; %mM/min
P.IDH_Km_icitm=2.55; %mM
P.IDH_Km_NADm=7.52;  %mM
P.IDH_Ki_NADHm=0.11; %mM
P.IDH_Ki_akgm=10;  %mM

%%% ACO citm<->icitm

P.ACOm_Km_citm=2.65; %mM
P.ACOm_Km_icitm=2.5; %mM
P.ACOm_Keq=0.9; %mM
P.ACOm_Vm=1; %mM

%%% CS oaam + AcCoAm -> citm + CoAm
P.CSm_Vm=1; %mM/min
P.CSm_Km_oaam=0.01; %mM
P.CSm_Km_AcCoAm=0.91; %mM


%% Generación de piruvato

%%% FUMHY : Fumarate hydratase 
% fum<--> malato + H2O
P.FUMHY_Vm=1; %mM/min
P.FUMHY_Km_fum=0.7; %mM
P.FUMHY_Km_mal=75; %mM
P.FUMHY_Keq=5; 
P.FUMHY_Ki_ATP=2.5; %mM

%%% Phosphoenol piruvate carbokinase
% ATP + oaa -> ADP + phosphoenolpiruvato + CO2;

P.PECK_V1=1; 
P.PECK_V2=1; 
P.PECK_Keq=1; 
P.PECK_Ki_ATP=1;
P.PECK_Kb_oaa=0.066;
P.PECK_Ka_ATP=0.034; 
P.PECK_Kq_pep=1; 
P.PECK_Kp_ADP=0.04; 
P.PECK_Kiiq_pep=1; 
P.PECK_Kiip_ADP=1; 

%%% Piruvate kinase
% phosphoenolpiruvato + ADP -> ATP + piruvate

P.PK_Vm=1; %mM/min
P.PK_Km_ADP=0.16; %mM
P.PK_Km_pep=1.18; %mM


%%  Initial concentrations
%%% Citoplasm

P.etoh0=0; %mM
P.NAD0=20.7;%20.7; %mM
P.aca0=0; %mM
P.NADH0=10.33; %mM
P.acetate0=0;%mM
P.ATP0=20; %20; %mM
P.CoA0=10 ; %mM
P.AcetylCoa0=0;  %mM
P.Bif0=0.8;   %mM
P.AMP0=10; %mM

P.cit0=0.1; %mM
P.iscit0=0.0001; %mM
P.succ0=0; %mM
P.glio0=0.0001; %mM
P.mal0=0; %mM
P.oaa0=4.5; %mM
P.fum0=0.1; %mM

P.pep0=0.86; %mM
P.ADP0=0.874; %mM
P.pir0=2.92; %mM
%%% Peroxisome

P.acetatep0=0; %mM
P.CoAp0=20; %mM
P.ATPp0=5; %mM
P.Bifp0=0.1; %mM
P.AcetylCoap0=0; %mM
P.AMPp0=5; %mM
P.oaap0=5; %mM
P.citp0=0; %mM

%%% Mitochondria

P.fumm0=0.91; %mM
P.sucm0=2.6; %mM
P.oaam0=0.01; %mM
P.ATPm0=0.1; %mM
P.ADPm0=4.34; %mM
P.AMPm0=0.15; %mM
P.pirm0=4.67; %mM
P.NADHm0=6.2; %mM
P.NADm0=1.7;  %mM
P.CoAm0=2.1; %mM
P.AcCoAm0=0.027; %mM
P.malm0=1.2; %mM
P.ATPm0=5.75; %mM
P.succCoAm0=0.177; %mM
P.akgm0=4.58; %mM
P.icitm0=1.7; %mM
P.citm0=16; %mM

%%% Extra %mM
P.acetatex0=20; %mM


nombres=["etoh" %1
    "NAD" %
    "aca" %3
    "NADH" %4
    "acetate" %5 
    "acetatep" %6
    "CoAp" %7
    "ATPp" %8
    "Bifp" %9
    "AcetylCoap" %1
    "AMPp" %11
    "acetatex"% 12
    "ATP" %13
    "CoA" %14
    "AcetylCoa" %15
    "Bif" %16 
    "AMP"%17 
    "oaap" %18
    "citp" %19
    "cit" %2
    "iscit" %21
    "succ" %22
    "glio" %23
    "mal" %24
    "oaa" %25
    "fum" %26
"fumm" %27
"sucm" %28
"oaam" %29  
"pep" %3
"ADP" %31
"pir" %32
"ATPm" %33
"ADPm" %34
"AMPm" %35
"pirm"  %36
"NADHm" %37
"NADm" %38
"CoAm" %39
"AcCoAm" %4 
"malm" %41
"succCoAm" %42
"akgm" %43
"icitm" %44
"citm" %45
]; 

NOMBRE_flujos=["ADH"
    "ALDH"
    "Acet_diff"
    "Acet_diffx"
    "Acs1p"
    "Acs1"
    "ATPt"
    "AMPt"
    "CSp"
    "CSexp"
    "Acon"
    "iscly"
    "MS"
    "MSH"
    "oaa_t"
    "sucfum"
    "etoh_uptake"
    "FUMHY"
    "PEP"
    "PK"
    "oaam_t"
    "CAK"
    "CAKm"
    "MPC"
    "ATPtm"
    "ADPtm"
    "RES"
    "PDH"
    "MAE"
    "MDHm"
    "FUMm"
    "SDH"
    "SCS"
    "AK"
    "IDH"
    "ACOm"
    "CSM"
];
%% NO BORRAR!!!!
P.coef=1;

% Kinetic multip.
P.l=0.01;

P.citcrit=0.1750;
 %% Initial conditions vector
P.X0=[P.etoh0 %1
    P.NAD0 %2
    P.aca0 %3
    P.NADH0 %4
    P.acetate0 %5 
    P.acetatep0 %6
    P.CoAp0 %7
    P.ATPp0 %8
    P.Bifp0 %9
    P.AcetylCoap0 %10
    P.AMPp0 %11
    P.acetatex0% 12
    P.ATP0 %13
    P.CoA0 %14
    P.AcetylCoa0 %15
    P.Bif0 %16 
    P.AMP0%17 
    P.oaap0 %18
    P.citp0 %19
    P.cit0 %20
    P.iscit0 %21
    P.succ0 %22
    P.glio0 %23
    P.mal0 %24
    P.oaa0 %25
    P.fum0 %26
P.fumm0 %27
P.sucm0 %28
P.oaam0 %29  
P.pep0 %30
P.ADP0 %31
P.pir0 %32
P.ATPm0 %33
P.ADPm0 %34
P.AMPm0 %35
P.pirm0  %36
P.NADHm0 %37
P.NADm0 %38
P.CoAm0 %39
P.AcCoAm0 %40 
P.malm0 %41
P.succCoAm0 %42
P.akgm0 %43
P.icitm0 %44
P.citm0 %45
]; 
 
%% Tests
%%Las  acetyl COa syntsas son ahora irreversibles 
%% Test reference
P.CS_V_deletion=1; 
 P.MS_deletion=1;
 P.Acs1_deletion=1;
 P.etoh=0.083;
 [tiempo,mets,flujos]=S(P,W,Z);  
%  Test cit2p
          P.CS_V_deletion=0;
           P.MS_deletion=1;
           P.Acs1_deletion=1;
           P.etoh=0.031417115;
          [tiempo_mut,mets_mut,flujos_mut]=S(P,W,Z);
 %  %% Test second mutation mls
            P.etoh=0.051;
           P.MS_deletion=0;
            P.CS_V_deletion=1;
            P.Acs1_deletion=1;
            [tiempo_mut2,mets_mut2,flujos_mut2]=S(P,W,Z);
  %% Test third mutation
           P.etoh=(7.8-7)/(46/1000*40*60);
           P.MS_deletion=1;
           P.CS_V_deletion=1;
           P.Acs1_deletion=0;
           P.time=[1:1:Tend]; %min
           P.X0(12)=1.6;
           [tiempo_mut3,mets_mut3,flujos_mut3]=S(P,W,Z);
 %%Leer archivos
file1 = readtable("Ref.txt");
time = 60*file1.Var1-60*20;
Acetref = 16.93*file1.Var2;
file2 = readtable("Mut.txt");
time2 = 60*file2.Var1-60*20;
Acmut = 16.93*file2.Var2;

file3=readtable("Mls.txt");
time3 = 60*file3.Var1-60*20;
Acmls = 16.93*file3.Var2;

file4=readtable("Acs.txt");
time4 = 60*file4.Var1-60*20;
Acs = 16.93*file4.Var2;

writematrix(mets,"metsref.csv");
writematrix(mets_mut,"metsmut1.csv");
writematrix(mets_mut2,"metsmut2.csv");
writematrix(mets_mut3,"metsmut3.csv");

writematrix(tiempo,"tref.csv");
writematrix(tiempo_mut,"tmut1.csv");
writematrix(tiempo_mut2,"tmut2.csv");
writematrix(tiempo_mut3,"tmut3.csv");

close all
function dvars = sysODE(t,vars,P,Z,W)
%% Descripción de variables
etoh=vars(1); %mM 
NAD=vars(2); %mM
aca=vars(3); %mM
NADH=vars(4); %mM
Acetate=vars(5); %mM
Acetatep=vars(6); %mM
CoAp=vars(7); %mM
ATPp=vars(8); %mM
Bifp=vars(9); %mM
AcetylCoAp=vars(10); %mM
AMPp=vars(11); %mM
Acetatex=vars(12); %mM
ATP=vars(13); %mM
CoA=vars(14); %mM
AcetylCoA=vars(15); %mM
Bif=vars(16); %mM
AMP=vars(17); %mM
oaap=vars(18); %mM
citp=vars(19); %mM
cit=vars(20); %mM
iscit=vars(21); %mM
succ=vars(22); %mM
gliox=vars(23); %mM
mal=vars(24);%mM
oaa=vars(25);  %mM
fum=vars(26);%mM 
fumm=vars(27);%mM
succm=vars(28); %mM
oaam=vars(29); %mM
pep=vars(30); %mM
ADP=vars(31); %mM
pir=vars(32); %mM
ATPm=vars(33); %mM
ADPm=vars(34); %mM
AMPm=vars(35); %mM
pirm=vars(36); %mM
NADHm=vars(37); %mM
NADm=vars(38); %mM
CoAm=vars(39); %mM
AcCoAm=vars(40); %mM
malm=vars(41); %mM
succCoAm=vars(42); %mM
akgm=vars(43); %mM
icitm=vars(44); %mM
citm=vars(45); %mM

%% Cinéticas mM/min.Obs: Todas las cinéticas de transporte entre compartmientos estan expresadas en mMol 
%mientras que todo el resto esta en mMolar
V_ADH=(P.ADH_Vm/(P.ADH_Km_aca*P.ADH_Km_NADH))*(-etoh*NAD/P.ADH_Keq)/((1+aca/P.ADH_Km_aca)*(1+NADH/P.ADH_Km_NADH)+(1+etoh/P.ADH_Km_etoh)*(1+NAD/P.ADH_Km_NAD)-1);

V_ALDH=(P.ALDH_Vm*NAD*aca)/(P.ALDH_Kia*P.ALDH_Kb+P.ALDH_Ka*aca+NAD*aca+aca^2/P.ALDH_Ki_Ald+aca^2*NADH/P.ALDH_Ki_Ald_NADH+aca*NADH/P.ALDH_Ki_NADH);

V_Acet_diff=P.Dp*(Acetatep-Acetate);

P.f=(0.7*(0.12^P.n/(0.12^P.n+mal^P.n))+0.3*P.citcrit^P.n/(P.citcrit^P.n+cit^P.n));  
P.f2=0.05^P.n/(0.05^P.n+AcetylCoA^P.n);


V_Acet_diffx=(1-P.f2)*(-P.Difex*(Acetatex-P.Acet_coeff*Acetate)...
*P.Acet_coeff_2*Acetate/(P.Acet_coeff_3+P.Acet_coeff_4*Acetatex)*...
(1-P.f)+(P.f)*P.kd*Acetate/(P.Acet_coeff_5+Acetate))+0.03*P.f2*Acetate^1.1/((5+Acetate)); 

V_Acs1p=P.Acs1_deletion*Z.V1*(Acetatep*ATPp*CoAp-0*AcetylCoAp*AMPp*Bifp/(Z.Keq))/(1+Z.Ka*Acetatep+Z.Kb*ATPp);

V_Acs1=P.Acs1_deletion*W.V1*(Acetate*ATP*CoA-0*AcetylCoA*AMP*Bif/(W.Keq))/(1+W.Ka*Acetate+W.Kb*ATP);

V_ATPt=0.1*P.ATP_V*ATP/(P.ATP_Km+ATP);

V_AMPt=P.AMP_V*AMPp/(P.AMP_Km+AMPp);

V_CSp=P.CS_V_deletion*(P.CS_V*oaap*AcetylCoAp)/((P.CS_Km_oaa+oaap)*(P.CS_Km_acCoA+AcetylCoAp));

V_CSexp=P.Cexp_Vm*citp/(P.Cexp_Km+citp);

V_Acon=P.AC_V*(cit-iscit/P.AC_Keq)/(1+cit/P.AC_Km_cit+iscit/P.AC_Km_iscit);

V_iscly=0.01*P.iscly_Vf*(iscit)/(P.iscly_Kms_iscit+iscit*(1+succ/P.iscly_Kip_succ)+P.iscly_Vf/(P.iscly_Vr*P.iscly_Keq)*(P.iscly_Kmq_glio*succ+P.iscly_Kmp_succ*gliox+gliox*succ));

num3=10*P.MS_deletion*P.MS_V1*P.MS_V2*(gliox*(AcetylCoA)-0*CoA*mal/P.MS_Keq);
denom3=P.MS_V2*P.MS_Ki_AcCoA+...;
    P.MS_V2*AcetylCoA+...;
    P.MS_V2*AcetylCoA*gliox+...;
    P.MS_V1*P.MS_Kq_mal*CoA/P.MS_Keq+...;
    P.MS_V1*P.MS_Kp_CoA*mal/P.MS_Keq+...;
    P.MS_V1*CoA*mal/P.MS_Keq+...;
    P.MS_V2*P.MS_Ka_AcCoA*gliox*mal/P.MS_Kiiq_mal+...;
    P.MS_V2*AcetylCoA*CoA/(P.MS_Kiip_CoA) ;
V_MS=num3/denom3;

num4=P.MSH_V/(P.MSH_Km_mal*P.MSH_Km_NAD)*(mal*NAD-0*oaa*NADH/P.MSH_Keq);
denom4=1+mal/P.MSH_Km_mal+NAD/P.MSH_Km_NAD+mal*NAD/(P.MSH_Km_NAD*P.MSH_Km_mal)+...;
      oaa/P.MSH_Km_oaa+NADH/P.MSH_Km_NADH+oaa/P.MSH_Km_oaa*NADH/P.MSH_Km_NADH;
V_MSH=num4/denom4;

V_oaa_t=P.oaat_V*(oaa/P.oaat_K1-oaap/P.oaat_K2)/(1+oaa/P.oaat_K1+oaap/P.oaat_K2);

V_sucfum=P.sucfum_V*succ*fumm/(P.sucfum_Kia_fumm*P.sucfum_Kb_suc+fumm*P.sucfum_Kb_suc+succ*P.sucfum_Ka_fumm+fumm*succ);

V_etoh_uptake=P.etoh ;

V_FUMHY=P.FUMHY_Vm/P.FUMHY_Km_fum*(fum-mal/P.FUMHY_Keq)/(1+fum/P.FUMHY_Km_fum+mal/P.FUMHY_Km_mal)*P.FUMHY_Ki_ATP/(P.FUMHY_Ki_ATP+ATP);

num4=P.PECK_V1*P.PECK_Kb_oaa*(ATP*oaa-ADP*pep/P.PECK_Keq);
denom4=P.PECK_Kb_oaa*P.PECK_Ki_ATP*P.PECK_Kb_oaa+...;
    P.PECK_Kb_oaa*P.PECK_Kb_oaa*ATP+...;
    P.PECK_Kb_oaa*ATP*oaa+...;
    P.PECK_V1*P.PECK_Kq_pep*ADP/P.PECK_Keq+...;
    P.PECK_V1*P.PECK_Kp_ADP*pep/P.PECK_Keq+...;
    P.PECK_V1*ADP*pep/P.PECK_Keq+...;
    P.PECK_Kb_oaa*P.PECK_Ka_ATP*oaa*pep/P.PECK_Kiiq_pep+...;
    P.PECK_Kb_oaa*P.PECK_Kb_oaa*ATP*ADP/P.PECK_Kiip_ADP ;
V_PEP=num4/denom4;

V_PK=P.PK_Vm*ADP*pep*ATP/(P.PK_Km_ADP*P.PK_Km_pep+P.PK_Km_ADP*pep+ADP*P.PK_Km_pep);

  %% Cinéticas mitocondriales de KUMMER

V_oaam_t=P.OAC_Vm*(oaa/P.OAC_Km_oaa-oaam/P.OAC_Km_oaam)/(1+oaa/P.OAC_Km_oaa+oaam/P.OAC_Km_oaam); %canonica %mMol

V_CAK=P.CAK_Vm/(P.CAK_Km_ATP*P.CAK_Km_AMP)*(ATP*AMP-ADP*ADP/P.CAK_Keq)/((1+ATP/P.CAK_Km_ATP)*(1+AMP/P.CAK_Km_AMP)+(1+ADP/P.CAK_Km_ADP)^2-1);

V_CAKm=P.CAKm_Vm/(P.CAKm_Km_ATP*P.CAKm_Km_AMP)*(ATPm*AMPm-ADPm*ADPm/P.CAKm_Keq)/((1+ATPm/P.CAKm_Km_ATP)*(1+AMPm/P.CAKm_Km_AMP)+(1+ADPm/P.CAKm_Km_ADP)^2-1);

V_MPC=P.MPC_Vm*(pir/P.MPC_Km_pir-pirm/P.MPC_Km_pirm)/(1+pir/P.MPC_Km_pir+pirm/P.MPC_Km_pirm); %canonica

V_ATPtm=P.ATPtm_Vm*ATPm/(P.ATPtm_Km+ATPm);

V_ADPtm=P.ADPtm_Vm*ADP/(P.ADPtm_Km+ADP);

V_RES=(P.RES_Vm/(P.RES_Km_NADH*P.RES_Km_NADm)*(NADH*NADm-NAD*NADHm/P.RES_Keq))/((1+NADH/P.RES_Km_NADH)*(1+NADm/P.RES_Km_NADm)+(1+NAD/P.RES_Km_NAD)*(1+NADHm/(P.RES_Km_NADHm))-1);

V_PDH=(P.PDH_Vm*pirm*CoAm*NADm/(P.PDH_Km_pirm*P.PDH_Km_CoAm*P.PDH_Km_NADm))/((1+pirm/P.PDH_Km_pirm)*(1+CoAm/P.PDH_Km_CoAm)*(1+NADm/P.PDH_Km_NADm))*P.PDH_Ki_AcCoAm*P.PDH_Ki_NADHm/((P.PDH_Ki_AcCoAm+AcCoAm)*(P.PDH_Ki_NADHm+NADHm));

V_MAE=P.MAE_Vm*malm*NADm/(P.MAE_Km_malm*P.MAE_Km_NADm)/((1+malm/P.MAE_Km_malm)*(1+NADm/P.MAE_Km_NADm));

V_MDHm=P.MDHm_Vm/(P.MDHm_Km_malm*P.MDHm_Km_NADm)*(malm*NADm-oaam*NADHm/P.MDHm_Keq)/((1+malm/P.MDHm_Km_malm)*((1+NADm/P.MDHm_Km_NADm)+(oaam/P.MDHm_Km_oaam)*(1+NADHm/P.MDHm_Km_NADHm)-1));

V_FUMm=P.FUMHY_Vm/(P.FUMm_Km_fumm)*(fumm-malm/P.FUMm_Keq)/(1+fumm/P.FUMm_Km_fumm+malm/P.FUMm_Km_malm)*P.FUMm_Ki_ATPm/(P.FUMm_Ki_ATPm+ATPm);

V_SDH=P.SDH_Vm/P.SDH_Km_succ*(succm-fumm/P.SDH_Keq)/(1+succm/P.SDH_Km_succ+fumm/P.SDH_Km_fum);

V_SCS=P.SCS_Vm/(P.SCS_Km_SucCoA*P.SCS_Km_ADPm)*(succCoAm*ADPm-succm*CoAm*ATPm/P.SCS_Keq)/((1+succCoAm/P.SCS_Km_SucCoA)*(1+ADPm/P.SCS_Km_ADPm)+(1+succm/P.SCS_Km_succm)*(1+CoAm/P.SCS_Km_CoAm)*(1+ATPm/P.SCS_Km_ATPm)-1);

V_AK=P.AK_Vm*akgm*CoAm*NADm/(P.AK_Km_akg*P.AK_Km_CoAm*P.AK_Km_NADm)/((1+akgm/P.AK_Km_akg)*(1+CoAm/P.AK_Km_CoAm)*(1+NADm/P.AK_Km_NADm))*P.AK_Ki_NADHm/(P.AK_Ki_NADHm+NADHm)*(1+ADPm/P.AK_Ka_ADPm);

V_IDH=P.IDH_Vm*(icitm/P.IDH_Km_icitm)^4*(NADm/P.IDH_Km_NADm)^2/((1+(icitm/P.IDH_Km_icitm)^4)*(1+(NADm/P.IDH_Km_NADm)^2))*P.IDH_Ki_NADHm/(P.IDH_Ki_NADHm+NADHm)*P.IDH_Ki_akgm/(P.IDH_Ki_akgm+akgm);

V_ACOm=P.ACOm_Vm/P.ACOm_Km_citm*(citm-icitm/P.ACOm_Keq)/(1+cit/P.ACOm_Km_citm+icitm/P.ACOm_Km_icitm);

V_CSM=P.CSm_Vm*oaam*AcCoAm/(P.CSm_Km_oaam*P.CSm_Km_AcCoAm)/((1+oaam/P.CSm_Km_oaam)*(1+AcCoAm/P.CSm_Km_AcCoAm));


%% Kinetic matrix
m=1; %ADH
 p=1; % aldh
 s=1; % Acs1p
 ss=1 ;%Acs1
 o=1;  % CS
 t=1; % oaap_t
dvars(1)=V_etoh_uptake/P.l+m*V_ADH;                               % etanol citoplasmático
dvars(2)=(m*V_ADH-p*V_ALDH)-V_MSH+V_RES/((P.Vcels-P.vm));       % NAD citoplasmático
dvars(3)=-m*V_ADH-p*V_ALDH;                                     % acetaldheído citoplasmático
dvars(4)=-m*V_ADH+p*V_ALDH+V_MSH-V_RES/((P.Vcels-P.vm));        % NADH citoplasmático
dvars(5)=m*V_ALDH+V_Acet_diff/P.Vcels-V_Acet_diffx-ss*V_Acs1; % Acetato citoplasmático
dvars(6)=-V_Acet_diff/P.vp-s*V_Acs1p;                   % Acetato peroxisomal
dvars(7)=-s*V_Acs1p + o*V_CSp;                               %Coa peroxisomal
dvars(8)=-s*V_Acs1p+V_ATPt/(P.vp);                               %ATP peroxisomal
dvars(9)=s*V_Acs1p; %Bifosfato peroxisomal
dvars(10)= s*V_Acs1p-o*V_CSp; %Acetyl CoA Peroxisomal
dvars(11)=s*V_Acs1p-V_AMPt/(P.vp); %AMP peroxisomal
dvars(12)=V_Acet_diffx/(P.Vtot) ; % Acetato externo
dvars(13)=(-ss*V_Acs1-V_ATPt/(P.Vcels-P.vp)-V_PEP-V_CAK+V_ATPtm/(P.Vcels-P.vp)+V_PK)*1; %ATP citoplasmático
dvars(14)= -ss*V_Acs1+V_MS; %Coa citoplasmático
dvars(15)= ss*V_Acs1-V_MS; %Acetyl CoA  citoplasmático
dvars(16)= ss*V_Acs1; % Bifosfato citoplasmático
dvars(17)= ss*V_Acs1+P.coef*V_AMPt/(P.Vcels-P.vp)-V_CAK; %AMP citoplasmático
dvars(18)=-o*V_CSp+t*V_oaa_t/P.vp; %oaa peroxisomal
dvars(19)=o*V_CSp-V_CSexp/P.vp; %citrato peroxisomal
dvars(20)=V_CSexp/(P.Vcels-P.vp)-V_Acon; %citrato citoplasmático
dvars(21)=V_Acon-V_iscly; %isocitrato citoplasmático
dvars(22)=V_iscly-V_sucfum/(P.Vcels-P.vm); %succinato citoplasmático
dvars(23)=V_iscly-V_MS; %glioxilato citoplasmático
dvars(24)=V_MS-V_MSH+V_FUMHY; %malato citoplasmático
dvars(25)=V_MSH-t*V_oaa_t/(P.Vcels-P.vp)-V_oaam_t/(P.Vcels-P.vm)-V_PEP; %oxalacetato citoplasmático
dvars(26)=V_sucfum/(P.Vcels-P.vm)-V_FUMHY;%fumarato citoplasmático
dvars(27)=-V_sucfum/P.vm-V_FUMm+V_SDH;  %fumarato mitocondrial
dvars(28)=V_sucfum/P.vm-V_SDH+V_SCS; %succinato mitocondrial
dvars(29)=V_oaam_t/P.vm+V_MDHm-V_CSM; %oxalacetato mitocondrial
dvars(30)=V_PEP; % fosfoenol piruvato citoplasmático
dvars(31)=V_PEP-V_PK+2*V_CAK-V_ADPtm-V_SCS; % ADP citoplasmático
dvars(32)=V_PK-V_MPC; %piruvato citoplasmático DIVIDR POR VOLUMENES TRANCA EL PROGRAMA!!!(1)
dvars(33)=-V_CAKm-V_ATPtm/P.vm+V_SCS; %ATP mitocondrial
dvars(34)=2*V_CAKm+V_ADPtm;%ADP mitocondrial
dvars(35)=-V_CAKm;%AMP mitocondrial
dvars(36)=V_MPC-V_PDH+V_MAE; %piruvato mitocondrial DIVIDR POR VOLUMENES TRANCA EL PROGRAMA.!!!(1)
dvars(37)=V_RES/P.vm+V_PDH+V_MAE+V_MDHm+V_AK+V_IDH; %NADH mitocondrial
dvars(38)=-V_RES/P.vm-V_PDH-V_MAE-V_MDHm-V_AK-V_IDH; % NAD mitocondrial
dvars(39)=-V_PDH+V_SCS-V_AK+V_CSM;            %CoA mitocondrial
dvars(40)=V_PDH-V_CSM;             %AcCoa mitocondrial
dvars(41)=-V_MAE-V_MDHm+V_FUMm;           %malato mitocondiral 
dvars(42)=-V_SCS+V_AK;                    % SuccCoA mitocondrial
dvars(43)= -V_AK+V_IDH;       %Akg mitocondrial
dvars(44)=-V_IDH+V_ACOm;              %Isocitrato mitocondrial
dvars(45)=-V_ACOm+V_CSM;           % citrato mitcondrial 
dvars=P.l*dvars';

end
function [tiempo,variables,dvariables]=S(P,W,Z)
    
    handle = @(t,vars) sysODE(t,vars,P,Z,W);
    opciones= odeset("RelTol",1e-3,"AbsTol",1e-6);
    [t,metab]=ode23s(handle,P.time,P.X0,opciones);
    tiempo=t;
    etoh=metab(:,1);
    NAD=metab(:,2);
    aca=metab(:,3);
    NADH=metab(:,4);
    Acetate=metab(:,5);
    Acetatep=metab(:,6);
    CoAp=metab(:,7);
    ATPp=metab(:,8);
    Bifp=metab(:,9);
    AcetylCoAp=metab(:,10);
    AMPp=metab(:,11);
    Acetatex=metab(:,12);
    ATP=metab(:,13);
    CoA=metab(:,14);
    AcetylCoA=metab(:,15);
    Bif=metab(:,16);
    AMP=metab(:,17);
    oaap=metab(:,18);
    citp=metab(:,19);
    cit=metab(:,20);
    iscit=metab(:,21);
    succ=metab(:,22);
    gliox=metab(:,23);
    mal=metab(:,24);
    oaa=metab(:,25);
    fum=metab(:,26);
    fumm=metab(:,27);
    succm=metab(:,28);
    oaam=metab(:,29);
    pep=metab(:,30);
    ADP=metab(:,31);
    pir=metab(:,32);
    ATPm=metab(:,33);
    ADPm=metab(:,34);
    AMPm=metab(:,35);
    pirm=metab(:,36);
    NADHm=metab(:,37);
    NADm=metab(:,38);
    CoAm=metab(:,39); 
    AcCoAm=metab(:,40); 
    malm=metab(:,41); 
    succCoAm=metab(:,42);
    akgm=metab(:,43); 
    icitm=metab(:,44); 
    citm=metab(:,45); 
    variables=metab;
%%% Obtención de flujos
%%%%%%%%%%%%%%%%%%%%%%%
V_ADH=(P.ADH_Vm./(P.ADH_Km_aca.*P.ADH_Km_NADH)).*(-etoh.*NAD./P.ADH_Keq)./((1+aca./P.ADH_Km_aca).*(1+NADH./P.ADH_Km_NADH)+(1+etoh./P.ADH_Km_etoh).*(1+NAD./P.ADH_Km_NAD)-1);

V_ALDH=(P.ALDH_Vm.*NAD.*aca)./(P.ALDH_Kia.*P.ALDH_Kb+P.ALDH_Ka.*aca+NAD.*aca+aca.^2./P.ALDH_Ki_Ald+aca.^2.*NADH./P.ALDH_Ki_Ald_NADH+aca.*NADH./P.ALDH_Ki_NADH);

V_Acet_diff=P.Dp.*(Acetatep-Acetate);

P.f=(0.7.*(0.12.^P.n./(0.12.^P.n+mal.^P.n))+0.3.*P.citcrit.^P.n./(P.citcrit.^P.n+cit.^P.n));  % meh
P.f2=0.05.^P.n./(0.05.^P.n+AcetylCoA.^P.n);


V_Acet_diffx=(1-P.f2).*(-P.Difex.*(Acetatex-P.Acet_coeff.*Acetate).*P.Acet_coeff_2.*Acetate./(P.Acet_coeff_3+P.Acet_coeff_4.*Acetatex).*(1-P.f)+(P.f).*P.kd.*Acetate./(P.Acet_coeff_5+Acetate))+0.03.*P.f2.*Acetate.^1.1./((5+Acetate)); 

V_Acs1p=P.Acs1_deletion.*Z.V1.*(Acetatep.*ATPp.*CoAp-0.*AcetylCoAp.*AMPp.*Bifp./(Z.Keq))./(1+Z.Ka.*Acetatep+Z.Kb.*ATPp);

V_Acs1=P.Acs1_deletion.*W.V1.*(Acetate.*ATP.*CoA-0.*AcetylCoA.*AMP.*Bif./(W.Keq))./(1+W.Ka.*Acetate+W.Kb.*ATP);

V_ATPt=0.1.*P.ATP_V.*ATP./(P.ATP_Km+ATP);

V_AMPt=P.AMP_V.*AMPp./(P.AMP_Km+AMPp);

V_CSp=P.CS_V_deletion.*(P.CS_V.*oaap.*AcetylCoAp)./((P.CS_Km_oaa+oaap).*(P.CS_Km_acCoA+AcetylCoAp));

V_CSexp=P.Cexp_Vm.*citp./(P.Cexp_Km+citp);

V_Acon=P.AC_V.*(cit-iscit./P.AC_Keq)./(1+cit./P.AC_Km_cit+iscit./P.AC_Km_iscit);

V_iscly=0.01.*P.iscly_Vf.*(iscit)./(P.iscly_Kms_iscit+iscit.*(1+succ./P.iscly_Kip_succ)+P.iscly_Vf./(P.iscly_Vr.*P.iscly_Keq).*(P.iscly_Kmq_glio.*succ+P.iscly_Kmp_succ.*gliox+gliox.*succ));

num3=10.*P.MS_deletion.*P.MS_V1.*P.MS_V2.*(gliox.*(AcetylCoA)-0.*CoA.*mal./P.MS_Keq);
denom3=P.MS_V2.*P.MS_Ki_AcCoA.*P.PECK_Kb_oaa+...;
    P.MS_V2.*P.PECK_Kb_oaa.*AcetylCoA+...;
    P.MS_V2.*AcetylCoA.*gliox+...;
    P.MS_V1.*P.MS_Kq_mal.*CoA./P.MS_Keq+...;
    P.MS_V1.*P.MS_Kp_CoA.*mal./P.MS_Keq+...;
    P.MS_V1.*CoA.*mal./P.MS_Keq+...;
    P.MS_V2.*P.MS_Ka_AcCoA.*gliox.*mal./P.MS_Kiiq_mal+...;
    P.MS_V2.*P.PECK_Kb_oaa.*AcetylCoA.*CoA./P.MS_Kiip_CoA ;
V_MS=num3./denom3;

num4=P.MSH_V./(P.MSH_Km_mal.*P.MSH_Km_NAD).*(mal.*NAD-0.*oaa.*NADH./P.MSH_Keq);
denom4=1+mal./P.MSH_Km_mal+NAD./P.MSH_Km_NAD+mal.*NAD./(P.MSH_Km_NAD.*P.MSH_Km_mal)+...;
      oaa./P.MSH_Km_oaa+NADH./P.MSH_Km_NADH+oaa./P.MSH_Km_oaa.*NADH./P.MSH_Km_NADH;
V_MSH=num4./denom4;

V_oaa_t=P.oaat_V.*(oaa./P.oaat_K1-oaap./P.oaat_K2)./(1+oaa./P.oaat_K1+oaap./P.oaat_K2);

V_sucfum=P.sucfum_V.*succ.*fumm./(P.sucfum_Kia_fumm.*P.sucfum_Kb_suc+fumm.*P.sucfum_Kb_suc+succ.*P.sucfum_Ka_fumm+fumm.*succ);

V_etoh_uptake=P.etoh ;

V_FUMHY=P.FUMHY_Vm./P.FUMHY_Km_fum.*(fum-mal./P.FUMHY_Keq)./(1+fum./P.FUMHY_Km_fum+mal./P.FUMHY_Km_mal).*P.FUMHY_Ki_ATP./(P.FUMHY_Ki_ATP+ATP);

num4=P.PECK_V1.*P.PECK_Kb_oaa.*(ATP.*oaa-ADP.*pep./P.PECK_Keq);
denom4=P.PECK_Kb_oaa.*P.PECK_Ki_ATP.*P.PECK_Kb_oaa+...;
    P.PECK_Kb_oaa.*P.PECK_Kb_oaa.*ATP+...;
    P.PECK_Kb_oaa.*ATP.*oaa+...;
    P.PECK_V1.*P.PECK_Kq_pep.*ADP./P.PECK_Keq+...;
    P.PECK_V1.*P.PECK_Kp_ADP.*pep./P.PECK_Keq+...;
    P.PECK_V1.*ADP.*pep./P.PECK_Keq+...;
    P.PECK_Kb_oaa.*P.PECK_Ka_ATP.*oaa.*pep./P.PECK_Kiiq_pep+...;
    P.PECK_Kb_oaa.*P.PECK_Kb_oaa.*ATP.*ADP./P.PECK_Kiip_ADP ;
V_PEP=num4./denom4;

V_PK=P.PK_Vm.*ADP.*pep.*ATP./(P.PK_Km_ADP.*P.PK_Km_pep+P.PK_Km_ADP.*pep+ADP.*P.PK_Km_pep);

  %% Cinéticas mitocondriales de KUMMER

V_oaam_t=P.OAC_Vm.*(oaa./P.OAC_Km_oaa-oaam./P.OAC_Km_oaam)./(1+oaa./P.OAC_Km_oaa+oaam./P.OAC_Km_oaam); %canonica %mMol

V_CAK=P.CAK_Vm./(P.CAK_Km_ATP.*P.CAK_Km_AMP).*(ATP.*AMP-ADP.*ADP./P.CAK_Keq)./((1+ATP./P.CAK_Km_ATP).*(1+AMP./P.CAK_Km_AMP)+(1+ADP./P.CAK_Km_ADP).^2-1);

V_CAKm=P.CAKm_Vm./(P.CAKm_Km_ATP.*P.CAKm_Km_AMP).*(ATPm.*AMPm-ADPm.*ADPm./P.CAKm_Keq)./((1+ATPm./P.CAKm_Km_ATP).*(1+AMPm./P.CAKm_Km_AMP)+(1+ADPm./P.CAKm_Km_ADP).^2-1);

V_MPC=P.MPC_Vm.*(pir./P.MPC_Km_pir-pirm./P.MPC_Km_pirm)./(1+pir./P.MPC_Km_pir+pirm./P.MPC_Km_pirm); %canonica

V_ATPtm=P.ATPtm_Vm.*ATPm./(P.ATPtm_Km+ATPm);

V_ADPtm=P.ADPtm_Vm.*ADP./(P.ADPtm_Km+ADP);

V_RES=(P.RES_Vm./(P.RES_Km_NADH.*P.RES_Km_NADm).*(NADH.*NADm-NAD.*NADHm./P.RES_Keq))./((1+NADH./P.RES_Km_NADH).*(1+NADm./P.RES_Km_NADm)+(1+NAD./P.RES_Km_NAD).*(1+NADHm./(P.RES_Km_NADHm))-1);

V_PDH=(P.PDH_Vm.*pirm.*CoAm.*NADm./(P.PDH_Km_pirm.*P.PDH_Km_CoAm.*P.PDH_Km_NADm))./((1+pirm./P.PDH_Km_pirm).*(1+CoAm./P.PDH_Km_CoAm).*(1+NADm./P.PDH_Km_NADm)).*P.PDH_Ki_AcCoAm.*P.PDH_Ki_NADHm./((P.PDH_Ki_AcCoAm+AcCoAm).*(P.PDH_Ki_NADHm+NADHm));

V_MAE=P.MAE_Vm.*malm.*NADm./(P.MAE_Km_malm.*P.MAE_Km_NADm)./((1+malm./P.MAE_Km_malm).*(1+NADm./P.MAE_Km_NADm));

V_MDHm=P.MDHm_Vm./(P.MDHm_Km_malm.*P.MDHm_Km_NADm).*(malm.*NADm-oaam.*NADHm./P.MDHm_Keq)./((1+malm./P.MDHm_Km_malm).*((1+NADm./P.MDHm_Km_NADm)+(oaam./P.MDHm_Km_oaam).*(1+NADHm./P.MDHm_Km_NADHm)-1));

V_FUMm=P.FUMHY_Vm./(P.FUMm_Km_fumm).*(fumm-malm./P.FUMm_Keq)./(1+fumm./P.FUMm_Km_fumm+malm./P.FUMm_Km_malm).*P.FUMm_Ki_ATPm./(P.FUMm_Ki_ATPm+ATPm);

V_SDH=P.SDH_Vm./P.SDH_Km_succ.*(succm-fumm./P.SDH_Keq)./(1+succm./P.SDH_Km_succ+fumm./P.SDH_Km_fum);

V_SCS=P.SCS_Vm./(P.SCS_Km_SucCoA.*P.SCS_Km_ADPm).*(succCoAm.*ADPm-succm.*CoAm.*ATPm./P.SCS_Keq)./((1+succCoAm./P.SCS_Km_SucCoA).*(1+ADPm./P.SCS_Km_ADPm)+(1+succm./P.SCS_Km_succm).*(1+CoAm./P.SCS_Km_CoAm).*(1+ATPm./P.SCS_Km_ATPm)-1);

V_AK=P.AK_Vm.*akgm.*CoAm.*NADm./(P.AK_Km_akg.*P.AK_Km_CoAm.*P.AK_Km_NADm)./((1+akgm./P.AK_Km_akg).*(1+CoAm./P.AK_Km_CoAm).*(1+NADm./P.AK_Km_NADm)).*P.AK_Ki_NADHm./(P.AK_Ki_NADHm+NADHm).*(1+ADPm./P.AK_Ka_ADPm);

V_IDH=P.IDH_Vm.*(icitm./P.IDH_Km_icitm).^4.*(NADm./P.IDH_Km_NADm).^2./((1+(icitm./P.IDH_Km_icitm).^4).*(1+(NADm./P.IDH_Km_NADm).^2)).*P.IDH_Ki_NADHm./(P.IDH_Ki_NADHm+NADHm).*P.IDH_Ki_akgm./(P.IDH_Ki_akgm+akgm);

V_ACOm=P.ACOm_Vm./P.ACOm_Km_citm.*(citm-icitm./P.ACOm_Keq)./(1+cit./P.ACOm_Km_citm+icitm./P.ACOm_Km_icitm);

V_CSM=P.CSm_Vm.*oaam.*AcCoAm./(P.CSm_Km_oaam.*P.CSm_Km_AcCoAm)./((1+oaam./P.CSm_Km_oaam).*(1+AcCoAm./P.CSm_Km_AcCoAm));

flujos(:,1)=    V_ADH;
flujos(:,2)=    V_ALDH;
flujos(:,3)=    V_Acet_diff;
flujos(:,4)=   V_Acet_diffx;
flujos(:,5)=    V_Acs1p;
flujos(:,6)=    V_Acs1;
flujos(:,7)=    V_ATPt;
flujos(:,8)=    V_AMPt;
flujos(:,9)=   V_CSp;
flujos(:,10)=    V_CSexp;
flujos(:,11)=    V_Acon;
flujos(:,12)=    V_iscly;
flujos(:,13)=    V_MS;
flujos(:,14)=    V_MSH;
flujos(:,15)=    V_oaa_t;
flujos(:,16)=    V_sucfum;
flujos(:,17)=    V_etoh_uptake;
flujos(:,18)=    V_FUMHY;
flujos(:,19)=    V_PEP;
flujos(:,20)=    V_PK;
flujos(:,21)=    V_oaam_t;
flujos(:,22)=    V_CAK;
flujos(:,23)=    V_CAKm;
flujos(:,24)=    V_MPC;
flujos(:,25)=    V_ATPtm;
flujos(:,26)=    V_ADPtm;
flujos(:,27)=    V_RES;
flujos(:,28)=    V_PDH;
flujos(:,29)=    V_MAE;
flujos(:,30)=    V_MDHm;
flujos(:,31)=    V_FUMm;
flujos(:,32)=    V_SDH;
flujos(:,33)=    V_SCS;
flujos(:,34)=    V_AK;
flujos(:,35)=    V_IDH;
flujos(:,36)=    V_ACOm;
flujos(:,37)= V_CSM;
dvariables =  P.l*flujos;    
end
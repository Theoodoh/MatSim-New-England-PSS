%% this m-file perfrmes small-signal analysis
initDyn

% x = [27.7496      48.8779      23.3334       14.588      49.8646      32.9692       49.982      21.0934      10.8242     0.415166    0.0216259    0.0700703         0.02     0.246889    0.0212768    0.0854628    0.0204438   0.00102864    0.0202441    0.0323328    0.0200123     0.116839    0.0203005      0.22581    0.0672576     0.310696    0.0200447    0.0512819     0.170185      0.10437    0.0200754     0.210651    0.0274088     0.172856     0.299569     0.210283    0.0200893     0.159848         0.02    0.0926751    0.0200991   0.00100259    0.0448687     0.176143    0.0411551];
x = [23.4494      3.07216      20.9583      12.3951      12.8864      13.8002      28.0156      43.1885      14.7327     0.435078    0.0873338     0.411417     0.231506     0.558516      0.18848     0.743266    0.0716085     0.209208     0.053606     0.430979     0.233621     0.449469     0.111869     0.246019     0.228878     0.204096    0.0863517     0.151375     0.156983     0.220933    0.0667642     0.352767     0.171217     0.530037      0.42231     0.365416     0.249382      0.29222     0.535743     0.459654    0.0784239     0.344652     0.503143     0.187671     0.357315];

Tw = 10;
%% G1
KG1 = x(1);
T11 = x(10);
T12 = x(11);
T13 = x(12);
T14 = x(13);
Kpss1 = KG1*T11*T13/(T12*T14);

%% G3
KG3 = x(2);
T31 = x(14);
T32 = x(15);
T33 = x(16);
T34 = x(17);
Kpss3 =  KG3*T31*T33/(T32*T34);

%% G4
KG4 = x(3);
T41 = x(18);
T42 = x(19);
T43 = x(20);
T44 = x(21);
Kpss4 =  KG4*T41*T43/(T42*T44);

%% G5
KG5 = x(4);
T51 = x(22);
T52 = x(23);
T53 = x(24);
T54 = x(25);
Kpss5 =  KG5*T51*T53/(T52*T54);

%% G6
KG6 = x(5);
T61 = x(26);
T62 = x(27);
T63 = x(28);
T64 = x(29);
Kpss6 =  KG6*T61*T63/(T62*T64);

%% G7
KG7 = x(6);
T71 = x(30);
T72 = x(31);
T73 = x(32);
T74 = x(33);
Kpss7 =  KG7*T71*T73/(T72*T74);

%% G8
KG8 = x(7);
T81 = x(34);
T82 = x(35);
T83 = x(36);
T84 = x(37);
Kpss8 =  KG8*T81*T83/(T82*T84);

%% G9
KG9 = x(8);
T91 = x(38);
T92 = x(39);
T93 = x(40);
T94 = x(41);
Kpss9 =  KG9*T91*T93/(T92*T94);

%% G10
KG10 = x(9);
T101 = x(42);
T102 = x(43);
T103 = x(44);
T104 = x(45);
Kpss10 =  KG10*T101*T103/(T102*T104);

%% Linearize Power System
% f11=linmod('sys10m');
f11=linmod('sys10m_pss');

% dx/dt = A.x + B.u
% y = C.x + D.u

Asys = f11.a ;
Bsys = f11.b ;
Csys = f11.c ;
Dsys = f11.d ;

%% Calculate Eigenvalues
egs = eig(Asys)
Ns = length(egs);
Damp=-real(egs)./sqrt(real(egs).^2+imag(egs).^2)
freq=abs(imag(egs))/(2*pi)

%% calculae Participation Factors
[Vs,D_eig] = eig(Asys);
Ws=inv(Vs);
for i=1:Ns
    for k=1:Ns
        Pfact1(k,i)=abs(Vs(k,i))*abs(Ws(i,k));
    end
end

for i=1:Ns
     Pfact(i,:)=Pfact1(i,:)/sum(Pfact1(i,:));
end

for i=1:Ns
    [s_val s_idx]=sort(Pfact(:,i),'descend');
    mod_idx(i,:)=s_idx(1:4)';
    pf_fact(i,:)=s_val(1:4)';
end
mod_idx;
pf_fact;

idx_zero = find(abs(egs)<1e-6)
EG_new = egs;
EG_new(idx_zero)=[];

idx_EMs = find(max(real(EG_new)))
EMs = EG_new(idx_EMs)

function J = pss_objf(x,FunIndex,Dim)

load('sys_IO','f11')
As = f11.a;
Bs = f11.b;
Cs = f11.c;
Ds = f11.d;

Tw = 10;
KG1 = x(1);
T11 = x(10);
T12 = x(11);
T13 = x(12);
T14 = x(13);
Kpss1 = KG1*T11*T13/(T12*T14);

KG3 = x(2);
T31 = x(14);
T32 = x(15);
T33 = x(16);
T34 = x(17);
Kpss3 =  KG3*T31*T33/(T32*T34);

KG4 = x(3);
T41 = x(18);
T42 = x(19);
T43 = x(20);
T44 = x(21);
Kpss4 =  KG4*T41*T43/(T42*T44);

KG5 = x(4);
T51 = x(22);
T52 = x(23);
T53 = x(24);
T54 = x(25);
Kpss5 =  KG5*T51*T53/(T52*T54);

KG6 = x(5);
T61 = x(26);
T62 = x(27);
T63 = x(28);
T64 = x(29);
Kpss6 =  KG6*T61*T63/(T62*T64);

KG7 = x(6);
T71 = x(30);
T72 = x(31);
T73 = x(32);
T74 = x(33);
Kpss7 =  KG7*T71*T73/(T72*T74);

KG8 = x(7);
T81 = x(34);
T82 = x(35);
T83 = x(36);
T84 = x(37);
Kpss8 =  KG8*T81*T83/(T82*T84);

KG9 = x(8);
T91 = x(38);
T92 = x(39);
T93 = x(40);
T94 = x(41);
Kpss9 =  KG9*T91*T93/(T92*T94);

KG10 = x(9);
T101 = x(42);
T102 = x(43);
T103 = x(44);
T104 = x(45);
Kpss10 =  KG10*T101*T103/(T102*T104);

b1 = [KG1*T11*T13*Tw (KG1*T11*Tw + KG1*T13*Tw) KG1*Tw 0];
a1 = [T12*T14*Tw  (T12*T14 + T12*Tw + T14*Tw) (T12 + T14 + Tw) 1];

b3 = [KG3*T31*T33*Tw (KG3*T31*Tw + KG3*T33*Tw) KG3*Tw 0];
a3 = [T32*T34*Tw  (T32*T34 + T32*Tw + T34*Tw) (T32 + T34 + Tw) 1];

b4 = [KG4*T41*T43*Tw (KG4*T41*Tw + KG4*T43*Tw) KG4*Tw 0];
a4 = [T42*T44*Tw  (T42*T44 + T42*Tw + T44*Tw) (T42 + T44 + Tw) 1];

b5 = [KG5*T51*T53*Tw (KG5*T51*Tw + KG5*T53*Tw) KG5*Tw 0];
a5 = [T52*T54*Tw  (T52*T54 + T52*Tw + T54*Tw) (T52 + T54 + Tw) 1];

b6 = [KG6*T61*T63*Tw (KG6*T61*Tw + KG6*T63*Tw) KG6*Tw 0];
a6 = [T62*T64*Tw  (T62*T64 + T62*Tw + T64*Tw) (T62 + T64 + Tw) 1];

b7 = [KG7*T71*T73*Tw (KG7*T71*Tw + KG7*T73*Tw) KG7*Tw 0];
a7 = [T72*T74*Tw  (T72*T74 + T72*Tw + T74*Tw) (T72 + T74 + Tw) 1];

b8 = [KG8*T81*T83*Tw (KG8*T81*Tw + KG8*T83*Tw) KG8*Tw 0];
a8 = [T82*T84*Tw  (T82*T84 + T82*Tw + T84*Tw) (T82 + T84 + Tw) 1];

b9 = [KG9*T91*T93*Tw (KG9*T91*Tw + KG9*T93*Tw) KG9*Tw 0];
a9 = [T92*T94*Tw  (T92*T94 + T92*Tw + T94*Tw) (T92 + T94 + Tw) 1];

b10 = [KG10*T101*T103*Tw (KG10*T101*Tw + KG10*T103*Tw) KG10*Tw 0];
a10 = [T102*T104*Tw  (T102*T104 + T102*Tw + T104*Tw) (T102 + T104 + Tw) 1];

[A_1 B_1 C_1 D_1 ]= tf2ss(b1,a1);
[A_2 B_2 C_2 D_2 ]= tf2ss(b3,a3);
[A_3 B_3 C_3 D_3 ]= tf2ss(b4,a4);
[A_4 B_4 C_4 D_4 ]= tf2ss(b5,a5);
[A_5 B_5 C_5 D_5 ]= tf2ss(b6,a6);
[A_6 B_6 C_6 D_6 ]= tf2ss(b7,a7);
[A_7 B_7 C_7 D_7 ]= tf2ss(b8,a8);
[A_8 B_8 C_8 D_8 ]= tf2ss(b9,a9);
[A_9 B_9 C_9 D_9 ]= tf2ss(b10,a10);


Af = blkdiag(A_1,A_2,A_3,A_4,A_5,A_6,A_7,A_8,A_9);
Bf = blkdiag(B_1,B_2,B_3,B_4,B_5,B_6,B_7,B_8,B_9);
Cf = blkdiag(C_1,C_2,C_3,C_4,C_5,C_6,C_7,C_8,C_9);
Df = blkdiag(D_1,D_2,D_3,D_4,D_5,D_6,D_7,D_8,D_9);

Asys_1 = As + Bs*Df*Cs;
Asys_2 = Bs*Cf;
Asys_3 = Bf*Cs;
Asys_4 = Af + Bf*Ds*Cf;
Asys = [Asys_1 Asys_2;
    Asys_3 Asys_4];

egs = eig(Asys);

[z_val z_idx]=sort(abs(egs),'descend');
egs_new=egs;
egs_new(z_idx(end-1:end))=[];

%% unstable modes
ss_idx = find(real(egs_new)>0);
uss = egs_new(ss_idx);

%% EM modes
% Damp=-real(egs)./sqrt(real(egs).^2+imag(egs).^2)
freq = abs(imag(egs_new))/(2*pi);
em_idx = find(freq>0 & freq<3);

objf = max(real(egs_new(em_idx)))+sum(egs_new(ss_idx));

if isempty(objf)
    objf = max(real(egs_new));
end
J = objf;
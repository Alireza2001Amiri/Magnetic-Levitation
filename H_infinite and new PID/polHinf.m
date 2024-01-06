% DG1    %LMI hinf-pole placement   %Feasp  
%clear
clc
%model--------------------------------------------------------------------
model
P = gmdc;
B2_1=P.b;    %Bii
B1_1=P.b;    %Bw1
C1_1=P.c;
C2_1=P.c;        %Cii
A1 = [];

%LMI Code-----------------------------------------------------------------
gama =10; %  shoul be selected optimal 
beta=1 ;
n = 2;
teta = acos(.5); % ia selected in desired value
st = sin(teta); ct = cos(teta);
nin = 1 ; nout = 1;
setlmis([ ] ) ;
    % Define the solution variables
    R = lmivar(1, [n 1]);
    S = lmivar(1, [n 1]);
    Ahat = lmivar(1, [n 1]);
    Bhat = lmivar (2, [n nout]);
    Chat = lmivar(2, [nin n]);
    Dhat = lmivar(1, [nout nin]);
    lmiterm([-1 1 1 R], 1,1);
    lmiterm([-1 2 1 0 ] , 1 ) ;
    lmiterm([-1 2 2 S], 1,1);
%LMI 1--------------------------------------------------------------------
A1i = P.a; 
r=2;
%Bi = BB{i};
%Ci = CC{i};
%rr=minreal(ss(A1i,B2_1,C2_1,0),1e-4);%rank(ctrb(rr.a,rr.b),1e-8)   
%A1i = rr.a ; B2_1 = rr.b ; C2_1 = rr.c;
Ast = A1i * st; Act = A1i * ct; B2st = B2_1 * st; B2ct = B2_1 * ct; C2st = C2_1 * st; C2ct =C2_1 * ct;
% Set up the LMIs
lmiterm([r 1 1 R], Ast, 1, 's');
lmiterm([r 1 1 Chat], B2st, 1 , 's');
lmiterm([r 1 1 R], 1, beta);
lmiterm([r 2 1 Ahat], 1, st);
lmiterm([r 2 1 0], Ast'+beta*eye(n));   %%%*****
lmiterm([r 2 1 Dhat'], C2st', B2_1');      % eslahi
lmiterm([r 2 2 S], 1, Ast, 's');
lmiterm([r 2 2 Bhat], 1, C2st, 's');
lmiterm([r 2 2 S], 1, beta);
lmiterm([r 3 1 R], 1, Act');       %%%%   *****
lmiterm([r 3 1 -Chat], 1, B2ct');
lmiterm([r 3 1 R], Act, -1);
lmiterm([r 3 1 Chat], B2ct, -1);
lmiterm([r 3 2 -Ahat], 1, ct);
lmiterm([r 3 2 0], -Act);
lmiterm([r 3 3 R], Ast, 1, 's');
lmiterm([r 3 3 Chat], B2st, 1,'s');
lmiterm([r 3 3 R], 1, beta);
lmiterm([r 4 1 Ahat], 1, -ct);
lmiterm([r 4 1 0], Act');
lmiterm([r 4 1 -Dhat'], C2ct', B2_1');      % eslahi
lmiterm([r 4 2 S], Act', 1);
lmiterm([r 4 2 -Bhat], C2ct', 1);
lmiterm([r 4 2 S], 1, -Act);
lmiterm([r 4 2 Bhat],1 , -C2ct);
lmiterm([r 4 3 Ahat],1 , st ) ;
lmiterm([r 4 3 0], Ast'+beta*eye(n));
lmiterm([r 4 3 Dhat'], C2st', B2_1');      % eslahi
lmiterm([r 4 4 S], 1, Ast, 's');
lmiterm([r 4 4 Bhat], 1 , C2st,'s');
lmiterm([r 4 4 S], 1, beta);

lmiterm([r+1 1 1 R],A1i,1,'s');                     % LMI #1: A*R+R*A'
lmiterm([r+1 1 1 Chat],B2_1,1,'s');                  % LMI #1: B2*CKK+CKK'*B2'
lmiterm([r+1 2 1 0],B1_1');                         % LMI #1: B1'
lmiterm([r+1 2 2 0],-gama*1);                     % LMI #1: -gama*1
lmiterm([r+1 3 1 Ahat],1,1);                       % LMI #1: AKK
lmiterm([r+1 3 1 0],A1i');                          % LMI #1: A'
lmiterm([r+1 3 2 S],1,B1_1);                        % LMI #1: S*B1
lmiterm([r+1 3 3 S],A1i',1,'s');                    % LMI #1: A'*S+S*A
lmiterm([r+1 3 3 Bhat],1,C2_1,'s');                  % LMI #1: BKK*C2+C2'*BKK'
lmiterm([r+1 4 1 R],C1_1,1);                        % LMI #1: C1*R
lmiterm([r+1 4 3 0],C1_1);                          % LMI #1: C1
lmiterm([r+1 4 4 0],-gama*1);                     % LMI #1: -gama*1


% System of LMIs---------------------------------------------------------
syslmi = getlmis;

% Minimize the cost function---------------------------------------------
% objective = alphagamma2 + eps * Trace(R + S)
options=[0,500,1e9,0,0];
[tmin,xfeas] = feasp(syslmi,options);
% Retrieve the LMI variables form the solution
R=dec2mat(syslmi,xfeas,R);
S=dec2mat(syslmi,xfeas,S);
Ahat=dec2mat(syslmi,xfeas,Ahat);
Bhat=dec2mat(syslmi,xfeas,Bhat);
Chat=dec2mat(syslmi,xfeas,Chat);
Dhat=dec2mat(syslmi,xfeas,Dhat);
% Determine Mti(= inv(M')and Ni(= inv(N))from SVD of MN' = I-R*S----------
[u,sd,v]=svd(eye(n)-R*S);% factorize I-RS
M =u; N=v*sd'; Ni=inv(N); Mti=inv(M'); Mt=M';
% Retrieve the controller parameters from the transformed variables
Ck1=Chat*Mti-Dhat*C2_1*R*Mti;
Bk1=Ni*Bhat-Ni*S*B2_1*Dhat;
Ak1=Ni*(Ahat-N*Bk1*C2_1*R-S*B2_1*Ck1*Mt-S*(A1i+B2_1*Dhat*C2_1)*R)*Mti;
% Controller state space--------------------------------------------------
K=ss(Ak1,Bk1,Ck1,Dhat);
Ks=ss2tf(Ak1,Bk1,Ck1,Dhat);
[a,b]=ss2tf(Ak1,Bk1,Ck1,Dhat);
kk1 = -tf(a,b);  %Controller without D-part

kk2 = inv(1-kk1*P.d)*kk1;  % Controller with D-part

[kkb,gg] = balreal(kk2);
kkbr = modred(kkb,3:4,'MatchDC');
bodeplot(kk2,'-',kkbr,'x')
[a3,b3]=ss2tf(kkbr.a,kkbr.b,kkbr.c,kkbr.d); 
kk3 = tf(a3,b3);                   % Reduction of controoler to two-order

kkkbr = modred(kkb,2:4,'MatchDC');
bodeplot(kk2,'-',kkkbr,'x')
[a4,b4]=ss2tf(kkkbr.a,kkkbr.b,kkkbr.c,kkkbr.d);
kk4 = tf(a4,b4);


% function SolitonWall_04_10_10_Box(n)
%         % Try-catch expression that quits the Matlab session if your code crashes
%         try
%                 % Initialize the parallel pool
%                 c=parcluster();
%                 t=tempname();
%                 mkdir(t)
%                 c.JobStorageLocation=t;
%                 c.NumWorkers=24;
%                 parpool(c,n);
%                 % The actual program calls
%                 SolitonWall2D_BFGS_Para;
%                 delete(gcp('nocreate'));
%         catch error
%                 getReport(error)
%                 disp('Error occured');
%                 exit(0)
%         end
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%      Main Program Script     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%function SolitonWall2D_BFGS_Para

clear all


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

xiDt=1; %2*xiDt^2=xiD^2
Deltap=1; 
addpath('ClassCalculateA2lijInOneElement'); %% add library

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[nodes,elements,Np,Nelem]=TriPartition2DSW;
Np
Nelem

%[Rbu,Cbu,Rbd,Cbd,Rbiu,Cbiu,Rbid,Cbid,RbIv,CbIv,Rv,Cv]=ClassifyNodes(nodes);
 [RbW,CbW,Rbu,Cbu,Rbd,Cbd,Rbiu,Cbiu,Rbid,Cbid,RbIv,CbIv,Rv,Cv]=ClassifyNodes(nodes);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% Boundary Value %%%%%%%%%

 Deltav1=0.11*Deltap;
 gamma=Deltav1/Deltap;
gamma1=(1+gamma^2)/((1+gamma)^2); % Using Jaakko's parametrization, which is the lowest Fsoc paramentrization
gamma2=(gamma^2)/((1+gamma)^2);
G3=gamma/(1+gamma);

thetaBVwall=pi/2; % statonary point of 2/4 pi-soliton,
%thetaBVup=asin(G3/2); % value of theta in y>=0
thetaBVup=pi-asin(G3/2); % value of theta in y>=0, +|Deltav2|
thetaBVdown=asin(-G3/2); % value of theta in y<=0, -|Deltav2|

thetaBVBIup=pi-asin(G3/2); % value of theta in y>=0, BV for 2/4 pi-Soliton of inseparable soliton 
thetaBVBIdown=asin(-G3/2); % value of theta in y<=0, BV for 2/4 pi-Soliton part of inseparable soliton 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

A=sparse(Np,Np); %% the exact A Matrix corressponds to original nodes matrix+extra-intra-boundary point matrix Cbi

thetakn1t=[ones(length(Cv),1); ones(length(CbIv),1);]; %exact function point variable
thetabu=sparse(ones(length(Cbu),1)*thetaBVup);
thetabd=sparse(ones(length(Cbd),1)*thetaBVdown);
thetabw=sparse(ones(length(CbW),1)*thetaBVwall);

thetabiu=sparse(ones(length(Cbiu),1)*thetaBVBIup);
thetabid=sparse(ones(length(Cbid),1)*thetaBVBIdown);

% thetabiu=sparse(ones(length(Cbi),1)*thetaBVup);
% thetabid=sparse(ones(length(Cbi),1)*thetaBVdown);

%alpha=[alphabu; alphabd; alphabiu; alphabid; alphakn1t;]; %%easy to save and take the value
Thetakn1=zeros(Np,1); %% the ture function point matrix corresponding to nodes matrix+extra-intra-boundary point matrix Cbi 
%Thetakn1=TransferVtoTheta2DSolitonWallQN(thetabu,thetabd,thetabiu,thetabid,thetakn1t,Cbu,Cbd,Cbi,CbIv,Cv,Thetakn1,Np,1);
Thetakn1=TransferVtoTheta2DSolitonWallQN(thetabw,thetabu,thetabd,thetabiu,thetabid,thetakn1t,CbW,Cbu,Cbd,Cbiu,Cbid,CbIv,Cv,Thetakn1,Np,1);


Is=sparse(eye(length(Cv)+length(CbIv)));
Hkn1=Is;

tol=1e-9;
D=[];
o=1;
O=[];
c=0;
C=[];
Lambda=[];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% A Matrix %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


for l=1:Nelem
    
    % function object
     functionA2lij=clib.ClassCalculateA2lijInOneElement.CalculateA2lij;

    % call the member function to calculate
     functionA2lij.CalculateA2lijInOneElement(elements,nodes,gamma1,gamma2,l);

    % return result
     A2lij=reshape(functionA2lij.returnResult(9),3,3);
    %functionA2lij.returnResult(9)

    clear functionA2lij A2lij;
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   
      for u=1:3 %% The lth element locaizes inside Domain, us the original index 
  
          A(El(u),El(1))=A2lij(u,1)+A(El(u),El(1));
   
          A(El(u),El(2))=A2lij(u,2)+A(El(u),El(2));
       
          A(El(u),El(3))=A2lij(u,3)+A(El(u),El(3));
       

      end

end

A=sparse(A);

%%%%%%%%%%%%%%%% Trnasfer the alpha to Alpha %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

X=nodes(1,:);Y=nodes(2,:);

gkn1=gradient2DdSolitonWallQN(Cv,CbIv,A,Thetakn1,elements,nodes,Np,gamma1,gamma);
                             
while norm(gkn1)>tol
      c=c+1;
      C(c)=c;
      c
      
      O(c)=norm(gkn1);
      ABSgkn1=norm(gkn1)

      if norm(gkn1)<=1e-2
         save('SoliTonWall_PdBLattice_box25_5_gamma011.mat','X','Y','Thetakn1','gamma','Rbu','Cbu','Rbd','Cbd','Rbid','Cbid','Rbiu','Cbiu','RbIv','CbIv','Rv','Cv',...
             'nodes','elements','Np','Nelem','gamma','ABSgkn1');
      end      

      dkn1=-Hkn1*gkn1;
      lambdakn1=fminsearch(@(lambda) PiThetaLambda2DSolitonWallQN(lambda,A,Thetakn1,dkn1,Cv,CbIv,elements,Nelem,nodes,Np,gamma1,gamma),0);
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      CV=[Cv CbIv];
      thetakn1t=zeros(length(CV),1);
      for k=1:length(CV) 
          thetakn1t(k,1)=Thetakn1(CV(k));
      end
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
      thetakt=thetakn1t+lambdakn1.*dkn1;
      Thetak=TransferVtoTheta2DSolitonWallQN(thetabw,thetabu,thetabd,thetabiu,thetabid,thetakt,CbW,Cbu,Cbd,Cbiu,Cbid,CbIv,Cv,Thetakn1,Np,2);
                                          
      skn1=thetakt-thetakn1t;
      gk=gradient2DdSolitonWallQN(Cv,CbIv,A,Thetak,elements,nodes,Np,gamma1,gamma);
      ykn1=gk-gkn1;
      rhokn1=1./((ykn1')*skn1);
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
      Hk=InHassian2DSolitonQN(Is,rhokn1,skn1,ykn1,Hkn1);
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
      Hkn1=Hk;
      Thetakn1=Thetak;
      gkn1=gradient2DdSolitonWallQN(Cv,CbIv,A,Thetakn1,elements,nodes,Np,gamma1,gamma);
          
      Lambda(c)=lambdakn1; 
end      

%end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [nodes,elements,Np,Nelem]=TriPartition2DSW

model=createpde(1);
geometryFromEdges(model,@SolveDomain);
M=generateMesh(model,'JiggleIter',10,'Hgrad',1.8,'Hmax',2.5);
nodes=M.Nodes;
elements=M.Elements;
Np=length(nodes);
Nelem=length(elements);


end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [x,y]=SolveDomain(bs,s)
xiDt=1;nbs=15;S=25.5000;D1=12.6500;D2=0.1;H1=0.01;H2=0.1;
if nargin==0,
  x=nbs; % number of boundary segments
  return
end

d=[
  0 0 0 0 0 0 0 0 0 0 0 0 0 0 0% start parameter value
  S/2 S/2 S (S/2)-H1 (S/2)-H1 S D1 H2 D2 H2 S-(D1+D2) S-(D1+D2) H2 D2 H1% end parameter value
  1 2 2 2 1 1 1 1 1 1 1 2 2 2 2% left hand region
  0 0 0 0 0 0 2 0 0 0 0 0 0 0 0% right hand region
];

bs1=bs(:)';

if find(bs1<1 | bs1>nbs),
  error(message('pde:lshapeg:InvalidBs'))
end

if nargin==1,
  x=d(:,bs1);
  return
end

x=zeros(size(s));
y=zeros(size(s));
[m,n]=size(bs);
if m==1 && n==1,
  bs=bs*ones(size(s)); % expand bs
elseif m~=size(s,1) || n~=size(s,2),
  error(message('pde:lshapeg:SizeBs'));
end

if ~isempty(s),

% boundary segment 1
ii=find(bs==1);
if length(ii)
x(ii)=interp1([d(1,1),d(2,1)],[0 0],s(ii));
y(ii)=interp1([d(1,1),d(2,1)],[S/2 0],s(ii));
end

% boundary segment 2
ii=find(bs==2);
if length(ii)
x(ii)=interp1([d(1,2),d(2,2)],[0 0],s(ii));
y(ii)=interp1([d(1,2),d(2,2)],[0 -S/2],s(ii));
end

% boundary segment 3
ii=find(bs==3);
if length(ii)
x(ii)=interp1([d(1,3),d(2,3)],[0 S],s(ii));
y(ii)=interp1([d(1,3),d(2,3)],[-S/2 -S/2],s(ii));
end

% boundary segment 4
ii=find(bs==4);
if length(ii)
x(ii)=interp1([d(1,4),d(2,4)],[S S],s(ii));
y(ii)=interp1([d(1,4),d(2,4)],[-S/2 -H1],s(ii));
end

% boundary segment 5
ii=find(bs==5);
if length(ii)
x(ii)=interp1([d(1,5),d(2,5)],[S S],s(ii));
y(ii)=interp1([d(1,5),d(2,5)],[H1 S/2],s(ii));
end

% boundary segment 6
ii=find(bs==6);
if length(ii)
x(ii)=interp1([d(1,6),d(2,6)],[S 0],s(ii));
y(ii)=interp1([d(1,6),d(2,6)],[S/2 S/2],s(ii));
end

%boundary segment 7
ii=find(bs==7);
if length(ii)
x(ii)=interp1([d(1,7),d(2,7)],[0 D1],s(ii));
y(ii)=interp1([d(1,7),d(2,7)],[0 0],s(ii));
end

%boundary segment 8
ii=find(bs==8);
if length(ii)
x(ii)=interp1([d(1,8),d(2,8)],[D1 D1],s(ii));
y(ii)=interp1([d(1,8),d(2,8)],[0 H2],s(ii));
end

%boundary segment 9
ii=find(bs==9);
if length(ii)
x(ii)=interp1([d(1,9),d(2,9)],[D1 D1+D2],s(ii));
y(ii)=interp1([d(1,9),d(2,9)],[H2 H2],s(ii));
end

%boundary segment 10
ii=find(bs==10);
if length(ii)
x(ii)=interp1([d(1,10),d(2,10)],[D1+D2 D1+D2],s(ii));
y(ii)=interp1([d(1,10),d(2,10)],[H2 H1],s(ii));
end

%boundary segment 11
ii=find(bs==11);
if length(ii)
x(ii)=interp1([d(1,11),d(2,11)],[D1+D2 S],s(ii));
y(ii)=interp1([d(1,11),d(2,11)],[H1 H1],s(ii));
end

%boundary segment 12
ii=find(bs==12);
if length(ii)
x(ii)=interp1([d(1,12),d(2,12)],[S D1+D2],s(ii));
y(ii)=interp1([d(1,12),d(2,12)],[-H1 -H1],s(ii));
end

%boundary segment 13
ii=find(bs==13);
if length(ii)
x(ii)=interp1([d(1,13),d(2,13)],[D1+D2 D1+D2],s(ii));
y(ii)=interp1([d(1,13),d(2,13)],[-H1 -H2],s(ii));
end

%boundary segment 14
ii=find(bs==14);
if length(ii)
x(ii)=interp1([d(1,14),d(2,14)],[D1+D2 D1],s(ii));
y(ii)=interp1([d(1,14),d(2,14)],[-H2 -H2],s(ii));
end

%boundary segment 15
ii=find(bs==15);
if length(ii)
x(ii)=interp1([d(1,15),d(2,15)],[D1 D1],s(ii));
y(ii)=interp1([d(1,15),d(2,15)],[-H2 0],s(ii));
end

end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [RbW,CbW,Rbu,Cbu,Rbd,Cbd,Rbiu,Cbiu,Rbid,Cbid,RbIv,CbIv,Rv,Cv]=ClassifyNodes(nodes)
xiDt=1;nbs=15;S=25.5000;D1=12.6500;D2=0.1;H1=0.01;H2=0.1;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[Rbu,Cbu]=find(nodes(2,:)>0 & abs(nodes(2,:))>=S/2);

[Rbd,Cbd]=find(nodes(2,:)<0 & abs(nodes(2,:))>=S/2);

[Rbiu,Cbiu]=find(nodes(2,:)>0 & abs(nodes(2,:))==H1 & abs(nodes(1,:))>=D1+D2);

[Rbid,Cbid]=find(nodes(2,:)<0 & abs(nodes(2,:))==H1 & abs(nodes(1,:))>=D1+D2);

[RbW,CbW]=find(abs(nodes(2,:))==0 & abs(nodes(1,:))<D1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[RbIv,CbIv]=find( (nodes(1,:)==D1 & abs(nodes(2,:))<H2) |...
                  (abs(nodes(2,:))==H2 & nodes(1,:)>=D1 & nodes(1,:)<=D1+D2) |...
                  (nodes(1,:)==D1+D2 & abs(nodes(2,:))<H2 & abs(nodes(2,:))>H1));

[Rv,Cv]=find(nodes(2,:)~=0 & abs(nodes(2,:))<S/2 & ...
              ( (abs(nodes(2,:))>H2) |... 
                (nodes(1,:)<D1 & abs(nodes(2,:))<=H2) | ...
                (abs(nodes(2,:))>H1 & nodes(1,:)>D1+D2 & nodes(1,:)<=S & abs(nodes(2,:))<=H2) ));

for l=1:length(Cbid) %% align the Cbiu and Cbid
    a(:,l)=nodes(:,Cbid(l));
    b(:,l)=nodes(:,Cbiu(l));
end  
aa(1,:)=sort(a(1,:),'ascend');bb(1,:)=sort(b(1,:),'ascend');
for l=1:length(aa)
    aaa(l)=Cbid(find(a(1,:)==aa(l)));
    bbb(l)=Cbiu(find(b(1,:)==bb(l)));    
end   
Cbid=aaa;Cbiu=bbb;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Alpha=TransferVtoTheta2DSolitonWallQN(alphabW,alphabu,alphabd,alphabiu,alphabid,alphat,CbW,Cbu,Cbd,Cbiu,Cbid,CbIv,Cv,Alpha,Np,option)
                                               
if option==1 % Constrcut the total function value vector via all boundary values and variable values, usually be used at very begining
   
   for l=1:length(CbW)

       Alpha(CbW(l))=alphabW(l);   
   end   


   for l=1:length(Cbu)

       Alpha(Cbu(l))=alphabu(l);   
   end   
   
   for l=1:length(Cbd)
       Alpha(Cbd(l))=alphabd(l);   
   end  
    
   for l=1:length(Cbid)  
       Alpha(Cbid(l))=alphabid(l);   
   end  
   
   for l=1:length(Cbiu) %% the value of theta on Cbiu point
       Alpha(Cbiu(l))=alphabiu(l);  
   end
    
   CV=[Cv CbIv];
   for l=1:length(CV)
       Alpha(CV(l))=alphat(l);   
   end  
    
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
   
elseif option==2 % Refresh the total function values vector only via temporary variable values, usually be used in every BFGS iteration 
    
%        for l=1:length(Cbid)  
%            Alpha(Cbid(l))=alphat(length(Cv)+length(CbIv)+l);   
%        end  
%    
%        for l=1:length(Cbiu) %% the value of theta on Cbiu point
%            Alpha(Cbiu(l))=alphat(length(Cv)+length(CbIv)+l)+pi;   
%        end
    
       CV=[Cv CbIv];
       for l=1:length(CV)
           Alpha(CV(l))=alphat(l);   
       end  
        
end    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function A=AMatix2DSolitonWallBFGS(Nelem,elements,nodes,A,Np,gamma1,gamma2)
% 
% 
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function gkn1=gradient2DdSolitonWallQN(Cv,CbIv,A,Alpha,elements,nodes,Np,gamma1,gamma)

CV=[Cv CbIv];

Av=sparse(zeros(length(CV),Np));
%Av=zeros(length(CV),Np+length(Cbi));
for k=1:length(CV)
    
%     if k<=length(Cv)+length(CbIv)
       Av(k,:)=A(CV(k),:);
%     else
%         Av(k,:)=A(CV(k),:)+A(Cbiu(k-(length(Cv)+length(CbIv))),:);
%     end   
end  
dPidnuA=2*Av*Alpha;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dPidnuB=zeros(length(CV),1);
% ForLoopIndex2=length(Cbid);
parfor k=1:length(CV)
%     if k<=length(Cv)+length(CbIv) %% normal variable
       [x,y]=find(elements==CV(k));
       ForLoopIndex1=length(y);
    
       dPidnuBk=0;
       for l=1:ForLoopIndex1
        
            El=[]; n1=[]; n2=[]; n3=[]; x1=[]; y1=[]; x2=[]; y2=[]; x3=[]; y3=[]; xi1=[]; xi2=[]; xi3=[]; eta=[]; xi=[]; omega=[]; D=[]; Si=[];
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            El=elements(:,y(l));

            n1=nodes(:,El(1));x1=n1(1);y1=n1(2);
            n2=nodes(:,El(2));x2=n2(1);y2=n2(2);
            n3=nodes(:,El(3));x3=n3(1);y3=n3(2);
    
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%% gamma1 and gamma3 %%%%%%%%%%%%%
        
            gamma3=gamma3DSW2D(y1,y2,y3,gamma);
       
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
            xi=[x2-x3 x3-x1 x1-x2];
            eta=[y2-y3 y3-y1 y1-y2];
            omega=[x2*y3-x3*y2 x3*y1-x1*y3 x1*y2-x2*y1];
    
            Di=sum(omega);
            Si=abs(Di)/2;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       
            ALPHA=zeros(3,1);
        
%           if isempty(e) %% The triangle element does not contain extra-intra-boundary point, use the orignal Nodes Matrix index 
        
            ALPHA(x(l))=Alpha(El(x(l)));   
           
            [xx,yy]=find(El~=El(x(l)));   
            ALPHA(xx(1))=Alpha(El(xx(1))); 
            ALPHA(xx(2))=Alpha(El(xx(2))); 
                      
            dPidnuBk=dPidnuBk+TriEleIntegral2DSolitonWall(x1,x2,x3,y1,y2,y3,ALPHA,x(l),xi,eta,omega,Di,gamma1,gamma3);
          % ALPHA
         
       end
    
      
    
    dPidnuB(k,1)=dPidnuBk;
    
end  

gkn1=dPidnuA+dPidnuB;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dPidnuBkl=TriEleIntegral2DSolitonWall(x1,x2,x3,y1,y2,y3,ALPHA,xl,xi,eta,omega,Di,gamma1,gamma3)

xidt=1;

% integrant={@(x,y)(1/4)*sin(2.*((eta(1).*x-xi(1).*y+omega(1))./Di).*ALPHA(1)+2.*((eta(2).*x-xi(2).*y+omega(2))./Di).*ALPHA(2)+2.*((eta(3).*x-xi(3).*y+omega(3))./Di).*ALPHA(3)).*2.*((eta(1).*x-xi(1).*y+omega(1))./Di)...
%            @(x,y)(1/4)*sin(2.*((eta(1).*x-xi(1).*y+omega(1))./Di).*ALPHA(1)+2.*((eta(2).*x-xi(2).*y+omega(2))./Di).*ALPHA(2)+2.*((eta(3).*x-xi(3).*y+omega(3))./Di).*ALPHA(3)).*2.*((eta(2).*x-xi(2).*y+omega(2))./Di)...
%            @(x,y)(1/4)*sin(2.*((eta(1).*x-xi(1).*y+omega(1))./Di).*ALPHA(1)+2.*((eta(2).*x-xi(2).*y+omega(2))./Di).*ALPHA(2)+2.*((eta(3).*x-xi(3).*y+omega(3))./Di).*ALPHA(3)).*2.*((eta(3).*x-xi(3).*y+omega(3))./Di)};

integrant={@(x,y)((1/2)*1).*((eta(1).*x-xi(1).*y+omega(1))./Di).*(-gamma3.*cos(1.*((eta(1).*x-xi(1).*y+omega(1))./Di).*ALPHA(1)+1.*((eta(2).*x-xi(2).*y+omega(2))./Di).*ALPHA(2)+1.*((eta(3).*x-xi(3).*y+omega(3))./Di).*ALPHA(3))+...
                  sin(2.*((eta(1).*x-xi(1).*y+omega(1))./Di).*ALPHA(1)+2.*((eta(2).*x-xi(2).*y+omega(2))./Di).*ALPHA(2)+2.*((eta(3).*x-xi(3).*y+omega(3))./Di).*ALPHA(3)))...
           @(x,y)((1/2)*1).*((eta(2).*x-xi(2).*y+omega(2))./Di).*(-gamma3.*cos(1.*((eta(1).*x-xi(1).*y+omega(1))./Di).*ALPHA(1)+1.*((eta(2).*x-xi(2).*y+omega(2))./Di).*ALPHA(2)+1.*((eta(3).*x-xi(3).*y+omega(3))./Di).*ALPHA(3))+...
                  sin(2.*((eta(1).*x-xi(1).*y+omega(1))./Di).*ALPHA(1)+2.*((eta(2).*x-xi(2).*y+omega(2))./Di).*ALPHA(2)+2.*((eta(3).*x-xi(3).*y+omega(3))./Di).*ALPHA(3)))...
           @(x,y)((1/2)*1).*((eta(3).*x-xi(3).*y+omega(3))./Di).*(-gamma3.*cos(1.*((eta(1).*x-xi(1).*y+omega(1))./Di).*ALPHA(1)+1.*((eta(2).*x-xi(2).*y+omega(2))./Di).*ALPHA(2)+1.*((eta(3).*x-xi(3).*y+omega(3))./Di).*ALPHA(3))+...
                  sin(2.*((eta(1).*x-xi(1).*y+omega(1))./Di).*ALPHA(1)+2.*((eta(2).*x-xi(2).*y+omega(2))./Di).*ALPHA(2)+2.*((eta(3).*x-xi(3).*y+omega(3))./Di).*ALPHA(3)))};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Xe=[x1 x2 x3];Ye=[y1 y2 y3];

if (x1~=x2) && (x1~=x3) && (x2~=x3) %% every side have finite slope
   
    RPXe=sort(Xe);  %fix the middle point
    [rx,cx]=find(Xe==RPXe(2)); % find the position of middle point via x-cooordinate
    
    [ry,cy]=find(Ye==Ye(cx)); %check the position and number of y-coordinates equaling to middle point y-coordinate
    
    if length(cy)==1  %% no any other point has same y-coordinate with middle point
       
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        RPYe=sort(Ye); 
        [Ry,Cy]=find(RPYe==Ye(cx)); %% judge whether the middle point is largest or smallest
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if Cy==3 %% head upper, middle point maxinum
           Kother=find(Xe~=Xe(cx));
           if ((Ye(Kother(2))-Ye(Kother(1)))/((Xe(Kother(2))-Xe(Kother(1)))))==0
               Llow=Ye(Kother(1));
           else    
               Llow=@(x)(x-Xe(Kother(1))).*((Ye(Kother(2))-Ye(Kother(1)))/((Xe(Kother(2))-Xe(Kother(1)))))+Ye(Kother(1)); 
           end
           
           switch Xe(Kother(1))>Xe(Kother(2))
               case 1
                    Lup1=@(x)(x-Xe(cx)).*((Ye(Kother(2))-Ye(cy))/(Xe(Kother(2))-Xe(cx)))+Ye(cy);
                    Lup2=@(x)(x-Xe(cx)).*((Ye(Kother(1))-Ye(cy))/(Xe(Kother(1))-Xe(cx)))+Ye(cy);
                    xlow1=Xe(Kother(2));xup1=Xe(cx);
                    xlow2=Xe(cx);xup2=Xe(Kother(1));
               case 0                 
                    Lup1=@(x)(x-Xe(cx)).*((Ye(Kother(1))-Ye(cy))/(Xe(Kother(1))-Xe(cx)))+Ye(cy);
                    Lup2=@(x)(x-Xe(cx)).*((Ye(Kother(2))-Ye(cy))/(Xe(Kother(2))-Xe(cx)))+Ye(cy);
                    xlow1=Xe(Kother(1));xup1=Xe(cx);
                    xlow2=Xe(cx);xup2=Xe(Kother(2));
           end 
           dPidnuBkl=integral2(integrant{xl},xlow1,xup1,Llow,Lup1)+integral2(integrant{xl},xlow2,xup2,Llow,Lup2,'Reltol',1e-9);   
        end 
       
        if Cy==1 %% head down, middle point minium
           Kother=find(Xe~=Xe(cx));
           if ((Ye(Kother(2))-Ye(Kother(1)))/((Xe(Kother(2))-Xe(Kother(1)))))==0
               Lup=Ye(Kother(1));
           else    
               Lup=@(x)(x-Xe(Kother(1))).*((Ye(Kother(2))-Ye(Kother(1)))/((Xe(Kother(2))-Xe(Kother(1)))))+Ye(Kother(1)); 
           end
         
           switch Xe(Kother(1))>Xe(Kother(2))
               case 1
                    Llow1=@(x)(x-Xe(cx)).*((Ye(Kother(2))-Ye(cy))/(Xe(Kother(2))-Xe(cx)))+Ye(cy);
                    Llow2=@(x)(x-Xe(cx)).*((Ye(Kother(1))-Ye(cy))/(Xe(Kother(1))-Xe(cx)))+Ye(cy);
                    xlow1=Xe(Kother(2));xup1=Xe(cx);
                    xlow2=Xe(cx);xup2=Xe(Kother(1));
               case 0                 
                    Llow1=@(x)(x-Xe(cx)).*((Ye(Kother(1))-Ye(cy))/(Xe(Kother(1))-Xe(cx)))+Ye(cy);
                    Llow2=@(x)(x-Xe(cx)).*((Ye(Kother(2))-Ye(cy))/(Xe(Kother(2))-Xe(cx)))+Ye(cy);
                    xlow1=Xe(Kother(1));xup1=Xe(cx);
                    xlow2=Xe(cx);xup2=Xe(Kother(2));
           end 
           dPidnuBkl=integral2(integrant{xl},xlow1,xup1,Llow1,Lup)+integral2(integrant{xl},xlow2,xup2,Llow2,Lup,'Reltol',1e-9);   
        end
       
        if Cy==2 %% the middle point is the 2nd large
           [Rrx,Crx]=find(Xe~=Xe(cx));[Rry,Cry]=find(Ye~=Ye(cy)); 
           %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
           %%% Need an Equation stringnt line, judge the relative position of line and middle point   
           y=(Xe(cx)-Xe(Crx(2))).*((Ye(Cry(1))-Ye(Cry(2)))/(Xe(Crx(1))-Xe(Crx(2))))+Ye(Cry(2));
           Clarge=y>Ye(cy);Clow=y<Ye(cy);
           %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
           
           if Clarge==0 && Clow==1    %% middle point higher than line, The way of integral same Cy=3
              Kother=find(Xe~=Xe(cx));
              if ((Ye(Kother(2))-Ye(Kother(1)))/((Xe(Kother(2))-Xe(Kother(1)))))==0
                  Llow=Ye(Kother(1));
              else    
                  Llow=@(x)(x-Xe(Kother(1))).*((Ye(Kother(2))-Ye(Kother(1)))/((Xe(Kother(2))-Xe(Kother(1)))))+Ye(Kother(1)); 
              end
     
              switch Xe(Kother(1))>Xe(Kother(2))
                      case 1
                           Lup1=@(x)(x-Xe(cx)).*((Ye(Kother(2))-Ye(cy))/(Xe(Kother(2))-Xe(cx)))+Ye(cy);
                           Lup2=@(x)(x-Xe(cx)).*((Ye(Kother(1))-Ye(cy))/(Xe(Kother(1))-Xe(cx)))+Ye(cy);
                           xlow1=Xe(Kother(2));xup1=Xe(cx);
                           xlow2=Xe(cx);xup2=Xe(Kother(1));
                      case 0                 
                           Lup1=@(x)(x-Xe(cx)).*((Ye(Kother(1))-Ye(cy))/(Xe(Kother(1))-Xe(cx)))+Ye(cy);
                           Lup2=@(x)(x-Xe(cx)).*((Ye(Kother(2))-Ye(cy))/(Xe(Kother(2))-Xe(cx)))+Ye(cy);
                           xlow1=Xe(Kother(1));xup1=Xe(cx);
                           xlow2=Xe(cx);xup2=Xe(Kother(2));
               end 
               dPidnuBkl=integral2(integrant{xl},xlow1,xup1,Llow,Lup1)+integral2(integrant{xl},xlow2,xup2,Llow,Lup2,'Reltol',1e-9);  
                
            elseif Clarge==1 && Clow==0  %% middle point lower than line,  The way of integral same Cy=1
                   Kother=find(Xe~=Xe(cx));
                   if ((Ye(Kother(2))-Ye(Kother(1)))/((Xe(Kother(2))-Xe(Kother(1)))))==0
                      Lup=Ye(Kother(1));
                   else    
                      Lup=@(x)(x-Xe(Kother(1))).*((Ye(Kother(2))-Ye(Kother(1)))/((Xe(Kother(2))-Xe(Kother(1)))))+Ye(Kother(1)); 
                   end
               
                   switch Xe(Kother(1))>Xe(Kother(2))
                         case 1
                              Llow1=@(x)(x-Xe(cx)).*((Ye(Kother(2))-Ye(cy))/(Xe(Kother(2))-Xe(cx)))+Ye(cy);
                              Llow2=@(x)(x-Xe(cx)).*((Ye(Kother(1))-Ye(cy))/(Xe(Kother(1))-Xe(cx)))+Ye(cy);
                              xlow1=Xe(Kother(2));xup1=Xe(cx);
                              xlow2=Xe(cx);xup2=Xe(Kother(1));
                         case 0                 
                              Llow1=@(x)(x-Xe(cx)).*((Ye(Kother(1))-Ye(cy))/(Xe(Kother(1))-Xe(cx)))+Ye(cy);
                              Llow2=@(x)(x-Xe(cx)).*((Ye(Kother(2))-Ye(cy))/(Xe(Kother(2))-Xe(cx)))+Ye(cy);
                              xlow1=Xe(Kother(1));xup1=Xe(cx);
                              xlow2=Xe(cx);xup2=Xe(Kother(2));
                   end 
                  dPidnuBkl=integral2(integrant{xl},xlow1,xup1,Llow1,Lup)+integral2(integrant{xl},xlow2,xup2,Llow2,Lup,'Reltol',1e-9);    
            end
         end   
    elseif length(cy)==2
           Cother=find(cy~=cx);
           Xo=Xe(cy(Cother));Yo=Ye(cy(Cother)); %% find out the other same high with middle point
           %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
           %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
           Crx=find(Xe~=Xe(cx) & Xe~=Xo);Cry=find(Ye~=Ye(cx) & Ye~=Yo); 
           %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
           % %% Need an Equation stringnt line 
           y=(Xe(cx)-Xe(Crx)).*((Yo-Ye(Cry))/(Xo-Xe(Crx)))+Ye(Cry);
           Clarge=y>Ye(cx);Clow=y<Ye(cx);
           %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
           %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
           if Clarge==0 && Clow==1    %% middle point higher than line, The way of integral same Cy=3
              Kother=find(Xe~=Xe(cx));identy=0; %% 
                  if ((Ye(Kother(2))-Ye(Kother(1)))/((Xe(Kother(2))-Xe(Kother(1)))))==0
                      Llow=Ye(Kother(1));
                      identy=1;
                  else    
                      Llow=@(x)(x-Xe(Kother(1))).*((Ye(Kother(2))-Ye(Kother(1)))/((Xe(Kother(2))-Xe(Kother(1)))))+Ye(Kother(1)); 
                  end
                  
                  switch Xe(Kother(1))>Xe(Kother(2))  
                      case 1 %% Xe(Kother(1)) ahead of Xe(Kother(2))
                           if Ye(Kother(1))==Ye(cx) %% Xe(Kother(1)) same high with moddle point
                              Lup1=@(x)(x-Xe(cx)).*((Ye(Kother(2))-Ye(cx))/(Xe(Kother(2))-Xe(cx)))+Ye(cx);
                              Lup2=Ye(cx);
                              xlow1=Xe(Kother(2));xup1=Xe(cx);
                              xlow2=Xe(cx);xup2=Xe(Kother(1));
                              Identy1=0;Identy2=0;
                              
                           elseif Ye(Kother(2))==Ye(cx) %% Xe(Kother(2)) same high with moddle point
                                  Lup1=Ye(cx);
                                  Lup2=@(x)(x-Xe(cx)).*((Ye(Kother(1))-Ye(cx))/(Xe(Kother(1))-Xe(cx)))+Ye(cx);
                                  xlow1=Xe(Kother(2));xup1=Xe(cx);
                                  xlow2=Xe(cx);xup2=Xe(Kother(1));
                                  Identy1=0; Indety2=1;   
                           end 
                      case 0 %% Xe(Kother(2)) ahead of Xe(Kother(1))
                           if Ye(Kother(2))==Ye(cx)  %% Xe(Kother(2)) same high with moddle point
                          
                              Lup1=@(x)(x-Xe(cx)).*((Ye(Kother(1))-Ye(cx))/(Xe(Kother(1))-Xe(cx)))+Ye(cx);
                              Lup2=Ye(cx);
                              xlow1=Xe(Kother(1));xup1=Xe(cx);
                              xlow2=Xe(cx);xup2=Xe(Kother(2));
                              Identy1=1;Identy2=0;
                           
                           elseif Ye(Kother(1))==Ye(cx) %% Xe(Kother(1)) same high with moddle point
                               
                                  Lup1=Ye(cx);
                                  Lup2=@(x)(x-Xe(cx)).*((Ye(Kother(2))-Ye(cx))/(Xe(Kother(2))-Xe(cx)))+Ye(cx);
                                  xlow1=Xe(Kother(1));xup1=Xe(cx);
                                  xlow2=Xe(cx);xup2=Xe(Kother(2));
                                  Identy1=1;Identy2=1;
                           end
                  end  
                  dPidnuBkl=integral2(integrant{xl},xlow1,xup1,Llow,Lup1)+integral2(integrant{xl},xlow2,xup2,Llow,Lup2,'Reltol',1e-9);         
           elseif Clarge==1 && Clow==0  %% middle point higher than line
                  Kother=find(Xe~=Xe(cx));identy=0;
                      if ((Ye(Kother(2))-Ye(Kother(1)))/((Xe(Kother(2))-Xe(Kother(1)))))==0
                          Lup=Ye(Kother(1));
                          identy=1;
                      else    
                          Lup=@(x)(x-Xe(Kother(1))).*((Ye(Kother(2))-Ye(Kother(1)))/((Xe(Kother(2))-Xe(Kother(1)))))+Ye(Kother(1)); 
                      end
                 
                      switch Xe(Kother(1))>Xe(Kother(2))
                            case 1 %% Xe(Kother(1)) ahead of Xe(Kother(2))
                                  if Ye(Kother(1))==Ye(cx)
                                
                                     Llow1=@(x)(x-Xe(cx)).*((Ye(Kother(2))-Ye(cx))/(Xe(Kother(2))-Xe(cx)))+Ye(cx);
                                     Llow2=Ye(cx);
                                     xlow1=Xe(Kother(2));xup1=Xe(cx);
                                     xlow2=Xe(cx);xup2=Xe(Kother(1));
                                     Identy1=0;Identy2=0;
                                      
                                  elseif Ye(Kother(2))==Ye(cx)   
                                      
                                         Llow1=Ye(cx);
                                         Llow2=@(x)(x-Xe(cx)).*((Ye(Kother(1))-Ye(cx))/(Xe(Kother(1))-Xe(cx)))+Ye(cx);
                                         xlow1=Xe(Kother(2));xup1=Xe(cx);
                                         xlow2=Xe(cx);xup2=Xe(Kother(1));
                                         Identy1=0;Identy2=1;
                                  end    
                                  
                            case 0  %% Xe(Kother(2)) ahead of Xe(kother(1))              
                                  if Ye(Kother(2))==Ye(cx)
                                
                                     Llow1=@(x)(x-Xe(cx)).*((Ye(Kother(1))-Ye(cx))/(Xe(Kother(1))-Xe(cx)))+Ye(cx);
                                     Llow2=Ye(cx);
                                     xlow1=Xe(Kother(1));xup1=Xe(cx);
                                     xlow2=Xe(cx);xup2=Xe(Kother(2));
                                     Identy1=1;Identy2=0;
                                      
                                  elseif Ye(Kother(1))==Ye(cx)   
                                      
                                         Llow1=Ye(cx);
                                         Llow2=@(x)(x-Xe(cx)).*((Ye(Kother(2))-Ye(cx))/(Xe(Kother(2))-Xe(cx)))+Ye(cx);
                                         xlow1=Xe(Kother(1));xup1=Xe(cx);
                                         xlow2=Xe(cx);xup2=Xe(Kother(2));
                                         Identy1=1;Identy2=1;
                                  end    
                      end 
                  dPidnuBkl=integral2(integrant{xl},xlow1,xup1,Llow1,Lup)+integral2(integrant{xl},xlow2,xup2,Llow2,Lup,'Reltol',1e-9);    %% middle point lower than line  
           end
              
     end

elseif (x1==x2) || (x1==x3) || (x2==x3)
    
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       Who=[x1==x2 x1==x3 x2==x3;];
       switch find(Who==1)
           case 1
                 lo1=x1>x3;lo2=x1<x3;
                 if lo1==1 && lo2==0
                    Lup=@(x)(x-x3).*((y2-y3)./(x2-x3))+y3;
                    Llow=@(x)(x-x3).*((y1-y3)./(x1-x3))+y3;
                    xlow=x3;xup=x1;
                 elseif lo1==0 && lo2==1
                        Lup=@(x)(x-x3).*((y1-y3)./(x1-x3))+y3;
                        Llow=@(x)(x-x3).*((y2-y3)./(x2-x3))+y3;
                        xlow=x1;xup=x3;
                 end    
           case 2
                 lo1=x1>x2;lo2=x1<x2;
                 if lo1==0 && lo2==1
                    Lup=@(x)(x-x2).*((y3-y2)./(x3-x2))+y2;
                    Llow=@(x)(x-x2).*((y1-y2)./(x1-x2))+y2;
                    xlow=x1;xup=x2;
                 elseif lo1==1 && lo2==0
                         Lup=@(x)(x-x2).*((y1-y2)./(x1-x2))+y2;
                         Llow=@(x)(x-x2).*((y3-y2)./(x3-x2))+y2;
                         xlow=x2;xup=x1;
                 end   
           case 3
                  lo1=x1>x3;lo2=x1<x3;
                 if lo1==1 && lo2==0
                    Lup=@(x)(x-x1).*((y2-y1)./(x2-x1))+y1;
                    Llow=@(x)(x-x1).*((y3-y1)./(x3-x1))+y1;
                    xlow=x3;xup=x1;
                 elseif lo1==0 && lo2==1
                         Lup=@(x)(x-x1).*((y3-y1)./(x3-x1))+y1;
                         Llow=@(x)(x-x1).*((y2-y1)./(x2-x1))+y1;
                         xlow=x1;xup=x3;
                 end           
       end
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       dPidnuBkl=integral2(integrant{xl},xlow,xup,Llow,Lup,'Reltol',1e-9);    
end

end


function Pitheta=PiThetaLambda2DSolitonWallQN(lambda,A,Alpha,d,Cv,CbIv,elements,Nelem,nodes,Np,gamma1,gamma)

CV=[Cv CbIv];
alphakn1t=zeros(length(CV),1);

 for k=1:length(CV)
     alphakn1t(k,1)=Alpha(CV(k));     
 end

alphakn1t=alphakn1t+lambda.*d;

for l=1:length(CV)
%     if l<=(length(Cv)+length(CbIv))
       Alpha(CV(l))=alphakn1t(l);   
%     else
%          Alpha(CV(l))=alphakn1t(l);
%          Alpha(Cbiu(l-(length(Cv)+length(CbIv))))=alphakn1t(l)+pi;
%     end    
end  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Pitheta1=((Alpha')*A)*Alpha;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Pitheta2Matrix=zeros(Nelem,1);
parfor l=1:Nelem

    El=[]; n1=[]; n2=[]; n3=[]; x1=[]; y1=[]; x2=[]; y2=[]; x3=[]; y3=[]; xi1=[]; xi2=[]; xi3=[]; eta=[]; xi=[]; omega=[]; D=[]; Si=[];
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    El=elements(:,l);

    n1=nodes(:,El(1));x1=n1(1);y1=n1(2);
    n2=nodes(:,El(2));x2=n2(1);y2=n2(2);
    n3=nodes(:,El(3));x3=n3(1);y3=n3(2);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%% gamma1 and gamma3 %%%%%%%%%%%%%
        
        gamma3=gamma3DSW2D(y1,y2,y3,gamma);
       
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    xi=[x2-x3 x3-x1 x1-x2];
    eta=[y2-y3 y3-y1 y1-y2];
    omega=[x2*y3-x3*y2 x3*y1-x1*y3 x1*y2-x2*y1];
    
    Di=sum(omega);
    Si=abs(Di)/2;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    ALPHA=zeros(size(El));
    ALPHA(1)=Alpha(El(1));
    ALPHA(2)=Alpha(El(2));
    ALPHA(3)=Alpha(El(3));
           
       %Pialpha2=Pialpha2+TriElePi2Integral2DSoliton(x1,x2,x3,y1,y2,y3,ALPHA,eta,xi,omega,Di); 
    Pitheta2Matrix(l,1)=TriElePi2Integral2DSolitonWall(x1,x2,x3,y1,y2,y3,ALPHA,eta,xi,omega,Di,gamma1,gamma3); 
              
end   

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Pitheta2=sum(Pitheta2Matrix);
Pitheta=Pitheta1+Pitheta2;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Il=TriElePi2Integral2DSolitonWall(x1,x2,x3,y1,y2,y3,ALPHA,eta,xi,omega,Di,gamma1,gamma3)

% integrant=@(x,y)(1/4)*cos(2.*((eta(1).*x-xi(1).*y+omega(1))./Di).*ALPHA(1)+2.*((eta(2).*x-xi(2).*y+omega(2))./Di).*ALPHA(2)+2.*((eta(3).*x-xi(3).*y+omega(3))./Di).*ALPHA(3));

integrant=@(x,y)((1/2)*1).*1.*(-gamma3.*sin(1.*((eta(1).*x-xi(1).*y+omega(1))./Di).*ALPHA(1)+1.*((eta(2).*x-xi(2).*y+omega(2))./Di).*ALPHA(2)+1.*((eta(3).*x-xi(3).*y+omega(3))./Di).*ALPHA(3))+...1
                  -(1/2).*cos(2.*((eta(1).*x-xi(1).*y+omega(1))./Di).*ALPHA(1)+2.*((eta(2).*x-xi(2).*y+omega(2))./Di).*ALPHA(2)+2.*((eta(3).*x-xi(3).*y+omega(3))./Di).*ALPHA(3)));
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Xe=[x1 x2 x3];Ye=[y1 y2 y3];

if (x1~=x2) && (x1~=x3) && (x2~=x3) %% every side have finite slope
   
    RPXe=sort(Xe);  %fix the middle point
    [rx,cx]=find(Xe==RPXe(2)); % find the position of middle point via x-cooordinate
    
    [ry,cy]=find(Ye==Ye(cx)); %check the position and number of y-coordinates equaling to middle point y-coordinate
    
    if length(cy)==1  %% no any other point has same y-coordinate with middle point
       
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        RPYe=sort(Ye); 
        [Ry,Cy]=find(RPYe==Ye(cx)); %% judge whether the middle point is largest or smallest
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if Cy==3 %% head upper, middle point maxinum
           Kother=find(Xe~=Xe(cx));
           if ((Ye(Kother(2))-Ye(Kother(1)))/((Xe(Kother(2))-Xe(Kother(1)))))==0
               Llow=Ye(Kother(1));
           else    
               Llow=@(x)(x-Xe(Kother(1))).*((Ye(Kother(2))-Ye(Kother(1)))/((Xe(Kother(2))-Xe(Kother(1)))))+Ye(Kother(1)); 
           end
           
           switch Xe(Kother(1))>Xe(Kother(2))
               case 1
                    Lup1=@(x)(x-Xe(cx)).*((Ye(Kother(2))-Ye(cy))/(Xe(Kother(2))-Xe(cx)))+Ye(cy);
                    Lup2=@(x)(x-Xe(cx)).*((Ye(Kother(1))-Ye(cy))/(Xe(Kother(1))-Xe(cx)))+Ye(cy);
                    xlow1=Xe(Kother(2));xup1=Xe(cx);
                    xlow2=Xe(cx);xup2=Xe(Kother(1));
               case 0                 
                    Lup1=@(x)(x-Xe(cx)).*((Ye(Kother(1))-Ye(cy))/(Xe(Kother(1))-Xe(cx)))+Ye(cy);
                    Lup2=@(x)(x-Xe(cx)).*((Ye(Kother(2))-Ye(cy))/(Xe(Kother(2))-Xe(cx)))+Ye(cy);
                    xlow1=Xe(Kother(1));xup1=Xe(cx);
                    xlow2=Xe(cx);xup2=Xe(Kother(2));
           end 
           Il=integral2(integrant,xlow1,xup1,Llow,Lup1)+integral2(integrant,xlow2,xup2,Llow,Lup2,'Reltol',1e-9);   
        end 
       
        if Cy==1 %% head down, middle point minium
           Kother=find(Xe~=Xe(cx));
           if ((Ye(Kother(2))-Ye(Kother(1)))/((Xe(Kother(2))-Xe(Kother(1)))))==0
               Lup=Ye(Kother(1));
           else    
               Lup=@(x)(x-Xe(Kother(1))).*((Ye(Kother(2))-Ye(Kother(1)))/((Xe(Kother(2))-Xe(Kother(1)))))+Ye(Kother(1)); 
           end
         
           switch Xe(Kother(1))>Xe(Kother(2))
               case 1
                    Llow1=@(x)(x-Xe(cx)).*((Ye(Kother(2))-Ye(cy))/(Xe(Kother(2))-Xe(cx)))+Ye(cy);
                    Llow2=@(x)(x-Xe(cx)).*((Ye(Kother(1))-Ye(cy))/(Xe(Kother(1))-Xe(cx)))+Ye(cy);
                    xlow1=Xe(Kother(2));xup1=Xe(cx);
                    xlow2=Xe(cx);xup2=Xe(Kother(1));
               case 0                 
                    Llow1=@(x)(x-Xe(cx)).*((Ye(Kother(1))-Ye(cy))/(Xe(Kother(1))-Xe(cx)))+Ye(cy);
                    Llow2=@(x)(x-Xe(cx)).*((Ye(Kother(2))-Ye(cy))/(Xe(Kother(2))-Xe(cx)))+Ye(cy);
                    xlow1=Xe(Kother(1));xup1=Xe(cx);
                    xlow2=Xe(cx);xup2=Xe(Kother(2));
           end 
           Il=integral2(integrant,xlow1,xup1,Llow1,Lup)+integral2(integrant,xlow2,xup2,Llow2,Lup,'Reltol',1e-9);   
        end
       
        if Cy==2 %% the middle point is the 2nd large
           [Rrx,Crx]=find(Xe~=Xe(cx));[Rry,Cry]=find(Ye~=Ye(cy)); 
           %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
           %%% Need an Equation stringnt line, judge the relative position of line and middle point   
           y=(Xe(cx)-Xe(Crx(2))).*((Ye(Cry(1))-Ye(Cry(2)))/(Xe(Crx(1))-Xe(Crx(2))))+Ye(Cry(2));
           Clarge=y>Ye(cy);Clow=y<Ye(cy);
           %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
           
           if Clarge==0 && Clow==1    %% middle point higher than line, The way of integral same Cy=3
              Kother=find(Xe~=Xe(cx));
              if ((Ye(Kother(2))-Ye(Kother(1)))/((Xe(Kother(2))-Xe(Kother(1)))))==0
                  Llow=Ye(Kother(1));
              else    
                  Llow=@(x)(x-Xe(Kother(1))).*((Ye(Kother(2))-Ye(Kother(1)))/((Xe(Kother(2))-Xe(Kother(1)))))+Ye(Kother(1)); 
              end
     
              switch Xe(Kother(1))>Xe(Kother(2))
                      case 1
                           Lup1=@(x)(x-Xe(cx)).*((Ye(Kother(2))-Ye(cy))/(Xe(Kother(2))-Xe(cx)))+Ye(cy);
                           Lup2=@(x)(x-Xe(cx)).*((Ye(Kother(1))-Ye(cy))/(Xe(Kother(1))-Xe(cx)))+Ye(cy);
                           xlow1=Xe(Kother(2));xup1=Xe(cx);
                           xlow2=Xe(cx);xup2=Xe(Kother(1));
                      case 0                 
                           Lup1=@(x)(x-Xe(cx)).*((Ye(Kother(1))-Ye(cy))/(Xe(Kother(1))-Xe(cx)))+Ye(cy);
                           Lup2=@(x)(x-Xe(cx)).*((Ye(Kother(2))-Ye(cy))/(Xe(Kother(2))-Xe(cx)))+Ye(cy);
                           xlow1=Xe(Kother(1));xup1=Xe(cx);
                           xlow2=Xe(cx);xup2=Xe(Kother(2));
               end 
               Il=integral2(integrant,xlow1,xup1,Llow,Lup1)+integral2(integrant,xlow2,xup2,Llow,Lup2,'Reltol',1e-9);  
                
            elseif Clarge==1 && Clow==0  %% middle point lower than line,  The way of integral same Cy=1
                   Kother=find(Xe~=Xe(cx));
                   if ((Ye(Kother(2))-Ye(Kother(1)))/((Xe(Kother(2))-Xe(Kother(1)))))==0
                      Lup=Ye(Kother(1));
                   else    
                      Lup=@(x)(x-Xe(Kother(1))).*((Ye(Kother(2))-Ye(Kother(1)))/((Xe(Kother(2))-Xe(Kother(1)))))+Ye(Kother(1)); 
                   end
               
                   switch Xe(Kother(1))>Xe(Kother(2))
                         case 1
                              Llow1=@(x)(x-Xe(cx)).*((Ye(Kother(2))-Ye(cy))/(Xe(Kother(2))-Xe(cx)))+Ye(cy);
                              Llow2=@(x)(x-Xe(cx)).*((Ye(Kother(1))-Ye(cy))/(Xe(Kother(1))-Xe(cx)))+Ye(cy);
                              xlow1=Xe(Kother(2));xup1=Xe(cx);
                              xlow2=Xe(cx);xup2=Xe(Kother(1));
                         case 0                 
                              Llow1=@(x)(x-Xe(cx)).*((Ye(Kother(1))-Ye(cy))/(Xe(Kother(1))-Xe(cx)))+Ye(cy);
                              Llow2=@(x)(x-Xe(cx)).*((Ye(Kother(2))-Ye(cy))/(Xe(Kother(2))-Xe(cx)))+Ye(cy);
                              xlow1=Xe(Kother(1));xup1=Xe(cx);
                              xlow2=Xe(cx);xup2=Xe(Kother(2));
                   end 
                  Il=integral2(integrant,xlow1,xup1,Llow1,Lup)+integral2(integrant,xlow2,xup2,Llow2,Lup,'Reltol',1e-9);    
            end
         end   
    elseif length(cy)==2
           Cother=find(cy~=cx);
           Xo=Xe(cy(Cother));Yo=Ye(cy(Cother)); %% find out the other same high with middle point
           %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
           %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
           Crx=find(Xe~=Xe(cx) & Xe~=Xo);Cry=find(Ye~=Ye(cx) & Ye~=Yo); 
           %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
           % %% Need an Equation stringnt line 
           y=(Xe(cx)-Xe(Crx)).*((Yo-Ye(Cry))/(Xo-Xe(Crx)))+Ye(Cry);
           Clarge=y>Ye(cx);Clow=y<Ye(cx);
           %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
           %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
           if Clarge==0 && Clow==1    %% middle point higher than line, The way of integral same Cy=3
              Kother=find(Xe~=Xe(cx));identy=0;
              if ((Ye(Kother(2))-Ye(Kother(1)))/((Xe(Kother(2))-Xe(Kother(1)))))==0
                  Llow=Ye(Kother(1));
                  identy=1;
              else    
                  Llow=@(x)(x-Xe(Kother(1))).*((Ye(Kother(2))-Ye(Kother(1)))/((Xe(Kother(2))-Xe(Kother(1)))))+Ye(Kother(1)); 
              end
              
               switch Xe(Kother(1))>Xe(Kother(2))
                      case 1
                           if Ye(Kother(1))==Ye(cx) %% Xe(Kother(1)) same high with moddle point
                              Lup1=@(x)(x-Xe(cx)).*((Ye(Kother(2))-Ye(cx))/(Xe(Kother(2))-Xe(cx)))+Ye(cx);
                              Lup2=Ye(cx);
                              xlow1=Xe(Kother(2));xup1=Xe(cx);
                              xlow2=Xe(cx);xup2=Xe(Kother(1));
                              Identy1=0;Identy2=0;
                              
                           elseif Ye(Kother(2))==Ye(cx) %% Xe(Kother(2)) same high with moddle point
                                  Lup1=Ye(cx);
                                  Lup2=@(x)(x-Xe(cx)).*((Ye(Kother(1))-Ye(cx))/(Xe(Kother(1))-Xe(cx)))+Ye(cx);
                                  xlow1=Xe(Kother(2));xup1=Xe(cx);
                                  xlow2=Xe(cx);xup2=Xe(Kother(1));
                                  Identy1=0; Indety2=1;   
                           end 
                      case 0                 
                           if Ye(Kother(2))==Ye(cx)  %% Xe(Kother(2)) same high with moddle point
                          
                              Lup1=@(x)(x-Xe(cx)).*((Ye(Kother(1))-Ye(cx))/(Xe(Kother(1))-Xe(cx)))+Ye(cx);
                              Lup2=Ye(cx);
                              xlow1=Xe(Kother(1));xup1=Xe(cx);
                              xlow2=Xe(cx);xup2=Xe(Kother(2));
                              Identy1=1;Identy2=0;
                           
                           elseif Ye(Kother(1))==Ye(cx) %% Xe(Kother(1)) same high with moddle point
                               
                                  Lup1=Ye(cx);
                                  Lup2=@(x)(x-Xe(cx)).*((Ye(Kother(2))-Ye(cx))/(Xe(Kother(2))-Xe(cx)))+Ye(cx);
                                  xlow1=Xe(Kother(1));xup1=Xe(cx);
                                  xlow2=Xe(cx);xup2=Xe(Kother(2));
                                  Identy1=1;Identy2=1;
                           end
               end 
               Il=integral2(integrant,xlow1,xup1,Llow,Lup1)+integral2(integrant,xlow2,xup2,Llow,Lup2,'Reltol',1e-9);         
           elseif Clarge==1 && Clow==0  %% middle point higher than line
                  Kother=find(Xe~=Xe(cx));identy=0;
                  if ((Ye(Kother(2))-Ye(Kother(1)))/((Xe(Kother(2))-Xe(Kother(1)))))==0
                     Lup=Ye(Kother(1));
                     identy=1;
                  else    
                     Lup=@(x)(x-Xe(Kother(1))).*((Ye(Kother(2))-Ye(Kother(1)))/((Xe(Kother(2))-Xe(Kother(1)))))+Ye(Kother(1)); 
                  end
                 
                  switch Xe(Kother(1))>Xe(Kother(2))
                         case 1
                              if Ye(Kother(1))==Ye(cx)
                                
                                     Llow1=@(x)(x-Xe(cx)).*((Ye(Kother(2))-Ye(cx))/(Xe(Kother(2))-Xe(cx)))+Ye(cx);
                                     Llow2=Ye(cx);
                                     xlow1=Xe(Kother(2));xup1=Xe(cx);
                                     xlow2=Xe(cx);xup2=Xe(Kother(1));
                                     Identy1=0;Identy2=0;
                                      
                              elseif Ye(Kother(2))==Ye(cx)   
                                      
                                         Llow1=Ye(cx);
                                         Llow2=@(x)(x-Xe(cx)).*((Ye(Kother(1))-Ye(cx))/(Xe(Kother(1))-Xe(cx)))+Ye(cx);
                                         xlow1=Xe(Kother(2));xup1=Xe(cx);
                                         xlow2=Xe(cx);xup2=Xe(Kother(1));
                                         Identy1=0;Identy2=1;
                              end
                         case 0                 
                              if Ye(Kother(2))==Ye(cx)
                                
                                     Llow1=@(x)(x-Xe(cx)).*((Ye(Kother(1))-Ye(cx))/(Xe(Kother(1))-Xe(cx)))+Ye(cx);
                                     Llow2=Ye(cx);
                                     xlow1=Xe(Kother(1));xup1=Xe(cx);
                                     xlow2=Xe(cx);xup2=Xe(Kother(2));
                                     Identy1=1;Identy2=0;
                                      
                              elseif Ye(Kother(1))==Ye(cx)   
                                      
                                         Llow1=Ye(cx);
                                         Llow2=@(x)(x-Xe(cx)).*((Ye(Kother(2))-Ye(cx))/(Xe(Kother(2))-Xe(cx)))+Ye(cx);
                                         xlow1=Xe(Kother(1));xup1=Xe(cx);
                                         xlow2=Xe(cx);xup2=Xe(Kother(2));
                                         Identy1=1;Identy2=1;
                              end    
                  end
                  Il=integral2(integrant,xlow1,xup1,Llow1,Lup)+integral2(integrant,xlow2,xup2,Llow2,Lup,'Reltol',1e-9);    %% middle point lower than line  
           end
              
     end

elseif (x1==x2) || (x1==x3) || (x2==x3)
    
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       Who=[x1==x2 x1==x3 x2==x3;];
       switch find(Who==1)
           case 1
                 lo1=x1>x3;lo2=x1<x3;
                 if lo1==1 && lo2==0
                    Lup=@(x)(x-x3).*((y2-y3)./(x2-x3))+y3;
                    Llow=@(x)(x-x3).*((y1-y3)./(x1-x3))+y3;
                    xlow=x3;xup=x1;
                 elseif lo1==0 && lo2==1
                        Lup=@(x)(x-x3).*((y1-y3)./(x1-x3))+y3;
                        Llow=@(x)(x-x3).*((y2-y3)./(x2-x3))+y3;
                        xlow=x1;xup=x3;
                 end    
           case 2
                 lo1=x1>x2;lo2=x1<x2;
                 if lo1==0 && lo2==1
                    Lup=@(x)(x-x2).*((y3-y2)./(x3-x2))+y2;
                    Llow=@(x)(x-x2).*((y1-y2)./(x1-x2))+y2;
                    xlow=x1;xup=x2;
                 elseif lo1==1 && lo2==0
                         Lup=@(x)(x-x2).*((y1-y2)./(x1-x2))+y2;
                         Llow=@(x)(x-x2).*((y3-y2)./(x3-x2))+y2;
                         xlow=x2;xup=x1;
                 end   
           case 3
                  lo1=x1>x3;lo2=x1<x3;
                 if lo1==1 && lo2==0
                    Lup=@(x)(x-x1).*((y2-y1)./(x2-x1))+y1;
                    Llow=@(x)(x-x1).*((y3-y1)./(x3-x1))+y1;
                    xlow=x3;xup=x1;
                 elseif lo1==0 && lo2==1
                         Lup=@(x)(x-x1).*((y3-y1)./(x3-x1))+y1;
                         Llow=@(x)(x-x1).*((y2-y1)./(x2-x1))+y1;
                         xlow=x1;xup=x3;
                 end           
       end
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       Il=integral2(integrant,xlow,xup,Llow,Lup,'Reltol',1e-9);    
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Hk=InHassian2DSolitonQN(Is,rhokn1,skn1,ykn1,Hkn1)

left=Is-rhokn1*(skn1*(ykn1'));
right=Is-rhokn1*(ykn1*(skn1'));

InHassian1=(left*Hkn1)*right;
Inhassian2=rhokn1*(skn1*(skn1'));

Hk=InHassian1+Inhassian2;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% gamma3 for integral %%%%%%%%%%%%%%%%%%%
function gamma3=gamma3DSW2D(y1,y2,y3,gamma)
  

if y1<0 || y2<0 || y3<0
    omega=-gamma; % -y part is -|DeltaV2| domain
else
    omega=gamma;
end    
    gamma3=omega/(1+gamma);
end

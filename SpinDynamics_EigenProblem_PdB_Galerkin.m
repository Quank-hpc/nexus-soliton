function SpinDynamics_EigenProblem_PdB_Galerkin 
try

   clear all

   load('SoliTonWall_01_40_PdBLattice_gamma008_Box20_0.mat');

   Deltap=1;Deltav1=0.08*Deltap;%Deltav2=0.1*Deltap;
   gamma1=Deltav1/Deltap;%gamma2=Deltav2/Deltap;
   
%   [nodes,elements,Np,Nelem]=TriPartition2D;
   [nodes,elements,Np,Nelem]=TriPartition2DSW;
   Np
   Nelem

%   [Rbu,Cbu,Rbd,Cbd,Rbiu,Cbiu,Rbid,Cbid,RbIv,CbIv,Rv,Cv]=ClassifyNodes(nodes);
%[Rbu,Cbu,Rbd,Cbd,Rbm,Cbm,RbiuA,CbiuA,RbidA,CbidA,RbiuB,CbiuB,RbidB,CbidB,RbIv,CbIv,Rv,Cv]=ClassifyNodes(nodes);
% [RbW,CbW,Rbu,Cbu,Rbd,Cbd,Rbm,Cbm,RbiuA,CbiuA,RbidA,CbidA,RbiuB,CbiuB,RbidB,CbidB,RbIv,CbIv,Rv,Cv]=ClassifyNodes(nodes);
[RbW,CbW,Rbu,Cbu,Rbd,Cbd,Rbiu,Cbiu,Rbid,Cbid,RbIv,CbIv,Rv,Cv]=ClassifyNodes(nodes);


   alpha=griddata(X,Y,Thetakn1,nodes(1,:),nodes(2,:));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   A=zeros(Np,Np); %% the exact A Matrix corressponds to original nodes matrix+extra-intra-boundary point matrix Cbi
   B=zeros(Np,Np);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% A B Matrix %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   [A,B]=ABMatix2D_SpinWave_U(Nelem,elements,nodes,A,B,alpha,gamma1);

   % for ii=1:Np
   %     for jj=[Cbu Cbd Cbiu Cbid] % may be need for calculation of NonGaugeTrans case 
   %         A(ii,jj)=0;B(ii,jj)=0;
   %     end
   % end

   for ii=[Cbu Cbd CbW]
    
       if ii==1
          V=[1 zeros(1,Np-1)]; 
       elseif ii==Np
             V=[zeros(1,Np-1) 1]; 
       else
             V=[zeros(1,ii-1) 1 zeros(1,Np-ii)];
       end    
    
       A(ii,:)=V;B(ii,:)=zeros(1,Np);
    
   end
   

   disp('start to diagnolization, waiting...');

   [EVector,EValue] = eig(A,B);
%  [EVector,EValue] = eig(A(length(Cb)+1:Np,length(Cb)+1:Np),B(length(Cb)+1:Np,length(Cb)+1:Np));
% [EVector,EValue] = eigs(A(length(Cb)+1:Np,length(Cb)+1:Np),B(length(Cb)+1:Np,length(Cb)+1:Np),24,'smallestreal');

   [EigenValue,indx]=sort(diag(real(EValue)));
    EigenVector=EVector(:,indx);

   disp('diagonalization terminated, saving...');
   x=-20:0.01:20;y=x;
   %[X,Y]=meshgrid(x,y);
   EigenVector10=EigenVector(:,1:10);
   save('spinwave_PdB_Lattice_Box20_0_gamma008.mat','A','B','EigenValue','EigenVector10','nodes','elements','x','y','X','Y');
   disp('saving terminates'); 
catch error
      disp(getReport(error))
      exit(1)
end

end
   


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [nodes,elements,Np,Nelem]=TriPartition2DSW

model=createpde(1);
geometryFromEdges(model,@SolveDomain);
M=generateMesh(model,'JiggleIter',10,'Hgrad',1.8,'Hmax',0.50);
nodes=M.Nodes;
elements=M.Elements;
Np=length(nodes);
Nelem=length(elements);


end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [x,y]=SolveDomain(bs,s)
xiDt=1;nbs=15;S=20.0;D1=9.9;D2=0.1;H1=0.01;H2=0.1;
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
xiDt=1;nbs=15;S=20.0;D1=9.9;D2=0.1;H1=0.01;H2=0.1;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[Rbu,Cbu]=find(nodes(2,:)>0 & abs(nodes(2,:))>=S/2);

[Rbd,Cbd]=find(nodes(2,:)<0 & abs(nodes(2,:))>=S/2);

[Rbiu,Cbiu]=find(nodes(2,:)>0 & abs(nodes(2,:))==H1 & abs(nodes(1,:))>=D1+D2);

[Rbid,Cbid]=find(nodes(2,:)<0 & abs(nodes(2,:))==H1 & abs(nodes(1,:))>=D1+D2);

[RbW,CbW]=find(abs(nodes(2,:))==0 & abs(nodes(1,:))<D1);

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
function [A,B]=ABMatix2D_SpinWave_U(Nelem,elements,nodes,A,B,alpha,gamma1)

%Coef2=6*gamma2*gamma2+(gamma1*gamma1+1);
%Coef1=3*gamma1*gamma1+(2*gamma2*gamma2+1);


for k=1:Nelem
    clc
    Nelem
    length(nodes(1,:))
    k
    clear Ek n1 n2 n3 x1 y1 x2 y2 x3 y3 xi1 xi2 xi3 eta xi omega D Si a2xlij a2ylij L1 L2 L3 a2lij A2lij
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Ek=elements(:,k);

    n1=nodes(:,Ek(1));x1=n1(1);y1=n1(2);
    n2=nodes(:,Ek(2));x2=n2(1);y2=n2(2);
    n3=nodes(:,Ek(3));x3=n3(1);y3=n3(2);
    
    xi=[x2-x3 x3-x1 x1-x2];
    eta=[y2-y3 y3-y1 y1-y2];
    omega=[x2*y3-x3*y2 x3*y1-x1*y3 x1*y2-x2*y1];
    
    Di=sum(omega);
    Si=abs(Di)/2;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%% gamma2 for different domain %%%%
 
    gamma2=gamma2SD(y1,y2,y3,gamma1);

    Coef1=3*gamma1*gamma1+(2*gamma2*gamma2+1);
    Coef2=6*gamma2*gamma2+(gamma1*gamma1+1);    

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    akji=zeros(3,3);
    for o=1:3
        akji(o,1)=Coef1*(1/(1*(Di.^2)))*eta(o)*eta(1)*Si+Coef2*(1/(1*(Di.^2)))*xi(o)*xi(1)*Si+...
	    Integral2D_kvtu(o,1,Ek,x1,y1,x2,y2,x3,y3,xi,eta,omega,Di,alpha,gamma1)+...
                  Integral2D_kvUiu(o,1,Ek,x1,y1,x2,y2,x3,y3,xi,eta,omega,Di,alpha,gamma1,gamma2);

        akji(o,2)=Coef1*(1/(1*(Di.^2)))*eta(o)*eta(2)*Si+Coef2*(1/(1*(Di.^2)))*xi(o)*xi(2)*Si+...
	    Integral2D_kvtu(o,2,Ek,x1,y1,x2,y2,x3,y3,xi,eta,omega,Di,alpha,gamma1)+...
                  Integral2D_kvUiu(o,2,Ek,x1,y1,x2,y2,x3,y3,xi,eta,omega,Di,alpha,gamma1,gamma2);

        akji(o,3)=Coef1*(1/(1*(Di.^2)))*eta(o)*eta(3)*Si+Coef2*(1/(1*(Di.^2)))*xi(o)*xi(3)*Si+...
	    Integral2D_kvtu(o,3,Ek,x1,y1,x2,y2,x3,y3,xi,eta,omega,Di,alpha,gamma1)+...
                  Integral2D_kvUiu(o,3,Ek,x1,y1,x2,y2,x3,y3,xi,eta,omega,Di,alpha,gamma1,gamma2);

	
        % akji(o,2)=(1/(1*(Di.^2)))*eta(o)*eta(2)*Si+(1/(1*(Di.^2)))*xi(o)*xi(2)*Si+...
	%   Integral2D_kvtu(o,2,Ek,x1,y1,x2,y2,x3,y3,xi,eta,omega,Di,alpha)+...
        %           Integral2D_kvUiu(o,2,Ek,x1,y1,x2,y2,x3,y3,xi,eta,omega,Di,alpha);
              
        % akji(o,3)=(1/(1*(Di.^2)))*eta(o)*eta(3)*Si+(1/(1*(Di.^2)))*xi(o)*xi(3)*Si+...
	%   Integral2D_kvtu(o,3,Ek,x1,y1,x2,y2,x3,y3,xi,eta,omega,Di,alpha)+...
        %           Integral2D_kvUiu(o,3,Ek,x1,y1,x2,y2,x3,y3,xi,eta,omega,Di,alpha);
            
    end  
    
    bkji=zeros(3,3);
    for o=1:3
        
        for f=1:3
            bkji(o,f)=Integral2D_kvu(o,f,Ek,x1,y1,x2,y2,x3,y3,xi,eta,omega,Di);
        end            
    end
    
   %%%%%%%%%%%%matrix symmetrized%%%%%%%%%
   Ak=(1/2)*(akji+akji');
   Bk=(1/2)*(bkji+bkji');
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
      for u=1:3 %% The lth element locaizes inside Domain, us the original index 
  
          A(Ek(u),Ek(1))=Ak(u,1)+A(Ek(u),Ek(1));
          A(Ek(u),Ek(2))=Ak(u,2)+A(Ek(u),Ek(2));
          A(Ek(u),Ek(3))=Ak(u,3)+A(Ek(u),Ek(3));
          
          B(Ek(u),Ek(1))=Bk(u,1)+B(Ek(u),Ek(1));
          B(Ek(u),Ek(2))=Bk(u,2)+B(Ek(u),Ek(2));
          B(Ek(u),Ek(3))=Bk(u,3)+B(Ek(u),Ek(3));
      end

end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Ivu=Integral2D_kvu(o,f,Ek,x1,y1,x2,y2,x3,y3,xi,eta,omega,Di)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       integrant=@(x,y) ((eta(o).*x-xi(o).*y+omega(o))./Di).*((eta(f).*x-xi(f).*y+omega(f))./Di);
       Ivu=TriEleIntegral2D(integrant,x1,x2,x3,y1,y2,y3);
     
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Ivtu=Integral2D_kvtu(o,f,Ek,x1,y1,x2,y2,x3,y3,xi,eta,omega,Di,alpha,gamma1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
       integrant1=@(x,y) ((eta(o).*x-xi(o).*y+omega(o))./Di).*(1i.*(2)).*Ptheta1Pphi1(x,y,f,alpha,Ek,xi,eta,omega,Di,gamma1);
       Ivtu1=TriEleIntegral2D(integrant1,x1,x2,x3,y1,y2,y3);

       integrant2=@(x,y) ((eta(o).*x-xi(o).*y+omega(o))./Di).*(1i.*(2).*(1+gamma1*gamma1)).*Ptheta2Pphi2(f,alpha,Ek,xi,eta,omega,Di);
       Ivtu2=TriEleIntegral2D(integrant2,x1,x2,x3,y1,y2,y3);

       Ivtu=Ivtu1+Ivtu2;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function IvUu=Integral2D_kvUiu(o,f,Ek,x1,y1,x2,y2,x3,y3,xi,eta,omega,Di,alpha,gamma1,gamma2)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       integrant=@(x,y) ((eta(o).*x-xi(o).*y+omega(o))./Di).*Uiext(x,y,alpha,Ek,xi,eta,omega,Di,gamma1,gamma2).*((eta(f).*x-xi(f).*y+omega(f))./Di);

       IvUu=TriEleIntegral2D(integrant,x1,x2,x3,y1,y2,y3);
     
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function I=TriEleIntegral2D(integrant,x1,x2,x3,y1,y2,y3)

xidt=1;

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
           I=integral2(integrant,xlow1,xup1,Llow,Lup1,'Reltol',1e-9)+integral2(integrant,xlow2,xup2,Llow,Lup2,'Reltol',1e-9);   
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
           I=integral2(integrant,xlow1,xup1,Llow1,Lup,'Reltol',1e-9)+integral2(integrant,xlow2,xup2,Llow2,Lup,'Reltol',1e-9);   
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
               I=integral2(integrant,xlow1,xup1,Llow,Lup1,'Reltol',1e-9)+integral2(integrant,xlow2,xup2,Llow,Lup2,'Reltol',1e-9);  
                
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
                  I=integral2(integrant,xlow1,xup1,Llow1,Lup,'Reltol',1e-9)+integral2(integrant,xlow2,xup2,Llow2,Lup,'Reltol',1e-9);    
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
                  I=integral2(integrant,xlow1,xup1,Llow,Lup1,'Reltol',1e-9)+integral2(integrant,xlow2,xup2,Llow,Lup2,'Reltol',1e-9);         
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
                  I=integral2(integrant,xlow1,xup1,Llow1,Lup,'Reltol',1e-9)+integral2(integrant,xlow2,xup2,Llow2,Lup,'Reltol',1e-9);    %% middle point lower than line  
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
       I=integral2(integrant,xlow,xup,Llow,Lup,'Reltol',1e-9);    
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function U=Uiext(x,y,alpha,Ek,xi,eta,omega,Di,gamma1,gamma2)

alphak=((eta(1).*x-xi(1).*y+omega(1))./Di).*alpha(Ek(1))+...
       ((eta(2).*x-xi(2).*y+omega(2))./Di).*alpha(Ek(2))+...
       ((eta(3).*x-xi(3).*y+omega(3))./Di).*alpha(Ek(3));
% U1=-sin(alphak).*sin(alphak);

%U1=(1i).*(1/2).*(gamma2*(gamma1-1).*cos(alphak)+(1-gamma1)*(1-gamma1).*sin(2.*alphak));
U1=(1i).*(1/2).*(-gamma2*(1+gamma1).*cos(alphak)+(1+gamma1)*(1+gamma1).*sin(2.*alphak)); %% from Jaakko's parametrization

%DalphakDj=[(alpha(Ek(1)).*(eta(1)./Di))+(alpha(Ek(2)).*(eta(2)./Di))+(alpha(Ek(3)).*(eta(3)./Di));...
%           (alpha(Ek(1)).*(-xi(1)./Di))+(alpha(Ek(2)).*(-xi(2)./Di))+(alpha(Ek(3)).*(-xi(3)./Di));]; 
%DalphakDj2=sum((DalphakDj')*(DalphakDj));
%U2=-DalphakDj2;

% U2=(1i).*(1/2).*(sin(2.*(alphak)));
% U2=(-1/2).*((1-gamma1).*((gamma1-1).*cos(2.*alphak)-5*gamma2.*sin(alphak))+(1+gamma1*gamma1+4*gamma1*gamma1));
U2=(-1/2).*((1+gamma1).*(-(1+gamma1).*cos(2.*alphak)-5*gamma2.*sin(alphak))+(1+gamma1*gamma1+4*gamma1*gamma1)); %% from Jaakko's parametrization

U=U1+U2;

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function PtPp=PthetaPphi(f,alpha,Ek,xi,eta,omega,Di)

% DalphakDj=[(alpha(Ek(1)).*(eta(1)./Di))+(alpha(Ek(2)).*(eta(2)./Di))+(alpha(Ek(3)).*(eta(3)./Di));...
%            (alpha(Ek(1)).*(-xi(1)./Di))+(alpha(Ek(2)).*(-xi(2)./Di))+(alpha(Ek(3)).*(-xi(3)./Di));]; 
% DphiDj=(1/Di).*[eta(f); -xi(f);];      

% PtPp=(DalphakDj')*(DphiDj);
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function f1_vtu=Ptheta1Pphi1(x,y,f,alpha,Ek,xi,eta,omega,Di,gamma1)
alphak=((eta(1).*x-xi(1).*y+omega(1))./Di).*alpha(Ek(1))+...
       ((eta(2).*x-xi(2).*y+omega(2))./Di).*alpha(Ek(2))+...
       ((eta(3).*x-xi(3).*y+omega(3))./Di).*alpha(Ek(3));

DalphakD1=(alpha(Ek(1)).*(eta(1)./Di))+(alpha(Ek(2)).*(eta(2)./Di))+(alpha(Ek(3)).*(eta(3)./Di));

DphiD1=(1/Di)*eta(f);

f1_vtu=(1+3*(gamma1*gamma1).*cos(2.*alphak)).*(DalphakD1).*(DphiD1);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function f2_vtu=Ptheta2Pphi2(f,alpha,Ek,xi,eta,omega,Di)

DalphakD2=(alpha(Ek(1)).*(-xi(1)./Di))+(alpha(Ek(2)).*(-xi(2)./Di))+(alpha(Ek(3)).*(-xi(3)./Di));

DphiD2=(1/Di).*(-xi(f));

f2_vtu=DalphakD2.*DphiD2;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% gamma2 for different domain%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function gamma2=gamma2SD(y1,y2,y3,gamma1)


 if y1<0 || y2<0 || y3<0
    gamma2=-gamma1;
 else
    gamma2=gamma1;
 end    
   
end

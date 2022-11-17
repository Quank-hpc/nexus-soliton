function calculate_ratio_intensity
try

  EigenModes_Data=dir('SpinDynamicsResponse_PdB_Box*_*_gamma0*.mat');

  NoOfEigenModesMat=length(EigenModes_Data);

% DimensionlessCoefficient_Matrix=[];
  FrequencyShift_Matrix=[];
  DesityOfHQVsBox_Matrix=[];
  RatioIntensityPerUnitArea_Matrix=[];
  save('RatioIntensity_Matrix_box05_box25_5_gamma000_gamma020.mat','DesityOfHQVsBox_Matrix','RatioIntensityPerUnitArea_Matrix','FrequencyShift_Matrix');

  for ii=1:NoOfEigenModesMat
      ii 
      disp(EigenModes_Data(ii).name);

     load(EigenModes_Data(ii).name,'A','B','EigenValue','indx','EigenVector10','nodes','elements',...
       'Np','Nelem','gamma','gamma1','X','Y','Thetakn1','RbW','CbW','Rbu','Cbu','Rbd','Cbd','Rbiu','Cbiu','Rbid','Cbid','RbIv','CbIv','Rv','Cv');

     D=max(nodes(1,:)); % step length is 0.5
     DesityOfHQVsBox=1/(D^2); 
     
     q=gamma; % step length is 0.01;

     row_number=round(((D-5)/0.5)+1)

     colum_number=round((q/0.01)+1)

     load('RatioIntensity_Matrix_box05_box25_5_gamma000_gamma020.mat','DesityOfHQVsBox_Matrix','RatioIntensityPerUnitArea_Matrix','FrequencyShift_Matrix');
     
     if row_number==6

        continue; % for box 7.5

     else
         
         RatioIntensity=calculate_ration_intensity_of_eigen_mode(elements,nodes,EigenVector10,DesityOfHQVsBox);         
         
       % DimensionlessCoefficient_Matrix(row_number,colum_number)=DimensionlessCoefficient;

         FrequencyShift_Matrix(row_number,colum_number)=EigenValue(1,1);

         DesityOfHQVsBox_Matrix(row_number,colum_number)=DesityOfHQVsBox;

         RatioIntensityPerUnitArea_Matrix(row_number,colum_number)=RatioIntensity;

     end

     RatioIntensity
     save('RatioIntensity_Matrix_box05_box25_5_gamma000_gamma020.mat','DesityOfHQVsBox_Matrix','RatioIntensityPerUnitArea_Matrix','FrequencyShift_Matrix');
     disp(' successivly save ');
  end

catch error

      disp(getReport(error))
      exit(1)

end


end


function RatioIntensity=calculate_ration_intensity_of_eigen_mode(elements,nodes,EigenVector10,DesityOfHQVsBox)
%% for calculating Ratio instensity of PdB

%I_IS=0;I_SI=0;
I=[];
IN_n=[];
ID_n=[];
 
 Nelem=length(elements);
 
 N_EV=1;
 
for n=1:N_EV
  EigenVector=EigenVector10(:,n);

  I_IS=0;
  I_SI=0;
%  parpool(4);
  parfor k=1:Nelem
         
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      Ek=elements(:,k);

      n1=nodes(:,Ek(1));x1=n1(1);y1=n1(2);
      n2=nodes(:,Ek(2));x2=n2(1);y2=n2(2);
      n3=nodes(:,Ek(3));x3=n3(1);y3=n3(2);
    
      xi=[x2-x3 x3-x1 x1-x2];
      eta=[y2-y3 y3-y1 y1-y2];
      omega=[x2*y3-x3*y2 x3*y1-x1*y3 x1*y2-x2*y1];
    
      Di=sum(omega);
     
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
      I_IS_ek=Integral2TriElm_IS(x1,y1,x2,y2,x3,y3,eta,xi,omega,Di,EigenVector,Ek);
      I_SI_ek=Integral2TriElm_SI(x1,y1,x2,y2,x3,y3,eta,xi,omega,Di,EigenVector,Ek);

      I_IS(k)=I_IS_ek;
      I_SI(k)=I_SI_ek;
  end

  IN=abs(sum(I_IS)).*abs(sum(I_IS));ID=sum(I_SI);
  INID=IN/ID;
  RatioIntensity=DesityOfHQVsBox*INID;
    
end

end

function I_IS=Integral2TriElm_IS(x1,y1,x2,y2,x3,y3,eta,xi,omega,Di,EigenVector,Ek)

integrant=@ (x,y) EigenVector(Ek(1)).*((eta(1).*x-xi(1).*y+omega(1))./Di)+...
                           +EigenVector(Ek(2)).*((eta(2).*x-xi(2).*y+omega(2))./Di)+...
                           +EigenVector(Ek(3)).*((eta(3).*x-xi(3).*y+omega(3))./Di);

I_IS=TriEleIntegral2D(integrant,x1,x2,x3,y1,y2,y3);           
                                
end

function I_SI=Integral2TriElm_SI(x1,y1,x2,y2,x3,y3,eta,xi,omega,Di,EigenVector,Ek)

integrant=@ (x,y) (abs(EigenVector(Ek(1)).*((eta(1).*x-xi(1).*y+omega(1))./Di)+...
                           +EigenVector(Ek(2)).*((eta(2).*x-xi(2).*y+omega(2))./Di)+...
                           +EigenVector(Ek(3)).*((eta(3).*x-xi(3).*y+omega(3))./Di))).^2;
                       
I_SI=TriEleIntegral2D(integrant,x1,x2,x3,y1,y2,y3);     

end

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


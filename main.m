clear all;
NdofElement=6; %% gdl dell'elemento
global geom; %%alloco dello scope globale
global b_E;
global epsilon; epsilon=1e-5;

global RefiningOptions; %%rendiamo globale il tutto
    ntest=4; %%numero di divisioin del triangolo
    Area=zeros(1,ntest); %%Vettore che contiene l'area max dei triangoli.
    h=zeros(1,ntest); %%Vettore che contiene la max dimensione dei triangoli.
    PeTest=zeros(1,ntest);
    e_L2=zeros(1,ntest); %%Vettore che contiene l'errore in L2
    e_INF=zeros(1,ntest);
    e_H1=zeros(1,ntest);
    TheMesh=Mesh_info; %% Struttura
    ConditioningNumbers=zeros(1,ntest);
beta=[1;3]; %%campo della convezionegtt

for t=1:ntest

    RefiningOptions.AreaValue=0.1/3^t;
    Sample_Square_Dirichlet; %%AVVIO COSTRUZIONE DELLA MESH
    


    

    %P2;
    P2_Corretto; %%COSTRUZIONE DEI P2
    
    TheMesh.geom=geom; %%Acquisisco la struttura dati nella classe mesh
    
    [A,b_E,Ad,u_d]=Matrice_Stabilizzata(BC,@fun,beta,TheMesh);
    %[A,b_E,Ad,u_d]=Stiffness_Matrix(BC,@fun,beta,TheMesh);
    %ConditioningNumbers(t)=cond(A);
    %Assemblo Neumann%
    b_E=NeumannAssembly;

    %%Risolvo il sistema
    x=A\(b_E-Ad*u_d);

    %Costruzione della soluzione%
    u=built_sol(x,u_d);
    
    trisurf(geom.elements.triangles(:,1:3),geom.elements.coordinates(:,1),geom.elements.coordinates(:,2),u);
    %hold on;
    ERR.INF=Err(@u_re,@gradiente,u,1,NdofElement);
    ERR.L2=Err(@u_re,@gradiente,u,2,NdofElement);
    ERR.H1=Err(@u_re,@gradiente,u,3,NdofElement);
    %h(t)=TheMesh.HeightMax;
    Area(t)=TheMesh.AreaMax;
    h(t)=Area(t)^0.5;
    %%Salvo gli errori associati alla mesh
    TheMesh.Meshes(t).e_L2=ERR.L2;
    TheMesh.Meshes(t).e_H1=ERR.H1;
    TheMesh.Meshes(t).e_INF=ERR.INF;
%   TheMesh.Meshes(t).h=h(t);
    TheMesh.Meshes(t).Area=Area(t);
   e_L2(t)=ERR.L2;
   e_INF(t)=ERR.INF;
   e_H1(t)=ERR.H1;
   PeTest(t)=PeMax(TheMesh,beta,epsilon);
%    PeTest(t)
end
pL2=polyfit(log(h),log([TheMesh.Meshes(:).e_L2]),1)
pH1=polyfit(log(h),log([TheMesh.Meshes(:).e_H1]),1)
pINF=polyfit(log(h),log([TheMesh.Meshes(:).e_INF]),1);
yLinf=polyval(pINF,log(h));
yL2=polyval(pL2,log(h));
yH1=polyval(pH1,log(h));

% plot(log(h),yLinf,'g',LineWidth=2);
% hold on
% plot(log(h),yH1,'b',LineWidth=2);
% plot(log(h),yL2,'r',LineWidth=2);
% xlabel('log(h)');
% ylabel('Errore');
% legend('Errore in L_{inf}','Errore in H1','Errore in L_{2}');
% title(["Ordini di convergenza-P2-Numero di Peclet che va da 1240 circa fino a 240 nell'ultimo refinement"]);
% %title(['Ordini di convergenza-P2-Neummann con stabilizzazione,Neumann uscente']);
% hold off

Conditionscale=polyfit(log(h),log(ConditioningNumbers),1);
ScaleConditioning=polyval(Conditionscale,log(h));
hold on;
plot(log(h),ScaleConditioning,"r",LineWidth=2);


function u=built_sol(x,u_d)
global geom
    u=zeros(length(geom.elements.coordinates(:,1)),1);
    for j=1:length(geom.elements.coordinates(:,1))
        if(geom.pivot.pivot(j)>0)
            u(j,1)=x(geom.pivot.pivot(j),1);
        else
            u(j,1)=u_d(-geom.pivot.pivot(j),1);
        end
    end
end

function f=fun(CG,beta)
     global epsilon;
     gradi=gradiente(CG);
     f=32*epsilon*((CG(1)*(1-CG(1))+CG(2)*(1-CG(2))))+beta'*(gradi)';
end

function b_E=NeumannAssembly
global b_E;
global geom;global epsilon;
    phi1N=[1,0,0]; 
    phi2N=[0,1,0];
    phi3N=[0,0,1];

    for l=1:length(geom.pivot.Ne(:,1))
        V0=geom.elements.borders(geom.pivot.Ne(l,1),1);        
        V1=geom.elements.borders(geom.pivot.Ne(l,1),2);
        Vc=geom.elements.borders(geom.pivot.Ne(l,1),5);
        cord_V1=geom.elements.coordinates(V1,:);
        cord_V0=geom.elements.coordinates(V0,:);
        cord_Vc=geom.elements.coordinates(Vc,:);

        marker=geom.pivot.Ne(l,2); % marker del lato
        ii_0=geom.pivot.pivot(V0);
        if(ii_0>0)           
            %b_E(ii_0,1)=b_E(ii_0,1)+norm(cord_V1-cord_V0)*(1/3*Normal_Gradient(marker,geom.elements.coordinates(V0,:))+1/6*Normal_Gradient(marker,geom.elements.coordinates(V1,:)));
            b_E(ii_0,1)=b_E(ii_0,1)+1/6*norm(cord_V1-cord_V0)*(phi1N(1)*Normal_Gradient(marker,geom.elements.coordinates(V0,:))+...
                4*phi1N(2)*Normal_Gradient(marker,geom.elements.coordinates(Vc,:))+phi1N(3)*Normal_Gradient(marker,geom.elements.coordinates(V1,:)))*epsilon;
           
        end
        ii_1=geom.pivot.pivot(V1);
        if(ii_1>0)           
            %b_E(ii_1,1)=b_E(ii_1,1)+norm(cord_V1-cord_V0)*(1/3*Normal_Gradient(marker,geom.elements.coordinates(V1,:))+1/6*Normal_Gradient(marker,geom.elements.coordinates(V0,:)));
            b_E(ii_1,1)=b_E(ii_1,1)+1/6*norm(cord_V1-cord_V0)*(phi3N(1)*Normal_Gradient(marker,geom.elements.coordinates(V0,:))+...
                4*phi3N(2)*Normal_Gradient(marker,geom.elements.coordinates(Vc,:))+phi3N(3)*Normal_Gradient(marker,geom.elements.coordinates(V1,:)))*epsilon;
        end
        ii_C=geom.pivot.pivot(Vc);
          if(ii_C>0)         
            %b_E(ii_C,1)=b_E(ii_1,1)+norm(cord_V1-cord_V0)*(1/3*Normal_Gradient(marker,geom.elements.coordinates(V1,:))+1/6*Normal_Gradient(marker,geom.elements.coordinates(V0,:)));
             b_E(ii_C,1)=b_E(ii_C,1)+1/6*norm(cord_V1-cord_V0)*(phi2N(1)*Normal_Gradient(marker,geom.elements.coordinates(V0,:))+...
                4*phi2N(2)*Normal_Gradient(marker,geom.elements.coordinates(Vc,:))+phi2N(3)*Normal_Gradient(marker,geom.elements.coordinates(V1,:)))*epsilon;
        end
        
    end %for principale
    
end

%valutazione del gradiente nella direzione normale:
function Value=Normal_Gradient(marker,Vertex)
    %vertex: coordinate in cui valuto la funzione:
    Normal_Vect=zeros(1,2); 
    switch marker
        case 2
            Normal_Vect=[0,-1];
        case 4
            Normal_Vect=[1,0];
        case 6
            Normal_Vect=[0,1];
        case 8
            Normal_Vect=[-1,0];    
    end %fine dello switch
    grad_Value=gradiente(Vertex);
    grad_Value=grad_Value';
    Value=Normal_Vect*grad_Value;
end

function grad=gradiente(Vertex)
    grad=zeros(1,2);
    grad(1,1)=(16-32*Vertex(1,1))*Vertex(1,2)*(1-Vertex(1,2));
    grad(1,2)=(16-32*Vertex(1,2))*Vertex(1,1)*(1-Vertex(1,1));
end

function Ure=u_re(Vertex)
    %%Evaluate the real solution:
    Ure=16*Vertex(1,1)*(1-Vertex(1,1))*Vertex(1,2)*(1-Vertex(1,2));
end
function PEMAX=PeMax(TheMesh,beta,epsilon)
global geom;
PEMAX=0;
for e=1:geom.nelements.nTriangles
    PeE=TheMesh.Peclet(beta,epsilon,e);
    if PeE>PEMAX
        PEMAX=TheMesh.Peclet(beta,epsilon,e);
    end
end
end
function PlotError(el2,eh1,h)
    subplot(2,1,1);
    loglog(h,el2,'r',LineWidth=2);
    ylabel('log(e_{L_{2}})');
    xlabel('log(h)');
    title('Norma in L_{2}')

    subplot(2,1,2);
    loglog(h,eh1,'b',LineWidth=2);
    xlabel('log(h)');
    ylabel('log(e_{H{1}})');
    title('Norma in H_{1}');
    
end
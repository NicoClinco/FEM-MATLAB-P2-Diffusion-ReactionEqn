function [A,b_E,Ad,u_d] = Stiffness_Matrix(BC,forzante,beta,~)
global geom; %cerco nello scope globale
global epsilon;
%%Leggo innanzitutto il numero di elementi
A=spalloc(max(geom.pivot.pivot(:,1)),max(geom.pivot.pivot(:,1)),10*max(geom.pivot.pivot(:,1))); %Stiffness
Ad=spalloc(max(geom.pivot.pivot(:,1)),-min(geom.pivot.pivot(:,1)),10*max(geom.pivot.pivot(:,1))); %Vettore di Dirichlet
u_d=zeros(-min(geom.pivot.pivot(:,1)),1); %%Condizione di Dirichlet
b_E=zeros(max(geom.pivot.pivot(:,1)),1); %RHS termine noto

delta_x=zeros(1,3);
delta_y=zeros(1,3);

%%COSTRUZIONE DELLA PARTE DEL MAPPING%%%%
[N3,N1,N2,OMEGA]=int_nodes_weights(5);
Nquad=length(OMEGA);
[phi,gradphi,~]=mapping.map(N3,N1,N2,OMEGA,6);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for e=1:geom.nelements.nTriangles

    delta_x(1)=geom.elements.coordinates(geom.elements.triangles(e,3),1)-geom.elements.coordinates(geom.elements.triangles(e,2),1);
    delta_y(1)=geom.elements.coordinates(geom.elements.triangles(e,2),2)-geom.elements.coordinates(geom.elements.triangles(e,3),2);
    delta_x(2)=geom.elements.coordinates(geom.elements.triangles(e,1),1)-geom.elements.coordinates(geom.elements.triangles(e,3),1);
    delta_y(2)=geom.elements.coordinates(geom.elements.triangles(e,3),2)-geom.elements.coordinates(geom.elements.triangles(e,1),2);
    delta_x(3)=geom.elements.coordinates(geom.elements.triangles(e,2),1)-geom.elements.coordinates(geom.elements.triangles(e,1),1);
    delta_y(3)=geom.elements.coordinates(geom.elements.triangles(e,1),2)-geom.elements.coordinates(geom.elements.triangles(e,2),2);
    LocalArea=geom.support.TInfo(e).Area; %area del triangolo che sto considerando
    %Tau=TheMesh.Tau(beta,epsilon,e); %%fattore di correzione

    %Matrice del cambio di coordinate:
    B=mapping.b(geom.elements.triangles(e,:));

    for j=1:6
        jj=geom.pivot.pivot(geom.elements.triangles(e,j)); %dof
        if(jj>0)
            for k=1:6
                kk=geom.pivot.pivot(geom.elements.triangles(e,k)); %dof riferito a k
                if(kk>0)
                    %%diffusione-convezione-reazione
                    diff=0; conv=0; react=0;
                    for nq=1:Nquad
                        %diffusione
                        diff=diff+(inv(B')*[gradphi(k).x(nq);gradphi(k).y(nq)])'*(inv(B')*[gradphi(j).x(nq);gradphi(j).y(nq)])*2*LocalArea*epsilon*OMEGA(nq);
                            
                        %convezione
                        conv=conv+(beta'*(inv(B')*[gradphi(k).x(nq);gradphi(k).y(nq)]))*phi(j,nq)*2*LocalArea*OMEGA(nq);
                        %reazione
                        %
                    end
                    A(jj,kk)=A(jj,kk)+diff+conv+react;

                else
                    %Parte di Dirichlet
                    %%diffusione-convezione-reazione
                    diff1=0; conv1=0; react1=0;
                    for nq=1:Nquad
                        %diffusione
                        diff1=diff1+(inv(B')*[gradphi(k).x(nq);gradphi(k).y(nq)])'*(inv(B')*[gradphi(j).x(nq);gradphi(j).y(nq)])*...
                            2*LocalArea*epsilon*OMEGA(nq);
                        %convezione
                        conv1=conv1+(beta'*(inv(B')*[gradphi(k).x(nq);gradphi(k).y(nq)]))*phi(j,nq)*2*LocalArea*OMEGA(nq);
                        %reazione
                        %
                    end
                    Ad(jj,-kk)=Ad(jj,-kk)+diff1+conv1+react1;
                    %%i nodi sono in ordine
                    u_d(-kk,1)=BC.Values(geom.pivot.Di(-kk,2));

                end %% if(kk>0)

            end %% for k=1:6
            %quadratura sulla forzante
            f=0;
            V3=[geom.elements.coordinates(geom.elements.triangles(e,3),1),geom.elements.coordinates(geom.elements.triangles(e,3),2)];
            for nq=1:Nquad              
                f=f+OMEGA(nq)*2*LocalArea*forzante(mapping.FE(V3,B,N1(nq),N2(nq)),beta)*phi(j,nq);
            end
              b_E(jj,1)=b_E(jj,1)+f;
        end %if jj>0

    end %for j=1:6

end %for e=1:nTriangoli


end
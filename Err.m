function [Error] = Err(u,UGrad,u_h,num,NdofE)
global geom;
%%Ndof: numero gradi di libertà sull'elemento considerato
    Error=0;
    if num==1
        %%norma infinito
        max_err=0;
        for j=1:length(geom.elements.coordinates(:,1)) 
            err=abs(u_h(j,1)-u(geom.elements.coordinates(j,:)));
            if err>max_err
                max_err=err;
            end     
        end
        Error=max_err;
    elseif num==2 %%norma L^2
        %%%%%%%%%%%%%%%%%%%%%%%%MAPPING SUI TRIANGOLI %%%%%%%%%%%%%%%%%%%%
            %%INIZIALIZZO LA MATRICE PER SALVARE N1,N2,N3
            [N3,N1,N2,OMEGA]=int_nodes_weights(5); %% N3: 1-x-y  %N2: y  %N1 : x
            Nquad=length(OMEGA); %Numero dei nodi di quadratura
            [Phi_,~,~]=mapping.map(N3,N1,N2,OMEGA,NdofE);
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            for e=1:length(geom.elements.triangles)
                LocalError=0; %errore sul triangolo
                B=mapping.b(geom.elements.triangles(e,:));
                V3=[geom.elements.coordinates(geom.elements.triangles(e,3),1),geom.elements.coordinates(geom.elements.triangles(e,3),2)];  
                for k=1:Nquad
                    %Calcolo la FE innanzitutto:                 
                    sum=0; %%sommatoria valutata nel nodo
                    for j=1:NdofE
                        sum=sum+u_h(geom.elements.triangles(e,j))*Phi_(j,k);
                    end
                    LocalError=LocalError+2*(geom.support.TInfo(e).Area)*OMEGA(k)*(u(mapping.FE(V3,B,N1(k),N2(k)))-sum)^2;
                end %k=1:Nquad
                 Error=Error+LocalError;
            end %e=1:nTriangles
            Error=Error^0.5;
            
    elseif num==3 %%norma H1
        %%Errore del gradiente:
        %%%%%%%%%%%MAPPING SUI
        %%%%%%%%%%%TRIANGOLI%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [N3,N1,N2,OMEGA]=int_nodes_weights(5);
        Nquad=length(OMEGA); %Numero dei nodi di quadratura
        [Phi_,Gradphi,~]=mapping.map(N3,N1,N2,OMEGA,NdofE);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        for e=1:length(geom.elements.triangles)
            LocalError=0; %%Errore locale SUL triangolo
            
            B=mapping.b(geom.elements.triangles(e,:)); %cambio di coordinate
            InvBT=inv(B');                             %matrice inversa
            V3=[geom.elements.coordinates(geom.elements.triangles(e,3),1),geom.elements.coordinates(geom.elements.triangles(e,3),2)];
            for k=1:Nquad
                %Ciclo interno con j
                sum=zeros(2,1);
                for j=1:NdofE
                    sum=sum+u_h(geom.elements.triangles(e,j))*(InvBT*[Gradphi(j).x(k);Gradphi(j).y(k)]); %%Vettore riga                
                end             
                LocalError=LocalError+2*(geom.support.TInfo(e).Area)*OMEGA(k)*[UGrad(mapping.FE(V3,B,N1(k),N2(k)))-sum']*[UGrad(mapping.FE(V3,B,N1(k),N2(k)))-sum']';

            end %end k
            Error=Error+LocalError;

        end %1:nTriangoli
        NormL2=Err(u,UGrad,u_h,2,NdofE); %Prendo anche la norma L2 :|.|H1 = |.|L2 + |.'|L2
        Error=(NormL2^2+Error)^0.5;      
    end %end if
    %%%FUNZIONE DEL MAPPING SUI TRIANGOLI
end %end funzione
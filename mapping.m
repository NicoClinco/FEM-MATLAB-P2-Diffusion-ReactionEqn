classdef mapping

    properties
    end

    methods (Static)
        %costruzione della soluzione per mappare roba
        function [PHI,GRADPHI,HessPhi]=map(N3,N1,N2,OMEGA,Ndof)
            Nquad=length(OMEGA);
            phi=zeros(Ndof,Nquad);
            if Ndof==3
                %%Caso dei P1
                phi(1,:)=N1(:);phi(2,:)=N2(:);phi(3,:)=N3(:);
                PHI=phi;
                gradphi1=struct('x',1*ones(1,3),'y',0*ones(1,3)); % Gradphi(1).x(k)=1 Gradphi(1).y(k)=0; d(phi1(k))/dx , d(phi1(k))/dy
                gradphi2=struct('x',0*ones(1,3),'y',1*ones(1,3)); % Gradphi(2).x(k)=0 Gradphi(2).y(k)=1
                gradphi3=struct('x',-1*ones(1,3),'y',-1*ones(1,3));
                GRADPHI=[gradphi1,gradphi2,gradphi3];

            elseif Ndof==6

                phi(1,:)=2*N1(:)'.*(N1(:)'-1/2*(ones(1,Nquad)));
                phi(2,:)=2*N2(:)'.*(N2(:)'-1/2*(ones(1,Nquad)));
                phi(3,:)=2*N3(:)'.*(N3(:)'-1/2*(ones(1,Nquad)));
                phi(4,:)=4*N1(:)'.*(N2(:)');
                phi(5,:)=4*N3(:)'.*(N2(:)');
                phi(6,:)=4*N3(:)'.*(N1(:)');
                PHI=phi;

                gradphi1=struct('x',4*N1(:)'-1*ones(1,Nquad),'y',0*ones(1,Nquad));
                gradphi2=struct('x',0*ones(1,Nquad),'y',4*N2(:)'-1*ones(1,Nquad)); 
                gradphi3=struct('x',4*(N1(:)'+N2(:)')-3*ones(1,Nquad),'y',4*(N1(:)'+N2(:)')-3*ones(1,Nquad));
                gradphi6=struct('x',4*(ones(1,Nquad)-2*N1(:)'-N2(:)'),'y',-4*N1(:)');
                gradphi4=struct('x',4*N2(:)','y',4*N1(:)');
                gradphi5=struct('x',-4*N2(:)','y',4*(1*ones(1,Nquad)-2*N2(:)'-N1(:)'));
                GRADPHI=[gradphi1,gradphi2,gradphi3,gradphi4,gradphi5,gradphi6];
                Hess1=[4,0;0,0];
                Hess2=[0,0;0,4];
                Hess3=[4,4;4,4];
                Hess4=[0,4;4,0];
                Hess5=[0,-4;-4,-8];
                Hess6=[-8,-4;-4,0];
                HessPhi=Hess1;
                HessPhi(:,:,2)=Hess2;
                HessPhi(:,:,3)=Hess3;
                HessPhi(:,:,4)=Hess4;
                HessPhi(:,:,5)=Hess5;
                HessPhi(:,:,6)=Hess6;

            end


        end %%end della funzione che costruisce PHI
        
        function B=b(triangolo)
            global geom;
            %triangolo: Sarebbe un vettore colonna che contiene i vertici
            %del triangolo 
            V1=geom.elements.coordinates(triangolo(1),:);
            V2=geom.elements.coordinates(triangolo(2),:);
            V3=geom.elements.coordinates(triangolo(3),:);

            B=[V1(1,1)-V3(1,1),V2(1,1)-V3(1,1);
                V1(1,2)-V3(1,2),V2(1,2)-V3(1,2)];

        end
        function roto=M(V1,V2,V3)
            %matrice di rototraslazione
            B=mapping.b(V1,V2,V3);
            roto=inv(B)*inv((B'));
        end
        function f=FE(V3,B,xq,yq)
            f=V3+(B*[xq;yq])';
        
        end
    end
end

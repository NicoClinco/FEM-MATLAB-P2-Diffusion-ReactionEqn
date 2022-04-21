classdef Mesh_info
    properties
        geom; %%Memorizzo qui dentro la struttura dati intera delle mesh differenti;
        Meshes; %%vettore di strutture dati: 
        error; %%vettore in cui vado a salvare il tipo di errore
        errorID=["e_{L2}","e_{H1}","e_{INF}"];
    
    end
    methods
        function A_max=AreaMax(obj)
            %%troviamo il massimo:
            A_max=0;
            for i=1:obj.geom.nelements.nTriangles
                A=obj.geom.support.TInfo(i).Area;
                if A>A_max
                    A_max=A;
                end
            end
        end %end metodo per le aree% 
        function h_max=HeightMax(obj)
            h_max=0;
            for i=1:length(obj.geom.elements.borders)
                h=norm(obj.geom.elements.coordinates(obj.geom.elements.borders(i,1),:)-obj.geom.elements.coordinates(obj.geom.elements.borders(i,2),:));
                if h>h_max
                    h_max=h;
                end
            end

        end %end metodo per h_max %
        function hE=HeighTriangle(obj,e)
                %e: indice dell'elemento
                hE=(obj.geom.support.TInfo(e).Area)^0.5;
        end
        function PeE=Peclet(obj,Beta,Eps,e)
            %Calcolo del Peclet Locale.
            PeE=1/3*(norm(Beta))*obj.HeighTriangle(e)/(2*Eps);

        end
        function t=Tau(obj,Beta,Eps,e)
            PE=obj.Peclet(Beta,Eps,e); %%Calcolo il Peclet
            HE=obj.HeighTriangle(e);
            if PE>1
                t=HE/(2*norm(Beta));
            else
                t=(HE^2)/(12*Eps);
            end
        end

    end
   
end
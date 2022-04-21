
nnode = geom.nelements.nVertexes;

for l=1:geom.nelements.nBorders

  n(1)=geom.elements.borders(l,1);
  n(2)=geom.elements.borders(l,2);
  e(1)=geom.elements.borders(l,3);
  e(2)=geom.elements.borders(l,4);
  
  nnode = nnode + 1;

  geom.elements.coordinates(nnode,:) = (geom.elements.coordinates(n(1),:)+...
					geom.elements.coordinates(n(2),:))/2;

  geom.elements.borders(l,5)=nnode; %aggiungo il numero del nodo in borders
  
  idx = [1 2 3];

  for el=e %ciclo su e(1) ed e(2) 
    
    if(el ~= -1)
      acc = 0;
      acc = idx * ( geom.elements.triangles(el,1:3)==n(1) )'; 
      %geom.elements.tri(el,1:3)==n(1) restituisce un vettore che contiene
      %1 dove è presente 1 mentre il resto sono zeri.es: [0,1,0]
      acc = acc + idx * ( geom.elements.triangles(el,1:3)==n(2) )';

      switch acc
	case 3
	  geom.elements.triangles(el,4) = nnode;        
	case 4
	  geom.elements.triangles(el,6) = nnode;
	case 5
	  geom.elements.triangles(el,5) = nnode;
	otherwise
	  disp('sconoscuto');
      end % switch acc      
    end

  end % for el=e

  VertexValue=[0 0];
  Vertex=[0 0];
  D= [0 0];
  
  idxV = 1:length(geom.input.BC.InputVertexValues); %indice dei vertici

  % Se il lato e` di bordo
  if( any(e==-1) )
    % Lato di Dirichlet
    if( geom.pivot.nodelist(n(1))~=0 && geom.pivot.nodelist(n(2))~=0 )  %Entrambi i nodi sono di Dirichlet
 
      if( geom.pivot.nodelist(n(1)) == geom.pivot.nodelist(n(2)) )   %Stesso valore di Dirichlet
	    Di = geom.pivot.nodelist(n(1));
	    geom.pivot.nodelist(nnode)= Di;
	    geom.pivot.Di(end+1,:) = [nnode, Di]; %%aggiorno la lista in geom.pivot.Di
	    geom.pivot.pivot(nnode) = min(geom.pivot.pivot)-1;
      else
      	if( any(geom.input.BC.InputVertexValues==geom.pivot.nodelist(n(1))) ) %Se n(1) ha lo stesso valore di un vertice:
            
      	  VertexValue(1) = 1;% Vertice(1) con marker speciale
          % e` il vettore con un uno in corrispondenza del vertice
          % del lato con marker speciale

          % Marker del lato in BC.Boundary.Values che contiene n(1): PEZZO
          % DI CODICE PER RICAVARE IL VALORE DEL MARKER DEL LATO
    	  Vertex(1) = geom.input.BC.Boundary.Values*...
              (geom.input.BC.InputVertexValues == ...
              geom.pivot.nodelist(n(1)))'; %Posizione in InputVertexValues di n(1)) [0,0,0,1] ad es.

          % Marker del vertice del poligono iniziale pari a n1: PEZZO DI
          % CODICE PER RICAVARE IL VALORE DEL MARKER DEL VERTICE
    	  InputVertexValue = geom.input.BC.InputVertexValues*...
    	      (geom.input.BC.InputVertexValues == geom.pivot.nodelist(n(1)))'; %Posizione in InputVertexValues di n(1))
    	  D(1) = geom.pivot.nodelist(n(2));          
        end

	if( any(geom.input.BC.InputVertexValues==geom.pivot.nodelist(n(2))) )
	  VertexValue(2) = 1;

	  Vertex(2) = geom.input.BC.Boundary.Values*...
	  (geom.input.BC.InputVertexValues == geom.pivot.nodelist(n(2)))';

	  InputVertexValue = geom.input.BC.InputVertexValues*...
	      (geom.input.BC.InputVertexValues == geom.pivot.nodelist(n(2)))';

	  D(2) = geom.pivot.nodelist(n(1));

	end

	if( sum(VertexValue) ~= 2 )
	  Di = VertexValue*D'; % nodo con condizione di Dirichlet
	  
          geom.pivot.nodelist(nnode)= Di;
          geom.pivot.Di(end+1,:) = [nnode, Di];
          geom.pivot.pivot(nnode) = min(geom.pivot.pivot)-1;
	else
          % diamo al nuovo nodo il marker del lato
          % il lato che stiamo analizzando e` un lato del poligono
          % iniziale:
          
	  % l'indice del lato e` quello del nodo di inizio di quel lato
%%%%%%%%%%%%%%%%%%%%%commentata%%%%%%%%%%%%%%%%%%
% 	  if( max(Vertex)-min(Vertex)>1 ) % siamo sul lato di chiusura
% 	    Di = geom.input.BC.Boundary.Values(max(Vertex));
% 	  else % siamo sui lati 1->2->3->4->
% 	    Di = geom.input.BC.Boundary.Values(min(Vertex));
% 	  end
%%%%%%%%%%%%%%%%%%%commentata%%%%%%%%%%%%%%%%%%%%%%
          % check della condizione di Neumann aperta
          if( rem(Di,2)== 0 ) % nodo con grado di liberta`, lato di
              % Dirichlet aperto
      	    geom.pivot.nodelist(nnode)= 0; 
    	    geom.pivot.pivot(nnode) = max(geom.pivot.pivot)+1;
          else
              geom.pivot.nodelist(nnode)= Di;
              geom.pivot.Di(end+1,:) = [nnode, Di];
              geom.pivot.pivot(nnode) = min(geom.pivot.pivot)-1;
          end
	end
      end % else
    % Lato di Neumann
    else %Qui ho lato con estremi di Neumann
      geom.pivot.nodelist(nnode) = 0;
      geom.pivot.pivot(nnode) = max(geom.pivot.pivot)+1;
    end % if( geom.pivot.nodelist(n(1))~=0 & geom.pivot.nodelist(n(2))~=0 )

  else %%ho un lato interno alla triangolazione
      geom.pivot.nodelist(nnode) = 0;
      geom.pivot.pivot(nnode) = max(geom.pivot.pivot)+1;
    end %if( any(e==-1) )
  
end % for l=1:geom.nelements.nBorders

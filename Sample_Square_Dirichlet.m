if(~exist('bbtr30'))
     addpath('../bbtr30')
     disp('../bbtr30 added to the path')
end

%----------------------------------------------------------------------------
%
% Triangolazione di un dominio quadrata
% con condizioni di Dirichlet sul bordo
%
%----------------------------------------------------------------------------
%
%  Autore: Stefano Berrone
%  Politecnico di Torino
%
%----------------------------------------------------------------------------

%clc
%clear all

% -------------------------------
% Inserimento dei vertici
% -------------------------------

Domain.InputVertex = [ 0 0
                       1 0
                       1 1
                       0 1];


% ---------------------------------------------
% Definizione del dominio a partire dai Vertici
% ---------------------------------------------

% Dichiaro le variabili per delimitare il dominio
Domain.Boundary.Values = 1:4;
% lato di bordo 1 dal nodo 1 al nodo 2
% lato di bordo 2 dal nodo 2 al nodo 3
% lato di bordo 3 dal nodo 3 al nodo 4
% lato di bordo 4 dal nodo 4 al nodo 1

Domain.Holes.Hole = [];       % non ci sono buchi nel dominio
Domain.Segments.Segment = []; % non ci sono lati forzati nel dominio

% --------------------------------------------------
% Definizione delle condizioni al contorno a partire
% dai Vertici e dai lati di bordo
% --------------------------------------------------

% valori numerici per le condizioni al contorno
BC.Values = [0 0.0 0 0.0 0 0.0 0 0.0 0.0];

% marker delle condizioni al contorno sui bordi del dominio
% dispari -> Dirichlet; pari -> Neumann
BC.Boundary.Values = [1 1 6 1]; %Lato Neumann
% marker dei Vertici iniziali
BC.InputVertexValues = [1 1 1 1];
% Questi indici posso essere anche indici ai valori numerici
% contenuti nel vettore BC.Values

BC.Holes.Hole = [];
BC.Segments.Segment = [];



% --------------------------------------------
% Inserimento dei parametri di triangolazione
% --------------------------------------------

RefiningOptions.CheckArea  = 'Y';
RefiningOptions.CheckAngle = 'N';
%RefiningOptions.AreaValue  = 0.001;
RefiningOptions.AngleValue = [30];
RefiningOptions.Subregions = [];


% --------------------------------------------
% Creazione della triangolazione e plottaggio
% --------------------------------------------

[geom] = bbtr30(Domain,BC,RefiningOptions);
%draw_grid (geom,1);

% --------------------------------------------------
% --------------------------------------------------


% --------------------------------------------------
% Rielaborazione dei prodotti del triangolatore
% per un piu` agevole trattamento delle condizioni
% al contorno
% --------------------------------------------------

geom.elements.coordinates = geom.elements.coordinates(...
				1:geom.nelements.nVertexes,:);
geom.elements.triangles = geom.elements.triangles(...
				1:geom.nelements.nTriangles,:);
geom.elements.borders = geom.elements.borders(...
				1:geom.nelements.nBorders,:);
geom.elements.neighbourhood = geom.elements.neighbourhood(...
				1:geom.nelements.nTriangles,:);

% --------------------------------------------------

j  = 1;
Dj = 1;
for i=1:size(geom.pivot.nodelist)
     if geom.pivot.nodelist(i)==0
        geom.pivot.pivot(i)=j;
        j = j+1;
     else
        geom.pivot.pivot(i)=-Dj;
        Dj = Dj + 1;
     end
end

% --------------------------------------------------

geom.pivot.pivot = transpose(geom.pivot.pivot);

% --------------------------------------------------

% geom.pivot.Di dopo le operazioni seguenti contiene l`indice dei nodi
% di Dirichlet e il corrispondente marker

[X,I] = sort(geom.pivot.Di(:,1));
geom.pivot.Di = geom.pivot.Di(I,:);

clear X I;

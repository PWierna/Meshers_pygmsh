clear,clc
Connectivities_and_Coordinates_2D;
element = 103;

%Get nodal coordinates:
elemn   = MACRO_MODEL.Conectivity(element,2:end);
ncoords = MACRO_MODEL.Coordinates(elemn,:);

%Coordinates of the master element:
pos_nodal_master = [-1,-1 ; 0,-1 ; 1,-1;
                    1,0 ; 1,1, ; 0,1 ; 
                    -1,1 ; -1,0 ; 0,0 ];

%Plot:
figure()
subplot(1,3,2) %Master element
plot(pos_nodal_master(:,1),pos_nodal_master(:,2),'or','LineWidth',1.5);hold on;grid on
for n=1:length(elemn)
    text(pos_nodal_master(n,1),pos_nodal_master(n,2),['\leftarrow ',char(string(n-1))],'FontSize',15)
end
xlabel('X');ylabel('Y');zlabel('Z');title('Master element')
subplot(1,3,1) %Current ordering
plot(MACRO_MODEL.Coordinates(elemn,1),MACRO_MODEL.Coordinates(elemn,2),'or','LineWidth',1.5);hold on;grid on
for n=1:length(elemn)
    text(MACRO_MODEL.Coordinates(elemn(n),1),MACRO_MODEL.Coordinates(elemn(n),2),['\leftarrow ',char(string(n-1))],'FontSize',15)
end
xlabel('X');ylabel('Y');zlabel('Z');title('Current')

%INSERT HERE THE ORDERING OBSERVED--------------------------------------------------%
idx_gmsh2ref  = [0,4,1,5,2,6,3,7,8]+1; 
%-----------------------------------------------------------------------------------%

idx_master2ref = 1:9;
[~,idx_ref2master] = sort(idx_master2ref);

%Obtain resortind index:
idx_gmsh2master = idx_gmsh2ref(idx_ref2master)
idx_reordering = idx_gmsh2master-1;

%Plot resorted element:
elemn=elemn(idx_gmsh2master);
subplot(1,3,3) %Reordered
plot(MACRO_MODEL.Coordinates(elemn,1),MACRO_MODEL.Coordinates(elemn,2),'or','LineWidth',1.5);hold on;grid on
for n=1:length(elemn)
    text(MACRO_MODEL.Coordinates(elemn(n),1),MACRO_MODEL.Coordinates(elemn(n),2),['\leftarrow ',char(string(n-1))],'FontSize',15)
end
xlabel('X');ylabel('Y');zlabel('Z');title('Reordered')

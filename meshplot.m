function meshgen = meshplot(coord,lc,mesh_algo)
%converts the data points to geomtery and meshes in GMSH.
fileid = fopen('points.geo','w'); %create .geo file

l = length(coord);
fprintf(fileid,'Mesh.Algorithm = %d;\n',mesh_algo); %type of element 
fprintf(fileid,'Mesh.ElementOrder=2;\n');
fprintf(fileid,'lc = %d;\n',lc);

%loop for bspline
for i = 1:length(coord) 
    Point = [i,coord(i,1),coord(i,2),0,lc];
    fprintf(fileid,'Point(%d)={%d,%d,%d,%d};\n',Point(1,:));
end 
spline=(1:l);
fprintf(fileid,'Spline(%d) = {',l+1);
allOneString = fprintf(fileid,'%.0f,' , spline);
allOneString = allOneString(1:end);
fprintf(fileid,'1};\n');
fprintf(fileid,'Curve loop(1) = {%d};\n',l+1);
fprintf(fileid,'Physical Curve(1) = {1};\n');

%coordinates for the plate
fprintf(fileid,'Point(%d)={0,0,0,%d};\n',l+2,lc);
fprintf(fileid,'Point(%d)={1,0,0,%d};\n',l+3,lc);
fprintf(fileid,'Point(%d)={1,1,0,%d};\n',l+4,lc);
fprintf(fileid,'Point(%d)={0,1,0,%d};\n',l+5,lc);
fprintf(fileid,'Line(1)={%d,%d};\n',l+2,l+3);
fprintf(fileid,'Line(2)={%d,%d};\n',l+3,l+4);
fprintf(fileid,'Line(3)={%d,%d};\n',l+4,l+5);
fprintf(fileid,'Line(4)={%d,%d};\n',l+5,l+2);
fprintf(fileid,'Curve loop(2)={1,2,3,4};\n');
fprintf(fileid,'Physical Curve(2)={2};\n');
fprintf(fileid,'Plane Surface(1) = {2,1};\n'); %plate with spline as hole
fprintf(fileid,'Physical Surface(1) = {1};\n');

%enable below line to mesh in quads
% fprintf(fileid,'Recombine Surface{1};\n'); 
fprintf(fileid,'Mesh 2;\n');
fprintf(fileid,'Save "points.m";\n');
fclose(fileid);

%run GMSH
%Look for non interative command line to call gmsh https://gmsh.info/doc/texinfo/gmsh.html
command = ('gmsh points.geo -2 format string m');
system(command);

%plot mesh
points;
fopen("points.m",'r');
x=msh.POS(:,1);
y=msh.POS(:,2);
L = [msh.LINES3(:,1:3),msh.LINES3(:,1)];
T = msh.TRIANGLES6(:,1:3);
t2=msh.TRIANGLES6(:,4:6);
x1= msh.POS(t2,1);
y1=msh.POS(t2,2);
figure();
meshgen = trimesh(T,x,y);
hold on
meshgen = trimesh(L,x,y);
hold on
plot(x1,y1,'o','markersize',1);



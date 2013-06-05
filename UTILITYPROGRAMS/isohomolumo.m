clear;
homo0 = 0.01;
lumo0 = 0.01;

  MESH = [ 200 200 200 ];
LIMITS = [ 15.0 15.0 15.0 ];


[I1,I2,I3,RHOHOMO]=textread('LUMODENS.dat','%f %f %f %f',MESH(1)*MESH(2)*MESH(3));
[I1,I2,I3,RHOLUMO]=textread('HOMODENS.dat','%f %f %f %f',MESH(1)*MESH(2)*MESH(3));

for i=1:MESH(1)*MESH(2)*MESH(3)
	DENSHOMO(I1(i),I2(i),I3(i)) = RHOHOMO(i);
	DENSLUMO(I1(i),I2(i),I3(i)) = RHOLUMO(i);
	x(I1(i)) = -LIMITS(1) + 2*LIMITS(1)*(I1(i)-1)/(MESH(1)-1);
	y(I2(i)) = -LIMITS(2) + 2*LIMITS(2)*(I2(i)-1)/(MESH(2)-1);
	z(I3(i)) = -LIMITS(3) + 2*LIMITS(3)*(I3(i)-1)/(MESH(3)-1);
	n1(I1(i)) = I1(i);
	n2(I2(i)) = I2(i);
	n3(I3(i)) = I3(i);
end

for i=1:MESH(1)
	for j=1:MESH(2)
		fhomo(i,j) = DENSHOMO(i,j,100);
		flumo(i,j) = DENSLUMO(i,j,100);
	end
end

[f,v] = isosurface(x,y,z,DENSHOMO,homo0);

phomo = patch('Faces',f,'Vertices',v)

isonormals(n1,n2,n3,DENSHOMO,phomo);
set(phomo,'FaceColor','b','EdgeColor','none','NormalMode','auto');
alpha(0.50);set(gca,'Xcolor','w');set(gca,'Ycolor','w');set(gca,'Zcolor','w')
camlight;
lighting gouraud;
set(gca,'Xcolor','w');set(gca,'Ycolor','w');set(gca,'Zcolor','w')
hold on;


[f,v] = isosurface(x,y,z,DENSLUMO,lumo0);

plumo = patch('Faces',f,'Vertices',v)

isonormals(n1,n2,n3,DENSLUMO,plumo);
set(plumo,'FaceColor','r','EdgeColor','none','NormalMode','auto');
alpha(0.50);
camlight;
lighting gouraud;

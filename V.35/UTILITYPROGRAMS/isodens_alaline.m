clear;
rho0 = 0.050;
rho1 = 0.250;
rho2 = 0.500;


  MESH = [ 200 200 200 ];
LIMITS = [ 10.0 10.0 10.0 ];


[I1,I2,I3,RHO]=textread('CHARGEDENS.4.dat','%f %f %f %f',MESH(1)*MESH(2)*MESH(3));

for i=1:MESH(1)*MESH(2)*MESH(3)
	DENS(I2(i),I1(i),I3(i)) = RHO(i);
	y(I1(i)) = -LIMITS(1) + 2*LIMITS(1)*(I1(i)-1)/(MESH(1)-1);
	x(I2(i)) = -LIMITS(2) + 2*LIMITS(2)*(I2(i)-1)/(MESH(2)-1);
	z(I3(i)) = -LIMITS(3) + 2*LIMITS(3)*(I3(i)-1)/(MESH(3)-1);
	n1(I1(i)) = I1(i);
	n2(I2(i)) = I2(i);
	n3(I3(i)) = I3(i);
end

for i=1:MESH(1)
	for j=1:MESH(2)
        	fdens(i,j) = DENS(i,j,200);
	end
end

[f,v] = isosurface(x,y,z,DENS,rho0);
p = patch('Faces',f,'Vertices',v)
isonormals(x,y,z,DENS,p);
set(p,'FaceColor','cyan','EdgeColor','none','NormalMode','auto');
alpha(0.50);
camlight;set(gca,'Xcolor','w');set(gca,'Ycolor','w');set(gca,'Zcolor','w')
lighting gouraud;
set(gca,'Xcolor','w');set(gca,'Ycolor','w');set(gca,'Zcolor','w')
hold;

[f1,v1] = isosurface(x,y,z,DENS,rho1);
p1 = patch('Faces',f1,'Vertices',v1)
isonormals(x,y,z,DENS,p1);
set(p1,'FaceColor','g','EdgeColor','none','NormalMode','auto');
alpha(0.50);

[f2,v2] = isosurface(x,y,z,DENS,rho2);
p2 = patch('Faces',f2,'Vertices',v2)
isonormals(x,y,z,DENS,p2);
set(p2,'FaceColor','r','EdgeColor','none','NormalMode','auto');
alpha(0.50);

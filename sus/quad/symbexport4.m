qm = [...
kcn+cn^2*kn+(g*mn3)/(2*cn*ln)-(cn*g*mn3)/(2*ln) 0 0 0 0 0
 0 kcn+cn^2*kn+(g*mn3)/(2*cn*ln)-(cn*g*mn3)/(2*ln) 0 0 0 0
 0 0 c1^2*k1+kc1+(g*m13)/(2*c1*l1)-(c1*g*m13)/(2*l1) 0 0 0
 0 0 0 c1^2*k1+kc1+(g*m13)/(2*c1*l1)-(c1*g*m13)/(2*l1) 0 0
 0 0 0 0 c2^2*k2+kc2+(g*m23)/(2*c2*l2)-(c2*g*m23)/(2*l2) 0
 0 0 0 0 0 c2^2*k2+kc2+(g*m23)/(2*c2*l2)-(c2*g*m23)/(2*l2)
];
cqsm = [...
0 -((2*cn*kn*ln-g*mn3)*(nn0-nn1))/(2*ln^2) (2*cn^3*kn*ln+g*mn3-cn^2*g*mn3)/(2*cn*ln) 0 0 -(2*cn^3*kn*ln*nn0-cn^2*g*mn3*nn0+g*mn3*nn1)/(2*cn*ln)
 0 ((2*cn*kn*ln-g*mn3)*(nn0-nn1))/(2*ln^2) (2*cn^3*kn*ln+g*mn3-cn^2*g*mn3)/(2*cn*ln) 0 0 (2*cn^3*kn*ln*nn0-cn^2*g*mn3*nn0+g*mn3*nn1)/(2*cn*ln)
 0 0 0 0 0 0
 0 0 0 0 0 0
 0 0 0 0 0 0
 0 0 0 0 0 0
];
xmlp = [...
(g*(2*dm^2+2*cn*dm*ln)*mn3)/(2*cn*ln)+2*c1^2*k1*su^2-(c1*g*m13*su^2)/l1+(g*m13*(2*dn^2+2*c1*dn*l1+2*su^2))/(2*c1*l1) -((dn*g*m13)/(c1*l1))+(dm*g*mn3)/(cn*ln) -2*c1^2*k1*su^2+(c1*g*m13*su^2)/l1+(g*m13*(2*d0*dn-2*su^2))/(2*c1*l1) (dn*g*m13)/(c1*l1) 0 0 0 0
 -((dn*g*m13)/(c1*l1))+(dm*g*mn3)/(cn*ln) (g*m13)/(c1*l1)+(g*mn3)/(cn*ln) -((d0*g*m13)/(c1*l1)) -((g*m13)/(c1*l1)) 0 0 0 0
 -2*c1^2*k1*su^2+(c1*g*m13*su^2)/l1+(g*m13*(2*d0*dn-2*su^2))/(2*c1*l1) -((d0*g*m13)/(c1*l1)) 2*c2^2*k2*si^2-(c2*g*m23*si^2)/l2+(g*m23*(2*d1^2+2*c2*d1*l2+2*si^2))/(2*c2*l2)+2*c1^2*k1*su^2-(c1*g*m13*su^2)/l1+(g*m13*(2*d0^2+2*c1*d0*l1+2*su^2))/(2*c1*l1) (d0*g*m13)/(c1*l1)-(d1*g*m23)/(c2*l2) -2*c2^2*k2*si^2+(c2*g*m23*si^2)/l2+(g*m23*(2*d1*d2-2*si^2))/(2*c2*l2) (d1*g*m23)/(c2*l2) 0 0
 (dn*g*m13)/(c1*l1) -((g*m13)/(c1*l1)) (d0*g*m13)/(c1*l1)-(d1*g*m23)/(c2*l2) (g*m13)/(c1*l1)+(g*m23)/(c2*l2) -((d2*g*m23)/(c2*l2)) -((g*m23)/(c2*l2)) 0 0
 0 0 -2*c2^2*k2*si^2+(c2*g*m23*si^2)/l2+(g*m23*(2*d1*d2-2*si^2))/(2*c2*l2) -((d2*g*m23)/(c2*l2)) 2*c2^2*k2*si^2-(c2*g*m23*si^2)/l2+(g*m23*(2*d2^2+2*c2*d2*l2+2*si^2))/(2*c2*l2)+2*c3^2*k3*sl^2-(c3*g*m3*sl^2)/l3+(g*m3*(2*d3^2+2*c3*d3*l3+2*sl^2))/(2*c3*l3) (d2*g*m23)/(c2*l2)-(d3*g*m3)/(c3*l3) -2*c3^2*k3*sl^2+(c3*g*m3*sl^2)/l3+(g*m3*(2*d3*d4-2*sl^2))/(2*c3*l3) (d3*g*m3)/(c3*l3)
 0 0 (d1*g*m23)/(c2*l2) -((g*m23)/(c2*l2)) (d2*g*m23)/(c2*l2)-(d3*g*m3)/(c3*l3) (g*m23)/(c2*l2)+(g*m3)/(c3*l3) -((d4*g*m3)/(c3*l3)) -((g*m3)/(c3*l3))
 0 0 0 0 -2*c3^2*k3*sl^2+(c3*g*m3*sl^2)/l3+(g*m3*(2*d3*d4-2*sl^2))/(2*c3*l3) -((d4*g*m3)/(c3*l3)) 2*c3^2*k3*sl^2-(c3*g*m3*sl^2)/l3+(g*m3*(2*d4^2+2*c3*d4*l3+2*sl^2))/(2*c3*l3) (d4*g*m3)/(c3*l3)
 0 0 0 0 (d3*g*m3)/(c3*l3) -((g*m3)/(c3*l3)) (d4*g*m3)/(c3*l3) (g*m3)/(c3*l3)
];
kmlp = [...
Iny 0 0 0 0 0 0 0
 0 mn 0 0 0 0 0 0
 0 0 I1y 0 0 0 0 0
 0 0 0 m1 0 0 0 0
 0 0 0 0 I2y 0 0 0
 0 0 0 0 0 m2 0 0
 0 0 0 0 0 0 I3y 0
 0 0 0 0 0 0 0 m3
];
cqxmlp = [...
0   0   0   0   0   0   0   0;...
0   0   0   0   0   0   0   0;...
0   0   0   0   0   0   0   0;...
0   0   0   0   0   0   0   0;...
0   0   0   0   0   0   0   0;...
0   0   0   0   0   0   0   0;...
];
cxsmlp = [...
-((dm*g*mn3)/(cn*ln)) 0 0 0 0 0
 -((g*mn3)/(cn*ln)) 0 0 0 0 0
 0 0 0 0 0 0
 0 0 0 0 0 0
 0 0 0 0 0 0
 0 0 0 0 0 0
 0 0 0 0 0 0
 0 0 0 0 0 0
];
xmtr = [...
(2*c1*cn^3*kn*l1^3*ln^3*nn1^2-c1*cn^2*l1^3*ln^2*nn1*(4*dm*kn*(nn0-nn1)+g*mn3*nn1)+c1*g*l1^3*mn3*(dm^2*(ln^2-(nn0-nn1)^2)+ln^2*nn0*nn1)+cn*ln*(dn^2*ln^2*(g*m13*(l1^2-(n0-n1)^2)+2*c1*k1*l1*(n0-n1)^2)+c1*dn*l1*ln^2*(4*c1*k1*l1*n0*(n0-n1)+g*m13*(l1^2+2*n0*(-n0+n1)))+l1^2*(2*c1^3*k1*l1*ln^2*n0^2-c1^2*g*ln^2*m13*n0^2+g*ln^2*m13*n0*n1+c1*dm*l1*(2*dm*kn*(nn0-nn1)^2+g*mn3*(ln^2+2*(nn0-nn1)*nn1)))))/(c1*cn*l1^3*ln^3) (c1*dm*g*l1^3*mn3*(-ln^2+(nn0-nn1)^2)+2*c1*cn^2*kn*l1^3*ln^2*(nn0-nn1)*nn1+cn*ln*(dn*ln^2*(g*m13*(l1^2-(n0-n1)^2)+2*c1*k1*l1*(n0-n1)^2)+c1*l1*(2*c1*k1*l1*ln^2*n0*(n0-n1)-2*dm*kn*l1^2*(nn0-nn1)^2+g*(ln^2*m13*n0*(-n0+n1)+l1^2*mn3*nn1*(-nn0+nn1)))))/(c1*cn*l1^3*ln^3) (d0*(dn*(g*m13*(l1^2-(n0-n1)^2)+2*c1*k1*l1*(n0-n1)^2)+c1*l1*(2*c1*k1*l1-g*m13)*n0*(n0-n1))+l1*n1*(-2*c1^3*k1*l1^2*n0-g*l1*m13*n0+c1*dn*g*m13*(n0-n1)+c1^2*l1*(g*m13*n0+2*dn*k1*(-n0+n1))))/(c1*l1^3) (dn*(g*m13*(-l1^2+(n0-n1)^2)-2*c1*k1*l1*(n0-n1)^2)-c1*l1*(2*c1*k1*l1-g*m13)*n0*(n0-n1))/(c1*l1^3) 0 0 0 0
 (c1*dm*g*l1^3*mn3*(-ln^2+(nn0-nn1)^2)+2*c1*cn^2*kn*l1^3*ln^2*(nn0-nn1)*nn1+cn*ln*(dn*ln^2*(g*m13*(l1^2-(n0-n1)^2)+2*c1*k1*l1*(n0-n1)^2)+c1*l1*(2*c1*k1*l1*ln^2*n0*(n0-n1)-2*dm*kn*l1^2*(nn0-nn1)^2+g*(ln^2*m13*n0*(-n0+n1)+l1^2*mn3*nn1*(-nn0+nn1)))))/(c1*cn*l1^3*ln^3) (cn*ln*(g*ln^2*m13*(l1^2-(n0-n1)^2)+2*c1*l1*(k1*ln^2*(n0-n1)^2+kn*l1^2*(nn0-nn1)^2))+c1*g*l1^3*mn3*(ln^2-(nn0-nn1)^2))/(c1*cn*l1^3*ln^3) (d0*(g*m13*(l1^2-(n0-n1)^2)+2*c1*k1*l1*(n0-n1)^2)-c1*l1*(2*c1*k1*l1-g*m13)*(n0-n1)*n1)/(c1*l1^3) (g*m13*(-l1^2+(n0-n1)^2)-2*c1*k1*l1*(n0-n1)^2)/(c1*l1^3) 0 0 0 0
 (d0*(dn*(g*m13*(l1^2-(n0-n1)^2)+2*c1*k1*l1*(n0-n1)^2)+c1*l1*(2*c1*k1*l1-g*m13)*n0*(n0-n1))+l1*n1*(-2*c1^3*k1*l1^2*n0-g*l1*m13*n0+c1*dn*g*m13*(n0-n1)+c1^2*l1*(g*m13*n0+2*dn*k1*(-n0+n1))))/(c1*l1^3) (d0*(g*m13*(l1^2-(n0-n1)^2)+2*c1*k1*l1*(n0-n1)^2)-c1*l1*(2*c1*k1*l1-g*m13)*(n0-n1)*n1)/(c1*l1^3) (2*c1*c2^3*k2*l1^3*l2^3*n2^2+c1*c2^2*l1^3*l2^2*n2*(-(g*m23*n2)+4*d1*k2*(n2-n3))+c1*g*l1^3*m23*(d1^2*(l2^2-(n2-n3)^2)+l2^2*n2*n3)+c2*l2*(d0^2*l2^2*(g*m13*(l1^2-(n0-n1)^2)+2*c1*k1*l1*(n0-n1)^2)+c1*d0*l1*l2^2*(4*c1*k1*l1*n1*(-n0+n1)+g*m13*(l1^2+2*(n0-n1)*n1))+l1^2*(g*l2^2*m13*n0*n1+2*c1^3*k1*l1*l2^2*n1^2-c1^2*g*l2^2*m13*n1^2+c1*d1*l1*(2*d1*k2*(n2-n3)^2+g*m23*(l2^2+2*n2*(-n2+n3))))))/(c1*c2*l1^3*l2^3) (c1*d1*g*l1^3*m23*(l2^2-(n2-n3)^2)+2*c1*c2^2*k2*l1^3*l2^2*n2*(n2-n3)+c2*l2*(d0*l2^2*(g*m13*(-l1^2+(n0-n1)^2)-2*c1*k1*l1*(n0-n1)^2)+c1*l1*(2*c1*k1*l1*l2^2*(n0-n1)*n1+2*d1*k2*l1^2*(n2-n3)^2+g*(l2^2*m13*n1*(-n0+n1)+l1^2*m23*n2*(-n2+n3)))))/(c1*c2*l1^3*l2^3) (d1*(d2*(g*m23*(l2^2-(n2-n3)^2)+2*c2*k2*l2*(n2-n3)^2)-c2*l2*(2*c2*k2*l2-g*m23)*(n2-n3)*n3)+l2*n2*(-2*c2^3*k2*l2^2*n3-g*l2*m23*n3+c2*d2*g*m23*(-n2+n3)+c2^2*l2*(2*d2*k2*(n2-n3)+g*m23*n3)))/(c2*l2^3) (d1*(g*m23*(-l2^2+(n2-n3)^2)-2*c2*k2*l2*(n2-n3)^2)-c2*l2*(2*c2*k2*l2-g*m23)*n2*(n2-n3))/(c2*l2^3) 0 0
 (dn*(g*m13*(-l1^2+(n0-n1)^2)-2*c1*k1*l1*(n0-n1)^2)-c1*l1*(2*c1*k1*l1-g*m13)*n0*(n0-n1))/(c1*l1^3) (g*m13*(-l1^2+(n0-n1)^2)-2*c1*k1*l1*(n0-n1)^2)/(c1*l1^3) (c1*d1*g*l1^3*m23*(l2^2-(n2-n3)^2)+2*c1*c2^2*k2*l1^3*l2^2*n2*(n2-n3)+c2*l2*(d0*l2^2*(g*m13*(-l1^2+(n0-n1)^2)-2*c1*k1*l1*(n0-n1)^2)+c1*l1*(2*c1*k1*l1*l2^2*(n0-n1)*n1+2*d1*k2*l1^2*(n2-n3)^2+g*(l2^2*m13*n1*(-n0+n1)+l1^2*m23*n2*(-n2+n3)))))/(c1*c2*l1^3*l2^3) (c2*l2*(g*l2^2*m13*(l1^2-(n0-n1)^2)+2*c1*l1*(k1*l2^2*(n0-n1)^2+k2*l1^2*(n2-n3)^2))+c1*g*l1^3*m23*(l2^2-(n2-n3)^2))/(c1*c2*l1^3*l2^3) (d2*(g*m23*(l2^2-(n2-n3)^2)+2*c2*k2*l2*(n2-n3)^2)-c2*l2*(2*c2*k2*l2-g*m23)*(n2-n3)*n3)/(c2*l2^3) (g*m23*(-l2^2+(n2-n3)^2)-2*c2*k2*l2*(n2-n3)^2)/(c2*l2^3) 0 0
 0 0 (d1*(d2*(g*m23*(l2^2-(n2-n3)^2)+2*c2*k2*l2*(n2-n3)^2)-c2*l2*(2*c2*k2*l2-g*m23)*(n2-n3)*n3)+l2*n2*(-2*c2^3*k2*l2^2*n3-g*l2*m23*n3+c2*d2*g*m23*(-n2+n3)+c2^2*l2*(2*d2*k2*(n2-n3)+g*m23*n3)))/(c2*l2^3) (d2*(g*m23*(l2^2-(n2-n3)^2)+2*c2*k2*l2*(n2-n3)^2)-c2*l2*(2*c2*k2*l2-g*m23)*(n2-n3)*n3)/(c2*l2^3) (2*c2*c3^3*k3*l2^3*l3^3*n4^2+c2*c3^2*l2^3*l3^2*n4*(-(g*m3*n4)+4*d3*k3*(n4-n5))+c2*g*l2^3*m3*(d3^2*(l3^2-(n4-n5)^2)+l3^2*n4*n5)+c3*l3*(d2^2*l3^2*(g*m23*(l2^2-(n2-n3)^2)+2*c2*k2*l2*(n2-n3)^2)+c2*d2*l2*l3^2*(4*c2*k2*l2*n3*(-n2+n3)+g*m23*(l2^2+2*(n2-n3)*n3))+l2^2*(g*l3^2*m23*n2*n3+2*c2^3*k2*l2*l3^2*n3^2-c2^2*g*l3^2*m23*n3^2+c2*d3*l2*(2*d3*k3*(n4-n5)^2+g*m3*(l3^2+2*n4*(-n4+n5))))))/(c2*c3*l2^3*l3^3) (c2*d3*g*l2^3*m3*(l3^2-(n4-n5)^2)+2*c2*c3^2*k3*l2^3*l3^2*n4*(n4-n5)+c3*l3*(d2*l3^2*(g*m23*(-l2^2+(n2-n3)^2)-2*c2*k2*l2*(n2-n3)^2)+c2*l2*(2*c2*k2*l2*l3^2*(n2-n3)*n3+2*d3*k3*l2^2*(n4-n5)^2+g*(l3^2*m23*n3*(-n2+n3)+l2^2*m3*n4*(-n4+n5)))))/(c2*c3*l2^3*l3^3) (d3*(d4*(g*m3*(l3^2-(n4-n5)^2)+2*c3*k3*l3*(n4-n5)^2)-c3*l3*(2*c3*k3*l3-g*m3)*(n4-n5)*n5)+l3*n4*(-2*c3^3*k3*l3^2*n5-g*l3*m3*n5+c3*d4*g*m3*(-n4+n5)+c3^2*l3*(2*d4*k3*(n4-n5)+g*m3*n5)))/(c3*l3^3) (d3*(g*m3*(-l3^2+(n4-n5)^2)-2*c3*k3*l3*(n4-n5)^2)-c3*l3*(2*c3*k3*l3-g*m3)*n4*(n4-n5))/(c3*l3^3)
 0 0 (d1*(g*m23*(-l2^2+(n2-n3)^2)-2*c2*k2*l2*(n2-n3)^2)-c2*l2*(2*c2*k2*l2-g*m23)*n2*(n2-n3))/(c2*l2^3) (g*m23*(-l2^2+(n2-n3)^2)-2*c2*k2*l2*(n2-n3)^2)/(c2*l2^3) (c2*d3*g*l2^3*m3*(l3^2-(n4-n5)^2)+2*c2*c3^2*k3*l2^3*l3^2*n4*(n4-n5)+c3*l3*(d2*l3^2*(g*m23*(-l2^2+(n2-n3)^2)-2*c2*k2*l2*(n2-n3)^2)+c2*l2*(2*c2*k2*l2*l3^2*(n2-n3)*n3+2*d3*k3*l2^2*(n4-n5)^2+g*(l3^2*m23*n3*(-n2+n3)+l2^2*m3*n4*(-n4+n5)))))/(c2*c3*l2^3*l3^3) (c3*l3*(g*l3^2*m23*(l2^2-(n2-n3)^2)+2*c2*l2*(k2*l3^2*(n2-n3)^2+k3*l2^2*(n4-n5)^2))+c2*g*l2^3*m3*(l3^2-(n4-n5)^2))/(c2*c3*l2^3*l3^3) (d4*(g*m3*(l3^2-(n4-n5)^2)+2*c3*k3*l3*(n4-n5)^2)-c3*l3*(2*c3*k3*l3-g*m3)*(n4-n5)*n5)/(c3*l3^3) (g*m3*(-l3^2+(n4-n5)^2)-2*c3*k3*l3*(n4-n5)^2)/(c3*l3^3)
 0 0 0 0 (d3*(d4*(g*m3*(l3^2-(n4-n5)^2)+2*c3*k3*l3*(n4-n5)^2)-c3*l3*(2*c3*k3*l3-g*m3)*(n4-n5)*n5)+l3*n4*(-2*c3^3*k3*l3^2*n5-g*l3*m3*n5+c3*d4*g*m3*(-n4+n5)+c3^2*l3*(2*d4*k3*(n4-n5)+g*m3*n5)))/(c3*l3^3) (d4*(g*m3*(l3^2-(n4-n5)^2)+2*c3*k3*l3*(n4-n5)^2)-c3*l3*(2*c3*k3*l3-g*m3)*(n4-n5)*n5)/(c3*l3^3) (d4^2*(g*m3*(l3^2-(n4-n5)^2)+2*c3*k3*l3*(n4-n5)^2)+l3^2*n5*(2*c3^3*k3*l3*n5+g*m3*(n4-c3^2*n5))+c3*d4*l3*(4*c3*k3*l3*n5*(-n4+n5)+g*m3*(l3^2+2*(n4-n5)*n5)))/(c3*l3^3) (d4*(g*m3*(-l3^2+(n4-n5)^2)-2*c3*k3*l3*(n4-n5)^2)+c3*l3*(2*c3*k3*l3-g*m3)*(n4-n5)*n5)/(c3*l3^3)
 0 0 0 0 (d3*(g*m3*(-l3^2+(n4-n5)^2)-2*c3*k3*l3*(n4-n5)^2)-c3*l3*(2*c3*k3*l3-g*m3)*n4*(n4-n5))/(c3*l3^3) (g*m3*(-l3^2+(n4-n5)^2)-2*c3*k3*l3*(n4-n5)^2)/(c3*l3^3) (d4*(g*m3*(-l3^2+(n4-n5)^2)-2*c3*k3*l3*(n4-n5)^2)+c3*l3*(2*c3*k3*l3-g*m3)*(n4-n5)*n5)/(c3*l3^3) (g*m3*(l3^2-(n4-n5)^2)+2*c3*k3*l3*(n4-n5)^2)/(c3*l3^3)
];
kmtr = [...
Inx 0 0 0 0 0 0 0
 0 mn 0 0 0 0 0 0
 0 0 I1x 0 0 0 0 0
 0 0 0 m1 0 0 0 0
 0 0 0 0 I2x 0 0 0
 0 0 0 0 0 m2 0 0
 0 0 0 0 0 0 I3x 0
 0 0 0 0 0 0 0 m3
];
cqxmtr = [...
(cn*dm*g*mn3*(nn0-nn1)+2*cn^3*kn*ln^2*nn1+g*ln*mn3*nn1-cn^2*ln*(2*dm*kn*(nn0-nn1)+g*mn3*nn1))/(2*cn*ln^2) ((2*cn*kn*ln-g*mn3)*(nn0-nn1))/(2*ln^2) 0 0 0 0 0 0
 (-2*cn^3*kn*ln^2*nn1-g*ln*mn3*nn1+cn*dm*g*mn3*(-nn0+nn1)+cn^2*ln*(2*dm*kn*(nn0-nn1)+g*mn3*nn1))/(2*cn*ln^2) -((2*cn*kn*ln-g*mn3)*(nn0-nn1))/(2*ln^2) 0 0 0 0 0 0
 (-2*c1^3*k1*l1^2*n0+c1*dn*g*m13*(n0-n1)-g*l1*m13*n1+c1^2*l1*(g*m13*n0+2*dn*k1*(-n0+n1)))/(2*c1*l1^2) -((2*c1*k1*l1-g*m13)*(n0-n1))/(2*l1^2) (c1*d0*g*m13*(n0-n1)+2*c1^3*k1*l1^2*n1+g*l1*m13*n1-c1^2*l1*(2*d0*k1*(n0-n1)+g*m13*n1))/(2*c1*l1^2) ((2*c1*k1*l1-g*m13)*(n0-n1))/(2*l1^2) 0 0 0 0
 (2*c1^3*k1*l1^2*n0+c1^2*l1*(-(g*m13*n0)+2*dn*k1*(n0-n1))+g*l1*m13*n1+c1*dn*g*m13*(-n0+n1))/(2*c1*l1^2) ((2*c1*k1*l1-g*m13)*(n0-n1))/(2*l1^2) (-2*c1^3*k1*l1^2*n1-g*l1*m13*n1+c1*d0*g*m13*(-n0+n1)+c1^2*l1*(2*d0*k1*(n0-n1)+g*m13*n1))/(2*c1*l1^2) -((2*c1*k1*l1-g*m13)*(n0-n1))/(2*l1^2) 0 0 0 0
 0 0 (-2*c2^3*k2*l2^2*n2+c2*d1*g*m23*(n2-n3)-g*l2*m23*n3+c2^2*l2*(g*m23*n2+2*d1*k2*(-n2+n3)))/(2*c2*l2^2) -((2*c2*k2*l2-g*m23)*(n2-n3))/(2*l2^2) (c2*d2*g*m23*(n2-n3)+2*c2^3*k2*l2^2*n3+g*l2*m23*n3-c2^2*l2*(2*d2*k2*(n2-n3)+g*m23*n3))/(2*c2*l2^2) ((2*c2*k2*l2-g*m23)*(n2-n3))/(2*l2^2) 0 0
 0 0 (2*c2^3*k2*l2^2*n2+c2^2*l2*(-(g*m23*n2)+2*d1*k2*(n2-n3))+g*l2*m23*n3+c2*d1*g*m23*(-n2+n3))/(2*c2*l2^2) ((2*c2*k2*l2-g*m23)*(n2-n3))/(2*l2^2) (-2*c2^3*k2*l2^2*n3-g*l2*m23*n3+c2*d2*g*m23*(-n2+n3)+c2^2*l2*(2*d2*k2*(n2-n3)+g*m23*n3))/(2*c2*l2^2) -((2*c2*k2*l2-g*m23)*(n2-n3))/(2*l2^2) 0 0
];
cxsmtr = [...
0 (dm*(g*mn3*(ln^2-(nn0-nn1)^2)+2*cn*kn*ln*(nn0-nn1)^2)-cn*ln*(2*cn*kn*ln-g*mn3)*(nn0-nn1)*nn1)/(cn*ln^3) 0 0 0 (nn0*(-2*cn^3*kn*ln^2*nn1-g*ln*mn3*nn1+cn*dm*g*mn3*(-nn0+nn1)+cn^2*ln*(2*dm*kn*(nn0-nn1)+g*mn3*nn1)))/(cn*ln^2)
 0 (g*mn3*(-ln^2+(nn0-nn1)^2)-2*cn*kn*ln*(nn0-nn1)^2)/(cn*ln^3) 0 0 0 -(((2*cn*kn*ln-g*mn3)*nn0*(nn0-nn1))/ln^2)
 0 0 0 0 0 0
 0 0 0 0 0 0
 0 0 0 0 0 0
 0 0 0 0 0 0
 0 0 0 0 0 0
 0 0 0 0 0 0
];
xmy = [...
(g*mn3*(2*(nn0-nn1)*nn1+2*nn1^2))/(4*cn*ln)+(g*mn3*(2*nn1^2-2*nn1*(-nn0+nn1)))/(4*cn*ln)+(k1*(n0-n1)^2*su^2)/l1^2-(g*m13*(n0-n1)^2*su^2)/(2*c1*l1^3)+(k1*(-n0+n1)^2*su^2)/l1^2-(g*m13*(-n0+n1)^2*su^2)/(2*c1*l1^3)+(g*m13*(2*n0^2-2*n0*(n0-n1)+2*su^2))/(4*c1*l1)+(g*m13*(2*n0^2+2*n0*(-n0+n1)+2*su^2))/(4*c1*l1) -((k1*(n0-n1)^2*su^2)/l1^2)+(g*m13*(n0-n1)^2*su^2)/(2*c1*l1^3)-(k1*(-n0+n1)^2*su^2)/l1^2+(g*m13*(-n0+n1)^2*su^2)/(2*c1*l1^3)+(g*m13*(-2*n0*n1-2*su^2))/(2*c1*l1) 0 0
 -((k1*(n0-n1)^2*su^2)/l1^2)+(g*m13*(n0-n1)^2*su^2)/(2*c1*l1^3)-(k1*(-n0+n1)^2*su^2)/l1^2+(g*m13*(-n0+n1)^2*su^2)/(2*c1*l1^3)+(g*m13*(-2*n0*n1-2*su^2))/(2*c1*l1) (k2*(n2-n3)^2*si^2)/l2^2-(g*m23*(n2-n3)^2*si^2)/(2*c2*l2^3)+(k2*(-n2+n3)^2*si^2)/l2^2-(g*m23*(-n2+n3)^2*si^2)/(2*c2*l2^3)+(g*m23*(2*n2^2-2*n2*(n2-n3)+2*si^2))/(4*c2*l2)+(g*m23*(2*n2^2+2*n2*(-n2+n3)+2*si^2))/(4*c2*l2)+(k1*(n0-n1)^2*su^2)/l1^2-(g*m13*(n0-n1)^2*su^2)/(2*c1*l1^3)+(k1*(-n0+n1)^2*su^2)/l1^2-(g*m13*(-n0+n1)^2*su^2)/(2*c1*l1^3)+(g*m13*(2*(n0-n1)*n1+2*n1^2+2*su^2))/(4*c1*l1)+(g*m13*(2*n1^2-2*n1*(-n0+n1)+2*su^2))/(4*c1*l1) -((k2*(n2-n3)^2*si^2)/l2^2)+(g*m23*(n2-n3)^2*si^2)/(2*c2*l2^3)-(k2*(-n2+n3)^2*si^2)/l2^2+(g*m23*(-n2+n3)^2*si^2)/(2*c2*l2^3)+(g*m23*(-2*n2*n3-2*si^2))/(2*c2*l2) 0
 0 -((k2*(n2-n3)^2*si^2)/l2^2)+(g*m23*(n2-n3)^2*si^2)/(2*c2*l2^3)-(k2*(-n2+n3)^2*si^2)/l2^2+(g*m23*(-n2+n3)^2*si^2)/(2*c2*l2^3)+(g*m23*(-2*n2*n3-2*si^2))/(2*c2*l2) (k2*(n2-n3)^2*si^2)/l2^2-(g*m23*(n2-n3)^2*si^2)/(2*c2*l2^3)+(k2*(-n2+n3)^2*si^2)/l2^2-(g*m23*(-n2+n3)^2*si^2)/(2*c2*l2^3)+(g*m23*(2*(n2-n3)*n3+2*n3^2+2*si^2))/(4*c2*l2)+(g*m23*(2*n3^2-2*n3*(-n2+n3)+2*si^2))/(4*c2*l2)+(k3*(n4-n5)^2*sl^2)/l3^2-(g*m3*(n4-n5)^2*sl^2)/(2*c3*l3^3)+(k3*(-n4+n5)^2*sl^2)/l3^2-(g*m3*(-n4+n5)^2*sl^2)/(2*c3*l3^3)+(g*m3*(2*n4^2-2*n4*(n4-n5)+2*sl^2))/(4*c3*l3)+(g*m3*(2*n4^2+2*n4*(-n4+n5)+2*sl^2))/(4*c3*l3) -((k3*(n4-n5)^2*sl^2)/l3^2)+(g*m3*(n4-n5)^2*sl^2)/(2*c3*l3^3)-(k3*(-n4+n5)^2*sl^2)/l3^2+(g*m3*(-n4+n5)^2*sl^2)/(2*c3*l3^3)+(g*m3*(-2*n4*n5-2*sl^2))/(2*c3*l3)
 0 0 -((k3*(n4-n5)^2*sl^2)/l3^2)+(g*m3*(n4-n5)^2*sl^2)/(2*c3*l3^3)-(k3*(-n4+n5)^2*sl^2)/l3^2+(g*m3*(-n4+n5)^2*sl^2)/(2*c3*l3^3)+(g*m3*(-2*n4*n5-2*sl^2))/(2*c3*l3) (k3*(n4-n5)^2*sl^2)/l3^2-(g*m3*(n4-n5)^2*sl^2)/(2*c3*l3^3)+(k3*(-n4+n5)^2*sl^2)/l3^2-(g*m3*(-n4+n5)^2*sl^2)/(2*c3*l3^3)+(g*m3*(2*(n4-n5)*n5+2*n5^2+2*sl^2))/(4*c3*l3)+(g*m3*(2*n5^2-2*n5*(-n4+n5)+2*sl^2))/(4*c3*l3)
];
kmy = [...
Inz 0 0 0
 0 I1z 0 0
 0 0 I2z 0
 0 0 0 I3z
];
cqxmy = [...
0   0   0   0;...
0   0   0   0;...
0   0   0   0;...
0   0   0   0;...
0   0   0   0;...
0   0   0   0;...
];
cxsmy = [...
0 0 0 -((g*mn3*nn0*nn1)/(cn*ln)) 0 0
 0 0 0 0 0 0
 0 0 0 0 0 0
 0 0 0 0 0 0
];
xmv = [...
2*c1^2*k1+2*cn^2*kn+(g*m13)/(c1*l1)-(c1*g*m13)/l1+(g*mn3)/(cn*ln)-(cn*g*mn3)/ln -2*c1^2*k1-(g*m13)/(c1*l1)+(c1*g*m13)/l1 0 0
 -2*c1^2*k1-(g*m13)/(c1*l1)+(c1*g*m13)/l1 2*c1^2*k1+2*c2^2*k2+(g*m13)/(c1*l1)-(c1*g*m13)/l1+(g*m23)/(c2*l2)-(c2*g*m23)/l2 -2*c2^2*k2-(g*m23)/(c2*l2)+(c2*g*m23)/l2 0
 0 -2*c2^2*k2-(g*m23)/(c2*l2)+(c2*g*m23)/l2 2*c2^2*k2+2*c3^2*k3+(g*m23)/(c2*l2)-(c2*g*m23)/l2+(g*m3)/(c3*l3)-(c3*g*m3)/l3 -2*c3^2*k3-(g*m3)/(c3*l3)+(c3*g*m3)/l3
 0 0 -2*c3^2*k3-(g*m3)/(c3*l3)+(c3*g*m3)/l3 2*c3^2*k3+(g*m3)/(c3*l3)-(c3*g*m3)/l3
];
kmv = [...
mn 0 0 0
 0 m1 0 0
 0 0 m2 0
 0 0 0 m3
];
cqxmv = [...
-(cn^2*kn)-(g*mn3)/(2*cn*ln)+(cn*g*mn3)/(2*ln) 0 0 0
 -(cn^2*kn)-(g*mn3)/(2*cn*ln)+(cn*g*mn3)/(2*ln) 0 0 0
 c1^2*k1+(g*m13)/(2*c1*l1)-(c1*g*m13)/(2*l1) -(c1^2*k1)-(g*m13)/(2*c1*l1)+(c1*g*m13)/(2*l1) 0 0
 c1^2*k1+(g*m13)/(2*c1*l1)-(c1*g*m13)/(2*l1) -(c1^2*k1)-(g*m13)/(2*c1*l1)+(c1*g*m13)/(2*l1) 0 0
 0 c2^2*k2+(g*m23)/(2*c2*l2)-(c2*g*m23)/(2*l2) -(c2^2*k2)-(g*m23)/(2*c2*l2)+(c2*g*m23)/(2*l2) 0
 0 c2^2*k2+(g*m23)/(2*c2*l2)-(c2*g*m23)/(2*l2) -(c2^2*k2)-(g*m23)/(2*c2*l2)+(c2*g*m23)/(2*l2) 0
];
cxsmv = [...
0 0 (-2*cn^3*kn*ln-g*mn3+cn^2*g*mn3)/(cn*ln) 0 0 0
 0 0 0 0 0 0
 0 0 0 0 0 0
 0 0 0 0 0 0
];
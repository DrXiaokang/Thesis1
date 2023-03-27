clear all
close all
setlmis([]);

cf0=95000;cr0=85500;lr=1.67;lf=1.11;vxmin=25;vxmax=25; mmin=1530;mmax=1680;Izmin=4200;Izmax=4600;
tao=0.019;lbd=0.1;gama1=0.37;


A1=[-(cf0+cr0)/(mmin*vxmax) -vxmax+(cr0*lr-cf0*lf)/(mmin*vxmax);(cr0*lr-cf0*lf)/(Izmin*vxmax) -(cf0*lf*lf+cr0*lr*lr)/(Izmin*vxmax)];
A2=[-(cf0+cr0)/(mmin*vxmax) -vxmax+(cr0*lr-cf0*lf)/(mmin*vxmax);(cr0*lr-cf0*lf)/(Izmax*vxmax) -(cf0*lf*lf+cr0*lr*lr)/(Izmax*vxmax)];
A3=[-(cf0+cr0)/(mmax*vxmin) -vxmin+(cr0*lr-cf0*lf)/(mmax*vxmin);(cr0*lr-cf0*lf)/(Izmin*vxmin) -(cf0*lf*lf+cr0*lr*lr)/(Izmin*vxmin)];
A4=[-(cf0+cr0)/(mmax*vxmin) -vxmin+(cr0*lr-cf0*lf)/(mmax*vxmin);(cr0*lr-cf0*lf)/(Izmax*vxmin) -(cf0*lf*lf+cr0*lr*lr)/(Izmax*vxmin)];


B1=[cf0/mmin;(lf*cf0)/Izmin];B2=[cf0/mmin;(lf*cf0)/Izmax];B3=[cf0/mmax;(lf*cf0)/Izmin];B4=[cf0/mmax;(lf*cf0)/Izmax];

E1=[cf0/mmin;(lf*cf0)/Izmin];E2=[cf0/mmin;(lf*cf0)/Izmax];E3=[cf0/mmax;(lf*cf0)/Izmin];E4=[cf0/mmax;(lf*cf0)/Izmax];

C1=[0 1];

D=[1 0;0 1];

Ls1=[49.1209
   31.0100];
Ls2=[54.3644
   30.6819];
Ls3=[38.1701
   30.0822];
Ls4=[48.5045
   31.0602];
 
Ld1=0.9815;
Ld2=1.0689;
Ld3=0.9848;
Ld4=1.0649;

Mz1=lmivar(2,[2,2]);Mz2=lmivar(2,[1,1]);

P10=lmivar(1,[2,1]);P11=lmivar(1,[2,0]);P12=lmivar(1,[2,0]);P13=lmivar(1,[2,0]);P14=lmivar(1,[2,0]);
P20=lmivar(2,[2,2]);P30=lmivar(2,[2,1]);
P40=lmivar(1,[2,1]);P41=lmivar(1,[2,0]);P42=lmivar(1,[2,0]);P43=lmivar(1,[2,0]);P44=lmivar(1,[2,0]);
P50=lmivar(2,[2,1]);
P60=lmivar(1,[1,1]);P61=lmivar(1,[1,0]);P62=lmivar(1,[1,0]);P63=lmivar(1,[1,0]);P64=lmivar(1,[1,0]);

Q10=lmivar(1,[2,1]);Q11=lmivar(1,[2,0]);Q12=lmivar(1,[2,0]);Q13=lmivar(1,[2,0]);Q14=lmivar(1,[2,0]);
Q20=lmivar(2,[2,2]);Q30=lmivar(2,[2,1]);
Q40=lmivar(1,[2,1]);Q41=lmivar(1,[2,0]);Q42=lmivar(1,[2,0]);Q43=lmivar(1,[2,0]);Q44=lmivar(1,[2,0]);
Q50=lmivar(2,[2,1]);
Q60=lmivar(1,[1,1]);Q61=lmivar(1,[1,0]);Q62=lmivar(1,[1,0]);Q63=lmivar(1,[1,0]);Q64=lmivar(1,[1,0]);


R10=lmivar(1,[2,1]);R11=lmivar(1,[2,0]);R12=lmivar(1,[2,0]);R13=lmivar(1,[2,0]);R14=lmivar(1,[2,0]);
R20=lmivar(2,[2,2]);R30=lmivar(2,[2,1]);
R40=lmivar(1,[2,1]);R41=lmivar(1,[2,0]);R42=lmivar(1,[2,0]);R43=lmivar(1,[2,0]);R44=lmivar(1,[2,0]);
R50=lmivar(2,[2,1]);
R60=lmivar(1,[1,1]);R61=lmivar(1,[1,0]);R62=lmivar(1,[1,0]);R63=lmivar(1,[1,0]);R64=lmivar(1,[1,0]);


Zg1=lmivar(2,[1,2]);Zg2=lmivar(2,[1,2]);Zg3=lmivar(2,[1,2]);Zg4=lmivar(2,[1,2]);



lmiterm([1,1,1,Q10],1,1);lmiterm([1,1,1,Q11],1,1);
lmiterm([1,1,1,R10],tao/(lbd*lbd),1);lmiterm([1,1,1,R11],tao/(lbd*lbd),1);
lmiterm([1,1,1,P10],-2/lbd,1);lmiterm([1,1,1,P11],-2/lbd,1);
lmiterm([1,1,2,Q20],1,1);
lmiterm([1,1,2,R20],tao/(lbd*lbd),1);
lmiterm([1,1,2,P20],-2/lbd,1);
lmiterm([1,1,3,Q30],1,1);
lmiterm([1,1,3,R30],tao/(lbd*lbd),1);
lmiterm([1,1,3,P30],-2/lbd,1);
lmiterm([1,2,2,Q40],1,1);lmiterm([1,2,2,Q41],1,1);
lmiterm([1,2,2,R40],tao/(lbd*lbd),1);lmiterm([1,2,2,R41],tao/(lbd*lbd),1);
lmiterm([1,2,2,P40],-2/lbd,1);lmiterm([1,2,2,P41],-2/lbd,1);
lmiterm([1,2,3,Q50],1,1);
lmiterm([1,2,3,R50],tao/(lbd*lbd),1);
lmiterm([1,2,3,P50],-2/lbd,1);
lmiterm([1,3,3,Q60],1,1);lmiterm([1,3,3,Q61],1,1);
lmiterm([1,3,3,R60],tao/(lbd*lbd),1);lmiterm([1,3,3,R61],tao/(lbd*lbd),1);
lmiterm([1,3,3,P60],-2/lbd,1);lmiterm([1,3,3,P61],-2/lbd,1);
lmiterm([1,1,4,P10],1,1);lmiterm([1,1,4,P11],1,1);
lmiterm([1,1,4,R10],-tao/lbd,1);lmiterm([1,1,4,R11],-tao/lbd,1);
lmiterm([1,1,4,-Mz1],1,1);lmiterm([1,1,4,-Mz1],lbd,A1');lmiterm([1,1,4,-Zg1],lbd,B1');
lmiterm([1,1,5,P20],1,1); 
lmiterm([1,1,5,R20],-tao/lbd,1);
lmiterm([1,1,6,P30],1,1);
lmiterm([1,1,6,R30],-tao/lbd,1);
lmiterm([1,2,4,-P20],1,1); 
lmiterm([1,2,4,-R20],-tao/lbd,1); lmiterm([1,2,4,-Zg1],-lbd,B1');
lmiterm([1,2,5,P40],1,1);lmiterm([1,2,5,P41],1,1);
lmiterm([1,2,5,R40],-tao/lbd,1);lmiterm([1,2,5,R41],-tao/lbd,1);
lmiterm([1,2,5,-Mz1],1,1);lmiterm([1,2,5,-Mz1],lbd,A1');lmiterm([1,2,5,-Mz1],-lbd,C1'*Ls1');
lmiterm([1,2,6,P50],1,1);
lmiterm([1,2,6,R50],-tao/lbd,1);lmiterm([1,2,6,-Mz1],-lbd,A1'*C1'*Ld1');
lmiterm([1,3,4,-P30],1,1);
lmiterm([1,3,4,-R30],-tao/lbd,1);lmiterm([1,3,4,-Mz2],lbd,B1');
lmiterm([1,3,5,-P50],1,1);
lmiterm([1,3,5,-R50],-tao/lbd,1);lmiterm([1,3,5,-Mz2],lbd,B1');
lmiterm([1,3,6,P60],1,1);lmiterm([1,3,6,P61],1,1);
lmiterm([1,3,6,R60],-tao/lbd,1);lmiterm([1,3,6,R61],-tao/lbd,1);
lmiterm([1,3,6,-Mz2],1,1);lmiterm([1,3,6,-Mz2],-lbd,B1'*C1'*Ld1');

lmiterm([1,1,7,0],0);
lmiterm([1,1,8,0],0);
lmiterm([1,1,9,0],0);
lmiterm([1,2,7,0],0);
lmiterm([1,2,8,0],0);
lmiterm([1,2,9,0],0);
lmiterm([1,3,7,0],0);
lmiterm([1,3,8,0],0);
lmiterm([1,3,9,0],0);
lmiterm([1,1,10,0],0);
lmiterm([1,2,10,0],0);
lmiterm([1,3,10,0],0);
lmiterm([1,1,11,0],0);
lmiterm([1,1,12,0],0);
lmiterm([1,1,13,0],0);
lmiterm([1,2,11,0],0);
lmiterm([1,2,12,0],0);
lmiterm([1,2,13,0],0);
lmiterm([1,3,11,0],0);
lmiterm([1,3,12,0],0);
lmiterm([1,3,13,0],0);
lmiterm([1,1,14,-Mz1],1,1);
lmiterm([1,1,15,0],0);
lmiterm([1,1,16,0],0);
lmiterm([1,2,14,0],0);
lmiterm([1,2,15,0],0);
lmiterm([1,2,16,0],0);
lmiterm([1,3,14,0],0);
lmiterm([1,3,15,0],0);
lmiterm([1,3,16,0],0);


lmiterm([1,4,4,Mz1],-lbd,1,'s');
lmiterm([1,4,4,R10],tao,1);lmiterm([1,4,4,R11],tao,1);
lmiterm([1,4,5,R20],tao,1);
lmiterm([1,4,6,R30],tao,1);
lmiterm([1,5,5,Mz1],-lbd,1,'s');
lmiterm([1,5,5,R40],tao,1);lmiterm([1,5,5,R41],tao,1);
lmiterm([1,5,6,R50],tao,1);
lmiterm([1,6,6,Mz2],-lbd,1,'s');
lmiterm([1,6,6,R60],tao,1);lmiterm([1,6,6,R61],tao,1);

lmiterm([1,4,7,0],0);
lmiterm([1,4,8,0],0);
lmiterm([1,4,9,0],0);
lmiterm([1,5,7,0],0);
lmiterm([1,5,8,0],0);
lmiterm([1,5,9,0],0);
lmiterm([1,6,7,0],0);
lmiterm([1,6,8,0],0);
lmiterm([1,6,9,0],0);

lmiterm([1,4,10,0],0);
lmiterm([1,5,10,0],0);
lmiterm([1,6,10,0],lbd);
lmiterm([1,4,11,0],0);
lmiterm([1,4,12,0],0);
lmiterm([1,4,13,0],0);
lmiterm([1,5,11,0],0);lmiterm([1,5,12,-Mz1],-lbd,C1'*Ls1');lmiterm([1,5,13,-Mz1],-lbd,A1'*C1'*Ld1');
lmiterm([1,6,11,0],0);
lmiterm([1,6,12,0],0);lmiterm([1,6,13,-Mz2],-lbd,B1'*C1'*Ld1');

lmiterm([1,4,14,0],0);
lmiterm([1,4,15,0],0);
lmiterm([1,4,16,0],0);
lmiterm([1,5,14,0],0);
lmiterm([1,5,15,0],0);
lmiterm([1,5,16,0],0);
lmiterm([1,6,14,0],0);
lmiterm([1,6,15,0],0);
lmiterm([1,6,16,0],0);

lmiterm([1,7,7,Q10],-1,1);lmiterm([1,7,7,Q11],-1,1);
lmiterm([1,7,8,Q20],-1,1);
lmiterm([1,7,9,Q30],-1,1);
lmiterm([1,8,8,Q40],-1,1);lmiterm([1,8,8,Q41],-1,1);
lmiterm([1,8,9,Q50],-1,1);
lmiterm([1,9,9,Q60],-1,1);lmiterm([1,9,9,Q61],-1,1);
lmiterm([1,7,10,0],0);
lmiterm([1,8,10,0],0);
lmiterm([1,9,10,0],0);
lmiterm([1,7,11,0],0);
lmiterm([1,7,12,0],0);
lmiterm([1,7,13,0],0);
lmiterm([1,8,11,0],0);
lmiterm([1,8,12,0],0);
lmiterm([1,8,13,0],0);
lmiterm([1,9,11,0],0);
lmiterm([1,9,12,0],0);
lmiterm([1,9,13,0],0);

lmiterm([1,7,14,0],0);
lmiterm([1,7,15,0],0);
lmiterm([1,7,16,0],0);
lmiterm([1,8,14,0],0);
lmiterm([1,8,15,0],0);
lmiterm([1,8,16,0],0);
lmiterm([1,9,14,0],0);
lmiterm([1,9,15,0],0);
lmiterm([1,9,16,0],0);

lmiterm([1,10,10,0],-gama1*gama1);
lmiterm([1,10,11,0],0);
lmiterm([1,10,12,0],0);
lmiterm([1,10,13,0],0);
lmiterm([1,10,14,0],0);
lmiterm([1,10,15,0],0);
lmiterm([1,10,16,0],0);

lmiterm([1,11,11,R10],-1/tao,1);lmiterm([1,11,11,R11],-1/tao,1);
lmiterm([1,11,12,R20],-1/tao,1);
lmiterm([1,11,13,R30],-1/tao,1);
lmiterm([1,12,12,R40],-1/tao,1);lmiterm([1,12,12,R41],-1/tao,1);
lmiterm([1,12,13,R50],-1/tao,1);
lmiterm([1,13,13,R60],-1/tao,1);lmiterm([1,13,13,R61],-1/tao,1);

lmiterm([1,11,14,0],0);
lmiterm([1,11,15,0],0);
lmiterm([1,11,16,0],0);
lmiterm([1,12,14,0],0);
lmiterm([1,12,15,0],0);
lmiterm([1,12,16,0],0);
lmiterm([1,13,14,0],0);
lmiterm([1,13,15,0],0);
lmiterm([1,13,16,0],0);

lmiterm([1,14,14,0],-1);
lmiterm([1,14,15,0],0);
lmiterm([1,14,16,0],0);
lmiterm([1,15,15,0],-1);
lmiterm([1,15,16,0],0);
lmiterm([1,16,16,0],-1);




lmiterm([2,1,1,Q10],1,1);lmiterm([2,1,1,Q12],1,1);
lmiterm([2,1,1,R10],tao/(lbd*lbd),1);lmiterm([2,1,1,R12],tao/(lbd*lbd),1);
lmiterm([2,1,1,P10],-2/lbd,1);lmiterm([2,1,1,P12],-2/lbd,1);
lmiterm([2,1,2,Q20],1,1);
lmiterm([2,1,2,R20],tao/(lbd*lbd),1);
lmiterm([2,1,2,P20],-2/lbd,1);
lmiterm([2,1,3,Q30],1,1);
lmiterm([2,1,3,R30],tao/(lbd*lbd),1);
lmiterm([2,1,3,P30],-2/lbd,1);
lmiterm([2,2,2,Q40],1,1);lmiterm([2,2,2,Q42],1,1);
lmiterm([2,2,2,R40],tao/(lbd*lbd),1);lmiterm([2,2,2,R42],tao/(lbd*lbd),1);
lmiterm([2,2,2,P40],-2/lbd,1);lmiterm([2,2,2,P42],-2/lbd,1);
lmiterm([2,2,3,Q50],1,1);
lmiterm([2,2,3,R50],tao/(lbd*lbd),1);
lmiterm([2,2,3,P50],-2/lbd,1);
lmiterm([2,3,3,Q60],1,1);lmiterm([2,3,3,Q62],1,1);
lmiterm([2,3,3,R60],tao/(lbd*lbd),1);lmiterm([2,3,3,R62],tao/(lbd*lbd),1);
lmiterm([2,3,3,P60],-2/lbd,1);lmiterm([2,3,3,P62],-2/lbd,1);
lmiterm([2,1,4,P10],1,1);lmiterm([2,1,4,P12],1,1);
lmiterm([2,1,4,R10],-tao/lbd,1);lmiterm([2,1,4,R12],-tao/lbd,1);
lmiterm([2,1,4,-Mz1],1,1);lmiterm([2,1,4,-Mz1],lbd,A2');lmiterm([2,1,4,-Zg2],lbd,B2');
lmiterm([2,1,5,P20],1,1); 
lmiterm([2,1,5,R20],-tao/lbd,1);
lmiterm([2,1,6,P30],1,1);
lmiterm([2,1,6,R30],-tao/lbd,1);
lmiterm([2,2,4,-P20],1,1); 
lmiterm([2,2,4,-R20],-tao/lbd,1); lmiterm([2,2,4,-Zg2],-lbd,B2');
lmiterm([2,2,5,P40],1,1);lmiterm([2,2,5,P42],1,1);
lmiterm([2,2,5,R40],-tao/lbd,1);lmiterm([2,2,5,R42],-tao/lbd,1);
lmiterm([2,2,5,-Mz1],1,1);lmiterm([2,2,5,-Mz1],lbd,A2');lmiterm([2,2,5,-Mz1],-lbd,C1'*Ls2');
lmiterm([2,2,6,P50],1,1);
lmiterm([2,2,6,R50],-tao/lbd,1);lmiterm([2,2,6,-Mz1],-lbd,A2'*C1'*Ld2');
lmiterm([2,3,4,-P30],1,1);
lmiterm([2,3,4,-R30],-tao/lbd,1);lmiterm([2,3,4,-Mz2],lbd,B2');
lmiterm([2,3,5,-P50],1,1);
lmiterm([2,3,5,-R50],-tao/lbd,1);lmiterm([2,3,5,-Mz2],lbd,B2');
lmiterm([2,3,6,P60],1,1);lmiterm([2,3,6,P62],1,1);
lmiterm([2,3,6,R60],-tao/lbd,1);lmiterm([2,3,6,R62],-tao/lbd,1);
lmiterm([2,3,6,-Mz2],1,1);lmiterm([2,3,6,-Mz2],-lbd,B2'*C1'*Ld2');

lmiterm([2,1,7,0],0);
lmiterm([2,1,8,0],0);
lmiterm([2,1,9,0],0);
lmiterm([2,2,7,0],0);
lmiterm([2,2,8,0],0);
lmiterm([2,2,9,0],0);
lmiterm([2,3,7,0],0);
lmiterm([2,3,8,0],0);
lmiterm([2,3,9,0],0);
lmiterm([2,1,10,0],0);
lmiterm([2,2,10,0],0);
lmiterm([2,3,10,0],0);
lmiterm([2,1,11,0],0);
lmiterm([2,1,12,0],0);
lmiterm([2,1,13,0],0);
lmiterm([2,2,11,0],0);
lmiterm([2,2,12,0],0);
lmiterm([2,2,13,0],0);
lmiterm([2,3,11,0],0);
lmiterm([2,3,12,0],0);
lmiterm([2,3,13,0],0);
lmiterm([2,1,14,-Mz1],1,1);
lmiterm([2,1,15,0],0);
lmiterm([2,1,16,0],0);
lmiterm([2,2,14,0],0);
lmiterm([2,2,15,0],0);
lmiterm([2,2,16,0],0);
lmiterm([2,3,14,0],0);
lmiterm([2,3,15,0],0);
lmiterm([2,3,16,0],0);

lmiterm([2,4,4,Mz1],-lbd,1,'s');
lmiterm([2,4,4,R10],tao,1);lmiterm([2,4,4,R12],tao,1);
lmiterm([2,4,5,R20],tao,1);
lmiterm([2,4,6,R30],tao,1);
lmiterm([2,5,5,Mz1],-lbd,1,'s');
lmiterm([2,5,5,R40],tao,1);lmiterm([2,5,5,R42],tao,1);
lmiterm([2,5,6,R50],tao,1);
lmiterm([2,6,6,Mz2],-lbd,1,'s');
lmiterm([2,6,6,R60],tao,1);lmiterm([2,6,6,R62],tao,1);

lmiterm([2,4,7,0],0);
lmiterm([2,4,8,0],0);
lmiterm([2,4,9,0],0);
lmiterm([2,5,7,0],0);
lmiterm([2,5,8,0],0);
lmiterm([2,5,9,0],0);
lmiterm([2,6,7,0],0);
lmiterm([2,6,8,0],0);
lmiterm([2,6,9,0],0);

lmiterm([2,4,10,0],0);
lmiterm([2,5,10,0],0);
lmiterm([2,6,10,0],lbd);
lmiterm([2,4,11,0],0);
lmiterm([2,4,12,0],0);
lmiterm([2,4,13,0],0);
lmiterm([2,5,11,0],0);lmiterm([2,5,12,-Mz1],-lbd,C1'*Ls2');lmiterm([2,5,13,-Mz1],-lbd,A2'*C1'*Ld2');
lmiterm([2,6,11,0],0);
lmiterm([2,6,12,0],0);lmiterm([2,6,13,-Mz2],-lbd,B2'*C1'*Ld2');

lmiterm([2,4,14,0],0);
lmiterm([2,4,15,0],0);
lmiterm([2,4,16,0],0);
lmiterm([2,5,14,0],0);
lmiterm([2,5,15,0],0);
lmiterm([2,5,16,0],0);
lmiterm([2,6,14,0],0);
lmiterm([2,6,15,0],0);
lmiterm([2,6,16,0],0);

lmiterm([2,7,7,Q10],-1,1);lmiterm([2,7,7,Q12],-1,1);
lmiterm([2,7,8,Q20],-1,1);
lmiterm([2,7,9,Q30],-1,1);
lmiterm([2,8,8,Q40],-1,1);lmiterm([2,8,8,Q42],-1,1);
lmiterm([2,8,9,Q50],-1,1);
lmiterm([2,9,9,Q60],-1,1);lmiterm([2,9,9,Q62],-1,1);
lmiterm([2,7,10,0],0);
lmiterm([2,8,10,0],0);
lmiterm([2,9,10,0],0);
lmiterm([2,7,11,0],0);
lmiterm([2,7,12,0],0);
lmiterm([2,7,13,0],0);
lmiterm([2,8,11,0],0);
lmiterm([2,8,12,0],0);
lmiterm([2,8,13,0],0);
lmiterm([2,9,11,0],0);
lmiterm([2,9,12,0],0);
lmiterm([2,9,13,0],0);

lmiterm([2,7,14,0],0);
lmiterm([2,7,15,0],0);
lmiterm([2,7,16,0],0);
lmiterm([2,8,14,0],0);
lmiterm([2,8,15,0],0);
lmiterm([2,8,16,0],0);
lmiterm([2,9,14,0],0);
lmiterm([2,9,15,0],0);
lmiterm([2,9,16,0],0);

lmiterm([2,10,10,0],-gama1*gama1);
lmiterm([2,10,11,0],0);
lmiterm([2,10,12,0],0);
lmiterm([2,10,13,0],0);
lmiterm([2,10,14,0],0);
lmiterm([2,10,15,0],0);
lmiterm([2,10,16,0],0);

lmiterm([2,11,11,R10],-1/tao,1);lmiterm([2,11,11,R12],-1/tao,1);
lmiterm([2,11,12,R20],-1/tao,1);
lmiterm([2,11,13,R30],-1/tao,1);
lmiterm([2,12,12,R40],-1/tao,1);lmiterm([2,12,12,R42],-1/tao,1);
lmiterm([2,12,13,R50],-1/tao,1);
lmiterm([2,13,13,R60],-1/tao,1);lmiterm([2,13,13,R62],-1/tao,1);

lmiterm([2,11,14,0],0);
lmiterm([2,11,15,0],0);
lmiterm([2,11,16,0],0);
lmiterm([2,12,14,0],0);
lmiterm([2,12,15,0],0);
lmiterm([2,12,16,0],0);
lmiterm([2,13,14,0],0);
lmiterm([2,13,15,0],0);
lmiterm([2,13,16,0],0);

lmiterm([2,14,14,0],-1);
lmiterm([2,14,15,0],0);
lmiterm([2,14,16,0],0);
lmiterm([2,15,15,0],-1);
lmiterm([2,15,16,0],0);
lmiterm([2,16,16,0],-1);




lmiterm([3,1,1,Q10],1,1);lmiterm([3,1,1,Q13],1,1);
lmiterm([3,1,1,R10],tao/(lbd*lbd),1);lmiterm([3,1,1,R13],tao/(lbd*lbd),1);
lmiterm([3,1,1,P10],-2/lbd,1);lmiterm([3,1,1,P13],-2/lbd,1);
lmiterm([3,1,2,Q20],1,1);
lmiterm([3,1,2,R20],tao/(lbd*lbd),1);
lmiterm([3,1,2,P20],-2/lbd,1);
lmiterm([3,1,3,Q30],1,1);
lmiterm([3,1,3,R30],tao/(lbd*lbd),1);
lmiterm([3,1,3,P30],-2/lbd,1);
lmiterm([3,2,2,Q40],1,1);lmiterm([3,2,2,Q43],1,1);
lmiterm([3,2,2,R40],tao/(lbd*lbd),1);lmiterm([3,2,2,R43],tao/(lbd*lbd),1);
lmiterm([3,2,2,P40],-2/lbd,1);lmiterm([3,2,2,P43],-2/lbd,1);
lmiterm([3,2,3,Q50],1,1);
lmiterm([3,2,3,R50],tao/(lbd*lbd),1);
lmiterm([3,2,3,P50],-2/lbd,1);
lmiterm([3,3,3,Q60],1,1);lmiterm([3,3,3,Q63],1,1);
lmiterm([3,3,3,R60],tao/(lbd*lbd),1);lmiterm([3,3,3,R63],tao/(lbd*lbd),1);
lmiterm([3,3,3,P60],-2/lbd,1);lmiterm([3,3,3,P63],-2/lbd,1);
lmiterm([3,1,4,P10],1,1);lmiterm([3,1,4,P13],1,1);
lmiterm([3,1,4,R10],-tao/lbd,1);lmiterm([3,1,4,R13],-tao/lbd,1);
lmiterm([3,1,4,-Mz1],1,1);lmiterm([3,1,4,-Mz1],lbd,A3');lmiterm([3,1,4,-Zg3],lbd,B3');
lmiterm([3,1,5,P20],1,1); 
lmiterm([3,1,5,R20],-tao/lbd,1);
lmiterm([3,1,6,P30],1,1);
lmiterm([3,1,6,R30],-tao/lbd,1);
lmiterm([3,2,4,-P20],1,1); 
lmiterm([3,2,4,-R20],-tao/lbd,1); lmiterm([3,2,4,-Zg3],-lbd,B3');
lmiterm([3,2,5,P40],1,1);lmiterm([3,2,5,P43],1,1);
lmiterm([3,2,5,R40],-tao/lbd,1);lmiterm([3,2,5,R43],-tao/lbd,1);
lmiterm([3,2,5,-Mz1],1,1);lmiterm([3,2,5,-Mz1],lbd,A3');lmiterm([3,2,5,-Mz1],-lbd,C1'*Ls3');
lmiterm([3,2,6,P50],1,1);
lmiterm([3,2,6,R50],-tao/lbd,1);lmiterm([3,2,6,-Mz1],-lbd,A3'*C1'*Ld3');
lmiterm([3,3,4,-P30],1,1);
lmiterm([3,3,4,-R30],-tao/lbd,1);lmiterm([3,3,4,-Mz2],lbd,B3');
lmiterm([3,3,5,-P50],1,1);
lmiterm([3,3,5,-R50],-tao/lbd,1);lmiterm([3,3,5,-Mz2],lbd,B3');
lmiterm([3,3,6,P60],1,1);lmiterm([3,3,6,P63],1,1);
lmiterm([3,3,6,R60],-tao/lbd,1);lmiterm([3,3,6,R63],-tao/lbd,1);
lmiterm([3,3,6,-Mz2],1,1);lmiterm([3,3,6,-Mz2],-lbd,B3'*C1'*Ld3');

lmiterm([3,1,7,0],0);
lmiterm([3,1,8,0],0);
lmiterm([3,1,9,0],0);
lmiterm([3,2,7,0],0);
lmiterm([3,2,8,0],0);
lmiterm([3,2,9,0],0);
lmiterm([3,3,7,0],0);
lmiterm([3,3,8,0],0);
lmiterm([3,3,9,0],0);
lmiterm([3,1,10,0],0);
lmiterm([3,2,10,0],0);
lmiterm([3,3,10,0],0);
lmiterm([3,1,11,0],0);
lmiterm([3,1,12,0],0);
lmiterm([3,1,13,0],0);
lmiterm([3,2,11,0],0);
lmiterm([3,2,12,0],0);
lmiterm([3,2,13,0],0);
lmiterm([3,3,11,0],0);
lmiterm([3,3,12,0],0);
lmiterm([3,3,13,0],0);
lmiterm([3,1,14,-Mz1],1,1);
lmiterm([3,1,15,0],0);
lmiterm([3,1,16,0],0);
lmiterm([3,2,14,0],0);
lmiterm([3,2,15,0],0);
lmiterm([3,2,16,0],0);
lmiterm([3,3,14,0],0);
lmiterm([3,3,15,0],0);
lmiterm([3,3,16,0],0);


lmiterm([3,4,4,Mz1],-lbd,1,'s');
lmiterm([3,4,4,R10],tao,1);lmiterm([3,4,4,R13],tao,1);
lmiterm([3,4,5,R20],tao,1);
lmiterm([3,4,6,R30],tao,1);
lmiterm([3,5,5,Mz1],-lbd,1,'s');
lmiterm([3,5,5,R40],tao,1);lmiterm([3,5,5,R43],tao,1);
lmiterm([3,5,6,R50],tao,1);
lmiterm([3,6,6,Mz2],-lbd,1,'s');
lmiterm([3,6,6,R60],tao,1);lmiterm([3,6,6,R63],tao,1);

lmiterm([3,4,7,0],0);
lmiterm([3,4,8,0],0);
lmiterm([3,4,9,0],0);
lmiterm([3,5,7,0],0);
lmiterm([3,5,8,0],0);
lmiterm([3,5,9,0],0);
lmiterm([3,6,7,0],0);
lmiterm([3,6,8,0],0);
lmiterm([3,6,9,0],0);

lmiterm([3,4,10,0],0);
lmiterm([3,5,10,0],0);
lmiterm([3,6,10,0],lbd);
lmiterm([3,4,11,0],0);
lmiterm([3,4,12,0],0);
lmiterm([3,4,13,0],0);
lmiterm([3,5,11,0],0);lmiterm([3,5,12,-Mz1],-lbd,C1'*Ls3');lmiterm([3,5,13,-Mz1],-lbd,A3'*C1'*Ld3');
lmiterm([3,6,11,0],0);
lmiterm([3,6,12,0],0);lmiterm([3,6,13,-Mz2],-lbd,B3'*C1'*Ld3');

lmiterm([3,4,14,0],0);
lmiterm([3,4,15,0],0);
lmiterm([3,4,16,0],0);
lmiterm([3,5,14,0],0);
lmiterm([3,5,15,0],0);
lmiterm([3,5,16,0],0);
lmiterm([3,6,14,0],0);
lmiterm([3,6,15,0],0);
lmiterm([3,6,16,0],0);

lmiterm([3,7,7,Q10],-1,1);lmiterm([3,7,7,Q13],-1,1);
lmiterm([3,7,8,Q20],-1,1);
lmiterm([3,7,9,Q30],-1,1);
lmiterm([3,8,8,Q40],-1,1);lmiterm([3,8,8,Q43],-1,1);
lmiterm([3,8,9,Q50],-1,1);
lmiterm([3,9,9,Q60],-1,1);lmiterm([3,9,9,Q63],-1,1);
lmiterm([3,7,10,0],0);
lmiterm([3,8,10,0],0);
lmiterm([3,9,10,0],0);
lmiterm([3,7,11,0],0);
lmiterm([3,7,12,0],0);
lmiterm([3,7,13,0],0);
lmiterm([3,8,11,0],0);
lmiterm([3,8,12,0],0);
lmiterm([3,8,13,0],0);
lmiterm([3,9,11,0],0);
lmiterm([3,9,12,0],0);
lmiterm([3,9,13,0],0);

lmiterm([3,7,14,0],0);
lmiterm([3,7,15,0],0);
lmiterm([3,7,16,0],0);
lmiterm([3,8,14,0],0);
lmiterm([3,8,15,0],0);
lmiterm([3,8,16,0],0);
lmiterm([3,9,14,0],0);
lmiterm([3,9,15,0],0);
lmiterm([3,9,16,0],0);

lmiterm([3,10,10,0],-gama1*gama1);
lmiterm([3,10,11,0],0);
lmiterm([3,10,12,0],0);
lmiterm([3,10,13,0],0);
lmiterm([3,10,14,0],0);
lmiterm([3,10,15,0],0);
lmiterm([3,10,16,0],0);

lmiterm([3,11,11,R10],-1/tao,1);lmiterm([3,11,11,R13],-1/tao,1);
lmiterm([3,11,12,R20],-1/tao,1);
lmiterm([3,11,13,R30],-1/tao,1);
lmiterm([3,12,12,R40],-1/tao,1);lmiterm([3,12,12,R43],-1/tao,1);
lmiterm([3,12,13,R50],-1/tao,1);
lmiterm([3,13,13,R60],-1/tao,1);lmiterm([3,13,13,R63],-1/tao,1);

lmiterm([3,11,14,0],0);
lmiterm([3,11,15,0],0);
lmiterm([3,11,16,0],0);
lmiterm([3,12,14,0],0);
lmiterm([3,12,15,0],0);
lmiterm([3,12,16,0],0);
lmiterm([3,13,14,0],0);
lmiterm([3,13,15,0],0);
lmiterm([3,13,16,0],0);

lmiterm([3,14,14,0],-1);
lmiterm([3,14,15,0],0);
lmiterm([3,14,16,0],0);
lmiterm([3,15,15,0],-1);
lmiterm([3,15,16,0],0);
lmiterm([3,16,16,0],-1);



lmiterm([4,1,1,Q10],1,1);lmiterm([4,1,1,Q14],1,1);
lmiterm([4,1,1,R10],tao/(lbd*lbd),1);lmiterm([4,1,1,R14],tao/(lbd*lbd),1);
lmiterm([4,1,1,P10],-2/lbd,1);lmiterm([4,1,1,P14],-2/lbd,1);
lmiterm([4,1,2,Q20],1,1);
lmiterm([4,1,2,R20],tao/(lbd*lbd),1);
lmiterm([4,1,2,P20],-2/lbd,1);
lmiterm([4,1,3,Q30],1,1);
lmiterm([4,1,3,R30],tao/(lbd*lbd),1);
lmiterm([4,1,3,P30],-2/lbd,1);
lmiterm([4,2,2,Q40],1,1);lmiterm([4,2,2,Q44],1,1);
lmiterm([4,2,2,R40],tao/(lbd*lbd),1);lmiterm([4,2,2,R44],tao/(lbd*lbd),1);
lmiterm([4,2,2,P40],-2/lbd,1);lmiterm([4,2,2,P44],-2/lbd,1);
lmiterm([4,2,3,Q50],1,1);
lmiterm([4,2,3,R50],tao/(lbd*lbd),1);
lmiterm([4,2,3,P50],-2/lbd,1);
lmiterm([4,3,3,Q60],1,1);lmiterm([4,3,3,Q64],1,1);
lmiterm([4,3,3,R60],tao/(lbd*lbd),1);lmiterm([4,3,3,R64],tao/(lbd*lbd),1);
lmiterm([4,3,3,P60],-2/lbd,1);lmiterm([4,3,3,P64],-2/lbd,1);
lmiterm([4,1,4,P10],1,1);lmiterm([4,1,4,P14],1,1);
lmiterm([4,1,4,R10],-tao/lbd,1);lmiterm([4,1,4,R14],-tao/lbd,1);
lmiterm([4,1,4,-Mz1],1,1);lmiterm([4,1,4,-Mz1],lbd,A4');lmiterm([4,1,4,-Zg4],lbd,B4');
lmiterm([4,1,5,P20],1,1); 
lmiterm([4,1,5,R20],-tao/lbd,1);
lmiterm([4,1,6,P30],1,1);
lmiterm([4,1,6,R30],-tao/lbd,1);
lmiterm([4,2,4,-P20],1,1); 
lmiterm([4,2,4,-R20],-tao/lbd,1); lmiterm([4,2,4,-Zg4],-lbd,B4');
lmiterm([4,2,5,P40],1,1);lmiterm([4,2,5,P44],1,1);
lmiterm([4,2,5,R40],-tao/lbd,1);lmiterm([4,2,5,R44],-tao/lbd,1);
lmiterm([4,2,5,-Mz1],1,1);lmiterm([4,2,5,-Mz1],lbd,A4');lmiterm([4,2,5,-Mz1],-lbd,C1'*Ls4');
lmiterm([4,2,6,P50],1,1);
lmiterm([4,2,6,R50],-tao/lbd,1);lmiterm([4,2,6,-Mz1],-lbd,A4'*C1'*Ld4');
lmiterm([4,3,4,-P30],1,1);
lmiterm([4,3,4,-R30],-tao/lbd,1);lmiterm([4,3,4,-Mz2],lbd,B4');
lmiterm([4,3,5,-P50],1,1);
lmiterm([4,3,5,-R50],-tao/lbd,1);lmiterm([4,3,5,-Mz2],lbd,B4');
lmiterm([4,3,6,P60],1,1);lmiterm([4,3,6,P64],1,1);
lmiterm([4,3,6,R60],-tao/lbd,1);lmiterm([4,3,6,R64],-tao/lbd,1);
lmiterm([4,3,6,-Mz2],1,1);lmiterm([4,3,6,-Mz2],-lbd,B4'*C1'*Ld4');

lmiterm([4,1,7,0],0);
lmiterm([4,1,8,0],0);
lmiterm([4,1,9,0],0);
lmiterm([4,2,7,0],0);
lmiterm([4,2,8,0],0);
lmiterm([4,2,9,0],0);
lmiterm([4,3,7,0],0);
lmiterm([4,3,8,0],0);
lmiterm([4,3,9,0],0);
lmiterm([4,1,10,0],0);
lmiterm([4,2,10,0],0);
lmiterm([4,3,10,0],0);
lmiterm([4,1,11,0],0);
lmiterm([4,1,12,0],0);
lmiterm([4,1,13,0],0);
lmiterm([4,2,11,0],0);
lmiterm([4,2,12,0],0);
lmiterm([4,2,13,0],0);
lmiterm([4,3,11,0],0);
lmiterm([4,3,12,0],0);
lmiterm([4,3,13,0],0);
lmiterm([4,1,14,-Mz1],1,1);
lmiterm([4,1,15,0],0);
lmiterm([4,1,16,0],0);
lmiterm([4,2,14,0],0);
lmiterm([4,2,15,0],0);
lmiterm([4,2,16,0],0);
lmiterm([4,3,14,0],0);
lmiterm([4,3,15,0],0);
lmiterm([4,3,16,0],0);


lmiterm([4,4,4,Mz1],-lbd,1,'s');
lmiterm([4,4,4,R10],tao,1);lmiterm([4,4,4,R14],tao,1);
lmiterm([4,4,5,R20],tao,1);
lmiterm([4,4,6,R30],tao,1);
lmiterm([4,5,5,Mz1],-lbd,1,'s');
lmiterm([4,5,5,R40],tao,1);lmiterm([4,5,5,R44],tao,1);
lmiterm([4,5,6,R50],tao,1);
lmiterm([4,6,6,Mz2],-lbd,1,'s');
lmiterm([4,6,6,R60],tao,1);lmiterm([4,6,6,R64],tao,1);

lmiterm([4,4,7,0],0);
lmiterm([4,4,8,0],0);
lmiterm([4,4,9,0],0);
lmiterm([4,5,7,0],0);
lmiterm([4,5,8,0],0);
lmiterm([4,5,9,0],0);
lmiterm([4,6,7,0],0);
lmiterm([4,6,8,0],0);
lmiterm([4,6,9,0],0);

lmiterm([4,4,10,0],0);
lmiterm([4,5,10,0],0);
lmiterm([4,6,10,0],lbd);
lmiterm([4,4,11,0],0);
lmiterm([4,4,12,0],0);
lmiterm([4,4,13,0],0);
lmiterm([4,5,11,0],0);lmiterm([4,5,12,-Mz1],-lbd,C1'*Ls4');lmiterm([4,5,13,-Mz1],-lbd,A4'*C1'*Ld4');
lmiterm([4,6,11,0],0);
lmiterm([4,6,12,0],0);lmiterm([4,6,13,-Mz2],-lbd,B4'*C1'*Ld4');

lmiterm([4,4,14,0],0);
lmiterm([4,4,15,0],0);
lmiterm([4,4,16,0],0);
lmiterm([4,5,14,0],0);
lmiterm([4,5,15,0],0);
lmiterm([4,5,16,0],0);
lmiterm([4,6,14,0],0);
lmiterm([4,6,15,0],0);
lmiterm([4,6,16,0],0);

lmiterm([4,7,7,Q10],-1,1);lmiterm([4,7,7,Q14],-1,1);
lmiterm([4,7,8,Q20],-1,1);
lmiterm([4,7,9,Q30],-1,1);
lmiterm([4,8,8,Q40],-1,1);lmiterm([4,8,8,Q44],-1,1);
lmiterm([4,8,9,Q50],-1,1);
lmiterm([4,9,9,Q60],-1,1);lmiterm([4,9,9,Q64],-1,1);
lmiterm([4,7,10,0],0);
lmiterm([4,8,10,0],0);
lmiterm([4,9,10,0],0);
lmiterm([4,7,11,0],0);
lmiterm([4,7,12,0],0);
lmiterm([4,7,13,0],0);
lmiterm([4,8,11,0],0);
lmiterm([4,8,12,0],0);
lmiterm([4,8,13,0],0);
lmiterm([4,9,11,0],0);
lmiterm([4,9,12,0],0);
lmiterm([4,9,13,0],0);

lmiterm([4,7,14,0],0);
lmiterm([4,7,15,0],0);
lmiterm([4,7,16,0],0);
lmiterm([4,8,14,0],0);
lmiterm([4,8,15,0],0);
lmiterm([4,8,16,0],0);
lmiterm([4,9,14,0],0);
lmiterm([4,9,15,0],0);
lmiterm([4,9,16,0],0);

lmiterm([4,10,10,0],-gama1*gama1);
lmiterm([4,10,11,0],0);
lmiterm([4,10,12,0],0);
lmiterm([4,10,13,0],0);
lmiterm([4,10,14,0],0);
lmiterm([4,10,15,0],0);
lmiterm([4,10,16,0],0);

lmiterm([4,11,11,R10],-1/tao,1);lmiterm([4,11,11,R14],-1/tao,1);
lmiterm([4,11,12,R20],-1/tao,1);
lmiterm([4,11,13,R30],-1/tao,1);
lmiterm([4,12,12,R40],-1/tao,1);lmiterm([4,12,12,R44],-1/tao,1);
lmiterm([4,12,13,R50],-1/tao,1);
lmiterm([4,13,13,R60],-1/tao,1);lmiterm([4,13,13,R64],-1/tao,1);

lmiterm([4,11,14,0],0);
lmiterm([4,11,15,0],0);
lmiterm([4,11,16,0],0);
lmiterm([4,12,14,0],0);
lmiterm([4,12,15,0],0);
lmiterm([4,12,16,0],0);
lmiterm([4,13,14,0],0);
lmiterm([4,13,15,0],0);
lmiterm([4,13,16,0],0);

lmiterm([4,14,14,0],-1);
lmiterm([4,14,15,0],0);
lmiterm([4,14,16,0],0);
lmiterm([4,15,15,0],-1);
lmiterm([4,15,16,0],0);
lmiterm([4,16,16,0],-1);



lmiterm([5,1,1,Q10],4/3,1);lmiterm([5,1,1,Q11],5/6,1);lmiterm([5,1,1,Q12],5/6,1);
lmiterm([5,1,1,R10],tao/(lbd*lbd),4/3);lmiterm([5,1,1,R11],tao/(lbd*lbd),5/6);lmiterm([5,1,1,R12],tao/(lbd*lbd),5/6);
lmiterm([5,1,1,P10],-2/lbd,4/3);lmiterm([5,1,1,P11],-2/lbd,5/6);lmiterm([5,1,1,P12],-2/lbd,5/6);
lmiterm([5,1,2,Q20],4/3,1);
lmiterm([5,1,2,R20],tao/(lbd*lbd),4/3);
lmiterm([5,1,2,P20],-2/lbd,4/3);
lmiterm([5,1,3,Q30],4/3,1);
lmiterm([5,1,3,R30],tao/(lbd*lbd),4/3);
lmiterm([5,1,3,P30],-2/lbd,4/3);
lmiterm([5,2,2,Q40],4/3,1);lmiterm([5,2,2,Q41],5/6,1);lmiterm([5,2,2,Q42],5/6,1);
lmiterm([5,2,2,R40],tao/(lbd*lbd),4/3);lmiterm([5,2,2,R41],tao/(lbd*lbd),5/6);lmiterm([5,2,2,R42],tao/(lbd*lbd),5/6);
lmiterm([5,2,2,P40],-2/lbd,4/3);lmiterm([5,2,2,P41],-2/lbd,5/6);lmiterm([5,2,2,P42],-2/lbd,5/6);
lmiterm([5,2,3,Q50],4/3,1);
lmiterm([5,2,3,R50],tao/(lbd*lbd),4/3);
lmiterm([5,2,3,P50],-2/lbd,4/3);
lmiterm([5,3,3,Q60],4/3,1);lmiterm([5,3,3,Q61],5/6,1);lmiterm([5,3,3,Q62],5/6,1);
lmiterm([5,3,3,R60],tao/(lbd*lbd),4/3);lmiterm([5,3,3,R61],tao/(lbd*lbd),5/6);lmiterm([5,3,3,R62],tao/(lbd*lbd),5/6);
lmiterm([5,3,3,P60],-2/lbd,4/3);lmiterm([5,3,3,P61],-2/lbd,5/6);lmiterm([5,3,3,P62],-2/lbd,5/6);
lmiterm([5,1,4,P10],4/3,1);lmiterm([5,1,4,P11],5/6,1);lmiterm([5,1,4,P12],5/6,1);
lmiterm([5,1,4,R10],-tao/lbd,4/3);lmiterm([5,1,4,R11],-tao/lbd,5/6);lmiterm([5,1,4,R12],-tao/lbd,5/6);
lmiterm([5,1,4,-Mz1],4/3,1);lmiterm([5,1,4,-Mz1],(5/6)*lbd,A1');lmiterm([5,1,4,-Mz1],(5/6)*lbd,A2');lmiterm([5,1,4,-Zg1],(5/6)*lbd,B2');lmiterm([5,1,4,-Zg2],(5/6)*lbd,B1');
lmiterm([5,1,5,P20],4/3,1); 
lmiterm([5,1,5,R20],-tao/lbd,4/3);
lmiterm([5,1,6,P30],4/3,1);
lmiterm([5,1,6,R30],-tao/lbd,4/3);
lmiterm([5,2,4,-P20],4/3,1); 
lmiterm([5,2,4,-R20],-tao/lbd,4/3); lmiterm([5,2,4,-Zg1],-(5/6)*lbd,B2');lmiterm([5,2,4,-Zg2],-(5/6)*lbd,B1');
lmiterm([5,2,5,P40],4/3,1);lmiterm([5,2,5,P41],5/6,1);lmiterm([5,2,5,P42],5/6,1);
lmiterm([5,2,5,R40],-tao/lbd,4/3);lmiterm([5,2,5,R41],-tao/lbd,5/6);lmiterm([5,2,5,R42],-tao/lbd,5/6);
lmiterm([5,2,5,-Mz1],4/3,1);lmiterm([5,2,5,-Mz1],(5/6)*lbd,A1');lmiterm([5,2,5,-Mz1],(5/6)*lbd,A2');lmiterm([5,2,5,-Mz1],-(5/6)*lbd,C1'*Ls1');lmiterm([5,2,5,-Mz1],-(5/6)*lbd,C1'*Ls2');
lmiterm([5,2,6,P50],4/3,1);
lmiterm([5,2,6,R50],-tao/lbd,4/3);lmiterm([5,2,6,-Mz1],-(5/6)*lbd,A1'*C1'*Ld2');lmiterm([5,2,6,-Mz1],-(5/6)*lbd,A2'*C1'*Ld1');
lmiterm([5,3,4,-P30],4/3,1);
lmiterm([5,3,4,-R30],-tao/lbd,4/3);lmiterm([5,3,4,-Mz2],(5/6)*lbd,B1');lmiterm([5,3,4,-Mz2],(5/6)*lbd,B2');
lmiterm([5,3,5,-P50],4/3,1);
lmiterm([5,3,5,-R50],-tao/lbd,4/3);lmiterm([5,3,5,-Mz2],(5/6)*lbd,B1');lmiterm([5,3,5,-Mz2],(5/6)*lbd,B2');
lmiterm([5,3,6,P60],4/3,1);lmiterm([5,3,6,P61],5/6,1);lmiterm([5,3,6,P62],5/6,1);
lmiterm([5,3,6,R60],-tao/lbd,4/3);lmiterm([5,3,6,R61],-tao/lbd,5/6);lmiterm([5,3,6,R62],-tao/lbd,5/6);
lmiterm([5,3,6,-Mz2],4/3,1);lmiterm([5,3,6,-Mz2],-(5/6)*lbd,B1'*C1'*Ld2');lmiterm([5,3,6,-Mz2],-(5/6)*lbd,B2'*C1'*Ld1');

lmiterm([5,1,7,0],0);
lmiterm([5,1,8,0],0);
lmiterm([5,1,9,0],0);
lmiterm([5,2,7,0],0);
lmiterm([5,2,8,0],0);
lmiterm([5,2,9,0],0);
lmiterm([5,3,7,0],0);
lmiterm([5,3,8,0],0);
lmiterm([5,3,9,0],0);
lmiterm([5,1,10,0],0);
lmiterm([5,2,10,0],0);
lmiterm([5,3,10,0],0);
lmiterm([5,1,11,0],0);
lmiterm([5,1,12,0],0);
lmiterm([5,1,13,0],0);
lmiterm([5,2,11,0],0);
lmiterm([5,2,12,0],0);
lmiterm([5,2,13,0],0);
lmiterm([5,3,11,0],0);
lmiterm([5,3,12,0],0);
lmiterm([5,3,13,0],0);
lmiterm([5,1,14,-Mz1],4/3,1);
lmiterm([5,1,15,0],0);
lmiterm([5,1,16,0],0);
lmiterm([5,2,14,0],0);
lmiterm([5,2,15,0],0);
lmiterm([5,2,16,0],0);
lmiterm([5,3,14,0],0);
lmiterm([5,3,15,0],0);
lmiterm([5,3,16,0],0);

lmiterm([5,4,4,Mz1],-lbd,4/3,'s');
lmiterm([5,4,4,R10],tao,4/3);lmiterm([5,4,4,R11],tao,5/6);lmiterm([5,4,4,R12],tao,5/6);
lmiterm([5,4,5,R20],tao,4/3);
lmiterm([5,4,6,R30],tao,4/3);
lmiterm([5,5,5,Mz1],-lbd,4/3,'s');
lmiterm([5,5,5,R40],tao,4/3);lmiterm([5,5,5,R41],tao,5/6);lmiterm([5,5,5,R42],tao,5/6);
lmiterm([5,5,6,R50],tao,4/3);
lmiterm([5,6,6,Mz2],-lbd,4/3,'s');
lmiterm([5,6,6,R60],tao,4/3);lmiterm([5,6,6,R61],tao,5/6);lmiterm([5,6,6,R62],tao,5/6);

lmiterm([5,4,7,0],0);
lmiterm([5,4,8,0],0);
lmiterm([5,4,9,0],0);
lmiterm([5,5,7,0],0);
lmiterm([5,5,8,0],0);
lmiterm([5,5,9,0],0);
lmiterm([5,6,7,0],0);
lmiterm([5,6,8,0],0);
lmiterm([5,6,9,0],0);

lmiterm([5,4,10,0],0);
lmiterm([5,5,10,0],0);
lmiterm([5,6,10,0],(4/3)*lbd);
lmiterm([5,4,11,0],0);
lmiterm([5,4,12,0],0);
lmiterm([5,4,13,0],0);
lmiterm([5,5,11,0],0);lmiterm([5,5,12,-Mz1],-(5/6)*lbd,C1'*Ls1');lmiterm([5,5,12,-Mz1],-(5/6)*lbd,C1'*Ls2');lmiterm([5,5,13,-Mz1],-(5/6)*lbd,A1'*C1'*Ld2');lmiterm([5,5,13,-Mz1],-(5/6)*lbd,A2'*C1'*Ld1');
lmiterm([5,6,11,0],0);
lmiterm([5,6,12,0],0);lmiterm([5,6,13,-Mz2],-(5/6)*lbd,B1'*C1'*Ld2');lmiterm([5,6,13,-Mz2],-(5/6)*lbd,B2'*C1'*Ld1');

lmiterm([5,4,14,0],0);
lmiterm([5,4,15,0],0);
lmiterm([5,4,16,0],0);
lmiterm([5,5,14,0],0);
lmiterm([5,5,15,0],0);
lmiterm([5,5,16,0],0);
lmiterm([5,6,14,0],0);
lmiterm([5,6,15,0],0);
lmiterm([5,6,16,0],0);

lmiterm([5,7,7,Q10],-4/3,1);lmiterm([5,7,7,Q11],-5/6,1);lmiterm([5,7,7,Q12],-5/6,1);
lmiterm([5,7,8,Q20],-4/3,1);
lmiterm([5,7,9,Q30],-4/3,1);
lmiterm([5,8,8,Q40],-4/3,1);lmiterm([5,8,8,Q41],-5/6,1);lmiterm([5,8,8,Q42],-5/6,1);
lmiterm([5,8,9,Q50],-4/3,1);
lmiterm([5,9,9,Q60],-4/3,1);lmiterm([5,9,9,Q61],-5/6,1);lmiterm([5,9,9,Q62],-5/6,1);
lmiterm([5,7,10,0],0);
lmiterm([5,8,10,0],0);
lmiterm([5,9,10,0],0);
lmiterm([5,7,11,0],0);
lmiterm([5,7,12,0],0);
lmiterm([5,7,13,0],0);
lmiterm([5,8,11,0],0);
lmiterm([5,8,12,0],0);
lmiterm([5,8,13,0],0);
lmiterm([5,9,11,0],0);
lmiterm([5,9,12,0],0);
lmiterm([5,9,13,0],0);

lmiterm([5,7,14,0],0);
lmiterm([5,7,15,0],0);
lmiterm([5,7,16,0],0);
lmiterm([5,8,14,0],0);
lmiterm([5,8,15,0],0);
lmiterm([5,8,16,0],0);
lmiterm([5,9,14,0],0);
lmiterm([5,9,15,0],0);
lmiterm([5,9,16,0],0);

lmiterm([5,10,10,0],-(4/3)*gama1*gama1);
lmiterm([5,10,11,0],0);
lmiterm([5,10,12,0],0);
lmiterm([5,10,13,0],0);
lmiterm([5,10,14,0],0);
lmiterm([5,10,15,0],0);
lmiterm([5,10,16,0],0);

lmiterm([5,11,11,R10],-1/tao,4/3);lmiterm([5,11,11,R11],-1/tao,5/6);lmiterm([5,11,11,R12],-1/tao,5/6);
lmiterm([5,11,12,R20],-1/tao,4/3);
lmiterm([5,11,13,R30],-1/tao,4/3);
lmiterm([5,12,12,R40],-1/tao,4/3);lmiterm([5,12,12,R41],-1/tao,5/6);lmiterm([5,12,12,R42],-1/tao,5/6);
lmiterm([5,12,13,R50],-1/tao,4/3);
lmiterm([5,13,13,R60],-1/tao,4/3);lmiterm([5,13,13,R61],-1/tao,5/6);lmiterm([5,13,13,R62],-1/tao,5/6);

lmiterm([5,11,14,0],0);
lmiterm([5,11,15,0],0);
lmiterm([5,11,16,0],0);
lmiterm([5,12,14,0],0);
lmiterm([5,12,15,0],0);
lmiterm([5,12,16,0],0);
lmiterm([5,13,14,0],0);
lmiterm([5,13,15,0],0);
lmiterm([5,13,16,0],0);

lmiterm([5,14,14,0],-4/3);
lmiterm([5,14,15,0],0);
lmiterm([5,14,16,0],0);
lmiterm([5,15,15,0],-4/3);
lmiterm([5,15,16,0],0);
lmiterm([5,16,16,0],-4/3);




lmiterm([6,1,1,Q10],4/3,1);lmiterm([6,1,1,Q11],5/6,1);lmiterm([6,1,1,Q13],5/6,1);
lmiterm([6,1,1,R10],tao/(lbd*lbd),4/3);lmiterm([6,1,1,R11],tao/(lbd*lbd),5/6);lmiterm([6,1,1,R13],tao/(lbd*lbd),5/6);
lmiterm([6,1,1,P10],-2/lbd,4/3);lmiterm([6,1,1,P11],-2/lbd,5/6);lmiterm([6,1,1,P13],-2/lbd,5/6);
lmiterm([6,1,2,Q20],4/3,1);
lmiterm([6,1,2,R20],tao/(lbd*lbd),4/3);
lmiterm([6,1,2,P20],-2/lbd,4/3);
lmiterm([6,1,3,Q30],4/3,1);
lmiterm([6,1,3,R30],tao/(lbd*lbd),4/3);
lmiterm([6,1,3,P30],-2/lbd,4/3);
lmiterm([6,2,2,Q40],4/3,1);lmiterm([6,2,2,Q41],5/6,1);lmiterm([6,2,2,Q43],5/6,1);
lmiterm([6,2,2,R40],tao/(lbd*lbd),4/3);lmiterm([6,2,2,R41],tao/(lbd*lbd),5/6);lmiterm([6,2,2,R43],tao/(lbd*lbd),5/6);
lmiterm([6,2,2,P40],-2/lbd,4/3);lmiterm([6,2,2,P41],-2/lbd,5/6);lmiterm([6,2,2,P43],-2/lbd,5/6);
lmiterm([6,2,3,Q50],4/3,1);
lmiterm([6,2,3,R50],tao/(lbd*lbd),4/3);
lmiterm([6,2,3,P50],-2/lbd,4/3);
lmiterm([6,3,3,Q60],4/3,1);lmiterm([6,3,3,Q61],5/6,1);lmiterm([6,3,3,Q63],5/6,1);
lmiterm([6,3,3,R60],tao/(lbd*lbd),4/3);lmiterm([6,3,3,R61],tao/(lbd*lbd),5/6);lmiterm([6,3,3,R63],tao/(lbd*lbd),5/6);
lmiterm([6,3,3,P60],-2/lbd,4/3);lmiterm([6,3,3,P61],-2/lbd,5/6);lmiterm([6,3,3,P63],-2/lbd,5/6);
lmiterm([6,1,4,P10],4/3,1);lmiterm([6,1,4,P11],5/6,1);lmiterm([6,1,4,P13],5/6,1);
lmiterm([6,1,4,R10],-tao/lbd,4/3);lmiterm([6,1,4,R11],-tao/lbd,5/6);lmiterm([6,1,4,R13],-tao/lbd,5/6);
lmiterm([6,1,4,-Mz1],4/3,1);lmiterm([6,1,4,-Mz1],(5/6)*lbd,A1');lmiterm([6,1,4,-Mz1],(5/6)*lbd,A3');lmiterm([6,1,4,-Zg1],(5/6)*lbd,B3');lmiterm([6,1,4,-Zg3],(5/6)*lbd,B1');
lmiterm([6,1,5,P20],4/3,1); 
lmiterm([6,1,5,R20],-tao/lbd,4/3);
lmiterm([6,1,6,P30],4/3,1);
lmiterm([6,1,6,R30],-tao/lbd,4/3);
lmiterm([6,2,4,-P20],4/3,1); 
lmiterm([6,2,4,-R20],-tao/lbd,4/3); lmiterm([6,2,4,-Zg1],-(5/6)*lbd,B3');lmiterm([6,2,4,-Zg3],-(5/6)*lbd,B1');
lmiterm([6,2,5,P40],4/3,1);lmiterm([6,2,5,P41],5/6,1);lmiterm([6,2,5,P43],5/6,1);
lmiterm([6,2,5,R40],-tao/lbd,4/3);lmiterm([6,2,5,R41],-tao/lbd,5/6);lmiterm([6,2,5,R43],-tao/lbd,5/6);
lmiterm([6,2,5,-Mz1],4/3,1);lmiterm([6,2,5,-Mz1],(5/6)*lbd,A1');lmiterm([6,2,5,-Mz1],(5/6)*lbd,A3');lmiterm([6,2,5,-Mz1],-(5/6)*lbd,C1'*Ls1');lmiterm([6,2,5,-Mz1],-(5/6)*lbd,C1'*Ls3');
lmiterm([6,2,6,P50],4/3,1);
lmiterm([6,2,6,R50],-tao/lbd,4/3);lmiterm([6,2,6,-Mz1],-(5/6)*lbd,A1'*C1'*Ld3');lmiterm([6,2,6,-Mz1],-(5/6)*lbd,A3'*C1'*Ld1');
lmiterm([6,3,4,-P30],4/3,1);
lmiterm([6,3,4,-R30],-tao/lbd,4/3);lmiterm([6,3,4,-Mz2],(5/6)*lbd,B1');lmiterm([6,3,4,-Mz2],(5/6)*lbd,B3');
lmiterm([6,3,5,-P50],4/3,1);
lmiterm([6,3,5,-R50],-tao/lbd,4/3);lmiterm([6,3,5,-Mz2],(5/6)*lbd,B1');lmiterm([6,3,5,-Mz2],(5/6)*lbd,B3');
lmiterm([6,3,6,P60],4/3,1);lmiterm([6,3,6,P61],5/6,1);lmiterm([6,3,6,P63],5/6,1);
lmiterm([6,3,6,R60],-tao/lbd,4/3);lmiterm([6,3,6,R61],-tao/lbd,5/6);lmiterm([6,3,6,R63],-tao/lbd,5/6);
lmiterm([6,3,6,-Mz2],4/3,1);lmiterm([6,3,6,-Mz2],-(5/6)*lbd,B1'*C1'*Ld3');lmiterm([6,3,6,-Mz2],-(5/6)*lbd,B3'*C1'*Ld1');

lmiterm([6,1,7,0],0);
lmiterm([6,1,8,0],0);
lmiterm([6,1,9,0],0);
lmiterm([6,2,7,0],0);
lmiterm([6,2,8,0],0);
lmiterm([6,2,9,0],0);
lmiterm([6,3,7,0],0);
lmiterm([6,3,8,0],0);
lmiterm([6,3,9,0],0);
lmiterm([6,1,10,0],0);
lmiterm([6,2,10,0],0);
lmiterm([6,3,10,0],0);
lmiterm([6,1,11,0],0);
lmiterm([6,1,12,0],0);
lmiterm([6,1,13,0],0);
lmiterm([6,2,11,0],0);
lmiterm([6,2,12,0],0);
lmiterm([6,2,13,0],0);
lmiterm([6,3,11,0],0);
lmiterm([6,3,12,0],0);
lmiterm([6,3,13,0],0);
lmiterm([6,1,14,-Mz1],4/3,1);
lmiterm([6,1,15,0],0);
lmiterm([6,1,16,0],0);
lmiterm([6,2,14,0],0);
lmiterm([6,2,15,0],0);
lmiterm([6,2,16,0],0);
lmiterm([6,3,14,0],0);
lmiterm([6,3,15,0],0);
lmiterm([6,3,16,0],0);

lmiterm([6,4,4,Mz1],-lbd,4/3,'s');
lmiterm([6,4,4,R10],tao,4/3);lmiterm([6,4,4,R11],tao,5/6);lmiterm([6,4,4,R13],tao,5/6);
lmiterm([6,4,5,R20],tao,4/3);
lmiterm([6,4,6,R30],tao,4/3);
lmiterm([6,5,5,Mz1],-lbd,4/3,'s');
lmiterm([6,5,5,R40],tao,4/3);lmiterm([6,5,5,R41],tao,5/6);lmiterm([6,5,5,R43],tao,5/6);
lmiterm([6,5,6,R50],tao,4/3);
lmiterm([6,6,6,Mz2],-lbd,4/3,'s');
lmiterm([6,6,6,R60],tao,4/3);lmiterm([6,6,6,R61],tao,5/6);lmiterm([6,6,6,R63],tao,5/6);

lmiterm([6,4,7,0],0);
lmiterm([6,4,8,0],0);
lmiterm([6,4,9,0],0);
lmiterm([6,5,7,0],0);
lmiterm([6,5,8,0],0);
lmiterm([6,5,9,0],0);
lmiterm([6,6,7,0],0);
lmiterm([6,6,8,0],0);
lmiterm([6,6,9,0],0);

lmiterm([6,4,10,0],0);
lmiterm([6,5,10,0],0);
lmiterm([6,6,10,0],(4/3)*lbd);
lmiterm([6,4,11,0],0);
lmiterm([6,4,12,0],0);
lmiterm([6,4,13,0],0);
lmiterm([6,5,11,0],0);lmiterm([6,5,12,-Mz1],-(5/6)*lbd,C1'*Ls1');lmiterm([6,5,12,-Mz1],-(5/6)*lbd,C1'*Ls3');lmiterm([6,5,13,-Mz1],-(5/6)*lbd,A1'*C1'*Ld3');lmiterm([6,5,13,-Mz1],-(5/6)*lbd,A3'*C1'*Ld1');
lmiterm([6,6,11,0],0);
lmiterm([6,6,12,0],0);lmiterm([6,6,13,-Mz2],-(5/6)*lbd,B1'*C1'*Ld3');lmiterm([6,6,13,-Mz2],-(5/6)*lbd,B3'*C1'*Ld1');

lmiterm([6,4,14,0],0);
lmiterm([6,4,15,0],0);
lmiterm([6,4,16,0],0);
lmiterm([6,5,14,0],0);
lmiterm([6,5,15,0],0);
lmiterm([6,5,16,0],0);
lmiterm([6,6,14,0],0);
lmiterm([6,6,15,0],0);
lmiterm([6,6,16,0],0);

lmiterm([6,7,7,Q10],-4/3,1);lmiterm([6,7,7,Q11],-5/6,1);lmiterm([6,7,7,Q13],-5/6,1);
lmiterm([6,7,8,Q20],-4/3,1);
lmiterm([6,7,9,Q30],-4/3,1);
lmiterm([6,8,8,Q40],-4/3,1);lmiterm([6,8,8,Q41],-5/6,1);lmiterm([6,8,8,Q43],-5/6,1);
lmiterm([6,8,9,Q50],-4/3,1);
lmiterm([6,9,9,Q60],-4/3,1);lmiterm([6,9,9,Q61],-5/6,1);lmiterm([6,9,9,Q63],-5/6,1);
lmiterm([6,7,10,0],0);
lmiterm([6,8,10,0],0);
lmiterm([6,9,10,0],0);
lmiterm([6,7,11,0],0);
lmiterm([6,7,12,0],0);
lmiterm([6,7,13,0],0);
lmiterm([6,8,11,0],0);
lmiterm([6,8,12,0],0);
lmiterm([6,8,13,0],0);
lmiterm([6,9,11,0],0);
lmiterm([6,9,12,0],0);
lmiterm([6,9,13,0],0);

lmiterm([6,7,14,0],0);
lmiterm([6,7,15,0],0);
lmiterm([6,7,16,0],0);
lmiterm([6,8,14,0],0);
lmiterm([6,8,15,0],0);
lmiterm([6,8,16,0],0);
lmiterm([6,9,14,0],0);
lmiterm([6,9,15,0],0);
lmiterm([6,9,16,0],0);

lmiterm([6,10,10,0],-(4/3)*gama1*gama1);
lmiterm([6,10,11,0],0);
lmiterm([6,10,12,0],0);
lmiterm([6,10,13,0],0);
lmiterm([6,10,14,0],0);
lmiterm([6,10,15,0],0);
lmiterm([6,10,16,0],0);

lmiterm([6,11,11,R10],-1/tao,4/3);lmiterm([6,11,11,R11],-1/tao,5/6);lmiterm([6,11,11,R13],-1/tao,5/6);
lmiterm([6,11,12,R20],-1/tao,4/3);
lmiterm([6,11,13,R30],-1/tao,4/3);
lmiterm([6,12,12,R40],-1/tao,4/3);lmiterm([6,12,12,R41],-1/tao,5/6);lmiterm([6,12,12,R43],-1/tao,5/6);
lmiterm([6,12,13,R50],-1/tao,4/3);
lmiterm([6,13,13,R60],-1/tao,4/3);lmiterm([6,13,13,R61],-1/tao,5/6);lmiterm([6,13,13,R63],-1/tao,5/6);

lmiterm([6,11,14,0],0);
lmiterm([6,11,15,0],0);
lmiterm([6,11,16,0],0);
lmiterm([6,12,14,0],0);
lmiterm([6,12,15,0],0);
lmiterm([6,12,16,0],0);
lmiterm([6,13,14,0],0);
lmiterm([6,13,15,0],0);
lmiterm([6,13,16,0],0);

lmiterm([6,14,14,0],-4/3);
lmiterm([6,14,15,0],0);
lmiterm([6,14,16,0],0);
lmiterm([6,15,15,0],-4/3);
lmiterm([6,15,16,0],0);
lmiterm([6,16,16,0],-4/3);



lmiterm([7,1,1,Q10],4/3,1);lmiterm([7,1,1,Q11],5/6,1);lmiterm([7,1,1,Q14],5/6,1);
lmiterm([7,1,1,R10],tao/(lbd*lbd),4/3);lmiterm([7,1,1,R11],tao/(lbd*lbd),5/6);lmiterm([7,1,1,R14],tao/(lbd*lbd),5/6);
lmiterm([7,1,1,P10],-2/lbd,4/3);lmiterm([7,1,1,P11],-2/lbd,5/6);lmiterm([7,1,1,P14],-2/lbd,5/6);
lmiterm([7,1,2,Q20],4/3,1);
lmiterm([7,1,2,R20],tao/(lbd*lbd),4/3);
lmiterm([7,1,2,P20],-2/lbd,4/3);
lmiterm([7,1,3,Q30],4/3,1);
lmiterm([7,1,3,R30],tao/(lbd*lbd),4/3);
lmiterm([7,1,3,P30],-2/lbd,4/3);
lmiterm([7,2,2,Q40],4/3,1);lmiterm([7,2,2,Q41],5/6,1);lmiterm([7,2,2,Q44],5/6,1);
lmiterm([7,2,2,R40],tao/(lbd*lbd),4/3);lmiterm([7,2,2,R41],tao/(lbd*lbd),5/6);lmiterm([7,2,2,R44],tao/(lbd*lbd),5/6);
lmiterm([7,2,2,P40],-2/lbd,4/3);lmiterm([7,2,2,P41],-2/lbd,5/6);lmiterm([7,2,2,P44],-2/lbd,5/6);
lmiterm([7,2,3,Q50],4/3,1);
lmiterm([7,2,3,R50],tao/(lbd*lbd),4/3);
lmiterm([7,2,3,P50],-2/lbd,4/3);
lmiterm([7,3,3,Q60],4/3,1);lmiterm([7,3,3,Q61],5/6,1);lmiterm([7,3,3,Q64],5/6,1);
lmiterm([7,3,3,R60],tao/(lbd*lbd),4/3);lmiterm([7,3,3,R61],tao/(lbd*lbd),5/6);lmiterm([7,3,3,R64],tao/(lbd*lbd),5/6);
lmiterm([7,3,3,P60],-2/lbd,4/3);lmiterm([7,3,3,P61],-2/lbd,5/6);lmiterm([7,3,3,P64],-2/lbd,5/6);
lmiterm([7,1,4,P10],4/3,1);lmiterm([7,1,4,P11],5/6,1);lmiterm([7,1,4,P14],5/6,1);
lmiterm([7,1,4,R10],-tao/lbd,4/3);lmiterm([7,1,4,R11],-tao/lbd,5/6);lmiterm([7,1,4,R14],-tao/lbd,5/6);
lmiterm([7,1,4,-Mz1],4/3,1);lmiterm([7,1,4,-Mz1],(5/6)*lbd,A1');lmiterm([7,1,4,-Mz1],(5/6)*lbd,A4');lmiterm([7,1,4,-Zg1],(5/6)*lbd,B4');lmiterm([7,1,4,-Zg4],(5/6)*lbd,B1');
lmiterm([7,1,5,P20],4/3,1); 
lmiterm([7,1,5,R20],-tao/lbd,4/3);
lmiterm([7,1,6,P30],4/3,1);
lmiterm([7,1,6,R30],-tao/lbd,4/3);
lmiterm([7,2,4,-P20],4/3,1); 
lmiterm([7,2,4,-R20],-tao/lbd,4/3); lmiterm([7,2,4,-Zg1],-(5/6)*lbd,B4');lmiterm([7,2,4,-Zg4],-(5/6)*lbd,B1');
lmiterm([7,2,5,P40],4/3,1);lmiterm([7,2,5,P41],5/6,1);lmiterm([7,2,5,P44],5/6,1);
lmiterm([7,2,5,R40],-tao/lbd,4/3);lmiterm([7,2,5,R41],-tao/lbd,5/6);lmiterm([7,2,5,R44],-tao/lbd,5/6);
lmiterm([7,2,5,-Mz1],4/3,1);lmiterm([7,2,5,-Mz1],(5/6)*lbd,A1');lmiterm([7,2,5,-Mz1],(5/6)*lbd,A4');lmiterm([7,2,5,-Mz1],-(5/6)*lbd,C1'*Ls1');lmiterm([7,2,5,-Mz1],-(5/6)*lbd,C1'*Ls4');
lmiterm([7,2,6,P50],4/3,1);
lmiterm([7,2,6,R50],-tao/lbd,4/3);lmiterm([7,2,6,-Mz1],-(5/6)*lbd,A1'*C1'*Ld4');lmiterm([7,2,6,-Mz1],-(5/6)*lbd,A4'*C1'*Ld1');
lmiterm([7,3,4,-P30],4/3,1);
lmiterm([7,3,4,-R30],-tao/lbd,4/3);lmiterm([7,3,4,-Mz2],(5/6)*lbd,B1');lmiterm([7,3,4,-Mz2],(5/6)*lbd,B4');
lmiterm([7,3,5,-P50],4/3,1);
lmiterm([7,3,5,-R50],-tao/lbd,4/3);lmiterm([7,3,5,-Mz2],(5/6)*lbd,B1');lmiterm([7,3,5,-Mz2],(5/6)*lbd,B4');
lmiterm([7,3,6,P60],4/3,1);lmiterm([7,3,6,P61],5/6,1);lmiterm([7,3,6,P64],5/6,1);
lmiterm([7,3,6,R60],-tao/lbd,4/3);lmiterm([7,3,6,R61],-tao/lbd,5/6);lmiterm([7,3,6,R64],-tao/lbd,5/6);
lmiterm([7,3,6,-Mz2],4/3,1);lmiterm([7,3,6,-Mz2],-(5/6)*lbd,B1'*C1'*Ld4');lmiterm([7,3,6,-Mz2],-(5/6)*lbd,B4'*C1'*Ld1');

lmiterm([7,1,7,0],0);
lmiterm([7,1,8,0],0);
lmiterm([7,1,9,0],0);
lmiterm([7,2,7,0],0);
lmiterm([7,2,8,0],0);
lmiterm([7,2,9,0],0);
lmiterm([7,3,7,0],0);
lmiterm([7,3,8,0],0);
lmiterm([7,3,9,0],0);
lmiterm([7,1,10,0],0);
lmiterm([7,2,10,0],0);
lmiterm([7,3,10,0],0);
lmiterm([7,1,11,0],0);
lmiterm([7,1,12,0],0);
lmiterm([7,1,13,0],0);
lmiterm([7,2,11,0],0);
lmiterm([7,2,12,0],0);
lmiterm([7,2,13,0],0);
lmiterm([7,3,11,0],0);
lmiterm([7,3,12,0],0);
lmiterm([7,3,13,0],0);
lmiterm([7,1,14,-Mz1],4/3,1);
lmiterm([7,1,15,0],0);
lmiterm([7,1,16,0],0);
lmiterm([7,2,14,0],0);
lmiterm([7,2,15,0],0);
lmiterm([7,2,16,0],0);
lmiterm([7,3,14,0],0);
lmiterm([7,3,15,0],0);
lmiterm([7,3,16,0],0);

lmiterm([7,4,4,Mz1],-lbd,4/3,'s');
lmiterm([7,4,4,R10],tao,4/3);lmiterm([7,4,4,R11],tao,5/6);lmiterm([7,4,4,R14],tao,5/6);
lmiterm([7,4,5,R20],tao,4/3);
lmiterm([7,4,6,R30],tao,4/3);
lmiterm([7,5,5,Mz1],-lbd,4/3,'s');
lmiterm([7,5,5,R40],tao,4/3);lmiterm([7,5,5,R41],tao,5/6);lmiterm([7,5,5,R44],tao,5/6);
lmiterm([7,5,6,R50],tao,4/3);
lmiterm([7,6,6,Mz2],-lbd,4/3,'s');
lmiterm([7,6,6,R60],tao,4/3);lmiterm([7,6,6,R61],tao,5/6);lmiterm([7,6,6,R64],tao,5/6);

lmiterm([7,4,7,0],0);
lmiterm([7,4,8,0],0);
lmiterm([7,4,9,0],0);
lmiterm([7,5,7,0],0);
lmiterm([7,5,8,0],0);
lmiterm([7,5,9,0],0);
lmiterm([7,6,7,0],0);
lmiterm([7,6,8,0],0);
lmiterm([7,6,9,0],0);

lmiterm([7,4,10,0],0);
lmiterm([7,5,10,0],0);
lmiterm([7,6,10,0],(4/3)*lbd);
lmiterm([7,4,11,0],0);
lmiterm([7,4,12,0],0);
lmiterm([7,4,13,0],0);
lmiterm([7,5,11,0],0);lmiterm([7,5,12,-Mz1],-(5/6)*lbd,C1'*Ls1');lmiterm([7,5,12,-Mz1],-(5/6)*lbd,C1'*Ls4');lmiterm([7,5,13,-Mz1],-(5/6)*lbd,A1'*C1'*Ld4');lmiterm([7,5,13,-Mz1],-(5/6)*lbd,A4'*C1'*Ld1');
lmiterm([7,6,11,0],0);
lmiterm([7,6,12,0],0);lmiterm([7,6,13,-Mz2],-(5/6)*lbd,B1'*C1'*Ld4');lmiterm([7,6,13,-Mz2],-(5/6)*lbd,B4'*C1'*Ld1');

lmiterm([7,4,14,0],0);
lmiterm([7,4,15,0],0);
lmiterm([7,4,16,0],0);
lmiterm([7,5,14,0],0);
lmiterm([7,5,15,0],0);
lmiterm([7,5,16,0],0);
lmiterm([7,6,14,0],0);
lmiterm([7,6,15,0],0);
lmiterm([7,6,16,0],0);

lmiterm([7,7,7,Q10],-4/3,1);lmiterm([7,7,7,Q11],-5/6,1);lmiterm([7,7,7,Q14],-5/6,1);
lmiterm([7,7,8,Q20],-4/3,1);
lmiterm([7,7,9,Q30],-4/3,1);
lmiterm([7,8,8,Q40],-4/3,1);lmiterm([7,8,8,Q41],-5/6,1);lmiterm([7,8,8,Q44],-5/6,1);
lmiterm([7,8,9,Q50],-4/3,1);
lmiterm([7,9,9,Q60],-4/3,1);lmiterm([7,9,9,Q61],-5/6,1);lmiterm([7,9,9,Q64],-5/6,1);
lmiterm([7,7,10,0],0);
lmiterm([7,8,10,0],0);
lmiterm([7,9,10,0],0);
lmiterm([7,7,11,0],0);
lmiterm([7,7,12,0],0);
lmiterm([7,7,13,0],0);
lmiterm([7,8,11,0],0);
lmiterm([7,8,12,0],0);
lmiterm([7,8,13,0],0);
lmiterm([7,9,11,0],0);
lmiterm([7,9,12,0],0);
lmiterm([7,9,13,0],0);

lmiterm([7,7,14,0],0);
lmiterm([7,7,15,0],0);
lmiterm([7,7,16,0],0);
lmiterm([7,8,14,0],0);
lmiterm([7,8,15,0],0);
lmiterm([7,8,16,0],0);
lmiterm([7,9,14,0],0);
lmiterm([7,9,15,0],0);
lmiterm([7,9,16,0],0);

lmiterm([7,10,10,0],-(4/3)*gama1*gama1);
lmiterm([7,10,11,0],0);
lmiterm([7,10,12,0],0);
lmiterm([7,10,13,0],0);
lmiterm([7,10,14,0],0);
lmiterm([7,10,15,0],0);
lmiterm([7,10,16,0],0);

lmiterm([7,11,11,R10],-1/tao,4/3);lmiterm([7,11,11,R11],-1/tao,5/6);lmiterm([7,11,11,R14],-1/tao,5/6);
lmiterm([7,11,12,R20],-1/tao,4/3);
lmiterm([7,11,13,R30],-1/tao,4/3);
lmiterm([7,12,12,R40],-1/tao,4/3);lmiterm([7,12,12,R41],-1/tao,5/6);lmiterm([7,12,12,R44],-1/tao,5/6);
lmiterm([7,12,13,R50],-1/tao,4/3);
lmiterm([7,13,13,R60],-1/tao,4/3);lmiterm([7,13,13,R61],-1/tao,5/6);lmiterm([7,13,13,R64],-1/tao,5/6);

lmiterm([7,11,14,0],0);
lmiterm([7,11,15,0],0);
lmiterm([7,11,16,0],0);
lmiterm([7,12,14,0],0);
lmiterm([7,12,15,0],0);
lmiterm([7,12,16,0],0);
lmiterm([7,13,14,0],0);
lmiterm([7,13,15,0],0);
lmiterm([7,13,16,0],0);

lmiterm([7,14,14,0],-4/3);
lmiterm([7,14,15,0],0);
lmiterm([7,14,16,0],0);
lmiterm([7,15,15,0],-4/3);
lmiterm([7,15,16,0],0);
lmiterm([7,16,16,0],-4/3);




lmiterm([8,1,1,Q10],4/3,1);lmiterm([8,1,1,Q12],5/6,1);lmiterm([8,1,1,Q13],5/6,1);
lmiterm([8,1,1,R10],tao/(lbd*lbd),4/3);lmiterm([8,1,1,R12],tao/(lbd*lbd),5/6);lmiterm([8,1,1,R13],tao/(lbd*lbd),5/6);
lmiterm([8,1,1,P10],-2/lbd,4/3);lmiterm([8,1,1,P12],-2/lbd,5/6);lmiterm([8,1,1,P13],-2/lbd,5/6);
lmiterm([8,1,2,Q20],4/3,1);
lmiterm([8,1,2,R20],tao/(lbd*lbd),4/3);
lmiterm([8,1,2,P20],-2/lbd,4/3);
lmiterm([8,1,3,Q30],4/3,1);
lmiterm([8,1,3,R30],tao/(lbd*lbd),4/3);
lmiterm([8,1,3,P30],-2/lbd,4/3);
lmiterm([8,2,2,Q40],4/3,1);lmiterm([8,2,2,Q42],5/6,1);lmiterm([8,2,2,Q43],5/6,1);
lmiterm([8,2,2,R40],tao/(lbd*lbd),4/3);lmiterm([8,2,2,R42],tao/(lbd*lbd),5/6);lmiterm([8,2,2,R43],tao/(lbd*lbd),5/6);
lmiterm([8,2,2,P40],-2/lbd,4/3);lmiterm([8,2,2,P42],-2/lbd,5/6);lmiterm([8,2,2,P43],-2/lbd,5/6);
lmiterm([8,2,3,Q50],4/3,1);
lmiterm([8,2,3,R50],tao/(lbd*lbd),4/3);
lmiterm([8,2,3,P50],-2/lbd,4/3);
lmiterm([8,3,3,Q60],4/3,1);lmiterm([8,3,3,Q62],5/6,1);lmiterm([8,3,3,Q63],5/6,1);
lmiterm([8,3,3,R60],tao/(lbd*lbd),4/3);lmiterm([8,3,3,R62],tao/(lbd*lbd),5/6);lmiterm([8,3,3,R63],tao/(lbd*lbd),5/6);
lmiterm([8,3,3,P60],-2/lbd,4/3);lmiterm([8,3,3,P62],-2/lbd,5/6);lmiterm([8,3,3,P63],-2/lbd,5/6);
lmiterm([8,1,4,P10],4/3,1);lmiterm([8,1,4,P12],5/6,1);lmiterm([8,1,4,P13],5/6,1);
lmiterm([8,1,4,R10],-tao/lbd,4/3);lmiterm([8,1,4,R12],-tao/lbd,5/6);lmiterm([8,1,4,R13],-tao/lbd,5/6);
lmiterm([8,1,4,-Mz1],4/3,1);lmiterm([8,1,4,-Mz1],(5/6)*lbd,A2');lmiterm([8,1,4,-Mz1],(5/6)*lbd,A3');lmiterm([8,1,4,-Zg2],(5/6)*lbd,B3');lmiterm([8,1,4,-Zg3],(5/6)*lbd,B2');
lmiterm([8,1,5,P20],4/3,1); 
lmiterm([8,1,5,R20],-tao/lbd,4/3);
lmiterm([8,1,6,P30],4/3,1);
lmiterm([8,1,6,R30],-tao/lbd,4/3);
lmiterm([8,2,4,-P20],4/3,1); 
lmiterm([8,2,4,-R20],-tao/lbd,4/3); lmiterm([8,2,4,-Zg2],-(5/6)*lbd,B3');lmiterm([8,2,4,-Zg3],-(5/6)*lbd,B2');
lmiterm([8,2,5,P40],4/3,1);lmiterm([8,2,5,P42],5/6,1);lmiterm([8,2,5,P43],5/6,1);
lmiterm([8,2,5,R40],-tao/lbd,4/3);lmiterm([8,2,5,R42],-tao/lbd,5/6);lmiterm([8,2,5,R43],-tao/lbd,5/6);
lmiterm([8,2,5,-Mz1],4/3,1);lmiterm([8,2,5,-Mz1],(5/6)*lbd,A2');lmiterm([8,2,5,-Mz1],(5/6)*lbd,A3');lmiterm([8,2,5,-Mz1],-(5/6)*lbd,C1'*Ls2');lmiterm([8,2,5,-Mz1],-(5/6)*lbd,C1'*Ls3');
lmiterm([8,2,6,P50],4/3,1);
lmiterm([8,2,6,R50],-tao/lbd,4/3);lmiterm([8,2,6,-Mz1],-(5/6)*lbd,A2'*C1'*Ld3');lmiterm([8,2,6,-Mz1],-(5/6)*lbd,A3'*C1'*Ld2');
lmiterm([8,3,4,-P30],4/3,1);
lmiterm([8,3,4,-R30],-tao/lbd,4/3);lmiterm([8,3,4,-Mz2],(5/6)*lbd,B2');lmiterm([8,3,4,-Mz2],(5/6)*lbd,B3');
lmiterm([8,3,5,-P50],4/3,1);
lmiterm([8,3,5,-R50],-tao/lbd,4/3);lmiterm([8,3,5,-Mz2],(5/6)*lbd,B2');lmiterm([8,3,5,-Mz2],(5/6)*lbd,B3');
lmiterm([8,3,6,P60],4/3,1);lmiterm([8,3,6,P62],5/6,1);lmiterm([8,3,6,P63],5/6,1);
lmiterm([8,3,6,R60],-tao/lbd,4/3);lmiterm([8,3,6,R62],-tao/lbd,5/6);lmiterm([8,3,6,R63],-tao/lbd,5/6);
lmiterm([8,3,6,-Mz2],4/3,1);lmiterm([8,3,6,-Mz2],-(5/6)*lbd,B2'*C1'*Ld3');lmiterm([8,3,6,-Mz2],-(5/6)*lbd,B3'*C1'*Ld2');

lmiterm([8,1,7,0],0);
lmiterm([8,1,8,0],0);
lmiterm([8,1,9,0],0);
lmiterm([8,2,7,0],0);
lmiterm([8,2,8,0],0);
lmiterm([8,2,9,0],0);
lmiterm([8,3,7,0],0);
lmiterm([8,3,8,0],0);
lmiterm([8,3,9,0],0);
lmiterm([8,1,10,0],0);
lmiterm([8,2,10,0],0);
lmiterm([8,3,10,0],0);
lmiterm([8,1,11,0],0);
lmiterm([8,1,12,0],0);
lmiterm([8,1,13,0],0);
lmiterm([8,2,11,0],0);
lmiterm([8,2,12,0],0);
lmiterm([8,2,13,0],0);
lmiterm([8,3,11,0],0);
lmiterm([8,3,12,0],0);
lmiterm([8,3,13,0],0);
lmiterm([8,1,14,-Mz1],4/3,1);
lmiterm([8,1,15,0],0);
lmiterm([8,1,16,0],0);
lmiterm([8,2,14,0],0);
lmiterm([8,2,15,0],0);
lmiterm([8,2,16,0],0);
lmiterm([8,3,14,0],0);
lmiterm([8,3,15,0],0);
lmiterm([8,3,16,0],0);

lmiterm([8,4,4,Mz1],-lbd,4/3,'s');
lmiterm([8,4,4,R10],tao,4/3);lmiterm([8,4,4,R12],tao,5/6);lmiterm([8,4,4,R13],tao,5/6);
lmiterm([8,4,5,R20],tao,4/3);
lmiterm([8,4,6,R30],tao,4/3);
lmiterm([8,5,5,Mz1],-lbd,4/3,'s');
lmiterm([8,5,5,R40],tao,4/3);lmiterm([8,5,5,R42],tao,5/6);lmiterm([8,5,5,R43],tao,5/6);
lmiterm([8,5,6,R50],tao,4/3);
lmiterm([8,6,6,Mz2],-lbd,4/3,'s');
lmiterm([8,6,6,R60],tao,4/3);lmiterm([8,6,6,R62],tao,5/6);lmiterm([8,6,6,R63],tao,5/6);

lmiterm([8,4,7,0],0);
lmiterm([8,4,8,0],0);
lmiterm([8,4,9,0],0);
lmiterm([8,5,7,0],0);
lmiterm([8,5,8,0],0);
lmiterm([8,5,9,0],0);
lmiterm([8,6,7,0],0);
lmiterm([8,6,8,0],0);
lmiterm([8,6,9,0],0);

lmiterm([8,4,10,0],0);
lmiterm([8,5,10,0],0);
lmiterm([8,6,10,0],(4/3)*lbd);
lmiterm([8,4,11,0],0);
lmiterm([8,4,12,0],0);
lmiterm([8,4,13,0],0);
lmiterm([8,5,11,0],0);lmiterm([8,5,12,-Mz1],-(5/6)*lbd,C1'*Ls2');lmiterm([8,5,12,-Mz1],-(5/6)*lbd,C1'*Ls3');lmiterm([8,5,13,-Mz1],-(5/6)*lbd,A2'*C1'*Ld3');lmiterm([8,5,13,-Mz1],-(5/6)*lbd,A3'*C1'*Ld2');
lmiterm([8,6,11,0],0);
lmiterm([8,6,12,0],0);lmiterm([8,6,13,-Mz2],-(5/6)*lbd,B2'*C1'*Ld3');lmiterm([8,6,13,-Mz2],-(5/6)*lbd,B3'*C1'*Ld2');

lmiterm([8,4,14,0],0);
lmiterm([8,4,15,0],0);
lmiterm([8,4,16,0],0);
lmiterm([8,5,14,0],0);
lmiterm([8,5,15,0],0);
lmiterm([8,5,16,0],0);
lmiterm([8,6,14,0],0);
lmiterm([8,6,15,0],0);
lmiterm([8,6,16,0],0);

lmiterm([8,7,7,Q10],-4/3,1);lmiterm([8,7,7,Q12],-5/6,1);lmiterm([8,7,7,Q13],-5/6,1);
lmiterm([8,7,8,Q20],-4/3,1);
lmiterm([8,7,9,Q30],-4/3,1);
lmiterm([8,8,8,Q40],-4/3,1);lmiterm([8,8,8,Q42],-5/6,1);lmiterm([8,8,8,Q43],-5/6,1);
lmiterm([8,8,9,Q50],-4/3,1);
lmiterm([8,9,9,Q60],-4/3,1);lmiterm([8,9,9,Q62],-5/6,1);lmiterm([8,9,9,Q63],-5/6,1);
lmiterm([8,7,10,0],0);
lmiterm([8,8,10,0],0);
lmiterm([8,9,10,0],0);
lmiterm([8,7,11,0],0);
lmiterm([8,7,12,0],0);
lmiterm([8,7,13,0],0);
lmiterm([8,8,11,0],0);
lmiterm([8,8,12,0],0);
lmiterm([8,8,13,0],0);
lmiterm([8,9,11,0],0);
lmiterm([8,9,12,0],0);
lmiterm([8,9,13,0],0);

lmiterm([8,7,14,0],0);
lmiterm([8,7,15,0],0);
lmiterm([8,7,16,0],0);
lmiterm([8,8,14,0],0);
lmiterm([8,8,15,0],0);
lmiterm([8,8,16,0],0);
lmiterm([8,9,14,0],0);
lmiterm([8,9,15,0],0);
lmiterm([8,9,16,0],0);

lmiterm([8,10,10,0],-(4/3)*gama1*gama1);
lmiterm([8,10,11,0],0);
lmiterm([8,10,12,0],0);
lmiterm([8,10,13,0],0);
lmiterm([8,10,14,0],0);
lmiterm([8,10,15,0],0);
lmiterm([8,10,16,0],0);

lmiterm([8,11,11,R10],-1/tao,4/3);lmiterm([8,11,11,R12],-1/tao,5/6);lmiterm([8,11,11,R13],-1/tao,5/6);
lmiterm([8,11,12,R20],-1/tao,4/3);
lmiterm([8,11,13,R30],-1/tao,4/3);
lmiterm([8,12,12,R40],-1/tao,4/3);lmiterm([8,12,12,R42],-1/tao,5/6);lmiterm([8,12,12,R43],-1/tao,5/6);
lmiterm([8,12,13,R50],-1/tao,4/3);
lmiterm([8,13,13,R60],-1/tao,4/3);lmiterm([8,13,13,R62],-1/tao,5/6);lmiterm([8,13,13,R63],-1/tao,5/6);

lmiterm([8,11,14,0],0);
lmiterm([8,11,15,0],0);
lmiterm([8,11,16,0],0);
lmiterm([8,12,14,0],0);
lmiterm([8,12,15,0],0);
lmiterm([8,12,16,0],0);
lmiterm([8,13,14,0],0);
lmiterm([8,13,15,0],0);
lmiterm([8,13,16,0],0);

lmiterm([8,14,14,0],-4/3);
lmiterm([8,14,15,0],0);
lmiterm([8,14,16,0],0);
lmiterm([8,15,15,0],-4/3);
lmiterm([8,15,16,0],0);
lmiterm([8,16,16,0],-4/3);




lmiterm([9,1,1,Q10],4/3,1);lmiterm([9,1,1,Q12],5/6,1);lmiterm([9,1,1,Q14],5/6,1);
lmiterm([9,1,1,R10],tao/(lbd*lbd),4/3);lmiterm([9,1,1,R12],tao/(lbd*lbd),5/6);lmiterm([9,1,1,R14],tao/(lbd*lbd),5/6);
lmiterm([9,1,1,P10],-2/lbd,4/3);lmiterm([9,1,1,P12],-2/lbd,5/6);lmiterm([9,1,1,P14],-2/lbd,5/6);
lmiterm([9,1,2,Q20],4/3,1);
lmiterm([9,1,2,R20],tao/(lbd*lbd),4/3);
lmiterm([9,1,2,P20],-2/lbd,4/3);
lmiterm([9,1,3,Q30],4/3,1);
lmiterm([9,1,3,R30],tao/(lbd*lbd),4/3);
lmiterm([9,1,3,P30],-2/lbd,4/3);
lmiterm([9,2,2,Q40],4/3,1);lmiterm([9,2,2,Q42],5/6,1);lmiterm([9,2,2,Q44],5/6,1);
lmiterm([9,2,2,R40],tao/(lbd*lbd),4/3);lmiterm([9,2,2,R42],tao/(lbd*lbd),5/6);lmiterm([9,2,2,R44],tao/(lbd*lbd),5/6);
lmiterm([9,2,2,P40],-2/lbd,4/3);lmiterm([9,2,2,P42],-2/lbd,5/6);lmiterm([9,2,2,P44],-2/lbd,5/6);
lmiterm([9,2,3,Q50],4/3,1);
lmiterm([9,2,3,R50],tao/(lbd*lbd),4/3);
lmiterm([9,2,3,P50],-2/lbd,4/3);
lmiterm([9,3,3,Q60],4/3,1);lmiterm([9,3,3,Q62],5/6,1);lmiterm([9,3,3,Q64],5/6,1);
lmiterm([9,3,3,R60],tao/(lbd*lbd),4/3);lmiterm([9,3,3,R62],tao/(lbd*lbd),5/6);lmiterm([9,3,3,R64],tao/(lbd*lbd),5/6);
lmiterm([9,3,3,P60],-2/lbd,4/3);lmiterm([9,3,3,P62],-2/lbd,5/6);lmiterm([9,3,3,P64],-2/lbd,5/6);
lmiterm([9,1,4,P10],4/3,1);lmiterm([9,1,4,P12],5/6,1);lmiterm([9,1,4,P14],5/6,1);
lmiterm([9,1,4,R10],-tao/lbd,4/3);lmiterm([9,1,4,R12],-tao/lbd,5/6);lmiterm([9,1,4,R14],-tao/lbd,5/6);
lmiterm([9,1,4,-Mz1],4/3,1);lmiterm([9,1,4,-Mz1],(5/6)*lbd,A2');lmiterm([9,1,4,-Mz1],(5/6)*lbd,A4');lmiterm([9,1,4,-Zg2],(5/6)*lbd,B4');lmiterm([9,1,4,-Zg4],(5/6)*lbd,B2');
lmiterm([9,1,5,P20],4/3,1); 
lmiterm([9,1,5,R20],-tao/lbd,4/3);
lmiterm([9,1,6,P30],4/3,1);
lmiterm([9,1,6,R30],-tao/lbd,4/3);
lmiterm([9,2,4,-P20],4/3,1); 
lmiterm([9,2,4,-R20],-tao/lbd,4/3); lmiterm([9,2,4,-Zg2],-(5/6)*lbd,B4');lmiterm([9,2,4,-Zg4],-(5/6)*lbd,B2');
lmiterm([9,2,5,P40],4/3,1);lmiterm([9,2,5,P42],5/6,1);lmiterm([9,2,5,P44],5/6,1);
lmiterm([9,2,5,R40],-tao/lbd,4/3);lmiterm([9,2,5,R42],-tao/lbd,5/6);lmiterm([9,2,5,R44],-tao/lbd,5/6);
lmiterm([9,2,5,-Mz1],4/3,1);lmiterm([9,2,5,-Mz1],(5/6)*lbd,A2');lmiterm([9,2,5,-Mz1],(5/6)*lbd,A4');lmiterm([9,2,5,-Mz1],-(5/6)*lbd,C1'*Ls2');lmiterm([9,2,5,-Mz1],-(5/6)*lbd,C1'*Ls4');
lmiterm([9,2,6,P50],4/3,1);
lmiterm([9,2,6,R50],-tao/lbd,4/3);lmiterm([9,2,6,-Mz1],-(5/6)*lbd,A2'*C1'*Ld4');lmiterm([9,2,6,-Mz1],-(5/6)*lbd,A4'*C1'*Ld2');
lmiterm([9,3,4,-P30],4/3,1);
lmiterm([9,3,4,-R30],-tao/lbd,4/3);lmiterm([9,3,4,-Mz2],(5/6)*lbd,B2');lmiterm([9,3,4,-Mz2],(5/6)*lbd,B4');
lmiterm([9,3,5,-P50],4/3,1);
lmiterm([9,3,5,-R50],-tao/lbd,4/3);lmiterm([9,3,5,-Mz2],(5/6)*lbd,B2');lmiterm([9,3,5,-Mz2],(5/6)*lbd,B4');
lmiterm([9,3,6,P60],4/3,1);lmiterm([9,3,6,P62],5/6,1);lmiterm([9,3,6,P64],5/6,1);
lmiterm([9,3,6,R60],-tao/lbd,4/3);lmiterm([9,3,6,R62],-tao/lbd,5/6);lmiterm([9,3,6,R64],-tao/lbd,5/6);
lmiterm([9,3,6,-Mz2],4/3,1);lmiterm([9,3,6,-Mz2],-(5/6)*lbd,B2'*C1'*Ld4');lmiterm([9,3,6,-Mz2],-(5/6)*lbd,B4'*C1'*Ld2');

lmiterm([9,1,7,0],0);
lmiterm([9,1,8,0],0);
lmiterm([9,1,9,0],0);
lmiterm([9,2,7,0],0);
lmiterm([9,2,8,0],0);
lmiterm([9,2,9,0],0);
lmiterm([9,3,7,0],0);
lmiterm([9,3,8,0],0);
lmiterm([9,3,9,0],0);
lmiterm([9,1,10,0],0);
lmiterm([9,2,10,0],0);
lmiterm([9,3,10,0],0);
lmiterm([9,1,11,0],0);
lmiterm([9,1,12,0],0);
lmiterm([9,1,13,0],0);
lmiterm([9,2,11,0],0);
lmiterm([9,2,12,0],0);
lmiterm([9,2,13,0],0);
lmiterm([9,3,11,0],0);
lmiterm([9,3,12,0],0);
lmiterm([9,3,13,0],0);
lmiterm([9,1,14,-Mz1],4/3,1);
lmiterm([9,1,15,0],0);
lmiterm([9,1,16,0],0);
lmiterm([9,2,14,0],0);
lmiterm([9,2,15,0],0);
lmiterm([9,2,16,0],0);
lmiterm([9,3,14,0],0);
lmiterm([9,3,15,0],0);
lmiterm([9,3,16,0],0);

lmiterm([9,4,4,Mz1],-lbd,4/3,'s');
lmiterm([9,4,4,R10],tao,4/3);lmiterm([9,4,4,R12],tao,5/6);lmiterm([9,4,4,R14],tao,5/6);
lmiterm([9,4,5,R20],tao,4/3);
lmiterm([9,4,6,R30],tao,4/3);
lmiterm([9,5,5,Mz1],-lbd,4/3,'s');
lmiterm([9,5,5,R40],tao,4/3);lmiterm([9,5,5,R42],tao,5/6);lmiterm([9,5,5,R44],tao,5/6);
lmiterm([9,5,6,R50],tao,4/3);
lmiterm([9,6,6,Mz2],-lbd,4/3,'s');
lmiterm([9,6,6,R60],tao,4/3);lmiterm([9,6,6,R62],tao,5/6);lmiterm([9,6,6,R64],tao,5/6);

lmiterm([9,4,7,0],0);
lmiterm([9,4,8,0],0);
lmiterm([9,4,9,0],0);
lmiterm([9,5,7,0],0);
lmiterm([9,5,8,0],0);
lmiterm([9,5,9,0],0);
lmiterm([9,6,7,0],0);
lmiterm([9,6,8,0],0);
lmiterm([9,6,9,0],0);

lmiterm([9,4,10,0],0);
lmiterm([9,5,10,0],0);
lmiterm([9,6,10,0],(4/3)*lbd);
lmiterm([9,4,11,0],0);
lmiterm([9,4,12,0],0);
lmiterm([9,4,13,0],0);
lmiterm([9,5,11,0],0);lmiterm([9,5,12,-Mz1],-(5/6)*lbd,C1'*Ls2');lmiterm([9,5,12,-Mz1],-(5/6)*lbd,C1'*Ls4');lmiterm([9,5,13,-Mz1],-(5/6)*lbd,A2'*C1'*Ld4');lmiterm([9,5,13,-Mz1],-(5/6)*lbd,A4'*C1'*Ld2');
lmiterm([9,6,11,0],0);
lmiterm([9,6,12,0],0);lmiterm([9,6,13,-Mz2],-(5/6)*lbd,B2'*C1'*Ld4');lmiterm([9,6,13,-Mz2],-(5/6)*lbd,B4'*C1'*Ld2');

lmiterm([9,4,14,0],0);
lmiterm([9,4,15,0],0);
lmiterm([9,4,16,0],0);
lmiterm([9,5,14,0],0);
lmiterm([9,5,15,0],0);
lmiterm([9,5,16,0],0);
lmiterm([9,6,14,0],0);
lmiterm([9,6,15,0],0);
lmiterm([9,6,16,0],0);

lmiterm([9,7,7,Q10],-4/3,1);lmiterm([9,7,7,Q12],-5/6,1);lmiterm([9,7,7,Q14],-5/6,1);
lmiterm([9,7,8,Q20],-4/3,1);
lmiterm([9,7,9,Q30],-4/3,1);
lmiterm([9,8,8,Q40],-4/3,1);lmiterm([9,8,8,Q42],-5/6,1);lmiterm([9,8,8,Q44],-5/6,1);
lmiterm([9,8,9,Q50],-4/3,1);
lmiterm([9,9,9,Q60],-4/3,1);lmiterm([9,9,9,Q62],-5/6,1);lmiterm([9,9,9,Q64],-5/6,1);
lmiterm([9,7,10,0],0);
lmiterm([9,8,10,0],0);
lmiterm([9,9,10,0],0);
lmiterm([9,7,11,0],0);
lmiterm([9,7,12,0],0);
lmiterm([9,7,13,0],0);
lmiterm([9,8,11,0],0);
lmiterm([9,8,12,0],0);
lmiterm([9,8,13,0],0);
lmiterm([9,9,11,0],0);
lmiterm([9,9,12,0],0);
lmiterm([9,9,13,0],0);

lmiterm([9,7,14,0],0);
lmiterm([9,7,15,0],0);
lmiterm([9,7,16,0],0);
lmiterm([9,8,14,0],0);
lmiterm([9,8,15,0],0);
lmiterm([9,8,16,0],0);
lmiterm([9,9,14,0],0);
lmiterm([9,9,15,0],0);
lmiterm([9,9,16,0],0);

lmiterm([9,10,10,0],-(4/3)*gama1*gama1);
lmiterm([9,10,11,0],0);
lmiterm([9,10,12,0],0);
lmiterm([9,10,13,0],0);
lmiterm([9,10,14,0],0);
lmiterm([9,10,15,0],0);
lmiterm([9,10,16,0],0);

lmiterm([9,11,11,R10],-1/tao,4/3);lmiterm([9,11,11,R12],-1/tao,5/6);lmiterm([9,11,11,R14],-1/tao,5/6);
lmiterm([9,11,12,R20],-1/tao,4/3);
lmiterm([9,11,13,R30],-1/tao,4/3);
lmiterm([9,12,12,R40],-1/tao,4/3);lmiterm([9,12,12,R42],-1/tao,5/6);lmiterm([9,12,12,R44],-1/tao,5/6);
lmiterm([9,12,13,R50],-1/tao,4/3);
lmiterm([9,13,13,R60],-1/tao,4/3);lmiterm([9,13,13,R62],-1/tao,5/6);lmiterm([9,13,13,R64],-1/tao,5/6);

lmiterm([9,11,14,0],0);
lmiterm([9,11,15,0],0);
lmiterm([9,11,16,0],0);
lmiterm([9,12,14,0],0);
lmiterm([9,12,15,0],0);
lmiterm([9,12,16,0],0);
lmiterm([9,13,14,0],0);
lmiterm([9,13,15,0],0);
lmiterm([9,13,16,0],0);

lmiterm([9,14,14,0],-4/3);
lmiterm([9,14,15,0],0);
lmiterm([9,14,16,0],0);
lmiterm([9,15,15,0],-4/3);
lmiterm([9,15,16,0],0);
lmiterm([9,16,16,0],-4/3);




lmiterm([10,1,1,Q10],4/3,1);lmiterm([10,1,1,Q13],5/6,1);lmiterm([10,1,1,Q14],5/6,1);
lmiterm([10,1,1,R10],tao/(lbd*lbd),4/3);lmiterm([10,1,1,R13],tao/(lbd*lbd),5/6);lmiterm([10,1,1,R14],tao/(lbd*lbd),5/6);
lmiterm([10,1,1,P10],-2/lbd,4/3);lmiterm([10,1,1,P13],-2/lbd,5/6);lmiterm([10,1,1,P14],-2/lbd,5/6);
lmiterm([10,1,2,Q20],4/3,1);
lmiterm([10,1,2,R20],tao/(lbd*lbd),4/3);
lmiterm([10,1,2,P20],-2/lbd,4/3);
lmiterm([10,1,3,Q30],4/3,1);
lmiterm([10,1,3,R30],tao/(lbd*lbd),4/3);
lmiterm([10,1,3,P30],-2/lbd,4/3);
lmiterm([10,2,2,Q40],4/3,1);lmiterm([10,2,2,Q43],5/6,1);lmiterm([10,2,2,Q44],5/6,1);
lmiterm([10,2,2,R40],tao/(lbd*lbd),4/3);lmiterm([10,2,2,R43],tao/(lbd*lbd),5/6);lmiterm([10,2,2,R44],tao/(lbd*lbd),5/6);
lmiterm([10,2,2,P40],-2/lbd,4/3);lmiterm([10,2,2,P43],-2/lbd,5/6);lmiterm([10,2,2,P44],-2/lbd,5/6);
lmiterm([10,2,3,Q50],4/3,1);
lmiterm([10,2,3,R50],tao/(lbd*lbd),4/3);
lmiterm([10,2,3,P50],-2/lbd,4/3);
lmiterm([10,3,3,Q60],4/3,1);lmiterm([10,3,3,Q63],5/6,1);lmiterm([10,3,3,Q64],5/6,1);
lmiterm([10,3,3,R60],tao/(lbd*lbd),4/3);lmiterm([10,3,3,R63],tao/(lbd*lbd),5/6);lmiterm([10,3,3,R64],tao/(lbd*lbd),5/6);
lmiterm([10,3,3,P60],-2/lbd,4/3);lmiterm([10,3,3,P63],-2/lbd,5/6);lmiterm([10,3,3,P64],-2/lbd,5/6);
lmiterm([10,1,4,P10],4/3,1);lmiterm([10,1,4,P13],5/6,1);lmiterm([10,1,4,P14],5/6,1);
lmiterm([10,1,4,R10],-tao/lbd,4/3);lmiterm([10,1,4,R13],-tao/lbd,5/6);lmiterm([10,1,4,R14],-tao/lbd,5/6);
lmiterm([10,1,4,-Mz1],4/3,1);lmiterm([10,1,4,-Mz1],(5/6)*lbd,A3');lmiterm([10,1,4,-Mz1],(5/6)*lbd,A4');lmiterm([10,1,4,-Zg3],(5/6)*lbd,B4');lmiterm([10,1,4,-Zg4],(5/6)*lbd,B3');
lmiterm([10,1,5,P20],4/3,1); 
lmiterm([10,1,5,R20],-tao/lbd,4/3);
lmiterm([10,1,6,P30],4/3,1);
lmiterm([10,1,6,R30],-tao/lbd,4/3);
lmiterm([10,2,4,-P20],4/3,1); 
lmiterm([10,2,4,-R20],-tao/lbd,4/3); lmiterm([10,2,4,-Zg3],-(5/6)*lbd,B4');lmiterm([10,2,4,-Zg4],-(5/6)*lbd,B3');
lmiterm([10,2,5,P40],4/3,1);lmiterm([10,2,5,P43],5/6,1);lmiterm([10,2,5,P44],5/6,1);
lmiterm([10,2,5,R40],-tao/lbd,4/3);lmiterm([10,2,5,R43],-tao/lbd,5/6);lmiterm([10,2,5,R44],-tao/lbd,5/6);
lmiterm([10,2,5,-Mz1],4/3,1);lmiterm([10,2,5,-Mz1],(5/6)*lbd,A3');lmiterm([10,2,5,-Mz1],(5/6)*lbd,A4');lmiterm([10,2,5,-Mz1],-(5/6)*lbd,C1'*Ls3');lmiterm([10,2,5,-Mz1],-(5/6)*lbd,C1'*Ls4');
lmiterm([10,2,6,P50],4/3,1);
lmiterm([10,2,6,R50],-tao/lbd,4/3);lmiterm([10,2,6,-Mz1],-(5/6)*lbd,A3'*C1'*Ld4');lmiterm([10,2,6,-Mz1],-(5/6)*lbd,A4'*C1'*Ld3');
lmiterm([10,3,4,-P30],4/3,1);
lmiterm([10,3,4,-R30],-tao/lbd,4/3);lmiterm([10,3,4,-Mz2],(5/6)*lbd,B3');lmiterm([10,3,4,-Mz2],(5/6)*lbd,B4');
lmiterm([10,3,5,-P50],4/3,1);
lmiterm([10,3,5,-R50],-tao/lbd,4/3);lmiterm([10,3,5,-Mz2],(5/6)*lbd,B3');lmiterm([10,3,5,-Mz2],(5/6)*lbd,B4');
lmiterm([10,3,6,P60],4/3,1);lmiterm([10,3,6,P63],5/6,1);lmiterm([10,3,6,P64],5/6,1);
lmiterm([10,3,6,R60],-tao/lbd,4/3);lmiterm([10,3,6,R63],-tao/lbd,5/6);lmiterm([10,3,6,R64],-tao/lbd,5/6);
lmiterm([10,3,6,-Mz2],4/3,1);lmiterm([10,3,6,-Mz2],-(5/6)*lbd,B3'*C1'*Ld4');lmiterm([10,3,6,-Mz2],-(5/6)*lbd,B4'*C1'*Ld3');

lmiterm([10,1,7,0],0);
lmiterm([10,1,8,0],0);
lmiterm([10,1,9,0],0);
lmiterm([10,2,7,0],0);
lmiterm([10,2,8,0],0);
lmiterm([10,2,9,0],0);
lmiterm([10,3,7,0],0);
lmiterm([10,3,8,0],0);
lmiterm([10,3,9,0],0);
lmiterm([10,1,10,0],0);
lmiterm([10,2,10,0],0);
lmiterm([10,3,10,0],0);
lmiterm([10,1,11,0],0);
lmiterm([10,1,12,0],0);
lmiterm([10,1,13,0],0);
lmiterm([10,2,11,0],0);
lmiterm([10,2,12,0],0);
lmiterm([10,2,13,0],0);
lmiterm([10,3,11,0],0);
lmiterm([10,3,12,0],0);
lmiterm([10,3,13,0],0);
lmiterm([10,1,14,-Mz1],4/3,1);
lmiterm([10,1,15,0],0);
lmiterm([10,1,16,0],0);
lmiterm([10,2,14,0],0);
lmiterm([10,2,15,0],0);
lmiterm([10,2,16,0],0);
lmiterm([10,3,14,0],0);
lmiterm([10,3,15,0],0);
lmiterm([10,3,16,0],0);

lmiterm([10,4,4,Mz1],-lbd,4/3,'s');
lmiterm([10,4,4,R10],tao,4/3);lmiterm([10,4,4,R13],tao,5/6);lmiterm([10,4,4,R14],tao,5/6);
lmiterm([10,4,5,R20],tao,4/3);
lmiterm([10,4,6,R30],tao,4/3);
lmiterm([10,5,5,Mz1],-lbd,4/3,'s');
lmiterm([10,5,5,R40],tao,4/3);lmiterm([10,5,5,R43],tao,5/6);lmiterm([10,5,5,R44],tao,5/6);
lmiterm([10,5,6,R50],tao,4/3);
lmiterm([10,6,6,Mz2],-lbd,4/3,'s');
lmiterm([10,6,6,R60],tao,4/3);lmiterm([10,6,6,R63],tao,5/6);lmiterm([10,6,6,R64],tao,5/6);

lmiterm([10,4,7,0],0);
lmiterm([10,4,8,0],0);
lmiterm([10,4,9,0],0);
lmiterm([10,5,7,0],0);
lmiterm([10,5,8,0],0);
lmiterm([10,5,9,0],0);
lmiterm([10,6,7,0],0);
lmiterm([10,6,8,0],0);
lmiterm([10,6,9,0],0);

lmiterm([10,4,10,0],0);
lmiterm([10,5,10,0],0);
lmiterm([10,6,10,0],(4/3)*lbd);
lmiterm([10,4,11,0],0);
lmiterm([10,4,12,0],0);
lmiterm([10,4,13,0],0);
lmiterm([10,5,11,0],0);lmiterm([10,5,12,-Mz1],-(5/6)*lbd,C1'*Ls3');lmiterm([10,5,12,-Mz1],-(5/6)*lbd,C1'*Ls4');lmiterm([10,5,13,-Mz1],-(5/6)*lbd,A3'*C1'*Ld4');lmiterm([10,5,13,-Mz1],-(5/6)*lbd,A4'*C1'*Ld3');
lmiterm([10,6,11,0],0);
lmiterm([10,6,12,0],0);lmiterm([10,6,13,-Mz2],-(5/6)*lbd,B3'*C1'*Ld4');lmiterm([10,6,13,-Mz2],-(5/6)*lbd,B4'*C1'*Ld3');

lmiterm([10,4,14,0],0);
lmiterm([10,4,15,0],0);
lmiterm([10,4,16,0],0);
lmiterm([10,5,14,0],0);
lmiterm([10,5,15,0],0);
lmiterm([10,5,16,0],0);
lmiterm([10,6,14,0],0);
lmiterm([10,6,15,0],0);
lmiterm([10,6,16,0],0);

lmiterm([10,7,7,Q10],-4/3,1);lmiterm([10,7,7,Q13],-5/6,1);lmiterm([10,7,7,Q14],-5/6,1);
lmiterm([10,7,8,Q20],-4/3,1);
lmiterm([10,7,9,Q30],-4/3,1);
lmiterm([10,8,8,Q40],-4/3,1);lmiterm([10,8,8,Q43],-5/6,1);lmiterm([10,8,8,Q44],-5/6,1);
lmiterm([10,8,9,Q50],-4/3,1);
lmiterm([10,9,9,Q60],-4/3,1);lmiterm([10,9,9,Q63],-5/6,1);lmiterm([10,9,9,Q64],-5/6,1);
lmiterm([10,7,10,0],0);
lmiterm([10,8,10,0],0);
lmiterm([10,9,10,0],0);
lmiterm([10,7,11,0],0);
lmiterm([10,7,12,0],0);
lmiterm([10,7,13,0],0);
lmiterm([10,8,11,0],0);
lmiterm([10,8,12,0],0);
lmiterm([10,8,13,0],0);
lmiterm([10,9,11,0],0);
lmiterm([10,9,12,0],0);
lmiterm([10,9,13,0],0);

lmiterm([10,7,14,0],0);
lmiterm([10,7,15,0],0);
lmiterm([10,7,16,0],0);
lmiterm([10,8,14,0],0);
lmiterm([10,8,15,0],0);
lmiterm([10,8,16,0],0);
lmiterm([10,9,14,0],0);
lmiterm([10,9,15,0],0);
lmiterm([10,9,16,0],0);

lmiterm([10,10,10,0],-(4/3)*gama1*gama1);
lmiterm([10,10,11,0],0);
lmiterm([10,10,12,0],0);
lmiterm([10,10,13,0],0);
lmiterm([10,10,14,0],0);
lmiterm([10,10,15,0],0);
lmiterm([10,10,16,0],0);

lmiterm([10,11,11,R10],-1/tao,4/3);lmiterm([10,11,11,R13],-1/tao,5/6);lmiterm([10,11,11,R14],-1/tao,5/6);
lmiterm([10,11,12,R20],-1/tao,4/3);
lmiterm([10,11,13,R30],-1/tao,4/3);
lmiterm([10,12,12,R40],-1/tao,4/3);lmiterm([10,12,12,R43],-1/tao,5/6);lmiterm([10,12,12,R44],-1/tao,5/6);
lmiterm([10,12,13,R50],-1/tao,4/3);
lmiterm([10,13,13,R60],-1/tao,4/3);lmiterm([10,13,13,R63],-1/tao,5/6);lmiterm([10,13,13,R64],-1/tao,5/6);

lmiterm([10,11,14,0],0);
lmiterm([10,11,15,0],0);
lmiterm([10,11,16,0],0);
lmiterm([10,12,14,0],0);
lmiterm([10,12,15,0],0);
lmiterm([10,12,16,0],0);
lmiterm([10,13,14,0],0);
lmiterm([10,13,15,0],0);
lmiterm([10,13,16,0],0);

lmiterm([10,14,14,0],-4/3);
lmiterm([10,14,15,0],0);
lmiterm([10,14,16,0],0);
lmiterm([10,15,15,0],-4/3);
lmiterm([10,15,16,0],0);
lmiterm([10,16,16,0],-4/3);



lmiterm([-11,1,1,P10],1,1);lmiterm([-11,1,1,P11],1,1);lmiterm([-11,1,2,P20],1,1);lmiterm([-11,1,3,P30],1,1);
lmiterm([-11,2,2,P40],1,1);lmiterm([-11,2,2,P41],1,1);lmiterm([-11,2,3,P50],1,1);
lmiterm([-11,3,3,P60],1,1);lmiterm([-11,3,3,P61],1,1);
lmiterm([-12,1,1,P10],1,1);lmiterm([-12,1,1,P12],1,1);lmiterm([-12,1,2,P20],1,1);lmiterm([-12,1,3,P30],1,1);
lmiterm([-12,2,2,P40],1,1);lmiterm([-12,2,2,P42],1,1);lmiterm([-12,2,3,P50],1,1);
lmiterm([-12,3,3,P60],1,1);lmiterm([-12,3,3,P62],1,1);
lmiterm([-13,1,1,P10],1,1);lmiterm([-13,1,1,P13],1,1);lmiterm([-13,1,2,P20],1,1);lmiterm([-13,1,3,P30],1,1);
lmiterm([-13,2,2,P40],1,1);lmiterm([-13,2,2,P43],1,1);lmiterm([-13,2,3,P50],1,1);
lmiterm([-13,3,3,P60],1,1);lmiterm([-13,3,3,P63],1,1);
lmiterm([-14,1,1,P10],1,1);lmiterm([-14,1,1,P14],1,1);lmiterm([-14,1,2,P20],1,1);lmiterm([-14,1,3,P30],1,1);
lmiterm([-14,2,2,P40],1,1);lmiterm([-14,2,2,P44],1,1);lmiterm([-14,2,3,P50],1,1);
lmiterm([-14,3,3,P60],1,1);lmiterm([-14,3,3,P64],1,1);

lmiterm([-15,1,1,Q10],1,1);lmiterm([-15,1,1,Q11],1,1);lmiterm([-15,1,2,Q20],1,1);lmiterm([-15,1,3,Q30],1,1);
lmiterm([-15,2,2,Q40],1,1);lmiterm([-15,2,2,Q41],1,1);lmiterm([-15,2,3,Q50],1,1);
lmiterm([-15,3,3,Q60],1,1);lmiterm([-15,3,3,Q61],1,1);
lmiterm([-16,1,1,Q10],1,1);lmiterm([-16,1,1,Q12],1,1);lmiterm([-16,1,2,Q20],1,1);lmiterm([-16,1,3,Q30],1,1);
lmiterm([-16,2,2,Q40],1,1);lmiterm([-16,2,2,Q42],1,1);lmiterm([-16,2,3,Q50],1,1);
lmiterm([-16,3,3,Q60],1,1);lmiterm([-16,3,3,Q62],1,1);
lmiterm([-17,1,1,Q10],1,1);lmiterm([-17,1,1,Q13],1,1);lmiterm([-17,1,2,Q20],1,1);lmiterm([-17,1,3,Q30],1,1);
lmiterm([-17,2,2,Q40],1,1);lmiterm([-17,2,2,Q43],1,1);lmiterm([-17,2,3,Q50],1,1);
lmiterm([-17,3,3,Q60],1,1);lmiterm([-17,3,3,Q63],1,1);
lmiterm([-18,1,1,Q10],1,1);lmiterm([-18,1,1,Q14],1,1);lmiterm([-18,1,2,Q20],1,1);lmiterm([-18,1,3,Q30],1,1);
lmiterm([-18,2,2,Q40],1,1);lmiterm([-18,2,2,Q44],1,1);lmiterm([-18,2,3,Q50],1,1);
lmiterm([-18,3,3,Q60],1,1);lmiterm([-18,3,3,Q64],1,1);


lmiterm([-19,1,1,R10],1,1);lmiterm([-19,1,1,R11],1,1);lmiterm([-19,1,2,R20],1,1);lmiterm([-19,1,3,R30],1,1);
lmiterm([-19,2,2,R40],1,1);lmiterm([-19,2,2,R41],1,1);lmiterm([-19,2,3,R50],1,1);
lmiterm([-19,3,3,R60],1,1);lmiterm([-19,3,3,R61],1,1);
lmiterm([-20,1,1,R10],1,1);lmiterm([-20,1,1,R12],1,1);lmiterm([-20,1,2,R20],1,1);lmiterm([-20,1,3,R30],1,1);
lmiterm([-20,2,2,R40],1,1);lmiterm([-20,2,2,R42],1,1);lmiterm([-20,2,3,R50],1,1);
lmiterm([-20,3,3,R60],1,1);lmiterm([-20,3,3,R62],1,1);
lmiterm([-21,1,1,R10],1,1);lmiterm([-21,1,1,R13],1,1);lmiterm([-21,1,2,R20],1,1);lmiterm([-21,1,3,R30],1,1);
lmiterm([-21,2,2,R40],1,1);lmiterm([-21,2,2,R43],1,1);lmiterm([-21,2,3,R50],1,1);
lmiterm([-21,3,3,R60],1,1);lmiterm([-21,3,3,R63],1,1);
lmiterm([-22,1,1,R10],1,1);lmiterm([-22,1,1,R14],1,1);lmiterm([-22,1,2,R20],1,1);lmiterm([-22,1,3,R30],1,1);
lmiterm([-22,2,2,R40],1,1);lmiterm([-22,2,2,R44],1,1);lmiterm([-22,2,3,R50],1,1);
lmiterm([-22,3,3,R60],1,1);lmiterm([-22,3,3,R64],1,1);












lmis=getlmis;   
[tmin,xfeas]=feasp(lmis);
tmin

Mc=dec2mat(lmis,xfeas,Mz1);
Zc1=dec2mat(lmis,xfeas,Zg1);
Zc2=dec2mat(lmis,xfeas,Zg2);
Zc3=dec2mat(lmis,xfeas,Zg3);
Zc4=dec2mat(lmis,xfeas,Zg4);


K{1}=Zc1*inv(Mc);
K{2}=Zc2*inv(Mc);
K{3}=Zc3*inv(Mc);
K{4}=Zc4*inv(Mc);

K1=Zc1*inv(Mc)
K2=Zc2*inv(Mc)
K3=Zc3*inv(Mc)
K4=Zc4*inv(Mc)

save controller K



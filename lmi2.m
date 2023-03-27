
clear all
close all
setlmis([]);

cf0=95000;cr0=85500;lr=1.67;lf=1.11;vxmin=25;vxmax=25; mmin=1530;mmax=1680;Izmin=4200;Izmax=4600;
tao=0.019;lbd=0.1;gama1=0.061;


A1=[-(cf0+cr0)/(mmin*vxmax) -vxmax+(cr0*lr-cf0*lf)/(mmin*vxmax);(cr0*lr-cf0*lf)/(Izmin*vxmax) -(cf0*lf*lf+cr0*lr*lr)/(Izmin*vxmax)];
A2=[-(cf0+cr0)/(mmin*vxmax) -vxmax+(cr0*lr-cf0*lf)/(mmin*vxmax);(cr0*lr-cf0*lf)/(Izmax*vxmax) -(cf0*lf*lf+cr0*lr*lr)/(Izmax*vxmax)];
A3=[-(cf0+cr0)/(mmax*vxmin) -vxmin+(cr0*lr-cf0*lf)/(mmax*vxmin);(cr0*lr-cf0*lf)/(Izmin*vxmin) -(cf0*lf*lf+cr0*lr*lr)/(Izmin*vxmin)];
A4=[-(cf0+cr0)/(mmax*vxmin) -vxmin+(cr0*lr-cf0*lf)/(mmax*vxmin);(cr0*lr-cf0*lf)/(Izmax*vxmin) -(cf0*lf*lf+cr0*lr*lr)/(Izmax*vxmin)];


B1=[cf0/mmin;(lf*cf0)/Izmin];B2=[cf0/mmin;(lf*cf0)/Izmax];B3=[cf0/mmax;(lf*cf0)/Izmin];B4=[cf0/mmax;(lf*cf0)/Izmax];

E1=[cf0/mmin;(lf*cf0)/Izmin];E2=[cf0/mmin;(lf*cf0)/Izmax];E3=[cf0/mmax;(lf*cf0)/Izmin];E4=[cf0/mmax;(lf*cf0)/Izmax];


C1=[0 1];
D=eye(2);

Mz1=lmivar(2,[2,2]);Mz2=lmivar(2,[1,1]);

P10=lmivar(1,[2,1]);P11=lmivar(1,[2,0]);P12=lmivar(1,[2,0]);P13=lmivar(1,[2,0]);P14=lmivar(1,[2,0]);
P20=lmivar(2,[2,1]);
P30=lmivar(1,[1,1]);P31=lmivar(1,[1,0]);P32=lmivar(1,[1,0]);P33=lmivar(1,[1,0]);P34=lmivar(1,[1,0]);

Q10=lmivar(1,[2,1]);Q11=lmivar(1,[2,0]);Q12=lmivar(1,[2,0]);Q13=lmivar(1,[2,0]);Q14=lmivar(1,[2,0]);
Q30=lmivar(1,[1,1]);Q31=lmivar(1,[1,0]);Q32=lmivar(1,[1,0]);Q33=lmivar(1,[1,0]);Q34=lmivar(1,[1,0]);
Q20=lmivar(2,[2,1]);

R10=lmivar(1,[2,1]);R11=lmivar(1,[2,0]);R12=lmivar(1,[2,0]);R13=lmivar(1,[2,0]);R14=lmivar(1,[2,0]);
R30=lmivar(1,[1,1]);R31=lmivar(1,[1,0]);R32=lmivar(1,[1,0]);R33=lmivar(1,[1,0]);R34=lmivar(1,[1,0]);
R20=lmivar(2,[2,1]);

Zs1=lmivar(2,[2,1]);Zs2=lmivar(2,[2,1]);Zs3=lmivar(2,[2,1]);Zs4=lmivar(2,[2,1]);
Zd1=lmivar(2,[1,1]);Zd2=lmivar(2,[1,1]);Zd3=lmivar(2,[1,1]);Zd4=lmivar(2,[1,1]);


lmiterm([1,1,1,Q10],1,1);lmiterm([1,1,1,Q11],1,1);
lmiterm([1,1,1,P10],-2/lbd,1);lmiterm([1,1,1,P11],-2/lbd,1);
lmiterm([1,1,1,R10],tao/(lbd*lbd),1);lmiterm([1,1,1,R11],tao/(lbd*lbd),1);
lmiterm([1,1,2,Q20],1,1);lmiterm([1,1,2,P20],-2/lbd,1);lmiterm([1,1,2,R20],tao/(lbd*lbd),1);
lmiterm([1,2,2,Q30],1,1);lmiterm([1,2,2,Q31],1,1);lmiterm([1,2,2,0],1);
lmiterm([1,2,2,P30],-2/lbd,1);lmiterm([1,2,2,P31],-2/lbd,1);
lmiterm([1,2,2,R30],tao/(lbd*lbd),1);lmiterm([1,2,2,R31],tao/(lbd*lbd),1);
lmiterm([1,1,3,P10],1,1);lmiterm([1,1,3,P11],1,1);lmiterm([1,1,3,Mz1],1,1);lmiterm([1,1,3,Mz1],A1',lbd);lmiterm([1,1,3,-Zs1],-C1',lbd);
lmiterm([1,1,3,R10],-tao/lbd,1);lmiterm([1,1,3,R11],-tao/lbd,1);

lmiterm([1,1,4,P20],1,1);lmiterm([1,1,4,R20],-tao/lbd,1);
lmiterm([1,1,4,-Zd1],-A1'*C1',lbd);lmiterm([1,2,3,-P20],1,1);lmiterm([1,2,3,-R20],-tao/lbd,1);
lmiterm([1,2,3,Mz1],B1',lbd);
lmiterm([1,2,4,P30],1,1);lmiterm([1,2,4,P31],1,1);
lmiterm([1,2,4,-Zd1],-B1'*C1',lbd);lmiterm([1,2,4,Mz2],1,1);
lmiterm([1,2,4,R30],-tao/lbd,1);lmiterm([1,2,4,R31],-tao/lbd,1);
lmiterm([1,1,5,0],0);lmiterm([1,1,6,0],0);
lmiterm([1,2,5,0],0);lmiterm([1,2,6,0],0);
lmiterm([1,1,7,0],0);
lmiterm([1,2,7,0],0);
lmiterm([1,1,8,0],0);lmiterm([1,1,9,0],0);
lmiterm([1,2,8,0],0);lmiterm([1,2,9,0],0);


lmiterm([1,3,3,Mz1],-lbd,1,'s');lmiterm([1,3,3,R10],tao,1);lmiterm([1,3,3,R11],tao,1);lmiterm([1,3,4,R20],tao,1);
lmiterm([1,4,4,Mz2],-lbd,1,'s');lmiterm([1,4,4,R30],tao,1);lmiterm([1,4,4,R31],tao,1);
lmiterm([1,3,5,0],0);lmiterm([1,3,6,0],0);
lmiterm([1,4,5,0],0);lmiterm([1,4,6,0],0);
lmiterm([1,3,7,0],0);
lmiterm([1,4,7,-Mz2],lbd,1);
lmiterm([1,3,8,Zs1],-lbd,C1);lmiterm([1,3,9,0],0);
lmiterm([1,4,8,Zd1],-lbd,C1*A1);lmiterm([1,4,9,Zd1],-lbd,C1*B1);



lmiterm([1,5,5,Q10],-1,1);lmiterm([1,5,5,Q11],-1,1);
lmiterm([1,5,6,Q20],-1,1);
lmiterm([1,6,6,Q30],-1,1);lmiterm([1,6,6,Q31],-1,1);
lmiterm([1,5,7,0],0);
lmiterm([1,6,7,0],0);
lmiterm([1,5,8,0],0);lmiterm([1,5,9,0],0);
lmiterm([1,6,8,0],0);lmiterm([1,6,9,0],0);


lmiterm([1,7,7,0],-gama1*gama1);
lmiterm([1,7,8,0],0);lmiterm([1,7,9,0],0);


lmiterm([1,8,8,R10],-1/tao,1);lmiterm([1,8,8,R11],-1/tao,1);lmiterm([1,8,9,R20],-1/tao,1);
lmiterm([1,9,9,R30],-1/tao,1);lmiterm([1,9,9,R31],-1/tao,1);



lmiterm([2,1,1,Q10],1,1);lmiterm([2,1,1,Q12],1,1);
lmiterm([2,1,1,P10],-2/lbd,1);lmiterm([2,1,1,P12],-2/lbd,1);
lmiterm([2,1,1,R10],tao/(lbd*lbd),1);lmiterm([2,1,1,R12],tao/(lbd*lbd),1);
lmiterm([2,1,2,Q20],1,1);lmiterm([2,1,2,P20],-2/lbd,1);lmiterm([2,1,2,R20],tao/(lbd*lbd),1);
lmiterm([2,2,2,Q30],1,1);lmiterm([2,2,2,Q32],1,1);lmiterm([2,2,2,0],1);
lmiterm([2,2,2,P30],-2/lbd,1);lmiterm([2,2,2,P32],-2/lbd,1);
lmiterm([2,2,2,R30],tao/(lbd*lbd),1);lmiterm([2,2,2,R32],tao/(lbd*lbd),1);
lmiterm([2,1,3,P10],1,1);lmiterm([2,1,3,P12],1,1);lmiterm([2,1,3,Mz1],1,1);lmiterm([2,1,3,Mz1],A2',lbd);lmiterm([2,1,3,-Zs2],-C1',lbd);
lmiterm([2,1,3,R10],-tao/lbd,1);lmiterm([2,1,3,R12],-tao/lbd,1);

lmiterm([2,1,4,P20],1,1);lmiterm([2,1,4,R20],-tao/lbd,1);
lmiterm([2,1,4,-Zd2],-A2'*C1',lbd);lmiterm([2,2,3,-P20],1,1);lmiterm([2,2,3,-R20],-tao/lbd,1);
lmiterm([2,2,3,Mz1],B2',lbd);
lmiterm([2,2,4,P30],1,1);lmiterm([2,2,4,P32],1,1);
lmiterm([2,2,4,-Zd2],-B2'*C1',lbd);lmiterm([2,2,4,Mz2],1,1);
lmiterm([2,2,4,R30],-tao/lbd,1);lmiterm([2,2,4,R32],-tao/lbd,1);
lmiterm([2,1,5,0],0);
lmiterm([2,1,6,0],0);
lmiterm([2,2,5,0],0);
lmiterm([2,2,6,0],0);
lmiterm([2,1,7,0],0);
lmiterm([2,2,7,0],0);
lmiterm([2,1,8,0],0);
lmiterm([2,1,9,0],0);
lmiterm([2,2,8,0],0);
lmiterm([2,2,9,0],0);

lmiterm([2,3,3,Mz1],-lbd,1,'s');
lmiterm([2,3,3,R10],tao,1);lmiterm([2,3,3,R12],tao,1);
lmiterm([2,3,4,R20],tao,1);
lmiterm([2,4,4,Mz2],-lbd,1,'s');
lmiterm([2,4,4,R30],tao,1);lmiterm([2,4,4,R32],tao,1);
lmiterm([2,3,5,0],0);
lmiterm([2,3,6,0],0);
lmiterm([2,4,5,0],0);
lmiterm([2,4,6,0],0);
lmiterm([2,3,7,0],0);
lmiterm([2,4,7,-Mz2],lbd,1);lmiterm([2,3,8,Zs2],-lbd,C1);
lmiterm([2,3,9,0],0);lmiterm([2,4,8,Zd2],-lbd,C1*A2);lmiterm([2,4,9,Zd2],-lbd,C1*B2);

lmiterm([2,5,5,Q10],-1,1);lmiterm([2,5,5,Q12],-1,1);
lmiterm([2,5,6,Q20],-1,1);
lmiterm([2,6,6,Q30],-1,1);lmiterm([2,6,6,Q32],-1,1);
lmiterm([2,5,7,0],0);
lmiterm([2,6,7,0],0);
lmiterm([2,5,8,0],0);lmiterm([2,5,9,0],0);
lmiterm([2,6,8,0],0);lmiterm([2,6,9,0],0);

lmiterm([2,7,7,0],-gama1*gama1);
lmiterm([2,7,8,0],0);
lmiterm([2,7,9,0],0);

lmiterm([2,8,8,R10],-1/tao,1);lmiterm([2,8,8,R12],-1/tao,1);
lmiterm([2,8,9,R20],-1/tao,1);
lmiterm([2,9,9,R30],-1/tao,1);lmiterm([2,9,9,R32],-1/tao,1);



lmiterm([3,1,1,Q10],1,1);lmiterm([3,1,1,Q13],1,1);
lmiterm([3,1,1,P10],-2/lbd,1);lmiterm([3,1,1,P13],-2/lbd,1);
lmiterm([3,1,1,R10],tao/(lbd*lbd),1);lmiterm([3,1,1,R13],tao/(lbd*lbd),1);
lmiterm([3,1,2,Q20],1,1);lmiterm([3,1,2,P20],-2/lbd,1);lmiterm([3,1,2,R20],tao/(lbd*lbd),1);
lmiterm([3,2,2,Q30],1,1);lmiterm([3,2,2,Q33],1,1);lmiterm([3,2,2,0],1);
lmiterm([3,2,2,P30],-2/lbd,1);lmiterm([3,2,2,P33],-2/lbd,1);
lmiterm([3,2,2,R30],tao/(lbd*lbd),1);lmiterm([3,2,2,R33],tao/(lbd*lbd),1);
lmiterm([3,1,3,P10],1,1);lmiterm([3,1,3,P13],1,1);lmiterm([3,1,3,Mz1],1,1);lmiterm([3,1,3,Mz1],A3',lbd);lmiterm([3,1,3,-Zs3],-C1',lbd);
lmiterm([3,1,3,R10],-tao/lbd,1);lmiterm([3,1,3,R13],-tao/lbd,1);

lmiterm([3,1,4,P20],1,1);lmiterm([3,1,4,R20],-tao/lbd,1);
lmiterm([3,1,4,-Zd3],-A3'*C1',lbd);lmiterm([3,2,3,-P20],1,1);lmiterm([3,2,3,-R20],-tao/lbd,1);
lmiterm([3,2,3,Mz1],B3',lbd);
lmiterm([3,2,4,P30],1,1);lmiterm([3,2,4,P33],1,1);
lmiterm([3,2,4,-Zd3],-B3'*C1',lbd);lmiterm([3,2,4,Mz2],1,1);
lmiterm([3,2,4,R30],-tao/lbd,1);lmiterm([3,2,4,R33],-tao/lbd,1);
lmiterm([3,1,5,0],0);lmiterm([3,1,6,0],0);
lmiterm([3,2,5,0],0);lmiterm([3,2,6,0],0);
lmiterm([3,1,7,0],0);
lmiterm([3,2,7,0],0);
lmiterm([3,1,8,0],0);lmiterm([3,1,9,0],0);
lmiterm([3,2,8,0],0);lmiterm([3,2,9,0],0);

lmiterm([3,3,3,Mz1],-lbd,1,'s');lmiterm([3,3,3,R10],tao,1);lmiterm([3,3,3,R13],tao,1);lmiterm([3,3,4,R20],tao,1);
lmiterm([3,4,4,Mz2],-lbd,1,'s');lmiterm([3,4,4,R30],tao,1);lmiterm([3,4,4,R33],tao,1);
lmiterm([3,3,5,0],0);lmiterm([3,3,6,0],0);
lmiterm([3,4,5,0],0);lmiterm([3,4,6,0],0);
lmiterm([3,3,7,0],0);
lmiterm([3,4,7,-Mz2],lbd,1);
lmiterm([3,3,8,Zs3],-lbd,C1);lmiterm([3,3,9,0],0);
lmiterm([3,4,8,Zd3],-lbd,C1*A3);lmiterm([3,4,9,Zd3],-lbd,C1*B3);

lmiterm([3,5,5,Q10],-1,1);lmiterm([3,5,5,Q13],-1,1);
lmiterm([3,5,6,Q20],-1,1);
lmiterm([3,6,6,Q30],-1,1);lmiterm([3,6,6,Q33],-1,1);
lmiterm([3,5,7,0],0);
lmiterm([3,6,7,0],0);
lmiterm([3,5,8,0],0);lmiterm([3,5,9,0],0);
lmiterm([3,6,8,0],0);lmiterm([3,6,9,0],0);

lmiterm([3,7,7,0],-gama1*gama1);
lmiterm([3,7,8,0],0);lmiterm([3,7,9,0],0);

lmiterm([3,8,8,R10],-1/tao,1);lmiterm([3,8,8,R13],-1/tao,1);lmiterm([3,8,9,R20],-1/tao,1);
lmiterm([3,9,9,R30],-1/tao,1);lmiterm([3,9,9,R33],-1/tao,1);



lmiterm([4,1,1,Q10],1,1);lmiterm([4,1,1,Q14],1,1);
lmiterm([4,1,1,P10],-2/lbd,1);lmiterm([4,1,1,P14],-2/lbd,1);
lmiterm([4,1,1,R10],tao/(lbd*lbd),1);lmiterm([4,1,1,R14],tao/(lbd*lbd),1);
lmiterm([4,1,2,Q20],1,1);lmiterm([4,1,2,P20],-2/lbd,1);lmiterm([4,1,2,R20],tao/(lbd*lbd),1);
lmiterm([4,2,2,Q30],1,1);lmiterm([4,2,2,Q34],1,1);lmiterm([4,2,2,0],1);
lmiterm([4,2,2,P30],-2/lbd,1);lmiterm([4,2,2,P34],-2/lbd,1);
lmiterm([4,2,2,R30],tao/(lbd*lbd),1);lmiterm([4,2,2,R34],tao/(lbd*lbd),1);
lmiterm([4,1,3,P10],1,1);lmiterm([4,1,3,P14],1,1);lmiterm([4,1,3,Mz1],1,1);lmiterm([4,1,3,Mz1],A4',lbd);lmiterm([4,1,3,-Zs4],-C1',lbd);
lmiterm([4,1,3,R10],-tao/lbd,1);lmiterm([4,1,3,R14],-tao/lbd,1);

lmiterm([4,1,4,P20],1,1);lmiterm([4,1,4,R20],-tao/lbd,1);
lmiterm([4,1,4,-Zd4],-A4'*C1',lbd);lmiterm([4,2,3,-P20],1,1);lmiterm([4,2,3,-R20],-tao/lbd,1);
lmiterm([4,2,3,Mz1],B4',lbd);
lmiterm([4,2,4,P30],1,1);lmiterm([4,2,4,P34],1,1);
lmiterm([4,2,4,-Zd4],-B4'*C1',lbd);lmiterm([4,2,4,Mz2],1,1);
lmiterm([4,2,4,R30],-tao/lbd,1);lmiterm([4,2,4,R34],-tao/lbd,1);
lmiterm([4,1,5,0],0);lmiterm([4,1,6,0],0);
lmiterm([4,2,5,0],0);lmiterm([4,2,6,0],0);
lmiterm([4,1,7,0],0);
lmiterm([4,2,7,0],0);
lmiterm([4,1,8,0],0);lmiterm([4,1,9,0],0);
lmiterm([4,2,8,0],0);lmiterm([4,2,9,0],0);

lmiterm([4,3,3,Mz1],-lbd,1,'s');lmiterm([4,3,3,R10],tao,1);lmiterm([4,3,3,R14],tao,1);lmiterm([4,3,4,R20],tao,1);
lmiterm([4,4,4,Mz2],-lbd,1,'s');lmiterm([4,4,4,R30],tao,1);lmiterm([4,4,4,R34],tao,1);
lmiterm([4,3,5,0],0);lmiterm([4,3,6,0],0);
lmiterm([4,4,5,0],0);lmiterm([4,4,6,0],0);
lmiterm([4,3,7,0],0);
lmiterm([4,4,7,-Mz2],lbd,1);
lmiterm([4,3,8,Zs4],-lbd,C1);lmiterm([4,3,9,0],0);
lmiterm([4,4,8,Zd4],-lbd,C1*A4);lmiterm([4,4,9,Zd4],-lbd,C1*B4);

lmiterm([4,5,5,Q10],-1,1);lmiterm([4,5,5,Q14],-1,1);
lmiterm([4,5,6,Q20],-1,1);
lmiterm([4,6,6,Q30],-1,1);lmiterm([4,6,6,Q34],-1,1);
lmiterm([4,5,7,0],0);
lmiterm([4,6,7,0],0);
lmiterm([4,5,8,0],0);lmiterm([4,5,9,0],0);
lmiterm([4,6,8,0],0);lmiterm([4,6,9,0],0);

lmiterm([4,7,7,0],-gama1*gama1);
lmiterm([4,7,8,0],0);lmiterm([4,7,9,0],0);
lmiterm([4,8,8,R10],-1/tao,1);lmiterm([4,8,8,R14],-1/tao,1);lmiterm([4,8,9,R20],-1/tao,1);
lmiterm([4,9,9,R30],-1/tao,1);lmiterm([4,9,9,R34],-1/tao,1);



lmiterm([5,1,1,Q10],4/3,1);lmiterm([5,1,1,Q11],5/6,1);lmiterm([5,1,1,Q12],5/6,1);
lmiterm([5,1,1,P10],-2/lbd,4/3);lmiterm([5,1,1,P11],-2/lbd,5/6);lmiterm([5,1,1,P12],-2/lbd,5/6);
lmiterm([5,1,1,R10],tao/(lbd*lbd),4/3);lmiterm([5,1,1,R11],tao/(lbd*lbd),5/6);lmiterm([5,1,1,R12],tao/(lbd*lbd),5/6);
lmiterm([5,1,2,Q20],4/3,1);
lmiterm([5,1,2,P20],-2/lbd,4/3);
lmiterm([5,1,2,R20],tao/(lbd*lbd),4/3);
lmiterm([5,2,2,Q30],4/3,1);lmiterm([5,2,2,Q31],5/6,1);lmiterm([5,2,2,Q32],5/6,1);
lmiterm([5,2,2,0],4/3);
lmiterm([5,2,2,P30],-2/lbd,4/3);lmiterm([5,2,2,P31],-2/lbd,5/6);lmiterm([5,2,2,P32],-2/lbd,5/6);
lmiterm([5,2,2,R30],tao/(lbd*lbd),4/3);lmiterm([5,2,2,R31],tao/(lbd*lbd),5/6);lmiterm([5,2,2,R32],tao/(lbd*lbd),5/6);
lmiterm([5,1,3,P10],4/3,1);lmiterm([5,1,3,P11],5/6,1);lmiterm([5,1,3,P12],5/6,1);
lmiterm([5,1,3,Mz1],4/3,1);lmiterm([5,1,3,Mz1],A1',(5/6)*lbd);lmiterm([5,1,3,Mz1],A2',(5/6)*lbd);lmiterm([5,1,3,-Zs1],-C1',(5/6)*lbd);lmiterm([5,1,3,-Zs2],-C1',(5/6)*lbd);
lmiterm([5,1,3,R10],-tao/lbd,4/3);lmiterm([5,1,3,R11],-tao/lbd,5/6);lmiterm([5,1,3,R12],-tao/lbd,5/6);

lmiterm([5,1,4,P20],4/3,1);
lmiterm([5,1,4,R20],-tao/lbd,4/3);lmiterm([5,1,4,-Zd1],-A2'*C1',(5/6)*lbd);lmiterm([5,1,4,-Zd2],-A1'*C1',(5/6)*lbd);
lmiterm([5,2,3,-P20],4/3,1);
lmiterm([5,2,3,-R20],-tao/lbd,4/3);lmiterm([5,2,3,Mz1],B1',(5/6)*lbd);lmiterm([5,2,3,Mz1],B2',(5/6)*lbd);
lmiterm([5,2,4,P30],4/3,1);lmiterm([5,2,4,P31],5/6,1);lmiterm([5,2,4,P32],5/6,1);lmiterm([5,2,4,-Zd1],-B2'*C1',(5/6)*lbd);lmiterm([5,2,4,-Zd2],-B1'*C1',(5/6)*lbd);
lmiterm([5,2,4,Mz2],4/3,1);
lmiterm([5,2,4,R30],-tao/lbd,4/3);lmiterm([5,2,4,R31],-tao/lbd,5/6);lmiterm([5,2,4,R32],-tao/lbd,5/6);
lmiterm([5,1,5,0],0);
lmiterm([5,1,6,0],0);
lmiterm([5,2,5,0],0);
lmiterm([5,2,6,0],0);
lmiterm([5,1,7,0],0);
lmiterm([5,2,7,0],0);
lmiterm([5,1,8,0],0);
lmiterm([5,1,9,0],0);
lmiterm([5,2,8,0],0);
lmiterm([5,2,9,0],0);

lmiterm([5,3,3,Mz1],-lbd,4/3,'s');
lmiterm([5,3,3,R10],tao,4/3);lmiterm([5,3,3,R11],tao,5/6);lmiterm([5,3,3,R12],tao,5/6);
lmiterm([5,3,4,R20],tao,4/3);
lmiterm([5,4,4,Mz2],-lbd,4/3,'s');
lmiterm([5,4,4,R30],tao,4/3);lmiterm([5,4,4,R31],tao,5/6);lmiterm([5,4,4,R32],tao,5/6);
lmiterm([5,3,5,0],0);
lmiterm([5,3,6,0],0);
lmiterm([5,4,5,0],0);
lmiterm([5,4,6,0],0);
lmiterm([5,3,7,0],0);
lmiterm([5,4,7,-Mz2],4/3,lbd);lmiterm([5,3,8,Zs1],-(5/6)*lbd,C1);lmiterm([5,3,8,Zs2],-(5/6)*lbd,C1);
lmiterm([5,3,9,0],0);lmiterm([5,4,8,Zd1],-(5/6)*lbd,C1*A2);lmiterm([5,4,8,Zd2],-(5/6)*lbd,C1*A1);lmiterm([5,4,9,Zd1],-(5/6)*lbd,C1*B2);lmiterm([5,4,9,Zd2],-(5/6)*lbd,C1*B1);

lmiterm([5,5,5,Q10],-4/3,1);lmiterm([5,5,5,Q11],-5/6,1);lmiterm([5,5,5,Q12],-5/6,1);
lmiterm([5,5,6,Q20],-4/3,1);
lmiterm([5,6,6,Q30],-4/3,1);lmiterm([5,6,6,Q31],-5/6,1);lmiterm([5,6,6,Q32],-5/6,1);
lmiterm([5,5,7,0],0);
lmiterm([5,6,7,0],0);
lmiterm([5,5,8,0],0);lmiterm([5,5,9,0],0);
lmiterm([5,6,8,0],0);lmiterm([5,6,9,0],0);

lmiterm([5,7,7,0],-4/3*gama1*gama1);
lmiterm([5,7,8,0],0);
lmiterm([5,7,9,0],0);

lmiterm([5,8,8,R10],-1/tao,4/3);lmiterm([5,8,8,R11],-1/tao,5/6);lmiterm([5,8,8,R12],-1/tao,5/6);
lmiterm([5,8,9,R20],-1/tao,4/3);
lmiterm([5,9,9,R30],-1/tao,4/3);lmiterm([5,9,9,R31],-1/tao,5/6);lmiterm([5,9,9,R32],-1/tao,5/6);




lmiterm([6,1,1,Q10],4/3,1);lmiterm([6,1,1,Q11],5/6,1);lmiterm([6,1,1,Q13],5/6,1);
lmiterm([6,1,1,P10],-2/lbd,4/3);lmiterm([6,1,1,P11],-2/lbd,5/6);lmiterm([6,1,1,P13],-2/lbd,5/6);
lmiterm([6,1,1,R10],tao/(lbd*lbd),4/3);lmiterm([6,1,1,R11],tao/(lbd*lbd),5/6);lmiterm([6,1,1,R13],tao/(lbd*lbd),5/6);
lmiterm([6,1,2,Q20],4/3,1);
lmiterm([6,1,2,P20],-2/lbd,4/3);
lmiterm([6,1,2,R20],tao/(lbd*lbd),4/3);
lmiterm([6,2,2,Q30],4/3,1);lmiterm([6,2,2,Q31],5/6,1);lmiterm([6,2,2,Q33],5/6,1);
lmiterm([6,2,2,0],4/3);
lmiterm([6,2,2,P30],-2/lbd,4/3);lmiterm([6,2,2,P31],-2/lbd,5/6);lmiterm([6,2,2,P33],-2/lbd,5/6);
lmiterm([6,2,2,R30],tao/(lbd*lbd),4/3);lmiterm([6,2,2,R31],tao/(lbd*lbd),5/6);lmiterm([6,2,2,R33],tao/(lbd*lbd),5/6);
lmiterm([6,1,3,P10],4/3,1);lmiterm([6,1,3,P11],5/6,1);lmiterm([6,1,3,P13],5/6,1);
lmiterm([6,1,3,Mz1],4/3,1);lmiterm([6,1,3,Mz1],A1',(5/6)*lbd);lmiterm([6,1,3,Mz1],A3',(5/6)*lbd);lmiterm([6,1,3,-Zs1],-C1',(5/6)*lbd);lmiterm([6,1,3,-Zs3],-C1',(5/6)*lbd);
lmiterm([6,1,3,R10],-tao/lbd,4/3);lmiterm([6,1,3,R11],-tao/lbd,5/6);lmiterm([6,1,3,R13],-tao/lbd,5/6);

lmiterm([6,1,4,P20],4/3,1);
lmiterm([6,1,4,R20],-tao/lbd,4/3);lmiterm([6,1,4,-Zd1],-A3'*C1',(5/6)*lbd);lmiterm([6,1,4,-Zd3],-A1'*C1',(5/6)*lbd);
lmiterm([6,2,3,-P20],4/3,1);
lmiterm([6,2,3,-R20],-tao/lbd,4/3);lmiterm([6,2,3,Mz1],B1',(5/6)*lbd);lmiterm([6,2,3,Mz1],B3',(5/6)*lbd);
lmiterm([6,2,4,P30],4/3,1);lmiterm([6,2,4,P31],5/6,1);lmiterm([6,2,4,P33],5/6,1);lmiterm([6,2,4,-Zd1],-B3'*C1',(5/6)*lbd);lmiterm([6,2,4,-Zd3],-B1'*C1',(5/6)*lbd);
lmiterm([6,2,4,Mz2],4/3,1);
lmiterm([6,2,4,R30],-tao/lbd,4/3);lmiterm([6,2,4,R31],-tao/lbd,5/6);lmiterm([6,2,4,R33],-tao/lbd,5/6);
lmiterm([6,1,5,0],0);
lmiterm([6,1,6,0],0);
lmiterm([6,2,5,0],0);
lmiterm([6,2,6,0],0);
lmiterm([6,1,7,0],0);
lmiterm([6,2,7,0],0);
lmiterm([6,1,8,0],0);
lmiterm([6,1,9,0],0);
lmiterm([6,2,8,0],0);
lmiterm([6,2,9,0],0);

lmiterm([6,3,3,Mz1],-lbd,4/3,'s');
lmiterm([6,3,3,R10],tao,4/3);lmiterm([6,3,3,R11],tao,5/6);lmiterm([6,3,3,R13],tao,5/6);
lmiterm([6,3,4,R20],tao,4/3);
lmiterm([6,4,4,Mz2],-lbd,4/3,'s');
lmiterm([6,4,4,R30],tao,4/3);lmiterm([6,4,4,R31],tao,5/6);lmiterm([6,4,4,R33],tao,5/6);
lmiterm([6,3,5,0],0);
lmiterm([6,3,6,0],0);
lmiterm([6,4,5,0],0);
lmiterm([6,4,6,0],0);
lmiterm([6,3,7,0],0);
lmiterm([6,4,7,-Mz2],4/3,lbd);lmiterm([6,3,8,Zs1],-(5/6)*lbd,C1);lmiterm([6,3,8,Zs3],-(5/6)*lbd,C1);
lmiterm([6,3,9,0],0);lmiterm([6,4,8,Zd1],-(5/6)*lbd,C1*A3);lmiterm([6,4,8,Zd3],-(5/6)*lbd,C1*A1);lmiterm([6,4,9,Zd1],-(5/6)*lbd,C1*B3);lmiterm([6,4,9,Zd3],-(5/6)*lbd,C1*B1);

lmiterm([6,5,5,Q10],-4/3,1);lmiterm([6,5,5,Q11],-5/6,1);lmiterm([6,5,5,Q13],-5/6,1);
lmiterm([6,5,6,Q20],-4/3,1);
lmiterm([6,6,6,Q30],-4/3,1);lmiterm([6,6,6,Q31],-5/6,1);lmiterm([6,6,6,Q33],-5/6,1);
lmiterm([6,5,7,0],0);
lmiterm([6,6,7,0],0);
lmiterm([6,5,8,0],0);lmiterm([6,5,9,0],0);
lmiterm([6,6,8,0],0);lmiterm([6,6,9,0],0);

lmiterm([6,7,7,0],-4/3*gama1*gama1);
lmiterm([6,7,8,0],0);
lmiterm([6,7,9,0],0);

lmiterm([6,8,8,R10],-1/tao,4/3);lmiterm([6,8,8,R11],-1/tao,5/6);lmiterm([6,8,8,R13],-1/tao,5/6);
lmiterm([6,8,9,R20],-1/tao,4/3);
lmiterm([6,9,9,R30],-1/tao,4/3);lmiterm([6,9,9,R31],-1/tao,5/6);lmiterm([6,9,9,R33],-1/tao,5/6);



lmiterm([7,1,1,Q10],4/3,1);lmiterm([7,1,1,Q11],5/6,1);lmiterm([7,1,1,Q14],5/6,1);
lmiterm([7,1,1,P10],-2/lbd,4/3);lmiterm([7,1,1,P11],-2/lbd,5/6);lmiterm([7,1,1,P14],-2/lbd,5/6);
lmiterm([7,1,1,R10],tao/(lbd*lbd),4/3);lmiterm([7,1,1,R11],tao/(lbd*lbd),5/6);lmiterm([7,1,1,R14],tao/(lbd*lbd),5/6);
lmiterm([7,1,2,Q20],4/3,1);
lmiterm([7,1,2,P20],-2/lbd,4/3);
lmiterm([7,1,2,R20],tao/(lbd*lbd),4/3);
lmiterm([7,2,2,Q30],4/3,1);lmiterm([7,2,2,Q31],5/6,1);lmiterm([7,2,2,Q34],5/6,1);
lmiterm([7,2,2,0],4/3);
lmiterm([7,2,2,P30],-2/lbd,4/3);lmiterm([7,2,2,P31],-2/lbd,5/6);lmiterm([7,2,2,P34],-2/lbd,5/6);
lmiterm([7,2,2,R30],tao/(lbd*lbd),4/3);lmiterm([7,2,2,R31],tao/(lbd*lbd),5/6);lmiterm([7,2,2,R34],tao/(lbd*lbd),5/6);
lmiterm([7,1,3,P10],4/3,1);lmiterm([7,1,3,P11],5/6,1);lmiterm([7,1,3,P14],5/6,1);
lmiterm([7,1,3,Mz1],4/3,1);lmiterm([7,1,3,Mz1],A1',(5/6)*lbd);lmiterm([7,1,3,Mz1],A4',(5/6)*lbd);lmiterm([7,1,3,-Zs1],-C1',(5/6)*lbd);lmiterm([7,1,3,-Zs4],-C1',(5/6)*lbd);
lmiterm([7,1,3,R10],-tao/lbd,4/3);lmiterm([7,1,3,R11],-tao/lbd,5/6);lmiterm([7,1,3,R14],-tao/lbd,5/6);

lmiterm([7,1,4,P20],4/3,1);
lmiterm([7,1,4,R20],-tao/lbd,4/3);lmiterm([7,1,4,-Zd1],-A4'*C1',(5/6)*lbd);lmiterm([7,1,4,-Zd4],-A1'*C1',(5/6)*lbd);
lmiterm([7,2,3,-P20],4/3,1);
lmiterm([7,2,3,-R20],-tao/lbd,4/3);lmiterm([7,2,3,Mz1],B1',(5/6)*lbd);lmiterm([7,2,3,Mz1],B4',(5/6)*lbd);
lmiterm([7,2,4,P30],4/3,1);lmiterm([7,2,4,P31],5/6,1);lmiterm([7,2,4,P34],5/6,1);lmiterm([7,2,4,-Zd1],-B4'*C1',(5/6)*lbd);lmiterm([7,2,4,-Zd4],-B1'*C1',(5/6)*lbd);
lmiterm([7,2,4,Mz2],4/3,1);
lmiterm([7,2,4,R30],-tao/lbd,4/3);lmiterm([7,2,4,R31],-tao/lbd,5/6);lmiterm([7,2,4,R34],-tao/lbd,5/6);
lmiterm([7,1,5,0],0);
lmiterm([7,1,6,0],0);
lmiterm([7,2,5,0],0);
lmiterm([7,2,6,0],0);
lmiterm([7,1,7,0],0);
lmiterm([7,2,7,0],0);
lmiterm([7,1,8,0],0);
lmiterm([7,1,9,0],0);
lmiterm([7,2,8,0],0);
lmiterm([7,2,9,0],0);

lmiterm([7,3,3,Mz1],-lbd,4/3,'s');
lmiterm([7,3,3,R10],tao,4/3);lmiterm([7,3,3,R11],tao,5/6);lmiterm([7,3,3,R14],tao,5/6);
lmiterm([7,3,4,R20],tao,4/3);
lmiterm([7,4,4,Mz2],-lbd,4/3,'s');
lmiterm([7,4,4,R30],tao,4/3);lmiterm([7,4,4,R31],tao,5/6);lmiterm([7,4,4,R34],tao,5/6);
lmiterm([7,3,5,0],0);
lmiterm([7,3,6,0],0);
lmiterm([7,4,5,0],0);
lmiterm([7,4,6,0],0);
lmiterm([7,3,7,0],0);
lmiterm([7,4,7,-Mz2],4/3,lbd);lmiterm([7,3,8,Zs1],-(5/6)*lbd,C1);lmiterm([7,3,8,Zs4],-(5/6)*lbd,C1);
lmiterm([7,3,9,0],0);lmiterm([7,4,8,Zd1],-(5/6)*lbd,C1*A4);lmiterm([7,4,8,Zd4],-(5/6)*lbd,C1*A1);lmiterm([7,4,9,Zd1],-(5/6)*lbd,C1*B4);lmiterm([7,4,9,Zd4],-(5/6)*lbd,C1*B1);

lmiterm([7,5,5,Q10],-4/3,1);lmiterm([7,5,5,Q11],-5/6,1);lmiterm([7,5,5,Q14],-5/6,1);
lmiterm([7,5,6,Q20],-4/3,1);
lmiterm([7,6,6,Q30],-4/3,1);lmiterm([7,6,6,Q31],-5/6,1);lmiterm([7,6,6,Q34],-5/6,1);
lmiterm([7,5,7,0],0);
lmiterm([7,6,7,0],0);
lmiterm([7,5,8,0],0);lmiterm([7,5,9,0],0);
lmiterm([7,6,8,0],0);lmiterm([7,6,9,0],0);

lmiterm([7,7,7,0],-4/3*gama1*gama1);
lmiterm([7,7,8,0],0);
lmiterm([7,7,9,0],0);

lmiterm([7,8,8,R10],-1/tao,4/3);lmiterm([7,8,8,R11],-1/tao,5/6);lmiterm([7,8,8,R14],-1/tao,5/6);
lmiterm([7,8,9,R20],-1/tao,4/3);
lmiterm([7,9,9,R30],-1/tao,4/3);lmiterm([7,9,9,R31],-1/tao,5/6);lmiterm([7,9,9,R34],-1/tao,5/6);


lmiterm([8,1,1,Q10],4/3,1);lmiterm([8,1,1,Q12],5/6,1);lmiterm([8,1,1,Q13],5/6,1);
lmiterm([8,1,1,P10],-2/lbd,4/3);lmiterm([8,1,1,P12],-2/lbd,5/6);lmiterm([8,1,1,P13],-2/lbd,5/6);
lmiterm([8,1,1,R10],tao/(lbd*lbd),4/3);lmiterm([8,1,1,R12],tao/(lbd*lbd),5/6);lmiterm([8,1,1,R13],tao/(lbd*lbd),5/6);
lmiterm([8,1,2,Q20],4/3,1);
lmiterm([8,1,2,P20],-2/lbd,4/3);
lmiterm([8,1,2,R20],tao/(lbd*lbd),4/3);
lmiterm([8,2,2,Q30],4/3,1);lmiterm([8,2,2,Q32],5/6,1);lmiterm([8,2,2,Q33],5/6,1);
lmiterm([8,2,2,0],4/3);
lmiterm([8,2,2,P30],-2/lbd,4/3);lmiterm([8,2,2,P32],-2/lbd,5/6);lmiterm([8,2,2,P33],-2/lbd,5/6);
lmiterm([8,2,2,R30],tao/(lbd*lbd),4/3);lmiterm([8,2,2,R32],tao/(lbd*lbd),5/6);lmiterm([8,2,2,R33],tao/(lbd*lbd),5/6);
lmiterm([8,1,3,P10],4/3,1);lmiterm([8,1,3,P12],5/6,1);lmiterm([8,1,3,P13],5/6,1);
lmiterm([8,1,3,Mz1],4/3,1);lmiterm([8,1,3,Mz1],A2',(5/6)*lbd);lmiterm([8,1,3,Mz1],A3',(5/6)*lbd);lmiterm([8,1,3,-Zs2],-C1',(5/6)*lbd);lmiterm([8,1,3,-Zs3],-C1',(5/6)*lbd);
lmiterm([8,1,3,R10],-tao/lbd,4/3);lmiterm([8,1,3,R12],-tao/lbd,5/6);lmiterm([8,1,3,R13],-tao/lbd,5/6);

lmiterm([8,1,4,P20],4/3,1);
lmiterm([8,1,4,R20],-tao/lbd,4/3);lmiterm([8,1,4,-Zd2],-A3'*C1',(5/6)*lbd);lmiterm([8,1,4,-Zd3],-A2'*C1',(5/6)*lbd);
lmiterm([8,2,3,-P20],4/3,1);
lmiterm([8,2,3,-R20],-tao/lbd,4/3);lmiterm([8,2,3,Mz1],B2',(5/6)*lbd);lmiterm([8,2,3,Mz1],B3',(5/6)*lbd);
lmiterm([8,2,4,P30],4/3,1);lmiterm([8,2,4,P32],5/6,1);lmiterm([8,2,4,P33],5/6,1);lmiterm([8,2,4,-Zd2],-B3'*C1',(5/6)*lbd);lmiterm([8,2,4,-Zd3],-B2'*C1',(5/6)*lbd);
lmiterm([8,2,4,Mz2],4/3,1);
lmiterm([8,2,4,R30],-tao/lbd,4/3);lmiterm([8,2,4,R32],-tao/lbd,5/6);lmiterm([8,2,4,R33],-tao/lbd,5/6);
lmiterm([8,1,5,0],0);
lmiterm([8,1,6,0],0);
lmiterm([8,2,5,0],0);
lmiterm([8,2,6,0],0);
lmiterm([8,1,7,0],0);
lmiterm([8,2,7,0],0);
lmiterm([8,1,8,0],0);
lmiterm([8,1,9,0],0);
lmiterm([8,2,8,0],0);
lmiterm([8,2,9,0],0);

lmiterm([8,3,3,Mz1],-lbd,4/3,'s');
lmiterm([8,3,3,R10],tao,4/3);lmiterm([8,3,3,R12],tao,5/6);lmiterm([8,3,3,R13],tao,5/6);
lmiterm([8,3,4,R20],tao,4/3);
lmiterm([8,4,4,Mz2],-lbd,4/3,'s');
lmiterm([8,4,4,R30],tao,4/3);lmiterm([8,4,4,R32],tao,5/6);lmiterm([8,4,4,R33],tao,5/6);
lmiterm([8,3,5,0],0);
lmiterm([8,3,6,0],0);
lmiterm([8,4,5,0],0);
lmiterm([8,4,6,0],0);
lmiterm([8,3,7,0],0);
lmiterm([8,4,7,-Mz2],4/3,lbd);lmiterm([8,3,8,Zs2],-(5/6)*lbd,C1);lmiterm([8,3,8,Zs3],-(5/6)*lbd,C1);
lmiterm([8,3,9,0],0);lmiterm([8,4,8,Zd2],-(5/6)*lbd,C1*A3);lmiterm([8,4,8,Zd3],-(5/6)*lbd,C1*A2);lmiterm([8,4,9,Zd2],-(5/6)*lbd,C1*B3);lmiterm([8,4,9,Zd3],-(5/6)*lbd,C1*B2);

lmiterm([8,5,5,Q10],-4/3,1);lmiterm([8,5,5,Q12],-5/6,1);lmiterm([8,5,5,Q13],-5/6,1);
lmiterm([8,5,6,Q20],-4/3,1);
lmiterm([8,6,6,Q30],-4/3,1);lmiterm([8,6,6,Q32],-5/6,1);lmiterm([8,6,6,Q33],-5/6,1);
lmiterm([8,5,7,0],0);
lmiterm([8,6,7,0],0);
lmiterm([8,5,8,0],0);lmiterm([8,5,9,0],0);
lmiterm([8,6,8,0],0);lmiterm([8,6,9,0],0);

lmiterm([8,7,7,0],-4/3*gama1*gama1);
lmiterm([8,7,8,0],0);
lmiterm([8,7,9,0],0);

lmiterm([8,8,8,R10],-1/tao,4/3);lmiterm([8,8,8,R12],-1/tao,5/6);lmiterm([8,8,8,R13],-1/tao,5/6);
lmiterm([8,8,9,R20],-1/tao,4/3);
lmiterm([8,9,9,R30],-1/tao,4/3);lmiterm([8,9,9,R32],-1/tao,5/6);lmiterm([8,9,9,R33],-1/tao,5/6);



lmiterm([9,1,1,Q10],4/3,1);lmiterm([9,1,1,Q12],5/6,1);lmiterm([9,1,1,Q14],5/6,1);
lmiterm([9,1,1,P10],-2/lbd,4/3);lmiterm([9,1,1,P12],-2/lbd,5/6);lmiterm([9,1,1,P14],-2/lbd,5/6);
lmiterm([9,1,1,R10],tao/(lbd*lbd),4/3);lmiterm([9,1,1,R12],tao/(lbd*lbd),5/6);lmiterm([9,1,1,R14],tao/(lbd*lbd),5/6);
lmiterm([9,1,2,Q20],4/3,1);
lmiterm([9,1,2,P20],-2/lbd,4/3);
lmiterm([9,1,2,R20],tao/(lbd*lbd),4/3);
lmiterm([9,2,2,Q30],4/3,1);lmiterm([9,2,2,Q32],5/6,1);lmiterm([9,2,2,Q34],5/6,1);
lmiterm([9,2,2,0],4/3);
lmiterm([9,2,2,P30],-2/lbd,4/3);lmiterm([9,2,2,P32],-2/lbd,5/6);lmiterm([9,2,2,P34],-2/lbd,5/6);
lmiterm([9,2,2,R30],tao/(lbd*lbd),4/3);lmiterm([9,2,2,R32],tao/(lbd*lbd),5/6);lmiterm([9,2,2,R34],tao/(lbd*lbd),5/6);
lmiterm([9,1,3,P10],4/3,1);lmiterm([9,1,3,P12],5/6,1);lmiterm([9,1,3,P14],5/6,1);
lmiterm([9,1,3,Mz1],4/3,1);lmiterm([9,1,3,Mz1],A2',(5/6)*lbd);lmiterm([9,1,3,Mz1],A4',(5/6)*lbd);lmiterm([9,1,3,-Zs2],-C1',(5/6)*lbd);lmiterm([9,1,3,-Zs4],-C1',(5/6)*lbd);
lmiterm([9,1,3,R10],-tao/lbd,4/3);lmiterm([9,1,3,R12],-tao/lbd,5/6);lmiterm([9,1,3,R14],-tao/lbd,5/6);

lmiterm([9,1,4,P20],4/3,1);
lmiterm([9,1,4,R20],-tao/lbd,4/3);lmiterm([9,1,4,-Zd2],-A4'*C1',(5/6)*lbd);lmiterm([9,1,4,-Zd4],-A2'*C1',(5/6)*lbd);
lmiterm([9,2,3,-P20],4/3,1);
lmiterm([9,2,3,-R20],-tao/lbd,4/3);lmiterm([9,2,3,Mz1],B2',(5/6)*lbd);lmiterm([9,2,3,Mz1],B4',(5/6)*lbd);
lmiterm([9,2,4,P30],4/3,1);lmiterm([9,2,4,P32],5/6,1);lmiterm([9,2,4,P34],5/6,1);lmiterm([9,2,4,-Zd2],-B4'*C1',(5/6)*lbd);lmiterm([9,2,4,-Zd4],-B2'*C1',(5/6)*lbd);
lmiterm([9,2,4,Mz2],4/3,1);
lmiterm([9,2,4,R30],-tao/lbd,4/3);lmiterm([9,2,4,R32],-tao/lbd,5/6);lmiterm([9,2,4,R34],-tao/lbd,5/6);
lmiterm([9,1,5,0],0);
lmiterm([9,1,6,0],0);
lmiterm([9,2,5,0],0);
lmiterm([9,2,6,0],0);
lmiterm([9,1,7,0],0);
lmiterm([9,2,7,0],0);
lmiterm([9,1,8,0],0);
lmiterm([9,1,9,0],0);
lmiterm([9,2,8,0],0);
lmiterm([9,2,9,0],0);

lmiterm([9,3,3,Mz1],-lbd,4/3,'s');
lmiterm([9,3,3,R10],tao,4/3);lmiterm([9,3,3,R12],tao,5/6);lmiterm([9,3,3,R14],tao,5/6);
lmiterm([9,3,4,R20],tao,4/3);
lmiterm([9,4,4,Mz2],-lbd,4/3,'s');
lmiterm([9,4,4,R30],tao,4/3);lmiterm([9,4,4,R32],tao,5/6);lmiterm([9,4,4,R34],tao,5/6);
lmiterm([9,3,5,0],0);
lmiterm([9,3,6,0],0);
lmiterm([9,4,5,0],0);
lmiterm([9,4,6,0],0);
lmiterm([9,3,7,0],0);
lmiterm([9,4,7,-Mz2],4/3,lbd);lmiterm([9,3,8,Zs2],-(5/6)*lbd,C1);lmiterm([9,3,8,Zs4],-(5/6)*lbd,C1);
lmiterm([9,3,9,0],0);lmiterm([9,4,8,Zd2],-(5/6)*lbd,C1*A4);lmiterm([9,4,8,Zd4],-(5/6)*lbd,C1*A2);lmiterm([9,4,9,Zd2],-(5/6)*lbd,C1*B4);lmiterm([9,4,9,Zd4],-(5/6)*lbd,C1*B2);

lmiterm([9,5,5,Q10],-4/3,1);lmiterm([9,5,5,Q12],-5/6,1);lmiterm([9,5,5,Q14],-5/6,1);
lmiterm([9,5,6,Q20],-4/3,1);
lmiterm([9,6,6,Q30],-4/3,1);lmiterm([9,6,6,Q32],-5/6,1);lmiterm([9,6,6,Q34],-5/6,1);
lmiterm([9,5,7,0],0);
lmiterm([9,6,7,0],0);
lmiterm([9,5,8,0],0);lmiterm([9,5,9,0],0);
lmiterm([9,6,8,0],0);lmiterm([9,6,9,0],0);

lmiterm([9,7,7,0],-4/3*gama1*gama1);
lmiterm([9,7,8,0],0);
lmiterm([9,7,9,0],0);

lmiterm([9,8,8,R10],-1/tao,4/3);lmiterm([9,8,8,R12],-1/tao,5/6);lmiterm([9,8,8,R14],-1/tao,5/6);
lmiterm([9,8,9,R20],-1/tao,4/3);
lmiterm([9,9,9,R30],-1/tao,4/3);lmiterm([9,9,9,R32],-1/tao,5/6);lmiterm([9,9,9,R34],-1/tao,5/6);




lmiterm([10,1,1,Q10],4/3,1);lmiterm([10,1,1,Q13],5/6,1);lmiterm([10,1,1,Q14],5/6,1);
lmiterm([10,1,1,P10],-2/lbd,4/3);lmiterm([10,1,1,P13],-2/lbd,5/6);lmiterm([10,1,1,P14],-2/lbd,5/6);
lmiterm([10,1,1,R10],tao/(lbd*lbd),4/3);lmiterm([10,1,1,R13],tao/(lbd*lbd),5/6);lmiterm([10,1,1,R14],tao/(lbd*lbd),5/6);
lmiterm([10,1,2,Q20],4/3,1);
lmiterm([10,1,2,P20],-2/lbd,4/3);
lmiterm([10,1,2,R20],tao/(lbd*lbd),4/3);
lmiterm([10,2,2,Q30],4/3,1);lmiterm([10,2,2,Q33],5/6,1);lmiterm([10,2,2,Q34],5/6,1);
lmiterm([10,2,2,0],4/3);
lmiterm([10,2,2,P30],-2/lbd,4/3);lmiterm([10,2,2,P33],-2/lbd,5/6);lmiterm([10,2,2,P34],-2/lbd,5/6);
lmiterm([10,2,2,R30],tao/(lbd*lbd),4/3);lmiterm([10,2,2,R33],tao/(lbd*lbd),5/6);lmiterm([10,2,2,R34],tao/(lbd*lbd),5/6);
lmiterm([10,1,3,P10],4/3,1);lmiterm([10,1,3,P13],5/6,1);lmiterm([10,1,3,P14],5/6,1);
lmiterm([10,1,3,Mz1],4/3,1);lmiterm([10,1,3,Mz1],A3',(5/6)*lbd);lmiterm([10,1,3,Mz1],A4',(5/6)*lbd);lmiterm([10,1,3,-Zs3],-C1',(5/6)*lbd);lmiterm([10,1,3,-Zs4],-C1',(5/6)*lbd);
lmiterm([10,1,3,R10],-tao/lbd,4/3);lmiterm([10,1,3,R13],-tao/lbd,5/6);lmiterm([10,1,3,R14],-tao/lbd,5/6);

lmiterm([10,1,4,P20],4/3,1);
lmiterm([10,1,4,R20],-tao/lbd,4/3);lmiterm([10,1,4,-Zd3],-A4'*C1',(5/6)*lbd);lmiterm([10,1,4,-Zd4],-A3'*C1',(5/6)*lbd);
lmiterm([10,2,3,-P20],4/3,1);
lmiterm([10,2,3,-R20],-tao/lbd,4/3);lmiterm([10,2,3,Mz1],B3',(5/6)*lbd);lmiterm([10,2,3,Mz1],B4',(5/6)*lbd);
lmiterm([10,2,4,P30],4/3,1);lmiterm([10,2,4,P33],5/6,1);lmiterm([10,2,4,P34],5/6,1);lmiterm([10,2,4,-Zd3],-B4'*C1',(5/6)*lbd);lmiterm([10,2,4,-Zd4],-B3'*C1',(5/6)*lbd);
lmiterm([10,2,4,Mz2],4/3,1);
lmiterm([10,2,4,R30],-tao/lbd,4/3);lmiterm([10,2,4,R33],-tao/lbd,5/6);lmiterm([10,2,4,R34],-tao/lbd,5/6);
lmiterm([10,1,5,0],0);
lmiterm([10,1,6,0],0);
lmiterm([10,2,5,0],0);
lmiterm([10,2,6,0],0);
lmiterm([10,1,7,0],0);
lmiterm([10,2,7,0],0);
lmiterm([10,1,8,0],0);
lmiterm([10,1,9,0],0);
lmiterm([10,2,8,0],0);
lmiterm([10,2,9,0],0);

lmiterm([10,3,3,Mz1],-lbd,4/3,'s');
lmiterm([10,3,3,R10],tao,4/3);lmiterm([10,3,3,R13],tao,5/6);lmiterm([10,3,3,R14],tao,5/6);
lmiterm([10,3,4,R20],tao,4/3);
lmiterm([10,4,4,Mz2],-lbd,4/3,'s');
lmiterm([10,4,4,R30],tao,4/3);lmiterm([10,4,4,R33],tao,5/6);lmiterm([10,4,4,R34],tao,5/6);
lmiterm([10,3,5,0],0);
lmiterm([10,3,6,0],0);
lmiterm([10,4,5,0],0);
lmiterm([10,4,6,0],0);
lmiterm([10,3,7,0],0);
lmiterm([10,4,7,-Mz2],4/3,lbd);lmiterm([10,3,8,Zs3],-(5/6)*lbd,C1);lmiterm([10,3,8,Zs4],-(5/6)*lbd,C1);
lmiterm([10,3,9,0],0);lmiterm([10,4,8,Zd3],-(5/6)*lbd,C1*A4);lmiterm([10,4,8,Zd4],-(5/6)*lbd,C1*A3);lmiterm([10,4,9,Zd3],-(5/6)*lbd,C1*B4);lmiterm([10,4,9,Zd4],-(5/6)*lbd,C1*B3);

lmiterm([10,5,5,Q10],-4/3,1);lmiterm([10,5,5,Q13],-5/6,1);lmiterm([10,5,5,Q14],-5/6,1);
lmiterm([10,5,6,Q20],-4/3,1);
lmiterm([10,6,6,Q30],-4/3,1);lmiterm([10,6,6,Q33],-5/6,1);lmiterm([10,6,6,Q34],-5/6,1);
lmiterm([10,5,7,0],0);
lmiterm([10,6,7,0],0);
lmiterm([10,5,8,0],0);lmiterm([10,5,9,0],0);
lmiterm([10,6,8,0],0);lmiterm([10,6,9,0],0);

lmiterm([10,7,7,0],-4/3*gama1*gama1);
lmiterm([10,7,8,0],0);
lmiterm([10,7,9,0],0);

lmiterm([10,8,8,R10],-1/tao,4/3);lmiterm([10,8,8,R13],-1/tao,5/6);lmiterm([10,8,8,R14],-1/tao,5/6);
lmiterm([10,8,9,R20],-1/tao,4/3);
lmiterm([10,9,9,R30],-1/tao,4/3);lmiterm([10,9,9,R33],-1/tao,5/6);lmiterm([10,9,9,R34],-1/tao,5/6);



lmiterm([-11,1,1,P10],1,1);lmiterm([-11,1,1,P11],1,1);lmiterm([-11,1,2,P20],1,1);lmiterm([-11,2,2,P30],1,1);lmiterm([-11,2,2,P31],1,1);
lmiterm([-12,1,1,P10],1,1);lmiterm([-12,1,1,P12],1,1);lmiterm([-12,1,2,P20],1,1);lmiterm([-12,2,2,P30],1,1);lmiterm([-12,2,2,P32],1,1);
lmiterm([-13,1,1,P10],1,1);lmiterm([-13,1,1,P13],1,1);lmiterm([-13,1,2,P20],1,1);lmiterm([-13,2,2,P30],1,1);lmiterm([-13,2,2,P33],1,1);
lmiterm([-14,1,1,P10],1,1);lmiterm([-14,1,1,P14],1,1);lmiterm([-14,1,2,P20],1,1);lmiterm([-14,2,2,P30],1,1);lmiterm([-14,2,2,P34],1,1);


lmiterm([-15,1,1,Q10],1,1);lmiterm([-15,1,1,Q11],1,1);lmiterm([-15,1,2,Q20],1,1);lmiterm([-15,2,2,Q30],1,1);lmiterm([-15,2,2,Q31],1,1);
lmiterm([-16,1,1,Q10],1,1);lmiterm([-16,1,1,Q12],1,1);lmiterm([-16,1,2,Q20],1,1);lmiterm([-16,2,2,Q30],1,1);lmiterm([-16,2,2,Q32],1,1);
lmiterm([-17,1,1,Q10],1,1);lmiterm([-17,1,1,Q13],1,1);lmiterm([-17,1,2,Q20],1,1);lmiterm([-17,2,2,Q30],1,1);lmiterm([-17,2,2,Q33],1,1);
lmiterm([-18,1,1,Q10],1,1);lmiterm([-18,1,1,Q14],1,1);lmiterm([-18,1,2,Q20],1,1);lmiterm([-18,2,2,Q30],1,1);lmiterm([-18,2,2,Q34],1,1);


lmiterm([-19,1,1,R10],1,1);lmiterm([-19,1,1,R11],1,1);lmiterm([-19,1,2,R20],1,1);lmiterm([-19,2,2,R30],1,1);lmiterm([-19,2,2,R31],1,1);
lmiterm([-20,1,1,R10],1,1);lmiterm([-20,1,1,R12],1,1);lmiterm([-20,1,2,R20],1,1);lmiterm([-20,2,2,R30],1,1);lmiterm([-20,2,2,R32],1,1);
lmiterm([-21,1,1,R10],1,1);lmiterm([-21,1,1,R13],1,1);lmiterm([-21,1,2,R20],1,1);lmiterm([-21,2,2,R30],1,1);lmiterm([-21,2,2,R33],1,1);
lmiterm([-22,1,1,R10],1,1);lmiterm([-22,1,1,R14],1,1);lmiterm([-22,1,2,R20],1,1);lmiterm([-22,2,2,R30],1,1);lmiterm([-22,2,2,R34],1,1);



lmis=getlmis;   
[tmin,xfeas]=feasp(lmis);
tmin

M1=dec2mat(lmis,xfeas,Mz1);
M2=dec2mat(lmis,xfeas,Mz2);

Kd1=dec2mat(lmis,xfeas,Zd1);
Kd2=dec2mat(lmis,xfeas,Zd2);
Kd3=dec2mat(lmis,xfeas,Zd3);
Kd4=dec2mat(lmis,xfeas,Zd4);


Ks1=dec2mat(lmis,xfeas,Zs1);
Ks2=dec2mat(lmis,xfeas,Zs2);
Ks3=dec2mat(lmis,xfeas,Zs3);
Ks4=dec2mat(lmis,xfeas,Zs4);

P10m=dec2mat(lmis,xfeas,P10);
P11m=dec2mat(lmis,xfeas,P11);

Ls1=inv(M1')*Ks1
Ls2=inv(M1')*Ks2
Ls3=inv(M1')*Ks3
Ls4=inv(M1')*Ks4


Ld1=inv(M2')*Kd1
Ld2=inv(M2')*Kd2
Ld3=inv(M2')*Kd3
Ld4=inv(M2')*Kd4
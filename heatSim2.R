#vypocita teplotu v bode za predpokladu dane teploty okolnich bodu, vykonu, teplotnich difuzivit, tepelne kapacity, dx a dt
library(Rcpp)

cppFunction(
"double calcTemp(double t1,double t2,double t3,double p2, double td2, double vhc2, double dx,double dt){
  return(t2+(td2*(t3-2*t2+t1)/(dx*dx)+p2/vhc2)*dt);
}"
)

cppFunction(
  "double calcTemp(double t1,double t2,double t3,double p2, double td2, double vhc2, double dx,double dt){
  double r;
  r=td2*dt/(dx*dx);  
  return(r*t3+(1-2*r)*t2+r*t1+dt*p2/vhc2);
  }"
)

generateNNextTemps<-function(t0,p,td,vhc,dx,dt,nt,nx){
    t<-array(data=0,dim=c(nt,nx)) #make array
    t[1,]<-t0 #start condition
    t[,1]<-t0[1] #border condition
    t[,nx]<-t0[nx] #border condition
    
    if(max(td*dt/(dx^2)>.5)){
      warning(paste(" Je-li hodnota alpha*dt/(dx^2) vetsi nez 0.5, pak FTCS bude nestabilni. Ted je",max(td*dt/(dx^2))))
    }
    
    for(i in 1:(nt-1)){#loop over time
        for(j in 2:(nx-1)){#loop over space
            t[i+1,j]<-calcTemp(t[i,j-1],t[i,j],t[i,j+1],
                p[i,j],
                td[j],
                vhc[j],
                dx,dt)
        }
    }
    t
}

generateNNextTempsBTCS<-function(t0,p,td,vhc,dx,dt,nt,nx){
    t<-array(data=0,dim=c(nt,nx)) #make array
    t[1,]<-t0 #start condition
    #t[,1]<-t0[1] #border condition
    #t[,nx]<-t0[nx] #border condition

    A<-matrix(0,nrow = nx,ncol = nx)
    B<-matrix(0,nrow = nx,ncol = 1)
    
    for(i in 2:(nrow(A)-1)){
        A[i,i-1]<- -td[i]/(dx*dx)
        A[i,i]<-1/dt+2*td[i]/(dx*dx)
        A[i,i+1]<- -td[i]/(dx*dx)
    }
    
    A[1,1]<-1
    A[1,2]<-0
    A[nrow(A),nrow(A)-1]<-0
    A[nrow(A),nrow(A)]<-1

    for(i in 2:nt){#loop over time
        B<-1/dt* t[i-1,]+p[i-1,]/vhc
        B[1]<-t0[1]#
        B[nx]<-t0[nx]#
        t[i,]<-solve(A,B)
        
        t[i,1]<-t0[1]#
        t[i,nx]<-t0[nx]
    }
    t
}


dt<-0.0001
dx<-0.1

tmax<-1
l<-2

nt<-round(tmax/dt)
nx<-round(l/dx)

t0<-c(rep(40,nx)) #vector of initial temperatures
# t0<-sin(rep(seq(from=0,to=pi,length=nx)))

td<-rep(8.8,nx)#vector of thermal diffusities (alpha; mm^2/s), si: 88
vhc<-c(rep(1.72E-3,nx))#volumetric heat capacity (specific heat capacity * density) (J/(mm^3*K)) Si: Michal:1,641E6, Ilja 1,72E6...

cyklus<-matrix(data=c(1000*c(11,0,11,0,12,0,12,0,9.2,0),#vykon v elementu tloustky dx
                      c(900,2080,920,4080,900,1100,900,3100,500,28000)),#pomerny cas proudu
               ncol=2)

active<-matrix(data=c(c(0,1,0),#pomer vykonu
                      c(10,1,10)),#pomerny rozmer
               ncol=2)

genPCyklus<-function(nt,nx,dx,cyklus,active){
    #prepare array of volumetric power densities (time x space)(q, J/(s*mm^3))
    
    k<-sum(cyklus[,2])
    tk<-ceiling(cyklus[,2]/k*nt)
    
    k2<-sum(active[,2])
    tk2<-ceiling(active[,2]/k2*nx)
    
    p2<-matrix(data=0,nrow=1,ncol=sum(tk))
    p2<-rep(cyklus[,1],tk)[1:nt]

    p3<-matrix(data=0,nrow=1,ncol=sum(tk2))
    p3<-rep(active[,1],tk2)[1:nx]
    
    p4<-matrix(data=p2,nrow=nt,ncol=nx)
    p5<-matrix(data=p3,nrow=nt,ncol=nx,byrow = T)
    
    p<-p4*p5
}

p<-genPCyklus(nt,nx,dx,cyklus,active)


Rprof(tmp <- tempfile())

tArr<-generateNNextTemps(t0, p,td,vhc,dx,dt,nt,nx)
tArr2<-generateNNextTempsBTCS(t0, p,td,vhc,dx,dt,nt,nx)

Rprof()
summaryRprof(tmp)

plot(tArr[,round(nx/2)],type="l")
lines(tArr2[,round(nx/2)],col="red")

lines(tArr[,round(nx/4)])
lines(tArr2[,round(nx/4)],col="red")

lines(tArr[,round(nx/10)])
lines(tArr2[,round(nx/10)],col="red")

max(tArr)

tArr[1:200,16:24]

plot(tArr[1,],type="l")
lines(tArr[900,])

lines(tArr[3000,])
lines(tArr[5000,])

persp((tArr[seq(1,nt,length.out = 100),]))
library(plotly)
plot_ly(z=~tArr[seq(1,nt,length.out = 100),]) %>% add_surface()


#-----
#http://www.nada.kth.se/~jjalap/numme/FDheat.pdf
# nx<-200
# nt<-1000
# L<-1
# dx<-L/nx
# dt<-0.05
# 
# td<-rep(0.1,nx)#vector of thermal diffusities (alpha; mm^2/s), si: 88
# vhc<-c(rep(1.72E-3,nx))#volumetric heat capacity (specific heat capacity * density) (J/(mm^3*K)) Si: Michal:1,641E6, Ilja 1,72E6...
# p<-array(data=0,dim=c(nt,nx))#volumetric power densities (q, J/(s*mm^3))
# t0<-sin(seq(from=0,to=pi,length=nx))
# 
# Rprof(tmp <- tempfile())
# tArr<-generateNNextTempsBTCS(t0,p,td,vhc,dx,dt,nt,nx)
# Rprof()
# summaryRprof(tmp)
# 
# Rprof(tmp <- tempfile())
# tArr2<-generateNNextTemps(t0,p,td,vhc,dx,dt,nt,nx)
# Rprof()
# summaryRprof(tmp)
# 
# plot(tArr[20,],type="l")
# lines(tArr2[20,],col="red")


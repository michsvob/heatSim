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

getVoltage<-function(current,temperature){
    v0RT<-0.843;rdRT<-0.026;v0HT<-0.672;rdHT<-0.029
    t<-c(25, 150)
    v<-c(v0RT+rdRT*current,v0HT+rdHT*current)
    f<-approxfun(t,v,rule=2)#interpolace napeti dle teploty. 
    #hodnoty mimo interval se priradi extremum v intervalu (aby napeti uplne neulitavalo)
    f(temperature)
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
        
        p[i]<-p[i]*mapply(getVoltage,current=p[i],temperature=t[i])
    }
    t
}

genPCyklus<-function(nt,nx,dx,cyklus,active){
    #prepare array of volumetric power densities (time x space)(q, J/(s*mm^3))
    #given is - how many time and space elements there should be
    #         - values in time and space segments
    #         - proportions of time and space
    
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

genTD<-function(nx,td_setting){
    k<-sum(td_setting[,2])
    tk<-ceiling(td_setting[,2]/k*nx)
    td<-matrix(data=0,nrow=1,ncol=sum(tk))
    td<-rep(td_setting[,1],tk)[1:nx]
}

#funkce, ktera danemu proudu a teplote priradi napeti


dt<-0.001
dx<-0.01

tmax<-.1
l<-4

nt<-round(tmax/dt)
nx<-round(l/dx)

t0<-c(rep(20,nx)) #vector of initial temperatures
# t0<-sin(rep(seq(from=0,to=pi,length=nx)))


vhc<-c(rep(1.72E-3,nx))#volumetric heat capacity (specific heat capacity * density) (J/(mm^3*K)) Si: Michal:1,641E6, Ilja 1,72E6...
#vhc je potreba jen u vrstev, kde se generuje teplo (v heat eq. je nasobena nulou)

cyklus<-matrix(data=c(1000*c(11,0,11,0,12,0,12,0,9.2,0),#vykon v elementu tloustky dx
                      c(900,2080,920,4080,900,1100,900,3100,500,28000)),#pomerny cas 
               ncol=2)



active<-matrix(data=c(c(0,  0, 1, 0, 0,  0)/1548,#pomer proudu deleny plochou Si
                      c(15,3.4,0.25,2,3.4,15)),#pomerny rozmer
               ncol=2)
td_setting<-matrix(data=c(c(111,111,88  ,54.3,111,111),#thermal diffusivity in segment (alpha; mm^2/s), si: 88
                          c(15 ,3.4,0.25,2   ,3.4,15)),#pomerny rozmer
                   ncol=2)

plocha_setting<-matrix(data=c(c(2026,2026,1548,2026,2026,2026)/1548,#plocha segmentu relativne k Si
                            c(15 ,3.4,0.25,2  ,3.4,15)),#pomerny rozmer
                   ncol=2)

plocha<-genTD(nx,plocha_setting)   
td<-genTD(nx,td_setting)
td_adjusted<-td*plocha


p<-genPCyklus(nt,nx,dx,cyklus,active)


Rprof(tmp <- tempfile())

#tArr<-generateNNextTemps(t0, p,td,vhc,dx,dt,nt,nx)
tArr2<-generateNNextTempsBTCS(t0, p,td_adjusted,vhc,dx,dt,nt,nx)

Rprof()
summaryRprof(tmp)

plot(tArr2[,round(nx/2)],col="red",type="l")
lines(tArr2[,round(nx/4)],col="red")
lines(tArr2[,round(nx/10)],col="red")

max(tArr2)


persp((tArr2[seq(1,nt,length.out = 100),]))
library(plotly)
plot_ly(z=~tArr2[seq(1,nt,length.out = 100),]) %>% add_surface()


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


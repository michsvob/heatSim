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

generateNNextTemps<-function(t0,p,td,vhc,simT,dx,dt,len,m){
    t<-array(data=0,dim=c(len,m)) #make array
    t[1,]<-t0 #start condition
    t[,1]<-t0[1] #border condition
    t[,m]<-t0[m] #border condition
    
    if(max(td*dt/(dx^2)>.5)){
      warning(paste(" Je-li hodnota alpha*dt/(dx^2) vetsi nez 0.5, pak FTCS bude nestabilni. Ted je",max(td*dt/(dx^2))))
    }
    
    for(i in 1:(len-1)){#loop over time
        for(j in 2:(m-1)){#loop over space
            t[i+1,j]<-calcTemp(t[i,j-1],t[i,j],t[i,j+1],
                p[i,j],
                td[j],
                vhc[j],
                dx,dt)
        }
    }
    t
}

generateNNextTempsBTCS<-function(t0,p,td,vhc,simT,dx,dt,len,m){
    t<-array(data=0,dim=c(len,m)) #make array
    t[1,]<-t0 #start condition
    t[,1]<-t0[1] #border condition
    t[,m]<-t0[m] #border condition

    A<-matrix(0,nrow = m,ncol = m)
    B<-matrix(0,nrow = m,ncol = 1)
    
    for(i in 2:(nrow(A)-1)){
        A[i,i-1]<- -td[i]/(dx*dx)
        A[i,i]<-1/dt+2*td[i]/(dx*dx)
        A[i,i+1]<- -td[i]/(dx*dx)
    }
    
    A[1,1]<-1
    A[1,2]<-0
    A[nrow(A),nrow(A)-1]<-0
    A[nrow(A),nrow(A)]<-1

    for(i in 2:len){#loop over time
        B<-1/dt* t[i-1,]+p[i-1]*dt/vhc[i-1]
        B[1]<-t0[1]#
        B[m]<-t0[m]#
        t[i,]<-solve(a=A,b=B)
    }
    t
}



# i_cyklus<-matrix(data=c(1000*c(11,0,11,0,12,0,12,0,9.2,0),
#                         0.001*c(900,2080,920,4080,900,1100,900,3100,500,28000)),ncol=2)
i_cyklus<-matrix(data=c(c(11),
                        c(100)),ncol=2)
# i_cyklus<-matrix(data=c(c(0),
#                         c(900)),ncol=2)

l<-3 #celkova delka v mm
dx<-.1 #vzdalenost mezi elementarnimi bunkami v mm
m<-ceiling(l/dx) #pocet bunek

simT<-sum(i_cyklus[,2]) #doba simulace v s
dt<-1 #dt (časový element simulace v s)
len<-ceiling(simT/dt) #pocet casovych useku


t0<-c(rep(40,m)) #vector of initial temperatures
td<-rep(88,m)#vector of thermal diffusities (alpha; mm^2/s), si: 88
vhc<-c(rep(1.72E-3,m))#volumetric heat capacity (specific heat capacity * density) (J/(mm^3*K)) Si: Michal:1,641E6, Ilja 1,72E6...

print(max(td*dt/(dx^2)))

u<-2 #ubytek (V)
tloustkaPN<-.05 #tloustka pn prechodu v mm
objemPN<-0.05*3000 #objem pn prechodu v mm^3
k<-tloustkaPN/dx #koeficient, kterym se nasobi objemovy vykon, 
# aby celkova energie mohla byt prirazena jedne elementarni jednotce (?)
  
p_cyklus<-k * i_cyklus[,1]*(u/objemPN) #vykon na mm^3, zde proud*napeti/objem PN prechodu

p<-array(data=0,dim=c(len,m))#volumetric power densities (q, J/(s*mm^3))

#naplneni vektoru p
a<-c()
for(i in 1:length(i_cyklus[,2])){
  a<-c(a,rep(p_cyklus[i],i_cyklus[i,2]/(1000*dt)))
}
p[,round(m/2)]<-a

Rprof(tmp <- tempfile())
tArr<-generateNNextTempsBTCS(t0,p,td,vhc,simT,dx,dt,len,m)
Rprof()
summaryRprof(tmp)

plot(tArr[,round(m/2)],type="l")
lines(tArr[,round(m/4)])
lines(tArr[,round(m/10)])
lines(tArr[,2])
max(tArr)

tArr[1:200,16:24]

plot(tArr[1,],type="l")
lines(tArr[900,])

lines(tArr[3000,])
lines(tArr[5000,])

persp((tArr[seq(1,len,length.out = 100),]))
library(plotly)
plot_ly(z=~tArr[seq(1,len,length.out = 100),]) %>% add_surface()


#-----
#http://www.nada.kth.se/~jjalap/numme/FDheat.pdf
m<-200
len<-1000
L<-1
dx<-L/m
dt<-0.05
simT<-len*dt
td<-rep(0.1,m)#vector of thermal diffusities (alpha; mm^2/s), si: 88
vhc<-c(rep(1.72E-3,m))#volumetric heat capacity (specific heat capacity * density) (J/(mm^3*K)) Si: Michal:1,641E6, Ilja 1,72E6...
p<-array(data=0,dim=c(len,m))#volumetric power densities (q, J/(s*mm^3))
t0<-sin(seq(from=0,to=pi,length=m))

Rprof(tmp <- tempfile())
tArr<-generateNNextTempsBTCS(t0,p,td,vhc,simT,dx,dt,len,m)
Rprof()
summaryRprof(tmp)

Rprof(tmp <- tempfile())
tArr2<-generateNNextTemps(t0,p,td,vhc,simT,dx,dt,len,m)
Rprof()
summaryRprof(tmp)

plot(tArr[20,],type="l")
lines(tArr2[20,],col="red")


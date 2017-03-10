#vypocita teplotu v bode za predpokladu dane teploty okolnich bodu, vykonu, teplotnich difuzivit, tepelne kapacity, dx a dt
library(Rcpp)

cppFunction("
    double calcTemp(double t1,double t2,double t3,double p2, double td2, double vhc2, double dx,double dt){
            return(t2+(td2*(t3-2*t2+t1)/(dx*dx)+p2/vhc2)*dt);
            }        
            ")

#vypocita teplotu v celem poli bodu
generateNextTemps<-function(t,p,td,vhc,dx,dt){
    t2<-double(length(t))
    t2[1]<-t[1] #border condition
    t2[length(t)]<-t[length(t)] #border condition
    
    for(i in 2:(length(t)-1)){
        t2[i]<-calcTemp(t[i-1],t[i],t[i+1],
                        p[i],
                        td[i],
                        vhc[i],
                        dx,dt)
    }
    t2
}

generateNNextTemps<-function(t0,p,td,vhc,simT,dx,dt){
    n<-simT/dt
    t<-array(dim=c(n,length(t0)))
    t[1,]<-t0
    
    for(i in 1:(n-1)){
        t[i+1,]<-generateNextTemps(t[i,],p[i,],td,vhc,dx,dt)
    }
    t
}

t0<-c(rep(20,40)) #vector of initial temperatures
td<-rep(88,40)#vector of thermal diffusities (alpha; mm^2/s)
vhc<-c(rep(1.72E-3,40))#volumetric heat capacity (specific heat capacity * density) (J/(mm^3*K)) Si: Michal:1,641E6, Ilja 1,72E6...


dx<-1 #vzdalenost mezi elementarnimi bunkami v mm
dt<-.0001 #dt (časový element simulace v s)

u<-2 #ubytek (V)

i_cyklus<-matrix(data=c(c(11,0,11,0,12,0,12,0,9.2,0),
                     c(900,2080,920,4080,900,1100,900,3100,500,28000)),ncol=2)
i_cyklus<-matrix(data=c(c(11),
                        c(900)),ncol=2)

i_cyklus[,2]<-i_cyklus[,2]/sum(i_cyklus[,2])
#sloupec 1 je proud, sloupec 2 je proporce casu

p_cyklus<-i_cyklus[,1]*u*dx

simT<-3 #doba simulace 
p<-array(data=0,dim=c(simT/dt,length(t0)))#volumetric power densities (q, J/(s*mm^3))

#naplneni vektoru p
j<-1
len<-length(p[,20])
k<-i_cyklus[1,2]*len
for(i in 1:len){
    if(i>k){
        j<-j+1
        k<-k+i_cyklus[j,2]*len
        }
    p[i,20]<-p_cyklus[j]
}


Rprof(tmp <- tempfile())
tArr<-generateNNextTemps(t0,p,td,vhc,simT,dx,dt)
Rprof()
summaryRprof(tmp)

plot(tArr[,20],type="l")
lines(tArr[,21])
lines(tArr[,10])
lines(tArr[,2])
max(tArr)
tArr[1:200,16:24]





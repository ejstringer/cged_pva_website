# https://rpubs.com/jrut/853787

#install.packages('pedigreemm')
library(pedigreemm)
#library(doBy)
library(pedigreemm)
ped<- data.frame(id=c('F', 'G', 'H', 'I', 'J'), 
                 p1=c(NA, 'F', 'F', 'G', 'H'),
                 p2= c(NA, NA, NA, 'G', 'H'))
ped

ped2<- editPed(ped$p1, ped$p2, ped$id)
ped2


ped3<- pedigree(ped2$sire, ped2$dam, ped2$label)
ped3
A<- getA(ped3)
A<- as.matrix(A) #converts the A matrix to a regular matrix
A


nbase<- 1000
base<- data.frame(id=c(1:nbase), p1=NA, p2=NA, allele1=NA, allele2=NA, gen=0)
base$allele1<- paste('A', c(1:nbase), "-1", sep="")
base$allele2<- paste('A', c(1:nbase), "-2", sep="")

allpop<- base
N<- 10
ngen<- 75

pop<- base[sample(base$id, size=N),]
pop

#loop for each generation
for(j in 1:ngen){ 
  
  ##create the empty data table to record information on new individuals created 
  popnew<- data.frame(id=c(nrow(allpop)+1):c(nrow(allpop)+N), p1=NA, p2=NA, allele1=NA, allele2=NA, gen=j)
  
  ##loop to create each new individual through random mating
  for(i in 1:nrow(popnew)){ 
    
    ##randomly select two parents and one allele per parent
    ixp1<- sample(1:N,1) #row index for parent 1
    ixp2<- sample(1:N,1) #row index for parent 2
    rowp1<- pop[ixp1,] #data on parent 1
    rowp2<- pop[ixp2,] #data in parent 2
    alsp1<- rowp1[,c('allele1', 'allele2')] #alleles p1
    alsp2<- rowp2[,c('allele1', 'allele2')] #alleles p2
    al1samp<- sample(alsp1,1) #randomly select one allele from parent 1
    al2samp<- sample(alsp2,1) #randomly select one allele from parent 2
    
    ##record information
    popnew[i,c('p1', 'p2')]<-  pop[c(ixp1, ixp2),'id']  #parent information 
    popnew[i,c('allele1','allele2')]<- c(al1samp,al2samp) #allele information
  }#new individuals 'popnew' created
  
  ##add new individuals to the information data.frame
  allpop<- rbind(allpop, popnew)
  
  ##replace the population
  pop<- popnew
}#ngen generations completed

allpop[allpop$gen == 2,]
allpop$allele1 %>% table

ped<- pedigree(allpop$p1, allpop$p2, allpop$id)
I<- inbreeding(ped)
allpop<- data.frame(allpop, I)

mnI<- allpop %>% group_by(gen) %>% summarise(I = mean(I))
head(mnI)

View(allpop)


ped4<- pedigree(allpop$p1, allpop$p2,
                allpop$id)

A<- getA(ped4)
A<- as.matrix(A) #converts the A matrix to a regular matrix
A

A[3:4,]

plot(mnI$gen, mnI$I)



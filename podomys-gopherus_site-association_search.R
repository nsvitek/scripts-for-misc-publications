# Analyses that form the basis of manuscript written with J. Bloch and J. Austin 

# Dependencies -----
library(dplyr)
library(viridis)
library(gridBase)
library(maps)
library(ggplot2)
library(ggrepel)
library(ggmap)

datadir<-"path/to/data/raw_data"

setwd(datadir)

# Pick Colors for bar plot -----
colorindex<-c(viridis(6))
colorindextol<-c("#332288", "#88CCEE", "#117733", "#DDCC77", "#CC6677","#AA4499")
colorindex3<-c("#4477AA", "#DDCC5F", "#CC6677") #also Tol colors

# podomys section -----
# load in query of FLMNH vp database for all occurences of Podomys
mouse<-read.csv("vp-query-result_floridanus_annotated.csv")
# str(mouse)

# read in series of queries to pull all records from Podomys sites (or at least, in case of Reddick 1A, all relevant burrowing animals)
Reddick1A<-10449 # number of specimens from site, csv wouldn't download
# Leisey1A<-23660 #same problem as Reddick1A, but Leisey mouse is "Podomys sp."
sites_of_interest<-read.csv("fauna_podomys_no-gopherus.csv")
sites1<-read.csv("vp-query-podomys_sites_1.csv")
sites2<-read.csv("vp-query-podomys_sites_2.csv")
sites3<-read.csv("vp-query-podomys-lit-sites.csv") #Coleman 2A, Haile 21A, and Arredondo 1A in literature but not catalogue
sites4<-read.csv("vp-query-reddick-burrowers.csv") #reddick had too many to download all
sites5<-read.csv("vp-query-monkey_jungle_hammock.csv") #MJH 2 and MJH considered separate sites, but misleading. Add in MJH
# sites<-rbind(sites_of_interest,sites1,sites2,sites3,sites4,sites5) #all specimens from all sites with podomys
sites<-rbind(sites_of_interest,sites1,sites2,sites3,sites4,sites5) #all specimens from all sites with podomys

# Make list of sites with Podomys
mousesites<-unique(mouse$Site) %>% .[-which(.==" Recent")] %>% as.character %>% 
  append(.,c(" Coleman 2A", " Haile 21A", " Arredondo 1A" 
             #, " Leisey Shell Pit 1A"
             )) %>%  
  .[order(.)] #Podomys sites
# vector of sites with Podomys. "Recent" doesn't count. Remove. Additions from literature not in catalog.
# mousesites_alt<-distinct(sites,Site) %>% select(.,County,Epoch,Site,Land.Mammal.Age) %>%
#   arrange(.,toupper(Site)) #for species account...doesn't work any more. 
mousesites_alt<-select(sites,Site,Land.Mammal.Age) %>% distinct(.,Site,Land.Mammal.Age) %>% 
  filter(.,Land.Mammal.Age!=" Clarendonian",Land.Mammal.Age!=" Hemphillian",Land.Mammal.Age!=" Barstovian") %>%
  filter(.,(Site!=" Vero"|Land.Mammal.Age!=" ")) %>%
  filter(.,(Site!=" Monkey Jungle Hammock 2"|Land.Mammal.Age!=" ")) #not elegant, but it works. 
mousesites[which(mousesites==" Monkey Jungle Hammock 2")]<-" Monkey Jungle Hammock" #standardize MJH for easier comparison
sites$Site[which(sites$Site==" Monkey Jungle Hammock 2")]<-" Monkey Jungle Hammock" #standardize MJH for easier comparison

# cbind(as.character(mousesites_alt$Site),as.character(mousesites_alt$Land.Mammal.Age))
# Check: The two provide the same information
mousesites[-which(mousesites %in% mousesites_alt$Site)]

# Make lists of sites with candidate commensal species, as well as a list of sites with turtles
# Sites with no turtles may have a general collection bias precluding Gopherus discovery 
tortoisesites<-sites %>% filter(.,Genus==" Gopherus") %>% select(.,Site) %>%unique %>% 
  unlist %>% as.character %>% .[order(.)] #vector of sites that have Gopherus and Podomys
# tortoisesites<-sites %>% filter(.,Species==" polyphemus") %>% select(.,Site) %>%unique %>% 
#   unlist %>% as.character %>% .[order(.)] #vector of sites that have Gopherus and Podomys
geomyssites<-sites %>% filter(.,Genus==" Geomys") %>% select(.,Site) %>%unique %>% 
  unlist %>% as.character %>% .[order(.)] #vector of sites that have Geomys
sigmodonsites<-sites %>% filter(.,Genus==" Sigmodon") %>% select(.,Site) %>%unique %>% 
  unlist %>% as.character %>% .[order(.)] #vector of sites that have Sigmodon
# sigmodonsites<-sites %>% filter(.,Species==" hispidus") %>% select(.,Site) %>%unique %>% 
#   unlist %>% as.character %>% .[order(.)] #vector of sites that have Sigmodon
dasypussites<-sites %>% filter(.,Genus==" Dasypus") %>% select(.,Site) %>%unique %>% 
  unlist %>% as.character %>% .[order(.)] #vector of sites that have Dasypus
polionotussites<-sites %>% filter(.,Species==" polionotus") %>% select(.,Site) %>%unique %>% 
  unlist %>% as.character %>% .[order(.)] #vector of sites that have Peromyscus polionotus
turtlesite<-sites %>% filter(.,Order==" Testudines") %>% select(.,Site) %>%unique %>% 
  unlist %>% as.character %>% .[order(.)] #vector of sites that have any turtle fossils to control for collecting bias

# create a presence/absence table
geomys<-sigmodon<-dasypus<-polionotus<-any_turtle<-gopherus<-podomys<-rep("present",length(mousesites))
geomys[-which(mousesites %in% geomyssites)]<-"absent"
sigmodon[-which(mousesites %in% sigmodonsites)]<-"absent"
dasypus[-which(mousesites %in% dasypussites)]<-"absent"
polionotus[-which(mousesites %in% polionotussites)]<-"absent"
any_turtle[-which(mousesites %in% turtlesite)]<-"absent"
gopherus[-which(mousesites %in% tortoisesites)]<-"absent"

# create an n for each site
counts<-sites %>% group_by(Site) %>% summarise(n=length(Catalog.Number)) %>% 
  as.matrix %>% .[-which(.[,1]==" Reddick 1A"|.[,1]==" Leisey Shell Pit 1A"),] %>% 
  rbind(.,c(" Reddick 1A",Reddick1A)) %>% #rbind(.,c(" Leisey Shell Pit 1A",Leisey1A))  %>% 
  as.data.frame %>% arrange(.,Site) %>% 
  cbind(.,podomys,gopherus,dasypus,geomys,polionotus,sigmodon,any_turtle,mousesites_alt$Land.Mammal.Age)
counts$n<-as.numeric(as.character(counts$n))

# remove ineligible sites
fauna<-counts[,] %>% .[-which(.[,1]==" Fort Meade Mine. #7 Dragline (Gardinier)"),]
# fauna<-counts[-which(counts$any_turtle=="absent"),] %>% .[-which(.[,1]==" Fort Meade Mine. #7 Dragline (Gardinier)"),]

# exclude poorly sampled sites (n < 50)
newfauna<-fauna[which(fauna$n>=50),]
candidates<-nrow(newfauna)
yesgopherus<-length(which(newfauna$gopherus=="present"))
yespolionotus<-length(which(newfauna$polionotus=="present"))
yessigmodon<-length(which(newfauna$sigmodon=="present"))
yesdasypus<-length(which(newfauna$dasypus=="present"))
yesgeomys<-length(which(newfauna$geomys=="present"))

# Test to see if Podomys and [taxon] are found together more often than expected by chance
# presence absence data are binomial, and probability(Podomys=present)==1, so
# can simplify to test where null is probability(Gopherus=present)==0.5
pfishg<-binom.test(x=yesgopherus,n=candidates,p=0.5,alternative="greater")
pfishs<-binom.test(x=yessigmodon,n=candidates,p=0.5,alternative="greater")
pfishd<-binom.test(x=yesdasypus,n=candidates,p=0.5,alternative="greater")
pfishgeo<-binom.test(x=yesgeomys,n=candidates,p=0.5,alternative="greater")
pfishpol<-binom.test(x=yespolionotus,n=candidates,p=0.5,alternative="greater")

# Make bar graph elements
newfauna<-droplevels(newfauna)
tmatrix<-NULL
for (age in 1:length(levels(newfauna$`mousesites_alt$Land.Mammal.Age`))){
  count1<-filter(newfauna,`mousesites_alt$Land.Mammal.Age`==levels(`mousesites_alt$Land.Mammal.Age`)[age]) %>%
    summarise(.,Podomys=length(which(podomys=="present")),
              `Peromyscus polionotus`=length(which(polionotus=="present")),
              Sigmodon=length(which(sigmodon=="present")),
              Geomys=length(which(geomys=="present")),
              Dasypus=length(which(dasypus=="present")),
              Gopherus=length(which(gopherus=="present")))
  tmatrix<-rbind(tmatrix,count1)
}

plotcount<-t(tmatrix) %>% apply(.,1,sum) %>% cbind(t(tmatrix),.)
colnames(plotcount)<-c("Irvingtonian","Rancholabrean","Total")

# #### Generate Report
# print("Number of Sites with Podomys:")
# length(mousesites)
# print("Number of Sites with Podomys + any turtle + enough specimens:")
# nrow(newfauna)
# print("Number of candidate sites without Gopherus:")
# length(which(newfauna$gopherus=="absent"))
# print("Percent of Podomys + turtle sites without Gopherus:")
# 1-length(which(newfauna$gopherus=="absent"))/nrow(newfauna)
# 
# write.csv(newfauna,"podomys_associates_report.csv")
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

# site section ----------- 
# Hypothesis: The two Podomys sites missing Gopherus are somehow environmentally weird or not conducive to Gopherus
# Create taxon lists for each of the two sites
taxa_Haile11B<-filter(sites, Site == " Haile 11B") %>%
  distinct(., Genus,Species) %>% #select(., Genus,Species,Site,Land.Mammal.Age) %>% #this part doesn't work
  arrange(.,Genus)
taxa_WMS<-filter(sites, Site == " Warm Mineral Springs") %>%
  distinct(., Genus,Species) %>% #select(., Genus,Species,Site,Land.Mammal.Age) %>% #this part doesn't work
  arrange(.,Genus)

# Save the data:
sink(file="gopherus-absent_site_taxon_list.txt",append=FALSE)
print("Taxon list for Haile 11B")
print(taxa_Haile11B)
print("Taxon list for Warm Mineral Springs")
print(taxa_WMS)
sink()

# gopherus section ------------ 
# Get a list of all sites with Gopherus, regardless of presence/absence of Podomys
# Additional from Franz & Quitmyer: Inglis 1F
sitelist_gopherus<-read.csv("vp-query-all_gopherus.csv") %>% #filter(.,State==" Florida") %>% #not necessary
  filter(.,Site!=" ",Site!=" Recent") %>% 
  select(.,Site) %>% unique %>% unlist %>% as.character %>% c(.," Inglis 1F") %>% .[order(.)] 

sitelist_gopherus_lma<-read.csv("vp-query-all_gopherus_latlon.csv") %>%
  filter(.,Site!=" ",Site!=" Recent") %>% 
  distinct(.,Site,Land.Mammal.Age,Latitude,Longitude) #%>% arrange(.,Site)
# sitelist_gopherus[-which(sitelist_gopherus %in% mousesites)] #check
# sitelist_gopherus[which(sitelist_gopherus %in% mousesites)] #double check

# will need to search AMNH for seminole field by hand :(
# To get the site info for those with both mice and turtles, filter out the sites object, then add the tables below (too many specimens to download all as one table)
# 1: Achan through Arredondo 1D
# 2: Auciulla River 3J mammalia (total too big)
# 3: Bald Knob to Forsberg
# 4: Fort Green to Haile 13E
# 5: Haile 13F to Horse Hill Low
# 6: Ichnetucknee
# 7: Inglis 1A mammalia (total too big)
# 8: Inglis 1C
# 9: Inglis 1D to Kingsford Mine
# 10: La Belle to Millenium Park [note Leisey 1A in separate file]
# 11: Moss Acres to Phosphoria Mine [note Palmetto Mine in separate file]
# 12: Palmetto Mine mammalia (total too big)
# 13: Ponte Vedra Beach to Reddick
# 14: Rock Springs to Santa Fe River 2
# 15: Santa Fe River 8 to Tri-Britton [note Thomas farm in separate file]
# 16: Thomas Farm Rodentia only, no records of Xenarthra
# 17: Turkey Foot High to Withlacoochie River 1A
# 18: The other two Withlacoochie River sites

Aucilla3J<-9018
Inglis1A<-9560
Leisey1A<-23660
PalmettoMine<-14857
ThomasFarm<-59138
sitesLeisey<-read.csv("vp-query-leisey1a.csv") #leisey had too many to download all
# get the site taxon data already read in for podomys search
sitesshared<-sites[which(sites$Site %in% sitelist_gopherus[which(sitelist_gopherus %in% mousesites)]),] %>% droplevels
# read in all taxa info for the rest of the gopherus sites
sitesnumericnames<-list.files(getwd(),pattern = "^vp-query-result_gopherus_sites[0-9]") 
sitesnumeric<-NULL
for (i in 1:length(sitesnumericnames)){
  temp1<-read.csv(sitesnumericnames[i])
  sitesnumeric<-rbind(sitesnumeric,temp1)
}
sitesnumeric<-select(sitesnumeric,-X) #get rid of useless unshared column

#consolidate into a single table
sites_gopherus<-rbind(sitesLeisey,sitesshared,sitesnumeric)
sites_gopherus$Site[which(sites_gopherus$Site==" Monkey Jungle Hammock 2")]<-" Monkey Jungle Hammock" #standardize MJH for easier comparison
sites_gopherus$Site[which(sites_gopherus$Site==" Leisey Shell Pit 1")]<-" Leisey Shell Pit 1A" #standardize Leisey for easier comparison

# Thought: Gopherus polyphemus considered only to go back to Late Pliocene by Franz and Quitmyer
# Further Thought: Comparison only makes sense if spatiotemporal extent kept constant. Keep to Irvingtonian-Rancholabrean
sites_gopherus<-filter(sites_gopherus,Land.Mammal.Age==" Irvingtonian"|Land.Mammal.Age==" Rancholabrean")
sites_gopherus<-droplevels(sites_gopherus)
# dim(sites_gopherus)

#Now run cleaning protocol similar to podomys. First, look for all relevant critters. Use "Rodentia" as broader comparison instead of Testudines
gsites_podomys<-sites_gopherus %>% filter(.,Genus==" Podomys") %>% select(.,Site) %>%unique %>% 
  unlist %>% as.character %>% .[order(.)] #vector of sites that have Gopherus and Podomys
gsites_geomys<-sites_gopherus %>% filter(.,Genus==" Geomys") %>% select(.,Site) %>%unique %>% 
  unlist %>% as.character %>% .[order(.)] #vector of sites that have Geomys
gsites_sigmodon<-sites_gopherus %>% filter(.,Genus==" Sigmodon") %>% select(.,Site) %>%unique %>% 
  unlist %>% as.character %>% .[order(.)] #vector of sites that have Sigmodon
gsites_dasypus<-sites_gopherus %>% filter(.,Genus==" Dasypus") %>% select(.,Site) %>%unique %>% 
  unlist %>% as.character %>% .[order(.)] #vector of sites that have Dasypus
gsites_polionotus<-sites_gopherus %>% filter(.,Species==" polionotus") %>% select(.,Site) %>%unique %>% 
  unlist %>% as.character %>% .[order(.)] #vector of sites that have Peromyscus polionotus
gsites_rodents<-sites_gopherus %>% filter(.,Order==" Rodentia") %>% select(.,Site) %>%unique %>% 
  unlist %>% as.character %>% .[order(.)] #vector of sites that have any turtle fossils to control for collecting bias

# get an n for each site 
gcounts<-sites_gopherus %>% group_by(Site) %>% summarise(n=length(Catalog.Number)) %>% 
  as.matrix %>% .[-which(.[,1]==" Reddick 1A"|
                           .[,1]==" Leisey Shell Pit 1A"|
                           .[,1]==" Aucilla River 3J (Page/ladson)"#|
                          # .[,1]==" Inglis 1A"|
                          # .[,1]==" Palmetto Mine (Agrico)"|
                          # .[,1]==" Thomas Farm"
                         ),] %>% 
  rbind(.,c(" Reddick 1A",Reddick1A)) %>% rbind(.,c(" Leisey Shell Pit 1A",Leisey1A))  %>%
  rbind(.,c(" Aucilla River 3J (Page/ladson)",Aucilla3J)) %>% #rbind(.,c(" Inglis 1A",Inglis1A))  %>%
  #rbind(.,c(" Palmetto Mine (Agrico)",PalmettoMine)) %>% #rbind(.,c(" Thomas Farm",ThomasFarm))  %>%
  as.data.frame %>% arrange(.,Site) 
gcounts$n<-as.numeric(as.character(gcounts$n))

# create a presence/absence table
ggeomys<-gsigmodon<-gdasypus<-gpolionotus<-any_rodent<-gpodomys<-ggopherus<-rep("present",length(gcounts$Site))
ggeomys[-which(gcounts$Site %in% gsites_geomys)]<-"absent"
gsigmodon[-which(gcounts$Site %in% gsites_sigmodon)]<-"absent"
gdasypus[-which(gcounts$Site %in% gsites_dasypus)]<-"absent"
gpolionotus[-which(gcounts$Site %in% gsites_polionotus)]<-"absent"
any_rodent[-which(gcounts$Site %in% gsites_rodents)]<-"absent"
gpodomys[-which(gcounts$Site %in% gsites_podomys)]<-"absent"

gcounts<-cbind(gcounts,ggopherus,gpodomys,gdasypus,ggeomys,gpolionotus,gsigmodon,any_rodent) 

# Make some manual fixes
litsites<-c(" Coleman 2A", " Haile 21A", " Arredondo 1A"," Monkey Jungle Hammock") #These three sites have podomys in literature, not in catalogue
gcounts$gpodomys[which(gcounts$Site %in% litsites)]<-"present"
# update seminole field info: it has geomys, dasypus, sigmodon in AMNH database.
gcounts[which(gcounts$Site==" Seminole Field"),
        which(colnames(gcounts)=="ggeomys"|
                colnames(gcounts)=="gdasypus"|
                colnames(gcounts)=="gsigmodon")]<-"present"

# info for ALL sites of same spatiotemporal extent as Podomys
gcounts_sameextent<-gcounts

# exclude poorly sampled sites (n < 50)
gcounts_50plus<-gcounts[which(gcounts$n>=50),]

#remove sites without rodents
gcounts_rodents<-filter(gcounts_50plus,any_rodent=="present") 

gfauna<-gcounts_rodents
nrow(gcounts_sameextent)-nrow(gcounts_50plus)
nrow(gcounts_50plus)
nrow(gcounts_rodents)

# Add in land mammal age data
head(sitelist_gopherus_lma)
lma<-lat<-lon<-rep(NA, nrow(gfauna))
for (i in 1:nrow(gfauna)){lma[i]<-as.character(sitelist_gopherus_lma$Land.Mammal.Age[which(sitelist_gopherus_lma$Site==as.character(gfauna$Site[i]))])}
for (i in 1:nrow(gfauna)){lat[i]<-as.character(sitelist_gopherus_lma$Latitude[which(sitelist_gopherus_lma$Site==as.character(gfauna$Site[i]))])}
for (i in 1:nrow(gfauna)){lon[i]<-as.character(sitelist_gopherus_lma$Longitude[which(sitelist_gopherus_lma$Site==as.character(gfauna$Site[i]))])}

gnewfauna<-cbind(gfauna,lma,lat,lon)
# Some manual fiddling: remove sites for which "Rodent" record is very poor 
# # (protocol: manually examine rodent records at sites missing sigmodon)
# # Aucilla River 3J: has Neotoma, so small stuff is there, but mostly Ondatra and Castor, with a little bit of Erithizon, Sciurids, and Rattus.
# # Haile 1A: 1 unidentified rodent fossil. Remove.
# # Bone Cave: 1 unidentified rodent limb bone. Remove.
# # Haile 21A: Has podomys and geomys. Also has sparse cataloged Rodent Record, found Sigmodon physically in collection
# # Hornsby Springs: 1 record of Castor
# # Monkey Jungle Hammock: should have Sigmodon based on Ober[& collection: MJH2], otherwise seems correct
# # Santa Fe River 1: 21 rodents, mostly Castoroides, Neochoerus, Hydrochoeurs, Ondatra, Castor, Neofiber.
# # Santa Fe River 2: 111 rodents: Castor, Castoroides [a lot], Neochoerus, Hydrochoerus, Ondatra, Geomys, Neofiber
# # Suwannee River: 1 record of Castor
# # Tri-Britton Site: 7 records, all of large rodents, Hydrocheridae and Castoroides
# # Wilson Quarry: 15 records: Neofiber and one murid.

gnewfauna$gsigmodon[which(gnewfauna$Site==" Haile 21A"|
                            gnewfauna$Site==" Monkey Jungle Hammock")]<-"present"
gnewfauna$ggeomys[which(gnewfauna$Site==" Haile 21A")]<-"present"


gnewfauna[which(gnewfauna$gpodomys=="absent"),]
gnewfauna<-gnewfauna[-which(#gnewfauna$Site==" Aucilla River 3J (Page/ladson)"|
                       gnewfauna$Site==" Haile 1A"|
                       gnewfauna$Site==" Bone Cave"|
                       gnewfauna$Site==" Hornsby Springs"|
                       gnewfauna$Site==" Santa Fe River 1"|
                       gnewfauna$Site==" Suwannee River"|
                       gnewfauna$Site==" Tri-Britton Site"|
                       gnewfauna$Site==" Wilson Quarry"
                       ),]

# Another point to keep in mind is that Peromyscus sp. is found at:
# # Reddick 1C, Haile 14A, Haile 8A, &  Haile 13A.
# # (not necessary search at Podomys sites. Only Peromyscus sp. records at Coleman 2A and Haile 21A...
# # Haile 21A searched physically, no other Peromyscus other than "Podomys n. sp." found, and Coleman 2A
# # discussed in published faunal list discussion: probably not polionotus)
gperomyscus_alt<-cbind(as.character(gnewfauna$Site),
                       as.character(gnewfauna$gpodomys),
                       as.character(gnewfauna$gpolionotus))
gperomyscus_alt[which(gperomyscus_alt[,1]==" Reddick 1C"|
                        gperomyscus_alt[,1]==" Haile 8A"|
                        gperomyscus_alt[,1]==" Haile 13A"|
                        gperomyscus_alt[,1]==" Haile 14A"),c(2:3)]<-"present"


# Get counts for statistical tests
gcandidates<-nrow(gnewfauna)
gyespodomys<-length(which(gnewfauna$gpodomys=="present"))
gyespolionotus<-length(which(gnewfauna$gpolionotus=="present"))
gyessigmodon<-length(which(gnewfauna$gsigmodon=="present"))
gyesdasypus<-length(which(gnewfauna$gdasypus=="present"))
gyesgeomys<-length(which(gnewfauna$ggeomys=="present"))
gyespolionotus_alt<-length(which(gperomyscus_alt[,3]=="present"))
gyespodomys_alt<-length(which(gperomyscus_alt[,2]=="present"))

# Test to see if Podomys and [taxon] are found together more often than expected by chance
# presence absence data are binomial, and probability(Podomys=present)==1, so
# can simplify to test where null is probability(Gopherus=present)==0.5
gpfishpg<-binom.test(x=gyespodomys,n=gcandidates,p=0.5,alternative="greater")
gpfishs<-binom.test(x=gyessigmodon,n=gcandidates,p=0.5,alternative="greater")
gpfishd<-binom.test(x=gyesdasypus,n=gcandidates,p=0.5,alternative="greater")
gpfishpol<-binom.test(x=gyespolionotus,n=gcandidates,p=0.5,alternative="greater")
gpfishgeo<-binom.test(x=gyesgeomys,n=gcandidates,p=0.5,alternative="greater")
gpfishpg_alt<-binom.test(x=gyespodomys_alt,n=gcandidates,p=0.5,alternative="greater")
gpfishpol_alt<-binom.test(x=gyespolionotus_alt,n=gcandidates,p=0.5,alternative="greater")

# Create elements for plots
gnewfauna<-gnewfauna[-which(gnewfauna$lma==" Santarosean"),] %>% droplevels(.)
gtmatrix<-NULL
for (age in 1:length(levels(gnewfauna$lma))){
  count1<-filter(gnewfauna,lma==levels(gnewfauna$lma)[age]) %>%
    summarise(.,Podomys=length(which(gpodomys=="present")),
              `Peromyscus polionotus`=length(which(gpolionotus=="present")),
              Sigmodon=length(which(gsigmodon=="present")),
              Geomys=length(which(ggeomys=="present")),
              Dasypus=length(which(gdasypus=="present")),
             Gopherus=length(which(ggopherus=="present")))
  gtmatrix<-rbind(gtmatrix,count1)
}

gplotcount<-t(gtmatrix) %>% apply(.,1,sum) %>% cbind(t(gtmatrix),.)
colnames(gplotcount)<-c("Irvingtonian","Rancholabrean","Total")
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
turtlemouse_stats<-print(c(pfishpol,
                           pfishs,
                           pfishgeo,
                           pfishd,
                           pfishg,
                           gpfishpg,
                           gpfishpol,
                           gpfishs,
                           gpfishgeo,
                           gpfishd,
                           gpfishpg_alt,
                           gpfishpol_alt)) %>% unlist %>% matrix(.,ncol=10,byrow=TRUE)
colnames(turtlemouse_stats) <- names(pfishg[c(1:4,4,5:9)])
bonferroni<-as.numeric(turtlemouse_stats[,3])*5
turtlemouse_stats<-cbind(turtlemouse_stats,bonferroni)
write.csv(turtlemouse_stats,file="podomys-gopherus_p-vals.csv")
# make plot --------------------------------------------------
pdf("podomys_barplot.pdf",width=3.14,height=6.28,pointsize=8)
par(mfrow=c(2,1))
par(mai=c(0.4,0.6,0.4,0))
barplot(plotcount, 
        xlab=NULL, col=colorindextol,ylim=c(0,max(plotcount)), #visualization can be better
        beside=TRUE,ylab="Number of Sites",cex.names=0.9,cex.lab=1.5,
        density = c(10,20,40,50,80,100)) #,xaxt='n')
# text(cex=1, x=c(6,14,21), y=-.5, colnames(plotcount), xpd=TRUE, srt=20, pos=2)
# savefont <- par(font=3)
# legend('topleft',legend=rownames(plotcount),pch=22,pt.bg=colorindex,pt.cex=1.5,bty="n",cex=0.9)
# par(savefont)
text(-3.5,15.5,labels="A",xpd=TRUE,cex=2)
segments(x0=seq(16.5,20.5,1),y0=as.numeric(turtlemouse_stats[c(1:5),4])*max(plotcount),
         x1=seq(16.5,20.5,1),y1=rep(max(plotcount),5))
segments(x0=15,x1=21,y0=max(plotcount)/2,y1=max(plotcount)/2,lwd=3,lty=6,col="black")
# pdf("gopherus_barplot.pdf",width=3.14,height=3.14,pointsize=8)
par(mai=c(0.6,0.6,0.2,0))
barplot(gplotcount, 
        xlab="Land Mammal Age", col=colorindextol,ylim=c(0,max(gplotcount)),
        beside=TRUE,ylab="Number of Sites", cex.names=0.9,cex.lab=1.5, density = c(10,20,40,50,80,100)) #,xaxt='n')
# text(cex=1, x=c(7,14,21,28), y=-1.25, colnames(gplotcount), xpd=TRUE, srt=20, pos=2)
savefont <- par(font=3)
legend('topleft',legend=rownames(plotcount),fill=colorindextol,pt.cex=1.5,
       bty="n",cex=1,density = c(10,20,40,50,80,100))
par(savefont)
text(-3.5,48,labels="B",xpd=TRUE,cex=2)
segments(x0=seq(15.5,19.5,1),y0=as.numeric(turtlemouse_stats[c(6:10),4])*max(gplotcount),
         x1=seq(15.5,19.5,1),y1=rep(max(gplotcount),5))
segments(x0=15,x1=21,y0=max(gplotcount)/2,y1=max(gplotcount)/2,lwd=3,lty=6,col="black")
dev.off()

# make map mouse -----------------------------------------------
# make labelled map of sites where fossils of Podomys were found (unpublished)
mousepoints<-read.csv("vp-query-result_podomys_site_latlon.csv",header=TRUE) %>%
  distinct(Site, Latitude, Longitude)
ylims<-range(mousepoints$Latitude)+c(-0.8,+0.8)
xlims<-range(mousepoints$Longitude)+c(-0.8,+0.8)

#old plot version
# map('state',region='Florida',lwd=2,fill=FALSE,add=TRUE)
# with(mousepoints, points(Latitude~Longitude,bg="red",pch=21,cex=2))
# with(mousepoints, text(Longitude,Latitude,Site,pos=4,adj=1))

fl_cn<-map_data('county',region='Florida')
fl_outline<-map_data('state',region='Florida')
# # Map of Podomys-only sites, with labels
# map1<-ggplot(data=mousepoints,aes(x=Longitude,y=Latitude))+
#   geom_polygon(data=fl_cn,aes(x=long,y=lat,group=group),fill=NA,color="gray")+
#   geom_polygon(data=fl_outline,aes(x=long,y=lat,group=group),fill=NA,color="black")+
#   coord_fixed(1.1)+
#   geom_point(data=mousepoints,aes(x=Longitude,y=Latitude),size=2)+
#   geom_label_repel(aes(x=Longitude,y=Latitude,label=Site),size=3,
#                   box.padding=unit(0.45,"lines"),label.padding=unit(0.05,"lines")) +
#   # coord_cartesian(xlim=xlims)+
#   theme_classic() +
#   theme(axis.line=element_blank(),axis.text.x=element_blank(),
#         axis.text.y=element_blank(),axis.ticks=element_blank(),
#         axis.title.x=element_blank(),axis.title.y=element_blank(),
#         text=element_text(family="serif",size=9))
# ggsave("Fig_podomys_map.eps",plot=map1,width=3.14,height=6.28,units="in")

# make map all --------
# Make map of all analyzed sites containing Podomys, Gopherus, or both (published)
mapdata<-gnewfauna[,c(3,4,11,12)] %>% as.data.frame
mapdata$lat<-as.numeric(levels(mapdata$lat))[mapdata$lat]
mapdata$lon<-as.numeric(levels(mapdata$lon))[mapdata$lon]
addtomapdata<-mousepoints[c(4,10),c(2,3)]
colnames(addtomapdata)<-c("lat","lon")
addtomapdata$ggopherus<-"absent"
addtomapdata$gpodomys<-"present"
mapdata<-rbind(mapdata,addtomapdata[,c(3,4,1,2)])
mapdata$combo<-paste("gopherus",mapdata$ggopherus,"podomys",mapdata$gpodomys,sep="_")

fl_cn<-map_data('county',region='Florida')
fl_outline<-map_data('state',region='Florida')
ylims<-range(mapdata$lat)+c(-0.8,+0.8)
xlims<-range(mapdata$lat)+c(-0.8,+0.8)


map2<-ggplot(data=mapdata,aes(x=lon,y=lat))+
  geom_polygon(data=fl_cn,aes(x=long,y=lat,group=group),fill=NA,color="gray")+
  geom_polygon(data=fl_outline,aes(x=long,y=lat,group=group),fill=NA,color="black")+
  coord_fixed(1.1)+
  geom_point(aes(shape=combo,fill=combo),size=2)+
  theme_classic() +
  scale_shape_manual(values=c(21,22,23),name="Species Occurrences",
                    labels=c("Gopherus - ; Podomys +","Gopherus + ; Podomys -","Gopherus + ; Podomys +"))  +
  scale_fill_manual(values= colorindex3,name="Species Occurrences",
        labels=c("Gopherus - ; Podomys +","Gopherus + ; Podomys -","Gopherus + ; Podomys +")) +
  theme(axis.line=element_blank(),axis.text.x=element_blank(),
        axis.text.y=element_blank(),axis.ticks=element_blank(),
        axis.title.x=element_blank(),axis.title.y=element_blank(),
        legend.position=c(0.35,0.2),text=element_text(size=10))

ggsave("Fig_podomys_map.pdf",plot=map2,width=3.14,height=4.28,units="in")
  
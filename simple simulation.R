
# ANCHOVY GROWTH THROUGH TIME INFO#######################################
# Anchovy hatch at 3mm, and reach 5mm in 7 days, and 10mm in 17 days. [10 days in 5-10mm interval]
# Anchovy reach 20mm 47 days [30 days in 10-20mm]
# Based on:
# Fukuhara, O., 1983. Development and growth of laboratory reared Engraulis 
#j aponica (Houttuyn) larvae. Journal of Fish Biology, 23(6), pp.641-652.


# ANCHOVY SURVIVAL/FEEDING BY SIZE & CRITICAL WINDOW #####################
# prior early labwork on feeding/swimming behavior of larval anchovy
# indicates that feeding success and volume of water searched for prey 
# per unit time grows exponentially across days/growth since hatch

# "laboratory experiments indicated that mortality due to starvation was
# higher just after yolk absorption than at any other point in life"

# "The density of rotifers and dinoflagellates required for larvae to meet metabolic needs
# was calculated from caloric and respiration data, and estimates of volume searched
# and feeding success. These calculations indicate that anchovy larvae just after yolk
# absorption require up to 37 times the food density of older larave"

# Based on:
# Hunter, J.R., 1972. Swimming and feeding behavior of larval anchovy 
# Engraulis mordax. Fish. Bull, 70(82), pp.l-834.

# SCENARIO ONE - CONSTANT MORTALITY ACROSS ALL AGES/SIZES#################
# mortality -.1 for all fish (10% die each day), start with 1000 larvae
m<- -.1
smalldayz<-10
largedayz<-30

small<-1000*exp(m*smalldayz)
large<-small*exp(m*largedayz)

big.small.ratio<-large/small
big.small.ratio

# SCENARIO TWO - EARLY MORTALITY HIGH, LATER MORTALITY LOW(ISH)##########
# mortality -.1 for small fish (10% die each day), start with 1000 larvae
# mortality -.05 for large fish (5% die each day)
m.s<- -.1
m.l<- -.05
smalldayz<-10
largedayz<-30

small<-1000*exp(m.s*smalldayz)
large<-small*exp(m.l*largedayz)

big.small.ratio2<-large/small
big.small.ratio2

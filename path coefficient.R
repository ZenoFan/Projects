library(lavaan)
library(semPlot)
library(Hmisc)
mydata=read.csv(file.choose())
pathdata=mydata[,-1]
colnames(pathdata)
library(car)
vif(lm(YIELD_ca_12_kgh=colnames(pathdata), data=pathdata))
aftervifdata=pathdata[,colnames(pathdata)!="Dry_wt_of_Fietile_tiller_Aplus7d_m2"]
vif(lm(YIELD_ca_12_kgh=colnames(pathdata), data=aftervifdata))

colnames(aftervifdata)
afterscale=aftervifdata
afterscale[,20]=scale(aftervifdata[,20],center=FALSE,scale=c(100))
afterscale[,4]=scale(aftervifdata[,4],center=FALSE,scale=c(100))
afterscale[,18]=scale(aftervifdata[,18],center=FALSE,scale=c(10))

"path analysis"

model_1=
'YIELD_ca_12_kgh~a1*Spike_partitioing_index+
a2*Lamina_partitioining_index+
a3*TrueStem_partitioing_index+
a4*Internode1_partitioing_index+
a5*Internode2_partitioing_index+
a6*Internode3_partitioing_index+
a7*X2ndplus3rd_partitioining_index+
a8*Internode1_Length+
a9*Internode2_Length+
a10*Internode3_Length+
a12*AGDM_maturity_m+
a13*AGDM_Aplus7d_m2+
a14*FE_Grains_g_CW_m.+
a15*Grain_partitioning_index_maturity+
a16*Straw_partitioning_index_maturity+
a17*Stem_partitioing_index_maturity

Spike_partitioing_index~~b1*Lamina_partitioining_index
Spike_partitioing_index~~b2*TrueStem_partitioing_index
Spike_partitioing_index~~b3*Internode1_partitioing_index
Spike_partitioing_index~~b4*Internode2_partitioing_index
Spike_partitioing_index~~b5*Internode3_partitioing_index
Spike_partitioing_index~~b6*X2ndplus3rd_partitioining_index
Spike_partitioing_index~~b7*Internode1_Length
Spike_partitioing_index~~b8*Internode2_Length
Spike_partitioing_index~~b9*Internode3_Length
Spike_partitioing_index~~b11*AGDM_maturity_m
Spike_partitioing_index~~b12*AGDM_Aplus7d_m2
Spike_partitioing_index~~b13*FE_Grains_g_CW_m.
Spike_partitioing_index~~b14*Grain_partitioning_index_maturity
Spike_partitioing_index~~b15*Straw_partitioning_index_maturity
Spike_partitioing_index~~b16*Stem_partitioing_index_maturity

Lamina_partitioining_index~~c2*TrueStem_partitioing_index
Lamina_partitioining_index~~c3*Internode1_partitioing_index
Lamina_partitioining_index~~c4*Internode2_partitioing_index
Lamina_partitioining_index~~c5*Internode3_partitioing_index
Lamina_partitioining_index~~c6*X2ndplus3rd_partitioining_index
Lamina_partitioining_index~~c7*Internode1_Length
Lamina_partitioining_index~~c8*Internode2_Length
Lamina_partitioining_index~~c9*Internode3_Length
Lamina_partitioining_index~~c11*AGDM_maturity_m
Lamina_partitioining_index~~c12*AGDM_Aplus7d_m2
Lamina_partitioining_index~~c13*FE_Grains_g_CW_m.
Lamina_partitioining_index~~c14*Grain_partitioning_index_maturity
Lamina_partitioining_index~~c15*Straw_partitioning_index_maturity
Lamina_partitioining_index~~c16*Stem_partitioing_index_maturity

TrueStem_partitioing_index~~d3*Internode1_partitioing_index
TrueStem_partitioing_index~~d4*Internode2_partitioing_index
TrueStem_partitioing_index~~d5*Internode3_partitioing_index
TrueStem_partitioing_index~~d6*X2ndplus3rd_partitioining_index
TrueStem_partitioing_index~~d7*Internode1_Length
TrueStem_partitioing_index~~d8*Internode2_Length
TrueStem_partitioing_index~~d9*Internode3_Length
TrueStem_partitioing_index~~d11*AGDM_maturity_m
TrueStem_partitioing_index~~d12*AGDM_Aplus7d_m2
TrueStem_partitioing_index~~d13*FE_Grains_g_CW_m.
TrueStem_partitioing_index~~d14*Grain_partitioning_index_maturity
TrueStem_partitioing_index~~d15*Straw_partitioning_index_maturity
TrueStem_partitioing_index~~d16*Stem_partitioing_index_maturity

Internode1_partitioing_index~~e4*Internode2_partitioing_index
Internode1_partitioing_index~~e5*Internode3_partitioing_index
Internode1_partitioing_index~~e6*X2ndplus3rd_partitioining_index
Internode1_partitioing_index~~e7*Internode1_Length
Internode1_partitioing_index~~e8*Internode2_Length
Internode1_partitioing_index~~e9*Internode3_Length
Internode1_partitioing_index~~e11*AGDM_maturity_m
Internode1_partitioing_index~~e12*AGDM_Aplus7d_m2
Internode1_partitioing_index~~e13*FE_Grains_g_CW_m.
Internode1_partitioing_index~~e14*Grain_partitioning_index_maturity
Internode1_partitioing_index~~e15*Straw_partitioning_index_maturity
Internode1_partitioing_index~~e16*Stem_partitioing_index_maturity

Internode2_partitioing_index~~f5*Internode3_partitioing_index
Internode2_partitioing_index~~f6*X2ndplus3rd_partitioining_index
Internode2_partitioing_index~~f7*Internode1_Length
Internode2_partitioing_index~~f8*Internode2_Length
Internode2_partitioing_index~~f9*Internode3_Length
Internode2_partitioing_index~~f11*AGDM_maturity_m
Internode2_partitioing_index~~f12*AGDM_Aplus7d_m2
Internode2_partitioing_index~~f13*FE_Grains_g_CW_m.
Internode2_partitioing_index~~f14*Grain_partitioning_index_maturity
Internode2_partitioing_index~~f15*Straw_partitioning_index_maturity
Internode2_partitioing_index~~f16*Stem_partitioing_index_maturity

Internode3_partitioing_index~~g6*X2ndplus3rd_partitioining_index
Internode3_partitioing_index~~g7*Internode1_Length
Internode3_partitioing_index~~g8*Internode2_Length
Internode3_partitioing_index~~g9*Internode3_Length
Internode3_partitioing_index~~g11*AGDM_maturity_m
Internode3_partitioing_index~~g12*AGDM_Aplus7d_m2
Internode3_partitioing_index~~g13*FE_Grains_g_CW_m.
Internode3_partitioing_index~~g14*Grain_partitioning_index_maturity
Internode3_partitioing_index~~g15*Straw_partitioning_index_maturity
Internode3_partitioing_index~~g16*Stem_partitioing_index_maturity

X2ndplus3rd_partitioining_index~~h7*Internode1_Length
X2ndplus3rd_partitioining_index~~h8*Internode2_Length
X2ndplus3rd_partitioining_index~~h9*Internode3_Length
X2ndplus3rd_partitioining_index~~h11*AGDM_maturity_m
X2ndplus3rd_partitioining_index~~h12*AGDM_Aplus7d_m2
X2ndplus3rd_partitioining_index~~h13*FE_Grains_g_CW_m.
X2ndplus3rd_partitioining_index~~h14*Grain_partitioning_index_maturity
X2ndplus3rd_partitioining_index~~h15*Straw_partitioning_index_maturity
X2ndplus3rd_partitioining_index~~h16*Stem_partitioing_index_maturity

Internode1_Length~~i8*Internode2_Length
Internode1_Length~~i9*Internode3_Length
Internode1_Length~~i11*AGDM_maturity_m
Internode1_Length~~i12*AGDM_Aplus7d_m2
Internode1_Length~~i13*FE_Grains_g_CW_m.
Internode1_Length~~i14*Grain_partitioning_index_maturity
Internode1_Length~~i15*Straw_partitioning_index_maturity
Internode1_Length~~i16*Stem_partitioing_index_maturity

Internode2_Length~~j9*Internode3_Length
Internode2_Length~~j11*AGDM_maturity_m
Internode2_Length~~j12*AGDM_Aplus7d_m2
Internode2_Length~~j13*FE_Grains_g_CW_m.
Internode2_Length~~j14*Grain_partitioning_index_maturity
Internode2_Length~~j15*Straw_partitioning_index_maturity
Internode2_Length~~j16*Stem_partitioing_index_maturity

Internode3_Length~~k11*AGDM_maturity_m
Internode3_Length~~k12*AGDM_Aplus7d_m2
Internode3_Length~~k13*FE_Grains_g_CW_m.
Internode3_Length~~k14*Grain_partitioning_index_maturity
Internode3_Length~~k15*Straw_partitioning_index_maturity
Internode3_Length~~k16*Stem_partitioing_index_maturity

AGDM_maturity_m~~m12*AGDM_Aplus7d_m2
AGDM_maturity_m~~m13*FE_Grains_g_CW_m.
AGDM_maturity_m~~m14*Grain_partitioning_index_maturity
AGDM_maturity_m~~m15*Straw_partitioning_index_maturity
AGDM_maturity_m~~m16*Stem_partitioing_index_maturity

AGDM_Aplus7d_m2~~n13*FE_Grains_g_CW_m.
AGDM_Aplus7d_m2~~n14*Grain_partitioning_index_maturity
AGDM_Aplus7d_m2~~n15*Straw_partitioning_index_maturity
AGDM_Aplus7d_m2~~n16*Stem_partitioing_index_maturity

FE_Grains_g_CW_m.~~o14*Grain_partitioning_index_maturity
FE_Grains_g_CW_m.~~o15*Straw_partitioning_index_maturity
FE_Grains_g_CW_m.~~o16*Stem_partitioing_index_maturity

Grain_partitioning_index_maturity~~p15*Straw_partitioning_index_maturity
Grain_partitioning_index_maturity~~p16*Stem_partitioing_index_maturity

Straw_partitioning_index_maturity~~q16*Stem_partitioing_index_maturity

'
fit <- lavaan(model_1,model.type='sem',estimator = "ML",fixed.x=F,orthogonal=F,data=afterscale,missing = "ml")

fit <- lavaan(model_1,model.type='sem',estimator = "ML",fixed.x=F,orthogonal=F,data=afterscale,auto.var=T,std.ov=F, bootstrap=1000,
              auto.cov.lv.x=F,missing = "ml", int.ov.free=T, int.lv.free=F,auto.cov.y=F,do.fit = TRUE)

varTable(fit)

summary(fit,standardized=F, rsq=TRUE)

coef(fit)
parameterEstimates(fit,standardize=F)
parameter=standardizedSolution(fit,type = "std.all")

write.table(parameter,file="pathway table.csv",sep=",",row.names=F)
fit2 <- sem(model_1, data=afterscale)
semPaths(fit,'std',
         edge.label.cex=0.5,
         intercepts=F,
         exoVar = F, 
         exoCov = T,curveAdjacent = T, style = "lisrel",layout='tree',edge.color='black',color.man='black',color.lat='black',asize=1)

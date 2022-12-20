### Data analysis for Santos et al. on SGD and OA buffering ###
### Created by Nyssa Silbiger #####
### Edited on 2021-03-12 ####


#### Load libraries ############
library(tidyverse)
library(here)
library(janitor)
library(lubridate)
library(ggrepel)
library(patchwork)
library(lutz)
library(suncalc)
library(brms)
library(tidybayes)
library(bayesplot)
library(stringr)
library(broom)

##### Read in the data #############
# SGDData<-read_csv(here("Data","Timeseries_edited.csv")) %>%
#   clean_names() %>%
#   mutate(datetime = dmy_hms(paste(date,time)),
#          date = dmy(date),
#          time = hms(time),
#          tz = tz_lookup_coords(lat = lat, lon = long, method = "accurate")) 
# 
# 
# ## Get sunrise and sunset times
# 
# # the timezone function is silly... it only takes one at a time so I need to split this up
# SiteInfo_OZ<-SGDData %>%
#   select(site, date,lat, "lon" = long, tz) %>%
#   distinct() %>%
#   filter(tz == "Australia/Brisbane")
# SunOz<-getSunlightTimes(data = SiteInfo_OZ, keep = c("sunrise", "sunset"), tz ="Australia/Brisbane" )
# 
# SiteInfo_ET<-SGDData %>%
#   select(site, date,lat, "lon" = long, tz) %>%
#   distinct() %>%
#   filter(tz == "Etc/GMT-11") 
# SunET<-  getSunlightTimes(data = SiteInfo_ET, keep = c("sunrise", "sunset"), tz ="Etc/GMT-11" )
# 
# SiteInfo_As<-SGDData %>%
#   select(site, date,lat, "lon" = long, tz) %>%
#   distinct() %>%
#   filter(tz == "Asia/Makassar")
# SunAs<-  getSunlightTimes(data = SiteInfo_As, keep = c("sunrise", "sunset"), tz ="Asia/Makassar" )
# 
# SiteInfo_Pac<-SGDData %>%
#   select(site, date,lat, "lon" = long, tz) %>%
#   distinct() %>%
#   filter(tz == "Pacific/Rarotonga")
# SunPac<-  getSunlightTimes(data = SiteInfo_Pac, keep = c("sunrise", "sunset"), tz ="Pacific/Rarotonga" )
# 
# SunTimes<-bind_rows(SunOz, SunET, SunAs, SunPac) %>%
#   select(date, sunrise, sunset) %>%
#   mutate(sunrise = sunrise - days(1),
#          sunset = sunset - days(1))  # they are all off by one day
#   
# 
# # join the datasets back and add a new column that says day or night
# SGDData_times<-SGDData %>%
#   left_join(SunTimes) %>%
#   mutate(sunrise = force_tz(sunrise, "UTC"),
#          sunset = force_tz(sunset, "UTC"),
#          datetime = force_tz(datetime, "UTC"))%>%
#   mutate(DayNight = ifelse(between(datetime, sunrise, sunset), "day", "night"))

#### Exploratory plots#############

### The dates are messaged up so I needed to do it manually ###
#SGDData <- read_csv(here("data","SGDData_times.csv"))
SGDData <- read_csv(here("data","Time_Space.csv"))


# standardize the data
SGDData <- SGDData %>%
 # mutate_at(.vars = c("rn_bq_m3", "nox_u_m", "salinity"), .funs = list(log = ~log(.))) %>%
  mutate_at(.vars = c("rn_bq_m3", "nox_u_m","aou_umol_l"  ,"excess_co2_umol_kg","temperature","p_h_in","advection_rate_cm_day", "water_level_m", "p_co_in_uatm", "salinity", "co_sst"), .funs = list(std = ~scale(.)))
  #rename( "rn_bq_m3_std" ="rn_bq_m3_log_std" ,"nox_u_m_std" ="nox_u_m_log_std", "Salinty_std"="salinity_log_std")


### site level sums ####

SGD_sum<- SGDData %>%
  select(site, advection_rate_cm_day, din, nox_u_m, excess_co2_umol_kg)%>%
  group_by(site) %>%
  summarise_if(is.numeric, sum) 

SGD_mean<- SGDData %>%
  select(site, advection_rate_cm_day, din, nox_u_m, excess_co2_umol_kg, p_co_in_uatm)%>%
  group_by(site) %>%
  summarise_if(is.numeric, mean) 

P_no_onetree<-SGD_sum %>%
  filter(site != "One Tree") %>% # one tree is crazy
  ggplot(aes(x = din, y = excess_co2_umol_kg))+
  geom_point(size = 2)+
  geom_smooth(method = "lm", se = FALSE, color = "black") +
  geom_label_repel(aes(label = site))+
  geom_hline(yintercept = 0, color = "grey", linetype = 2) +
  labs(x = expression(paste("Cumulative DIN ",mu, "mol L"^-1 )),
       y = expression(paste("Cumulative Excess CO"[2]," ", mu, "mol kg"^-1 )))+
  theme_bw()


P_no_onetree<-SGD_mean %>%
  filter(site != "One Tree") %>% # one tree is crazy
  ggplot(aes(x = din, y = excess_co2_umol_kg))+
  geom_point(size = 2)+
  geom_smooth(method = "lm", se = FALSE, color = "black") +
  geom_label_repel(aes(label = site))+
  scale_y_continuous(limits = c(-5, 35))+
  geom_hline(yintercept = 0, color = "grey", linetype = 2)+
  labs(x = expression(paste("Mean DIN ",mu, "mol L"^-1 )),
      # y = expression(paste("Mean pCO"[2]," ", mu, "atm")),
       y = expression(paste("Mean Excess CO"[2]," ", mu, "mol kg"^-1)))+
  theme_bw()


mod1<-lm(p_co_in_uatm ~ excess_co2_umol_kg, data = SGD_mean %>%
           filter(site != "One Tree"))

anova(mod1)


P_onetree<-SGD_sum %>%
  ggplot(aes(x = din, y = excess_co2_umol_kg))+
  geom_point(size = 2)+
 # geom_smooth(method = "lm", se = FALSE, color = "black") +
  geom_label_repel(aes(label = site))+
  geom_hline(yintercept = 0, color = "grey", linetype = 2)+
  labs(x = expression(paste("Cumulative DIN ",mu, "mol L"^-1 )),
       y = ""
       #y = expression(paste("Cumulative Excess CO"[2]," ", mu, "mol kg"^-1)) 
       )+
  theme_bw()

P_onetree<-SGD_mean %>%
  ggplot(aes(x = din, y =excess_co2_umol_kg))+
  geom_point(size = 2)+
  scale_y_continuous(limits = c(-5, 35))+
  # geom_smooth(method = "lm", se = FALSE, color = "black") +
  geom_label_repel(aes(label = site))+
  geom_hline(yintercept = 0, color = "grey", linetype = 2)+
  labs(x = expression(paste("Cumulative DIN ",mu, "mol L"^-1 )),
       y = ""
       #y = expression(paste("Cumulative Excess CO"[2]," ", mu, "mol kg"^-1)) 
  )+
  theme_bw()

P_no_onetree + P_onetree +plot_annotation(tag_levels = "A")
  ggsave(filename = "Output/MeanN_CO2.pdf", width = 8, height = 4)


### site level averages ####

SGD_mean<- SGDData %>%
  group_by(site) %>%
  summarise_if(is.numeric, .fun = mean) 

SGD_mean %>%
  filter(site != "One Tree") %>% 
  ggplot(aes(x = din, y = p_co_in_uatm))+
  geom_point(size = 2)+
   geom_smooth(method = "lm", se = FALSE, color = "black") +
  geom_label_repel(aes(label = site))+
 # geom_hline(yintercept = 0, color = "grey", linetype = 2)+
  labs(x = expression(paste("Mean DIN ",mu, "mol L"^-1 )),
        y = expression(paste("Mean pCO"[2]," ", mu, "atm")))+
  theme_bw()

### Bayes SEM analysis ###
# SGD ~ depth
SGDmod<-bf(rn_bq_m3_std~water_level_m_std)
# N ~ SGD*DayNight
#Nmod<-bf(nox_u_m_std~DayNight*rn_bq_m3_std)
Nmod<-bf(nox_u_m_std~rn_bq_m3_std)
# Temperature ~ SGD*DayNight
Tempmod<-bf(temperature_std ~ DayNight*rn_bq_m3_std)
# pH ~ DayNight*(AOU + SGD)
#pHmod<-bf(p_h_in_std~DayNight*(aou_umol_l_std+rn_bq_m3_std))
pHmod<-bf(p_h_in_std~aou_umol_l_std+rn_bq_m3_std)
pCO2mod<-bf(p_co_in_uatm_std~aou_umol_l_std+rn_bq_m3_std)

# AOU ~ DayNight*(N+Temperature)
#AOUmod<-bf(aou_umol_l_std~DayNight*(temperature_std+nox_u_m_std))

AOUmod<-bf(aou_umol_l_std~temperature_std+DayNight*nox_u_m_std)

# excess CO2 ~ DayNight*(SGD +AOU)
# CO2mod<-bf(excess_co2_umol_kg_std ~DayNight*(rn_bq_m3_std+aou_umol_l_std))

# without day night interaction
CO2mod<-bf(excess_co2_umol_kg_std ~rn_bq_m3_std+aou_umol_l_std)
## Later
# Delta TA ~ pH + temp

site <-unique(SGDData$site)[1]

# Function to run Bayesian SEM and make the posterior predictive checks and plot marginal effects
RunSEM<-function(site){
  

fit_brms <- brm(#SGDmod+
                    Nmod+
                    Tempmod+ 
                    #pHmod+
                    AOUmod +
                    pCO2mod+
                    set_rescor(FALSE),
                  data=SGDData[SGDData$site == site,]
                  ,cores=4, chains = 3)
# calculate LOO (leave one out) diagnostics
SGD_loo<-loo(fit_brms, reloo = TRUE) # looks good!

# p1<-pp_check(fit_brms, resp="rnbqm3std") +
#   scale_color_manual(values=c("red", "black"))+
#   ggtitle("Radon")
p2<-pp_check(fit_brms, resp="noxumstd") +
  scale_color_manual(values=c("red", "black"))+
  ggtitle("N+N")
p3<-pp_check(fit_brms, resp="temperaturestd") +
  scale_color_manual(values=c("red", "black"))+
  ggtitle("Temperature")
# p4<-pp_check(fit_brms, resp="phinstd") +
#   scale_color_manual(values=c("red", "black"))+
#   ggtitle("pH")
p4<-pp_check(fit_brms, resp="pcoinuatmstd") +
  scale_color_manual(values=c("red", "black"))+
  ggtitle("pCO2")

p5<-pp_check(fit_brms, resp="aouumollstd") +
  scale_color_manual(values=c("red", "black"))+
  ggtitle("AOU")
# p6<-pp_check(fit_brms, resp="excessco2umolkgstd") +
#   scale_color_manual(values=c("red", "black"))+
#   ggtitle("Excess CO2")

#p1+
  p2+p3+p4+p5+
    #p6+
    plot_layout(guides = "collect") +
  plot_annotation(title = 'Heron Island Posterior Predictive Checks', tag_levels = "A")
  ggsave(here("Output",paste(site,"Posteriorchecks.pdf")), width = 5, height = 5)

# plot the results
# Model 1
# SGD ~ depth

# R<-conditional_effects(fit_brms, "water_level_m_std", resp = "rnbqm3std", method = "predict", resolution = 1000)
# R1<-R$rnbqm3std.rnbqm3std_water_level_m_std %>% # back transform the scaled effects for the plot
#   mutate(estimate = estimate__*attr(SGDData$rn_bq_m3_std,"scaled:scale")+attr(SGDData$rn_bq_m3_std,"scaled:center"),
#          lower = lower__*attr(SGDData$rn_bq_m3_std,"scaled:scale")+attr(SGDData$rn_bq_m3_std,"scaled:center"),
#          upper = upper__*attr(SGDData$rn_bq_m3_std,"scaled:scale")+attr(SGDData$rn_bq_m3_std,"scaled:center")) %>%
#   mutate(WaterLevel = water_level_m_std*attr (SGDData$water_level_m_std,"scaled:scale")+attr(SGDData$water_level_m_std,"scaled:center")
#   )%>%
#   ggplot()+ # back trasform the log transformed data for better visual
#   geom_line(aes(x = WaterLevel, y = estimate), lwd = 1, color = 'grey')+
#   geom_ribbon(aes(x = WaterLevel,ymin=lower, ymax=upper), linetype=1.5, alpha=0.3, fill = "grey")+
#   geom_point(data = SGDData[SGDData$site==site,], aes(x = water_level_m, y = rn_bq_m3, color = DayNight)) +
#   xlab("Water level (m)")+
#   ylab(expression(atop("Radon", paste("(bq m"^3,")"))))+
#   ggtitle("Model 1")+
#  # scale_x_continuous(breaks = c(0.2,1,5,10,25))+
# #  scale_y_continuous(breaks = c(0,0.1,1,5,10,30))+
#   theme_minimal()+
#   theme(plot.title = element_text(hjust = 0.5, size = 14), legend.position = 'bottom')

## Model 2
# N ~ SGD*DayNight
R<-conditional_effects(fit_brms, "rn_bq_m3_std:DayNight", resp = "noxumstd", method = "predict", resolution = 1000)
R2<-R$`noxumstd.noxumstd_rn_bq_m3_std:DayNight` %>% # back transform the scaled effects for the plot
  mutate(estimate = estimate__*attr(SGDData$nox_u_m_std,"scaled:scale")+attr(SGDData$nox_u_m_std,"scaled:center"),
         lower = lower__*attr(SGDData$nox_u_m_std,"scaled:scale")+attr(SGDData$nox_u_m_std,"scaled:center"),
         upper = upper__*attr(SGDData$nox_u_m_std,"scaled:scale")+attr(SGDData$nox_u_m_std,"scaled:center")) %>%
  mutate(radon = rn_bq_m3_std*attr (SGDData$rn_bq_m3_std,"scaled:scale")+attr(SGDData$rn_bq_m3_std,"scaled:center")
  )%>%
  ggplot()+ # back trasform the log transformed data for better visual
  geom_line(aes(x = radon, y = estimate, color = DayNight), lwd = 1)+
  geom_ribbon(aes(x = radon,ymin=lower, ymax=upper, fill = DayNight), linetype=1.5, alpha=0.3)+
  geom_point(data = SGDData[SGDData$site==site,], aes(x = rn_bq_m3, y = nox_u_m, color = DayNight)) +
  ylab(expression(atop("Nitrate + Nitrite", paste("(",mu, "mol L"^-1,")"))))+
  xlab(expression(atop("Radon", paste("(bq m"^3,")"))))+
  ggtitle("Model 2")+
  # scale_x_continuous(breaks = c(0.2,1,5,10,25))+
  #  scale_y_continuous(breaks = c(0,0.1,1,5,10,30))+
  theme_minimal()+
  theme(plot.title = element_text(hjust = 0.5, size = 14), legend.position = 'bottom')

  # Model 3
  # Temperature ~ SGD*DayNight
  
R<-conditional_effects(fit_brms, "rn_bq_m3_std:DayNight", resp = "temperaturestd", method = "predict", resolution = 1000)
  R3<-R$`temperaturestd.temperaturestd_rn_bq_m3_std:DayNight` %>% # back transform the scaled effects for the plot
    mutate(estimate = estimate__*attr(SGDData$temperature_std,"scaled:scale")+attr(SGDData$temperature_std,"scaled:center"),
           lower = lower__*attr(SGDData$temperature_std,"scaled:scale")+attr(SGDData$temperature_std,"scaled:center"),
           upper = upper__*attr(SGDData$temperature_std,"scaled:scale")+attr(SGDData$temperature_std,"scaled:center")) %>%
    mutate(radon = rn_bq_m3_std*attr (SGDData$rn_bq_m3_std,"scaled:scale")+attr(SGDData$rn_bq_m3_std,"scaled:center")
    )%>%
    ggplot()+ # back trasform the log transformed data for better visual
    geom_line(aes(x = radon, y = estimate, color = DayNight), lwd = 1)+
    geom_ribbon(aes(x = radon,ymin=lower, ymax=upper, fill = DayNight), linetype=1.5, alpha=0.3)+
    geom_point(data = SGDData[SGDData$site==site,], aes(x = rn_bq_m3, y = temperature, color = DayNight)) +
    ylab(expression(atop("Temperature",paste("(", degree, "C)"))))+
  xlab(expression(atop("Radon", paste("(bq m"^3,")"))))+
    ggtitle("Model 3")+
    # scale_x_continuous(breaks = c(0.2,1,5,10,25))+
    #  scale_y_continuous(breaks = c(0,0.1,1,5,10,30))+
    theme_minimal()+
    theme(plot.title = element_text(hjust = 0.5, size = 14), legend.position = 'bottom')
  
  # # Model 4a
  # # pH ~ DayNight*(AOU + SGD)
  # R<-conditional_effects(fit_brms, "rn_bq_m3_std:DayNight", resp = "phinstd", method = "predict", resolution = 1000)
  # R4<-R$`phinstd.phinstd_rn_bq_m3_std:DayNight` %>% # back transform the scaled effects for the plot
  #   mutate(estimate = estimate__*attr(SGDData$p_h_in_std,"scaled:scale")+attr(SGDData$p_h_in_std,"scaled:center"),
  #          lower = lower__*attr(SGDData$p_h_in_std,"scaled:scale")+attr(SGDData$p_h_in_std,"scaled:center"),
  #          upper = upper__*attr(SGDData$p_h_in_std,"scaled:scale")+attr(SGDData$p_h_in_std,"scaled:center")) %>%
  #   mutate(radon = rn_bq_m3_std*attr (SGDData$rn_bq_m3_std,"scaled:scale")+attr(SGDData$rn_bq_m3_std,"scaled:center")
  #   )%>%
  #   ggplot()+ # back trasform the log transformed data for better visual
  #   geom_line(aes(x = radon, y = estimate, color = DayNight), lwd = 1)+
  #   geom_ribbon(aes(x = radon,ymin=lower, ymax=upper, fill = DayNight), linetype=1.5, alpha=0.3)+
  #   geom_point(data = SGDData[SGDData$site==site,], aes(x = rn_bq_m3, y = p_h_in, color = DayNight)) +
  #   ylab("pH")+
  #   xlab(expression(atop("Radon", paste("(bq m"^3,")"))))+
  #   ggtitle("Model 4")+
  #   # scale_x_continuous(breaks = c(0.2,1,5,10,25))+
  #   #  scale_y_continuous(breaks = c(0,0.1,1,5,10,30))+
  #   theme_minimal()+
  #   theme(plot.title = element_text(hjust = 0.5, size = 14), legend.position = 'bottom')
  # 
  # Model 4a
  # pH ~ DayNight*(AOU + SGD)
  # R<-conditional_effects(fit_brms, "rn_bq_m3_std", resp = "phinstd", method = "predict", resolution = 1000)
  # R4<-R$phinstd.phinstd_rn_bq_m3_std %>% # back transform the scaled effects for the plot
  #   mutate(estimate = estimate__*attr(SGDData$p_h_in_std,"scaled:scale")+attr(SGDData$p_h_in_std,"scaled:center"),
  #          lower = lower__*attr(SGDData$p_h_in_std,"scaled:scale")+attr(SGDData$p_h_in_std,"scaled:center"),
  #          upper = upper__*attr(SGDData$p_h_in_std,"scaled:scale")+attr(SGDData$p_h_in_std,"scaled:center")) %>%
  #   mutate(radon = rn_bq_m3_std*attr (SGDData$rn_bq_m3_std,"scaled:scale")+attr(SGDData$rn_bq_m3_std,"scaled:center")
  #   )%>%
  #   ggplot()+ # back trasform the log transformed data for better visual
  #   geom_line(aes(x = radon, y = estimate), lwd = 1)+
  #   geom_ribbon(aes(x = radon,ymin=lower, ymax=upper), linetype=1.5, alpha=0.3)+
  #   geom_point(data = SGDData[SGDData$site==site,], aes(x = rn_bq_m3, y = p_h_in)) +
  #   ylab("pH")+
  #   xlab(expression(atop("Radon", paste("(bq m"^3,")"))))+
  #   ggtitle("Model 4")+
  #   # scale_x_continuous(breaks = c(0.2,1,5,10,25))+
  #   #  scale_y_continuous(breaks = c(0,0.1,1,5,10,30))+
  #   theme_minimal()+
  #   theme(plot.title = element_text(hjust = 0.5, size = 14), legend.position = 'bottom')
  
  # Model 4a
  # pCO2 ~ DayNight*(AOU + SGD)
  R<-conditional_effects(fit_brms, "rn_bq_m3_std", resp = "pcoinuatmstd", method = "predict", resolution = 1000)
  R4<-R$pcoinuatmstd.pcoinuatmstd_rn_bq_m3_std %>% # back transform the scaled effects for the plot
    mutate(estimate = estimate__*attr(SGDData$p_co_in_uatm_std,"scaled:scale")+attr(SGDData$p_co_in_uatm_std,"scaled:center"),
           lower = lower__*attr(SGDData$p_co_in_uatm_std,"scaled:scale")+attr(SGDData$p_co_in_uatm_std,"scaled:center"),
           upper = upper__*attr(SGDData$p_co_in_uatm_std,"scaled:scale")+attr(SGDData$p_co_in_uatm_std,"scaled:center")) %>%
    mutate(radon = rn_bq_m3_std*attr (SGDData$rn_bq_m3_std,"scaled:scale")+attr(SGDData$rn_bq_m3_std,"scaled:center")
    )%>%
    ggplot()+ # back trasform the log transformed data for better visual
    geom_line(aes(x = radon, y = estimate), lwd = 1)+
    geom_ribbon(aes(x = radon,ymin=lower, ymax=upper), linetype=1.5, alpha=0.3)+
    geom_point(data = SGDData[SGDData$site==site,], aes(x = rn_bq_m3, y = p_co_in_uatm)) +
    ylab(expression(paste("pCO"[2]," ", mu, "atm")))+
    xlab(expression(atop("Radon", paste("(bq m"^3,")"))))+
    ggtitle("Model 4")+
    # scale_x_continuous(breaks = c(0.2,1,5,10,25))+
    #  scale_y_continuous(breaks = c(0,0.1,1,5,10,30))+
    theme_minimal()+
    theme(plot.title = element_text(hjust = 0.5, size = 14), legend.position = 'bottom')
  
  # Model 4b
  # pH ~ DayNight*(AOU + SGD)
  # R<-conditional_effects(fit_brms, "aou_umol_l_std:DayNight", resp = "phinstd", method = "predict", resolution = 1000)
  # R4b<-R$`phinstd.phinstd_aou_umol_l_std:DayNight` %>% # back transform the scaled effects for the plot
  #   mutate(estimate = estimate__*attr(SGDData$p_h_in_std,"scaled:scale")+attr(SGDData$p_h_in_std,"scaled:center"),
  #          lower = lower__*attr(SGDData$p_h_in_std,"scaled:scale")+attr(SGDData$p_h_in_std,"scaled:center"),
  #          upper = upper__*attr(SGDData$p_h_in_std,"scaled:scale")+attr(SGDData$p_h_in_std,"scaled:center")) %>%
  #   mutate(aou = aou_umol_l_std*attr (SGDData$aou_umol_l_std,"scaled:scale")+attr(SGDData$aou_umol_l_std,"scaled:center")
  #   )%>%
  #   ggplot()+ # back trasform the log transformed data for better visual
  #   geom_line(aes(x = aou, y = estimate, color = DayNight), lwd = 1)+
  #   geom_ribbon(aes(x = aou,ymin=lower, ymax=upper, fill = DayNight), linetype=1.5, alpha=0.3)+
  #   geom_point(data = SGDData[SGDData$site==site,], aes(x = aou_umol_l, y = p_h_in, color = DayNight)) +
  #   ylab("pH")+
  #   xlab(expression(atop("AOU", paste("(",mu, "mol L"^-1,")"))))+
  #   ggtitle("Model 4")+
  #   # scale_x_continuous(breaks = c(0.2,1,5,10,25))+
  #   #  scale_y_continuous(breaks = c(0,0.1,1,5,10,30))+
  #   theme_minimal()+
  #   theme(plot.title = element_text(hjust = 0.5, size = 14), legend.position = 'bottom')

  # pH ~ DayNight*(AOU + SGD)
  # R<-conditional_effects(fit_brms, "aou_umol_l_std", resp = "phinstd", method = "predict", resolution = 1000)
  # R4b<-R$phinstd.phinstd_aou_umol_l_std %>% # back transform the scaled effects for the plot
  #   mutate(estimate = estimate__*attr(SGDData$p_h_in_std,"scaled:scale")+attr(SGDData$p_h_in_std,"scaled:center"),
  #          lower = lower__*attr(SGDData$p_h_in_std,"scaled:scale")+attr(SGDData$p_h_in_std,"scaled:center"),
  #          upper = upper__*attr(SGDData$p_h_in_std,"scaled:scale")+attr(SGDData$p_h_in_std,"scaled:center")) %>%
  #   mutate(aou = aou_umol_l_std*attr (SGDData$aou_umol_l_std,"scaled:scale")+attr(SGDData$aou_umol_l_std,"scaled:center")
  #   )%>%
  #   ggplot()+ # back trasform the log transformed data for better visual
  #   geom_line(aes(x = aou, y = estimate), lwd = 1)+
  #   geom_ribbon(aes(x = aou,ymin=lower, ymax=upper), linetype=1.5, alpha=0.3)+
  #   geom_point(data = SGDData[SGDData$site==site,], aes(x = aou_umol_l, y = p_h_in)) +
  #   ylab("pH")+
  #   xlab(expression(atop("AOU", paste("(",mu, "mol L"^-1,")"))))+
  #   ggtitle("Model 4")+
  #   # scale_x_continuous(breaks = c(0.2,1,5,10,25))+
  #   #  scale_y_continuous(breaks = c(0,0.1,1,5,10,30))+
  #   theme_minimal()+
  #   theme(plot.title = element_text(hjust = 0.5, size = 14), legend.position = 'bottom')

  R<-conditional_effects(fit_brms, "aou_umol_l_std", resp = "pcoinuatmstd", method = "predict", resolution = 1000)
  R4b<-R$pcoinuatmstd.pcoinuatmstd_aou_umol_l_std %>% # back transform the scaled effects for the plot
    mutate(estimate = estimate__*attr(SGDData$p_co_in_uatm_std,"scaled:scale")+attr(SGDData$p_co_in_uatm_std,"scaled:center"),
           lower = lower__*attr(SGDData$p_co_in_uatm_std,"scaled:scale")+attr(SGDData$p_co_in_uatm_std,"scaled:center"),
           upper = upper__*attr(SGDData$p_co_in_uatm_std,"scaled:scale")+attr(SGDData$p_co_in_uatm_std,"scaled:center")) %>%
    mutate(aou = aou_umol_l_std*attr (SGDData$aou_umol_l_std,"scaled:scale")+attr(SGDData$aou_umol_l_std,"scaled:center")
    )%>%
    ggplot()+ # back trasform the log transformed data for better visual
    geom_line(aes(x = aou, y = estimate), lwd = 1)+
    geom_ribbon(aes(x = aou,ymin=lower, ymax=upper), linetype=1.5, alpha=0.3)+
    geom_point(data = SGDData[SGDData$site==site,], aes(x = aou_umol_l, y = p_co_in_uatm)) +
    ylab(expression(paste("pCO"[2]," ", mu, "atm")))+
    xlab(expression(atop("AOU", paste("(",mu, "mol L"^-1,")"))))+
    ggtitle("Model 4")+
    # scale_x_continuous(breaks = c(0.2,1,5,10,25))+
    #  scale_y_continuous(breaks = c(0,0.1,1,5,10,30))+
    theme_minimal()+
    theme(plot.title = element_text(hjust = 0.5, size = 14), legend.position = 'bottom')
  
  
  # Model 5
  # AOU ~ DayNight*(N+Temperature)
  R<-conditional_effects(fit_brms, "temperature_std", resp = "aouumollstd", method = "predict", resolution = 1000)
  R5<-R$aouumollstd.aouumollstd_temperature_std %>% # back transform the scaled effects for the plot
    mutate(estimate = estimate__*attr(SGDData$aou_umol_l_std,"scaled:scale")+attr(SGDData$aou_umol_l_std,"scaled:center"),
           lower = lower__*attr(SGDData$aou_umol_l_std,"scaled:scale")+attr(SGDData$aou_umol_l_std,"scaled:center"),
           upper = upper__*attr(SGDData$aou_umol_l_std,"scaled:scale")+attr(SGDData$aou_umol_l_std,"scaled:center")) %>%
    mutate(temperature = temperature_std*attr (SGDData$temperature_std,"scaled:scale")+attr(SGDData$temperature_std,"scaled:center")
    )%>%
    ggplot()+ # back trasform the log transformed data for better visual
    geom_line(aes(x = temperature, y = estimate), lwd = 1)+
    geom_ribbon(aes(x = temperature,ymin=lower, ymax=upper), linetype=1.5, alpha=0.3)+
    geom_point(data = SGDData[SGDData$site==site,], aes(x = temperature, y = aou_umol_l)) +
    xlab(expression(atop("Temperature",paste("(", degree, "C)"))))+
    ylab(expression(atop("AOU", paste("(",mu, "mol L"^-1,")"))))+
    ggtitle("Model 5")+
    # scale_x_continuous(breaks = c(0.2,1,5,10,25))+
    #  scale_y_continuous(breaks = c(0,0.1,1,5,10,30))+
    theme_minimal()+
    theme(plot.title = element_text(hjust = 0.5, size = 14), legend.position = 'bottom')
  
  # R<-conditional_effects(fit_brms, "temperature_std:DayNight", resp = "aouumollstd", method = "predict", resolution = 1000)
  # R5<-R$`aouumollstd.aouumollstd_temperature_std:DayNight` %>% # back transform the scaled effects for the plot
  #   mutate(estimate = estimate__*attr(SGDData$aou_umol_l_std,"scaled:scale")+attr(SGDData$aou_umol_l_std,"scaled:center"),
  #          lower = lower__*attr(SGDData$aou_umol_l_std,"scaled:scale")+attr(SGDData$aou_umol_l_std,"scaled:center"),
  #          upper = upper__*attr(SGDData$aou_umol_l_std,"scaled:scale")+attr(SGDData$aou_umol_l_std,"scaled:center")) %>%
  #   mutate(temperature = temperature_std*attr (SGDData$temperature_std,"scaled:scale")+attr(SGDData$temperature_std,"scaled:center")
  #   )%>%
  #   ggplot()+ # back trasform the log transformed data for better visual
  #   geom_line(aes(x = temperature, y = estimate, color = DayNight), lwd = 1)+
  #   geom_ribbon(aes(x = temperature,ymin=lower, ymax=upper, fill = DayNight), linetype=1.5, alpha=0.3)+
  #   geom_point(data = SGDData[SGDData$site==site,], aes(x = temperature, y = aou_umol_l, color = DayNight)) +
  #   xlab(expression(atop("Temperature",paste("(", degree, "C)"))))+
  #   ylab(expression(atop("AOU", paste("(",mu, "mol L"^-1,")"))))+
  #   ggtitle("Model 5")+
  #   # scale_x_continuous(breaks = c(0.2,1,5,10,25))+
  #   #  scale_y_continuous(breaks = c(0,0.1,1,5,10,30))+
  #   theme_minimal()+
  #   theme(plot.title = element_text(hjust = 0.5, size = 14), legend.position = 'bottom')
  # 
  # # Model 5b
  # # AOU ~ DayNight*(N+Temperature)
  # 
  R<-conditional_effects(fit_brms, "nox_u_m_std:DayNight", resp = "aouumollstd", method = "predict", resolution = 1000)
  R5b<-R$`aouumollstd.aouumollstd_nox_u_m_std:DayNight` %>% # back transform the scaled effects for the plot
    mutate(estimate = estimate__*attr(SGDData$aou_umol_l_std,"scaled:scale")+attr(SGDData$aou_umol_l_std,"scaled:center"),
           lower = lower__*attr(SGDData$aou_umol_l_std,"scaled:scale")+attr(SGDData$aou_umol_l_std,"scaled:center"),
           upper = upper__*attr(SGDData$aou_umol_l_std,"scaled:scale")+attr(SGDData$aou_umol_l_std,"scaled:center")) %>%
    mutate(NN = nox_u_m_std*attr (SGDData$nox_u_m_std,"scaled:scale")+attr(SGDData$nox_u_m_std,"scaled:center")
    )%>%
    ggplot()+ # back trasform the log transformed data for better visual
    geom_line(aes(x = NN, y = estimate, color = DayNight), lwd = 1)+
    geom_ribbon(aes(x = NN,ymin=lower, ymax=upper, fill = DayNight), linetype=1.5, alpha=0.3)+
    geom_point(data = SGDData[SGDData$site==site,], aes(x = nox_u_m, y = aou_umol_l, color = DayNight)) +
    xlab(expression(atop("Nitrate + Nitrite", paste("(",mu, "mol L"^-1,")"))))+
    ylab(expression(atop("AOU", paste("(",mu, "mol L"^-1,")"))))+
    ggtitle("Model 5")+
    # scale_x_continuous(breaks = c(0.2,1,5,10,25))+
    #  scale_y_continuous(breaks = c(0,0.1,1,5,10,30))+
    theme_minimal()+
    theme(plot.title = element_text(hjust = 0.5, size = 14), legend.position = 'bottom')
  
# Model 6
# excess CO2 ~ DayNight*(SGD +AOU)
  
# R<-conditional_effects(fit_brms, "rn_bq_m3_std:DayNight", resp = "excessco2umolkgstd", method = "predict", resolution = 1000)
#   R6<-R$`excessco2umolkgstd.excessco2umolkgstd_rn_bq_m3_std:DayNight` %>% # back transform the scaled effects for the plot
#     mutate(estimate = estimate__*attr(SGDData$excess_co2_umol_kg_std,"scaled:scale")+attr(SGDData$excess_co2_umol_kg_std,"scaled:center"),
#            lower = lower__*attr(SGDData$excess_co2_umol_kg_std,"scaled:scale")+attr(SGDData$excess_co2_umol_kg_std,"scaled:center"),
#            upper = upper__*attr(SGDData$excess_co2_umol_kg_std,"scaled:scale")+attr(SGDData$excess_co2_umol_kg_std,"scaled:center")) %>%
#     mutate(radon = rn_bq_m3_std*attr(SGDData$rn_bq_m3_std,"scaled:scale")+attr(SGDData$rn_bq_m3_std,"scaled:center")
#     )%>%
#     ggplot()+ # back trasform the log transformed data for better visual
#     geom_line(aes(x = radon, y = estimate, color = DayNight), lwd = 1)+
#     geom_ribbon(aes(x = radon,ymin=lower, ymax=upper, fill = DayNight), linetype=1.5, alpha=0.3)+
#     geom_point(data = SGDData[SGDData$site==site,], aes(x = rn_bq_m3, y = excess_co2_umol_kg, color = DayNight)) +
#     xlab(expression(atop("Radon", paste("(bq m"^3,")"))))+
#     ylab(expression(atop("Excess CO"[2], paste("(",mu, "mol kg"^-1,")"))))+
#     ggtitle("Model 6")+
#     # scale_x_continuous(breaks = c(0.2,1,5,10,25))+
#     #  scale_y_continuous(breaks = c(0,0.1,1,5,10,30))+
#     theme_minimal()+
#     theme(plot.title = element_text(hjust = 0.5, size = 14), legend.position = 'bottom')

  # R<-conditional_effects(fit_brms, "rn_bq_m3_std", resp = "excessco2umolkgstd", method = "predict", resolution = 1000)
  # R6<-R$excessco2umolkgstd.excessco2umolkgstd_rn_bq_m3_std %>% # back transform the scaled effects for the plot
  #   mutate(estimate = estimate__*attr(SGDData$excess_co2_umol_kg_std,"scaled:scale")+attr(SGDData$excess_co2_umol_kg_std,"scaled:center"),
  #          lower = lower__*attr(SGDData$excess_co2_umol_kg_std,"scaled:scale")+attr(SGDData$excess_co2_umol_kg_std,"scaled:center"),
  #          upper = upper__*attr(SGDData$excess_co2_umol_kg_std,"scaled:scale")+attr(SGDData$excess_co2_umol_kg_std,"scaled:center")) %>%
  #   mutate(radon = rn_bq_m3_std*attr(SGDData$rn_bq_m3_std,"scaled:scale")+attr(SGDData$rn_bq_m3_std,"scaled:center")
  #   )%>%
  #   ggplot()+ # back trasform the log transformed data for better visual
  #   geom_line(aes(x = radon, y = estimate), lwd = 1)+
  #   geom_ribbon(aes(x = radon,ymin=lower, ymax=upper), linetype=1.5, alpha=0.3)+
  #   geom_point(data = SGDData[SGDData$site==site,], aes(x = rn_bq_m3, y = excess_co2_umol_kg)) +
  #   xlab(expression(atop("Radon", paste("(bq m"^3,")"))))+
  #   ylab(expression(atop("Excess CO"[2], paste("(",mu, "mol kg"^-1,")"))))+
  #   ggtitle("Model 6")+
  #   # scale_x_continuous(breaks = c(0.2,1,5,10,25))+
  #   #  scale_y_continuous(breaks = c(0,0.1,1,5,10,30))+
  #   theme_minimal()+
  #   theme(plot.title = element_text(hjust = 0.5, size = 14), legend.position = 'bottom')
  
  # Model 6b
  # excess CO2 ~ DayNight*(SGD +AOU)
  
  # R<-conditional_effects(fit_brms, "aou_umol_l_std:DayNight", resp = "excessco2umolkgstd", method = "predict", resolution = 1000)
  # R6b<-R$`excessco2umolkgstd.excessco2umolkgstd_aou_umol_l_std:DayNight` %>% # back transform the scaled effects for the plot
  #   mutate(estimate = estimate__*attr(SGDData$excess_co2_umol_kg_std,"scaled:scale")+attr(SGDData$excess_co2_umol_kg_std,"scaled:center"),
  #          lower = lower__*attr(SGDData$excess_co2_umol_kg_std,"scaled:scale")+attr(SGDData$excess_co2_umol_kg_std,"scaled:center"),
  #          upper = upper__*attr(SGDData$excess_co2_umol_kg_std,"scaled:scale")+attr(SGDData$excess_co2_umol_kg_std,"scaled:center")) %>%
  #   mutate(aou = aou_umol_l_std*attr(SGDData$aou_umol_l_std,"scaled:scale")+attr(SGDData$aou_umol_l_std,"scaled:center")
  #   )%>%
  #   ggplot()+ # back trasform the log transformed data for better visual
  #   geom_line(aes(x = aou, y = estimate, color = DayNight), lwd = 1)+
  #   geom_ribbon(aes(x = aou,ymin=lower, ymax=upper, fill = DayNight), linetype=1.5, alpha=0.3)+
  #   geom_point(data = SGDData[SGDData$site==site,], aes(x = aou_umol_l, y = excess_co2_umol_kg, color = DayNight)) +
  #   xlab(expression(atop("AOU", paste("(",mu, "mol L"^-1,")"))))+
  #   ylab(expression(atop("Excess CO"[2], paste("(",mu, "mol kg"^-1,")"))))+
  #   ggtitle("Model 6")+
  #   # scale_x_continuous(breaks = c(0.2,1,5,10,25))+
  #   #  scale_y_continuous(breaks = c(0,0.1,1,5,10,30))+
  #   theme_minimal()+
  #   theme(plot.title = element_text(hjust = 0.5, size = 14), legend.position = 'bottom')
  # 
  # R<-conditional_effects(fit_brms, "aou_umol_l_std", resp = "excessco2umolkgstd", method = "predict", resolution = 1000)
  # R6b<-R$excessco2umolkgstd.excessco2umolkgstd_aou_umol_l_std %>% # back transform the scaled effects for the plot
  #   mutate(estimate = estimate__*attr(SGDData$excess_co2_umol_kg_std,"scaled:scale")+attr(SGDData$excess_co2_umol_kg_std,"scaled:center"),
  #          lower = lower__*attr(SGDData$excess_co2_umol_kg_std,"scaled:scale")+attr(SGDData$excess_co2_umol_kg_std,"scaled:center"),
  #          upper = upper__*attr(SGDData$excess_co2_umol_kg_std,"scaled:scale")+attr(SGDData$excess_co2_umol_kg_std,"scaled:center")) %>%
  #   mutate(aou = aou_umol_l_std*attr(SGDData$aou_umol_l_std,"scaled:scale")+attr(SGDData$aou_umol_l_std,"scaled:center")
  #   )%>%
  #   ggplot()+ # back trasform the log transformed data for better visual
  #   geom_line(aes(x = aou, y = estimate), lwd = 1)+
  #   geom_ribbon(aes(x = aou,ymin=lower, ymax=upper), linetype=1.5, alpha=0.3)+
  #   geom_point(data = SGDData[SGDData$site==site,], aes(x = aou_umol_l, y = excess_co2_umol_kg)) +
  #   xlab(expression(atop("AOU", paste("(",mu, "mol L"^-1,")"))))+
  #   ylab(expression(atop("Excess CO"[2], paste("(",mu, "mol kg"^-1,")"))))+
  #   ggtitle("Model 6")+
  #   # scale_x_continuous(breaks = c(0.2,1,5,10,25))+
  #   #  scale_y_continuous(breaks = c(0,0.1,1,5,10,30))+
  #   theme_minimal()+
  #   theme(plot.title = element_text(hjust = 0.5, size = 14), legend.position = 'bottom')
  
  ## bring them all together in patchwork
  R<-(R2|R3)/(R4|R4b)/(R5|R5b)+
  #/(R6|R6b)+
    plot_layout(guides = "collect")+
    plot_annotation(tag_levels = "A")&
    theme(axis.title.x = element_text(size = 14),
          axis.title.y = element_text(size = 14),
          axis.text = element_text(size = 11))
  
  ggsave(here("Output",paste(site,"marginaleffects.pdf")),R, width = 12, height = 18, useDingbats = FALSE)
  
  var_name <- paste(fit_brms, site, sep="_") # Construct the name
  var_name <- str_replace_all(var_name, " ", "")   # remove white space
  return(fit_brms)
  
}

# run the SEM for Heron Island
HI_fit<-RunSEM(site ="Heron Island")
CI_fit<-RunSEM(site ="Cook Islands")
LI_fit<-RunSEM(site ="Lizard Island")
OTI_fit<-RunSEM(site ="One Tree")
LHI_fit<-RunSEM(site ="Lord Howe Island")
NP_fit<-RunSEM(site ="Nusa Penida")

## Extract the posteriors

## get the posterior
get_post<-function(site, fit){

post <- posterior_samples(fit)

Cof<-post %>% 
  select(starts_with("b"),-ends_with("Intercept")) %>%
  gather() %>% 
  group_by(key)%>%
  median_hdci()%>%
  mutate(sig = ifelse(sign(.lower)==sign(.upper),'yes','no'))%>%# if not significant make it grey
  separate(col = key,into = c("b", "dependent", "independent"),sep = "_")%>% #loose the b and bring the values back together
  # mutate(key = paste(dependent, independent))%>%
  mutate(dependent = factor(dependent, levels = c("aouumollstd","pcoinuatmstd","noxumstd","temperaturestd")))%>%
  mutate(dependent  = recode(dependent, aouumollstd = "AOU",
                             pcoinuatmstd = "pCO2", 
                             noxumstd = "NOx", 
                           #  phinstd = "pH", 
                             temperaturestd  ="Temperature"), Site = site)

return(Cof)
}

HICoF<-get_post(site = "Heron Island", fit = HI_fit)
CICoF<-get_post(site = "Cook Islands", fit = CI_fit)
LICoF<-get_post(site = "Lizzard Island", fit = LI_fit)
OTICoF<-get_post(site = "One Tree", fit = OTI_fit)
LHICoF<-get_post(site = "Lord Howe Island", fit = LHI_fit)
NPCoF<-get_post(site = "Nusa Penida", fit = NP_fit)


AllCoefs<- HICoF %>%
  bind_rows(CICoF)%>%
  bind_rows(LICoF)%>%
  bind_rows(OTICoF)%>%
  bind_rows(LHICoF)%>%
  bind_rows(NPCoF)

#Make the plot
CoefPlot<-AllCoefs%>%
  # ggplot(aes(x = value, y = reorder(independent, value), shape = sig, color = Site)) +  # note how we used `reorder()` to arrange the coefficients
  ggplot(aes(x = value, y = independent, shape = sig, color = Site)) +  # note 
  geom_vline(xintercept = 0, alpha = 1/10, color = 'firebrick4') +
  geom_point(size = 2)+
  geom_errorbarh(aes(xmin = .lower, xmax = .upper), height = 0)+
  #scale_alpha_manual(values = c(0.2,1))+
  scale_shape_manual(values = c(1,16))+
  scale_color_brewer(palette = "Set2")+
#  scale_color_manual(values = c("firebrick4", "orange"), name = " ")+
  
  labs(#title = "Black Point",
    x = NULL, y = NULL) +
  theme_bw() +
  guides(shape = FALSE)+
  theme(legend.title = element_blank(),
    panel.grid = element_blank(),
        panel.grid.major.y = element_line(color = alpha("firebrick4", 1/4), linetype = 3),
        axis.text.y  = element_text(hjust = 0),
        axis.ticks.y = element_blank(),
        legend.position = "bottom",
        legend.text=element_text(size=14),
        strip.background = element_blank(),
        strip.text = element_text(size = 14, face = "bold")
  )+
  facet_grid(~dependent, scales = "free_y", space='free')
  ggplot2::ggsave("Output/coefficientsAll.pdf", width = 10, height = 5, useDingbats = FALSE)


IndivPlots<-function(site){
  AllCoefs%>%
    filter(Site == site) %>%
    ggplot(aes(x = value, y = independent, shape = sig)) +  # note how we used `reorder()` to arrange the coefficients
    geom_vline(xintercept = 0, alpha = 1/10, color = 'firebrick4') +
    geom_point(size = 2)+
    geom_errorbarh(aes(xmin = .lower, xmax = .upper), height = 0)+
    #scale_alpha_manual(values = c(0.2,1))+
    scale_shape_manual(values = c(1,16))+
    #scale_color_brewer(palette = "Set2")+
      scale_color_manual(values = c("firebrick4"), name = " ")+
    
    labs(title = site,
      x = NULL, y = NULL) +
    theme_bw() +
    guides(shape = FALSE)+
    theme(legend.title = element_blank(),
          panel.grid = element_blank(),
          panel.grid.major.y = element_line(color = alpha("firebrick4", 1/4), linetype = 3),
          axis.text.y  = element_text(hjust = 0),
          axis.ticks.y = element_blank(),
          legend.position = "bottom",
          legend.text=element_text(size=14),
          strip.background = element_blank(),
          strip.text = element_text(size = 14, face = "bold")
    )+
    facet_grid(~dependent, scales = "free_y", space='free')
}

HIplot<-IndivPlots(site = "Heron Island")
CIplot<-IndivPlots(site = "Cook Islands")
LIplot<-IndivPlots(site = "Lizzard Island")
OTIplot<-IndivPlots(site = "One Tree")
LHIplot<-IndivPlots(site = "Lord Howe Island")
NPplot<-IndivPlots(site = "Nusa Penida")

HIplot/CIplot/LIplot/OTIplot/LHIplot/NPplot
  ggsave("Output/CoefIndiv.pdf", width = 8, height = 20)
  
  ### try some residual plots
  ## plot the relationship between temp and aou, then take residuals of model to see if nutrients explain added variance
  
  SGDData_transect<-SGDData %>%
    filter(Time_Space == "Transect") %>%
    drop_na(salinity)
  
  SGDData_ts<-SGDData %>%
    filter(Time_Space == "Time Series") %>%
    drop_na(salinity)
  
  mod<-lm(aou_umol_l~temperature, data = SGDData_transect %>% filter(site == "Nusa Penida"))
  dat<-augment(mod, SGDData_transect%>% filter(site == "Nusa Penida"))

  dat %>%ggplot(aes(x = log(nox_u_m), y = .resid, color = DayNight))+
    geom_point()+
    geom_smooth(method = "lm")
  
  
 N_sum<- SGDData_ts%>% 
    group_by(site)%>%
    summarise(maxN = max(nox_u_m, na.rm = TRUE),
              meanN = mean(nox_u_m, na.rm = TRUE),
              covN = var(nox_u_m, na.rm = TRUE)/mean(nox_u_m, na.rm = TRUE))
    
dat2<-  SGDData_ts %>% 
  #  mutate(co_sst_std_log = scale(log(co_sst)))%>%
    group_by(site) %>% 
    nest() %>% 
    mutate(model = map(data, ~lm(p_h_in_std ~ aou_umol_l_std, data = .x) %>% 
                         tidy)) %>% 
    unnest(model) %>%
    filter(term == 'aou_umol_l_std') %>%
    select(site, estimate_aou = estimate, std.error_aou = std.error,p.value) %>% 
    mutate(sig = ifelse(p.value <= 0.05, "yes","no")) %>%
    left_join(N_sum) 

p1<-ggplot(dat2, aes(x = maxN, y = estimate_aou))+
  geom_point(size = 3, aes(color = sig))+
  geom_errorbar(aes(ymin = estimate_aou-std.error_aou, ymax = estimate_aou+std.error_aou), width = 0.1)+
  geom_smooth(method = "lm", formula = y~poly(x,2),se = FALSE, color = "black")+
  geom_label_repel(aes(label = site),nudge_x = 0.1)+
  # annotate(geom = "text", x = 2.7, y = 0.01, label = expression("net productivity did not effect pCO"[2]), size = 3, color = "blue")+
  # annotate(geom = "text", x = 2.7, y =- 0.2, label = expression("net productivity decreased pCO"[2]), size = 3, color = "blue")+
  # annotate(geom = "text", x = 2.7, y = 0.2, label = expression("net productivity increased pCO"[2]), size = 3, color = "blue")+
  scale_color_manual(values = c("grey","black"))+
  geom_hline(yintercept = 0, color = "grey", linetype = 2) +
  labs(y = expression("Standardized effect of productivity on pH"),
       x = expression(paste("Max Nitrate + Nitrite", " (",mu, "mol L"^-1,")"))
       )+
  theme_bw()+
  theme(legend.position = "null",
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14))


anova(lm(estimate_aou~maxN, data = dat2))


dat3<-  SGDData_ts %>% 
  #drop_na(co_sst_std) %>%
 # mutate(co_sst_std_log = scale(log(co_sst)))%>%
  group_by(site) %>% 
  nest() %>% 
  mutate(model = map(data, ~lm(p_h_in_std ~ rn_bq_m3_std, data = .x) %>% 
                       tidy)) %>% 
  unnest(model) %>%
  filter(term == 'rn_bq_m3_std') %>%
  select(site, estimate_sal = estimate, std.error_sal = std.error,p.value) %>% 
  mutate(sig = ifelse(p.value <= 0.05, "yes","no")) %>%
  left_join(N_sum)
  

p2<-ggplot(dat3, aes(x = maxN, y = estimate_sal))+
  geom_point(size = 3, aes(color = sig))+
  geom_errorbar(aes(ymin = estimate_sal-std.error_sal, ymax = estimate_sal+std.error_sal), width = 0.1)+
#  geom_smooth(method = "lm", se = FALSE, color = "black")+
  geom_label_repel(aes(label = site), nudge_x = 0.1)+
  scale_color_manual(values = c("grey","black"))+
  geom_hline(yintercept = 0, color = "grey", linetype = 2) +
  labs(y = expression("Standardized effect of SGD on pH"),
       x = expression(paste("Max Nitrate + Nitrite", " (",mu, "mol L"^-1,")"))
  )+
  # annotate(geom = "text", x = 2.7, y = 0.1, label = expression("SGD did not effect pCO"[2]), size = 3, color = "blue")+
  # annotate(geom = "text", x = 2.7, y = 1.5, label = expression("SGD decreased pCO"[2]), size = 3, color = "blue")+
  # annotate(geom = "text", x = 2.7, y = -1.5, label = expression("SGD increased pCO"[2]), size = 3, color = "blue")+
  
  theme_bw()+
  theme(legend.position = "null",
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14))


p1+p2+plot_annotation(tag_levels = "A")&theme(panel.grid = element_line(color = "white"))

ggsave(here("Output","maxN_effects.pdf"), width = 12, height = 5)

p_aou<-SGDData_ts %>%
  ggplot(aes(x = aou_umol_l, y = p_h_in))+
  geom_point()+
  geom_smooth(method = "lm",data = subset(SGDData_ts, site %in% c("Lord Howe Island","Heron Island","Nusa Penida","Cook Islands")))+
  facet_wrap(~site, scales = "free", ncol = 1)+
  theme_bw()+
  labs(x = expression(atop("Biological Productivity", paste("AOU (",mu, "mol L"^-1,")"))),
       y = "pH")+
  theme(legend.position = "null",
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14))
ggsave(here("Output","AOU_pH.pdf"), width = 8, height = 6)

p_radon<-SGDData_ts %>%
  ggplot(aes(x = rn_bq_m3, y = p_h_in))+
  geom_point()+
  geom_smooth(method = "lm", data = subset(SGDData_ts, site =="Nusa Penida"))+
  facet_wrap(~site, scales = "free", ncol = 1)+
  theme_bw()+
  labs(x = expression(paste("Radon", "(bq m"^3,")")),
       y = " ")+
  theme(legend.position = "null",
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14))

p_aou +p_radon

ggsave(here("Output","Radon_aou_pH.pdf"), width = 8, height = 16)

### transect
N_sum<- SGDData_transect %>% 
  group_by(site)%>%
  summarise(maxN = max(nox_u_m, na.rm = TRUE),
            meanN = mean(nox_u_m, na.rm = TRUE),
            covN = var(nox_u_m, na.rm = TRUE)/mean(nox_u_m, na.rm = TRUE))

dat2<-  SGDData_transect %>% 
  group_by(site) %>% 
  nest() %>% 
  mutate(model = map(data, ~lm(p_h_in_std ~ aou_umol_l_std, data = .x) %>% 
                       tidy)) %>% 
  unnest(model) %>%
  filter(term == 'aou_umol_l_std') %>%
  select(site, estimate_aou = estimate, std.error_aou = std.error,p.value) %>% 
  mutate(sig = ifelse(p.value <= 0.05, "yes","no")) %>%
  left_join(N_sum)


p1<-ggplot(dat2, aes(x = maxN, y = estimate_aou))+
  geom_point(size = 3, aes(color = sig))+
  geom_errorbar(aes(ymin = estimate_aou-std.error_aou, ymax = estimate_aou+std.error_aou), width = 0.1)+
  geom_smooth(method = "lm", formula = y~poly(x,2),se = FALSE, color = "black")+
  geom_label_repel(aes(label = site),nudge_x = 0.1)+
  # annotate(geom = "text", x = 2.7, y = 0.01, label = expression("net productivity did not effect pCO"[2]), size = 3, color = "blue")+
  # annotate(geom = "text", x = 2.7, y =- 0.2, label = expression("net productivity decreased pCO"[2]), size = 3, color = "blue")+
  # annotate(geom = "text", x = 2.7, y = 0.2, label = expression("net productivity increased pCO"[2]), size = 3, color = "blue")+
  scale_color_manual(values = c("grey","black"))+
  geom_hline(yintercept = 0, color = "grey", linetype = 2) +
  labs(y = expression("Standardized effect of productivity on pH"),
       x = expression(paste("Max Nitrate + Nitrite", " (",mu, "mol L"^-1,")"))
  )+
  theme_bw()+
  theme(legend.position = "null",
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14))

ggsave(here("Output","maxNeffects_transect.pdf"), width = 6, height = 5)

dat3<-  SGDData_transect %>% 
  #drop_na(co_sst_std) %>%
  # mutate(co_sst_std_log = scale(log(co_sst)))%>%
  group_by(site) %>% 
  nest() %>% 
  mutate(model = map(data, ~lm(p_h_in_std ~ rn_bq_m3_std, data = .x) %>% 
                       tidy)) %>% 
  unnest(model) %>%
  filter(term == 'rn_bq_m3_std') %>%
  select(site, estimate_sal = estimate, std.error_sal = std.error,p.value) %>% 
  mutate(sig = ifelse(p.value <= 0.05, "yes","no")) %>%
  left_join(N_sum)


p2<-ggplot(dat3, aes(x = maxN, y = estimate_sal))+
  geom_point(size = 3, aes(color = sig))+
  geom_errorbar(aes(ymin = estimate_sal-std.error_sal, ymax = estimate_sal+std.error_sal), width = 0.1)+
  #  geom_smooth(method = "lm", se = FALSE, color = "black")+
  geom_label_repel(aes(label = site), nudge_x = 0.1)+
  scale_color_manual(values = c("grey","black"))+
  geom_hline(yintercept = 0, color = "grey", linetype = 2) +
  labs(y = expression("Standardized effect of SGD on pH"),
       x = expression(paste("Max Nitrate + Nitrite", " (",mu, "mol L"^-1,")"))
  )+
  # annotate(geom = "text", x = 2.7, y = 0.1, label = expression("SGD did not effect pCO"[2]), size = 3, color = "blue")+
  # annotate(geom = "text", x = 2.7, y = 1.5, label = expression("SGD decreased pCO"[2]), size = 3, color = "blue")+
  # annotate(geom = "text", x = 2.7, y = -1.5, label = expression("SGD increased pCO"[2]), size = 3, color = "blue")+
  
  theme_bw()+
  theme(legend.position = "null",
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14))
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
SGDData <- read_csv(here("data","SGDData_times.csv"))

# standardize the data
SGDData <- SGDData %>%
  mutate_at(.vars = c("rn_bq_m3", "nox_u_m","aou_umol_l"  ,"excess_co2_umol_kg","temperature","p_h_in","advection_rate_cm_day", "water_level_m"), .funs = list(std = ~scale(.)))


### site level sums ####

SGD_sum<- SGDData %>%
  select(site, advection_rate_cm_day, din, nox_u_m, excess_co2_umol_kg)%>%
  group_by(site) %>%
  summarise_if(is.numeric, sum) 

P_no_onetree<-SGD_sum %>%
  filter(site != "One Tree") %>% # one tree is crazy
  ggplot(aes(x = din, y = excess_co2_umol_kg))+
  geom_point(size = 2)+
  geom_smooth(method = "lm", se = FALSE, color = "black") +
  geom_label_repel(aes(label = site))+
  geom_hline(yintercept = 0, color = "grey", linetype = 2)+
  labs(x = expression(paste("Cumulative DIN ",mu, "mol L"^-1 )),
       y = expression(paste("Cumulative Excess CO"[2]," ", mu, "mol kg"^-1 )))+
  theme_bw()


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

P_no_onetree + P_onetree +plot_annotation(tag_levels = "A")+
  ggsave(filename = "Output/cumulativeN_CO2.pdf", width = 8, height = 4)


### site level averages ####

SGD_mean<- SGDData %>%
  group_by(site) %>%
  summarise_if(is.numeric, .fun = mean) 

SGD_mean %>%
  filter(site != "One Tree") %>% 
  ggplot(aes(x = din, y = p_h_in))+
  geom_point(size = 2)+
   geom_smooth(method = "lm", se = FALSE, color = "black") +
  geom_label_repel(aes(label = site))+
 # geom_hline(yintercept = 0, color = "grey", linetype = 2)+
  labs(x = "Mean DIN umol/L",
       y = "Mean pH")+
  theme_bw()

### Bayes SEM analysis ###
# SGD ~ depth
SGDmod<-bf(rn_bq_m3_std~water_level_m_std)
# N ~ SGD*DayNight
Nmod<-bf(nox_u_m_std~DayNight*rn_bq_m3_std)
# Temperature ~ SGD*DayNight
Tempmod<-bf(temperature_std ~ DayNight*rn_bq_m3_std)
# pH ~ DayNight*(AOU + SGD)
pHmod<-bf(p_h_in_std~DayNight*(aou_umol_l_std+rn_bq_m3_std))
# AOU ~ DayNight*(N+Temperature)
AOUmod<-bf(aou_umol_l_std~DayNight*(temperature_std+nox_u_m_std))
# excess CO2 ~ DayNight*(SGD +AOU)
CO2mod<-bf(excess_co2_umol_kg_std ~DayNight*(rn_bq_m3_std+aou_umol_l_std))
## Later
# Delta TA ~ pH + temp

# run the SEM for Heron Island
heron_fit_brms <- brm(SGDmod+
                    Nmod+
                    Tempmod+ 
                    pHmod+
                    AOUmod +
                    CO2mod+
                    set_rescor(FALSE),
                  data=SGDData[SGDData$site == "Heron Island",]
                  ,cores=4, chains = 3)
# calculate LOO (leave one out) diagnostics
SGD_loo_heron<-loo(heron_fit_brms, reloo = TRUE) # looks good!

p1<-pp_check(heron_fit_brms, resp="rnbqm3std") +
  scale_color_manual(values=c("red", "black"))+
  ggtitle("Radon")
p2<-pp_check(heron_fit_brms, resp="noxumstd") +
  scale_color_manual(values=c("red", "black"))+
  ggtitle("N+N")
p3<-pp_check(heron_fit_brms, resp="temperaturestd") +
  scale_color_manual(values=c("red", "black"))+
  ggtitle("Temperature")
p4<-pp_check(heron_fit_brms, resp="phinstd") +
  scale_color_manual(values=c("red", "black"))+
  ggtitle("pH")
p5<-pp_check(heron_fit_brms, resp="aouumollstd") +
  scale_color_manual(values=c("red", "black"))+
  ggtitle("AOU")
p6<-pp_check(heron_fit_brms, resp="excessco2umolkgstd") +
  scale_color_manual(values=c("red", "black"))+
  ggtitle("Excess CO2")

p1+p2+p3+p4+p5+p6+plot_layout(guides = "collect") +
  plot_annotation(title = 'Heron Island Posterior Predictive Checks', tag_levels = "A")+
  ggsave("Output/Posteriorchecks_HeronIsland.pdf", width = 5, height = 5)

# plot the results
# Model 1
# SGD ~ depth

R<-conditional_effects(heron_fit_brms, "water_level_m_std", resp = "rnbqm3std", method = "predict", resolution = 1000)
R1<-R$rnbqm3std.rnbqm3std_water_level_m_std %>% # back transform the scaled effects for the plot
  mutate(estimate = estimate__*attr(SGDData$rn_bq_m3_std,"scaled:scale")+attr(SGDData$rn_bq_m3_std,"scaled:center"),
         lower = lower__*attr(SGDData$rn_bq_m3_std,"scaled:scale")+attr(SGDData$rn_bq_m3_std,"scaled:center"),
         upper = upper__*attr(SGDData$rn_bq_m3_std,"scaled:scale")+attr(SGDData$rn_bq_m3_std,"scaled:center")) %>%
  mutate(WaterLevel = water_level_m_std*attr (SGDData$water_level_m_std,"scaled:scale")+attr(SGDData$water_level_m_std,"scaled:center")
  )%>%
  ggplot()+ # back trasform the log transformed data for better visual
  geom_line(aes(x = WaterLevel, y = estimate), lwd = 1, color = 'grey')+
  geom_ribbon(aes(x = WaterLevel,ymin=lower, ymax=upper), linetype=1.5, alpha=0.3, fill = "grey")+
  geom_point(data = SGDData[SGDData$site=='Heron Island',], aes(x = water_level_m, y = rn_bq_m3, color = DayNight)) +
  xlab("Water level (m)")+
  ylab(expression(atop("Radon", paste("(bq m"^3,")"))))+
  ggtitle("Model 1")+
 # scale_x_continuous(breaks = c(0.2,1,5,10,25))+
#  scale_y_continuous(breaks = c(0,0.1,1,5,10,30))+
  theme_minimal()+
  theme(plot.title = element_text(hjust = 0.5, size = 14), legend.position = 'bottom')

## Model 2
# N ~ SGD*DayNight
R<-conditional_effects(heron_fit_brms, "rn_bq_m3_std:DayNight", resp = "noxumstd", method = "predict", resolution = 1000)
R2<-R$`noxumstd.noxumstd_rn_bq_m3_std:DayNight` %>% # back transform the scaled effects for the plot
  mutate(estimate = estimate__*attr(SGDData$nox_u_m_std,"scaled:scale")+attr(SGDData$nox_u_m_std,"scaled:center"),
         lower = lower__*attr(SGDData$nox_u_m_std,"scaled:scale")+attr(SGDData$nox_u_m_std,"scaled:center"),
         upper = upper__*attr(SGDData$nox_u_m_std,"scaled:scale")+attr(SGDData$nox_u_m_std,"scaled:center")) %>%
  mutate(radon = rn_bq_m3_std*attr (SGDData$rn_bq_m3_std,"scaled:scale")+attr(SGDData$rn_bq_m3_std,"scaled:center")
  )%>%
  ggplot()+ # back trasform the log transformed data for better visual
  geom_line(aes(x = radon, y = estimate, color = DayNight), lwd = 1)+
  geom_ribbon(aes(x = radon,ymin=lower, ymax=upper, fill = DayNight), linetype=1.5, alpha=0.3)+
  geom_point(data = SGDData[SGDData$site=='Heron Island',], aes(x = rn_bq_m3, y = nox_u_m, color = DayNight)) +
  ylab(expression(atop("Nitrate + Nitrite", paste("(",mu, "mol L"^-1,")"))))+
  xlab(expression(atop("Radon", paste("(bq m"^3,")"))))+
  ggtitle("Model 2")+
  # scale_x_continuous(breaks = c(0.2,1,5,10,25))+
  #  scale_y_continuous(breaks = c(0,0.1,1,5,10,30))+
  theme_minimal()+
  theme(plot.title = element_text(hjust = 0.5, size = 14), legend.position = 'bottom')

  # Model 3
  # Temperature ~ SGD*DayNight
  
R<-conditional_effects(heron_fit_brms, "rn_bq_m3_std:DayNight", resp = "temperaturestd", method = "predict", resolution = 1000)
  R3<-R$`temperaturestd.temperaturestd_rn_bq_m3_std:DayNight` %>% # back transform the scaled effects for the plot
    mutate(estimate = estimate__*attr(SGDData$temperature_std,"scaled:scale")+attr(SGDData$temperature_std,"scaled:center"),
           lower = lower__*attr(SGDData$temperature_std,"scaled:scale")+attr(SGDData$temperature_std,"scaled:center"),
           upper = upper__*attr(SGDData$temperature_std,"scaled:scale")+attr(SGDData$temperature_std,"scaled:center")) %>%
    mutate(radon = rn_bq_m3_std*attr (SGDData$rn_bq_m3_std,"scaled:scale")+attr(SGDData$rn_bq_m3_std,"scaled:center")
    )%>%
    ggplot()+ # back trasform the log transformed data for better visual
    geom_line(aes(x = radon, y = estimate, color = DayNight), lwd = 1)+
    geom_ribbon(aes(x = radon,ymin=lower, ymax=upper, fill = DayNight), linetype=1.5, alpha=0.3)+
    geom_point(data = SGDData[SGDData$site=='Heron Island',], aes(x = rn_bq_m3, y = temperature, color = DayNight)) +
    ylab(expression(atop("Temperature",paste("(", degree, "C)"))))+
  xlab(expression(atop("Radon", paste("(bq m"^3,")"))))+
    ggtitle("Model 3")+
    # scale_x_continuous(breaks = c(0.2,1,5,10,25))+
    #  scale_y_continuous(breaks = c(0,0.1,1,5,10,30))+
    theme_minimal()+
    theme(plot.title = element_text(hjust = 0.5, size = 14), legend.position = 'bottom')
  
  # Model 4a
  # pH ~ DayNight*(AOU + SGD)
  R<-conditional_effects(heron_fit_brms, "rn_bq_m3_std:DayNight", resp = "phinstd", method = "predict", resolution = 1000)
  R4<-R$`phinstd.phinstd_rn_bq_m3_std:DayNight` %>% # back transform the scaled effects for the plot
    mutate(estimate = estimate__*attr(SGDData$p_h_in_std,"scaled:scale")+attr(SGDData$p_h_in_std,"scaled:center"),
           lower = lower__*attr(SGDData$p_h_in_std,"scaled:scale")+attr(SGDData$p_h_in_std,"scaled:center"),
           upper = upper__*attr(SGDData$p_h_in_std,"scaled:scale")+attr(SGDData$p_h_in_std,"scaled:center")) %>%
    mutate(radon = rn_bq_m3_std*attr (SGDData$rn_bq_m3_std,"scaled:scale")+attr(SGDData$rn_bq_m3_std,"scaled:center")
    )%>%
    ggplot()+ # back trasform the log transformed data for better visual
    geom_line(aes(x = radon, y = estimate, color = DayNight), lwd = 1)+
    geom_ribbon(aes(x = radon,ymin=lower, ymax=upper, fill = DayNight), linetype=1.5, alpha=0.3)+
    geom_point(data = SGDData[SGDData$site=='Heron Island',], aes(x = rn_bq_m3, y = p_h_in, color = DayNight)) +
    ylab("pH")+
    xlab(expression(atop("Radon", paste("(bq m"^3,")"))))+
    ggtitle("Model 4")+
    # scale_x_continuous(breaks = c(0.2,1,5,10,25))+
    #  scale_y_continuous(breaks = c(0,0.1,1,5,10,30))+
    theme_minimal()+
    theme(plot.title = element_text(hjust = 0.5, size = 14), legend.position = 'bottom')
  
  # Model 4b
  # pH ~ DayNight*(AOU + SGD)
  R<-conditional_effects(heron_fit_brms, "aou_umol_l_std:DayNight", resp = "phinstd", method = "predict", resolution = 1000)
  R4b<-R$`phinstd.phinstd_aou_umol_l_std:DayNight` %>% # back transform the scaled effects for the plot
    mutate(estimate = estimate__*attr(SGDData$p_h_in_std,"scaled:scale")+attr(SGDData$p_h_in_std,"scaled:center"),
           lower = lower__*attr(SGDData$p_h_in_std,"scaled:scale")+attr(SGDData$p_h_in_std,"scaled:center"),
           upper = upper__*attr(SGDData$p_h_in_std,"scaled:scale")+attr(SGDData$p_h_in_std,"scaled:center")) %>%
    mutate(aou = aou_umol_l_std*attr (SGDData$aou_umol_l_std,"scaled:scale")+attr(SGDData$aou_umol_l_std,"scaled:center")
    )%>%
    ggplot()+ # back trasform the log transformed data for better visual
    geom_line(aes(x = aou, y = estimate, color = DayNight), lwd = 1)+
    geom_ribbon(aes(x = aou,ymin=lower, ymax=upper, fill = DayNight), linetype=1.5, alpha=0.3)+
    geom_point(data = SGDData[SGDData$site=='Heron Island',], aes(x = aou_umol_l, y = p_h_in, color = DayNight)) +
    ylab("pH")+
    xlab(expression(atop("AOU", paste("(",mu, "mol L"^-1,")"))))+
    ggtitle("Model 4")+
    # scale_x_continuous(breaks = c(0.2,1,5,10,25))+
    #  scale_y_continuous(breaks = c(0,0.1,1,5,10,30))+
    theme_minimal()+
    theme(plot.title = element_text(hjust = 0.5, size = 14), legend.position = 'bottom')
  
  # Model 5
  # AOU ~ DayNight*(N+Temperature)
  
  R<-conditional_effects(heron_fit_brms, "temperature_std:DayNight", resp = "aouumollstd", method = "predict", resolution = 1000)
  R5<-R$`aouumollstd.aouumollstd_temperature_std:DayNight` %>% # back transform the scaled effects for the plot
    mutate(estimate = estimate__*attr(SGDData$aou_umol_l_std,"scaled:scale")+attr(SGDData$aou_umol_l_std,"scaled:center"),
           lower = lower__*attr(SGDData$aou_umol_l_std,"scaled:scale")+attr(SGDData$aou_umol_l_std,"scaled:center"),
           upper = upper__*attr(SGDData$aou_umol_l_std,"scaled:scale")+attr(SGDData$p_h_in_std,"scaled:center")) %>%
    mutate(temperature = temperature_std*attr (SGDData$temperature_std,"scaled:scale")+attr(SGDData$temperature_std,"scaled:center")
    )%>%
    ggplot()+ # back trasform the log transformed data for better visual
    geom_line(aes(x = temperature, y = estimate, color = DayNight), lwd = 1)+
    geom_ribbon(aes(x = temperature,ymin=lower, ymax=upper, fill = DayNight), linetype=1.5, alpha=0.3)+
    geom_point(data = SGDData[SGDData$site=='Heron Island',], aes(x = temperature, y = aou_umol_l, color = DayNight)) +
    xlab(expression(atop("Temperature",paste("(", degree, "C)"))))+
    ylab(expression(atop("AOU", paste("(",mu, "mol L"^-1,")"))))+
    ggtitle("Model 5")+
    # scale_x_continuous(breaks = c(0.2,1,5,10,25))+
    #  scale_y_continuous(breaks = c(0,0.1,1,5,10,30))+
    theme_minimal()+
    theme(plot.title = element_text(hjust = 0.5, size = 14), legend.position = 'bottom')
  
  # Model 5b
  # AOU ~ DayNight*(N+Temperature)
  
  R<-conditional_effects(heron_fit_brms, "nox_u_m_std:DayNight", resp = "aouumollstd", method = "predict", resolution = 1000)
  R5b<-R$`aouumollstd.aouumollstd_nox_u_m_std:DayNight` %>% # back transform the scaled effects for the plot
    mutate(estimate = estimate__*attr(SGDData$aou_umol_l_std,"scaled:scale")+attr(SGDData$aou_umol_l_std,"scaled:center"),
           lower = lower__*attr(SGDData$aou_umol_l_std,"scaled:scale")+attr(SGDData$aou_umol_l_std,"scaled:center"),
           upper = upper__*attr(SGDData$aou_umol_l_std,"scaled:scale")+attr(SGDData$aou_umol_l_std,"scaled:center")) %>%
    mutate(NN = nox_u_m_std*attr (SGDData$nox_u_m_std,"scaled:scale")+attr(SGDData$nox_u_m_std,"scaled:center")
    )%>%
    ggplot()+ # back trasform the log transformed data for better visual
    geom_line(aes(x = NN, y = estimate, color = DayNight), lwd = 1)+
    geom_ribbon(aes(x = NN,ymin=lower, ymax=upper, fill = DayNight), linetype=1.5, alpha=0.3)+
    geom_point(data = SGDData[SGDData$site=='Heron Island',], aes(x = nox_u_m, y = aou_umol_l, color = DayNight)) +
    xlab(expression(atop("Nitrate + Nitrite", paste("(",mu, "mol L"^-1,")"))))+
    ylab(expression(atop("AOU", paste("(",mu, "mol L"^-1,")"))))+
    ggtitle("Model 5")+
    # scale_x_continuous(breaks = c(0.2,1,5,10,25))+
    #  scale_y_continuous(breaks = c(0,0.1,1,5,10,30))+
    theme_minimal()+
    theme(plot.title = element_text(hjust = 0.5, size = 14), legend.position = 'bottom')
  
# Model 6
# excess CO2 ~ DayNight*(SGD +AOU)
  
R<-conditional_effects(heron_fit_brms, "rn_bq_m3_std:DayNight", resp = "excessco2umolkgstd", method = "predict", resolution = 1000)
  R6<-R$`excessco2umolkgstd.excessco2umolkgstd_rn_bq_m3_std:DayNight` %>% # back transform the scaled effects for the plot
    mutate(estimate = estimate__*attr(SGDData$excess_co2_umol_kg_std,"scaled:scale")+attr(SGDData$excess_co2_umol_kg_std,"scaled:center"),
           lower = lower__*attr(SGDData$excess_co2_umol_kg_std,"scaled:scale")+attr(SGDData$excess_co2_umol_kg_std,"scaled:center"),
           upper = upper__*attr(SGDData$excess_co2_umol_kg_std,"scaled:scale")+attr(SGDData$excess_co2_umol_kg_std,"scaled:center")) %>%
    mutate(radon = rn_bq_m3_std*attr(SGDData$rn_bq_m3_std,"scaled:scale")+attr(SGDData$rn_bq_m3_std,"scaled:center")
    )%>%
    ggplot()+ # back trasform the log transformed data for better visual
    geom_line(aes(x = radon, y = estimate, color = DayNight), lwd = 1)+
    geom_ribbon(aes(x = radon,ymin=lower, ymax=upper, fill = DayNight), linetype=1.5, alpha=0.3)+
    geom_point(data = SGDData[SGDData$site=='Heron Island',], aes(x = rn_bq_m3, y = excess_co2_umol_kg, color = DayNight)) +
    xlab(expression(atop("Radon", paste("(bq m"^3,")"))))+
    ylab(expression(atop("Excess CO"[2], paste("(",mu, "mol kg"^-1,")"))))+
    ggtitle("Model 6")+
    # scale_x_continuous(breaks = c(0.2,1,5,10,25))+
    #  scale_y_continuous(breaks = c(0,0.1,1,5,10,30))+
    theme_minimal()+
    theme(plot.title = element_text(hjust = 0.5, size = 14), legend.position = 'bottom')

  # Model 6b
  # excess CO2 ~ DayNight*(SGD +AOU)
  
  R<-conditional_effects(heron_fit_brms, "aou_umol_l_std:DayNight", resp = "excessco2umolkgstd", method = "predict", resolution = 1000)
  R6b<-R$`excessco2umolkgstd.excessco2umolkgstd_aou_umol_l_std:DayNight` %>% # back transform the scaled effects for the plot
    mutate(estimate = estimate__*attr(SGDData$excess_co2_umol_kg_std,"scaled:scale")+attr(SGDData$excess_co2_umol_kg_std,"scaled:center"),
           lower = lower__*attr(SGDData$excess_co2_umol_kg_std,"scaled:scale")+attr(SGDData$excess_co2_umol_kg_std,"scaled:center"),
           upper = upper__*attr(SGDData$excess_co2_umol_kg_std,"scaled:scale")+attr(SGDData$excess_co2_umol_kg_std,"scaled:center")) %>%
    mutate(aou = aou_umol_l_std*attr(SGDData$aou_umol_l_std,"scaled:scale")+attr(SGDData$aou_umol_l_std,"scaled:center")
    )%>%
    ggplot()+ # back trasform the log transformed data for better visual
    geom_line(aes(x = aou, y = estimate, color = DayNight), lwd = 1)+
    geom_ribbon(aes(x = aou,ymin=lower, ymax=upper, fill = DayNight), linetype=1.5, alpha=0.3)+
    geom_point(data = SGDData[SGDData$site=='Heron Island',], aes(x = aou_umol_l, y = excess_co2_umol_kg, color = DayNight)) +
    xlab(expression(atop("AOU", paste("(",mu, "mol L"^-1,")"))))+
    ylab(expression(atop("Excess CO"[2], paste("(",mu, "mol kg"^-1,")"))))+
    ggtitle("Model 6")+
    # scale_x_continuous(breaks = c(0.2,1,5,10,25))+
    #  scale_y_continuous(breaks = c(0,0.1,1,5,10,30))+
    theme_minimal()+
    theme(plot.title = element_text(hjust = 0.5, size = 14), legend.position = 'bottom')
  
  ## bring them all together in patchwork
  R<-(R1|R2)/(R4|R4b)/(R5|R5b)/(R6|R6b)/(R3)+plot_layout(guides = "collect")+
    plot_annotation(tag_levels = "A")&
    theme(axis.title.x = element_text(size = 14),
          axis.title.y = element_text(size = 14),
          axis.text = element_text(size = 11))
  
  ggsave("Output/marginaleffects_heronisland.pdf",R, width = 10, height = 18, useDingbats = FALSE)
  
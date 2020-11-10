library(tidyhydat)
library(weathercan)
library(lubridate)
library(dplyr)
library(zoo)
library(ggplot2)
library(leaflet) 

# water year function
water_year <- function(date, day_start) {
  year_day <- yday(date)
  years <- year(date)
  water_year <- vector(length=length(date))
  for(i in 1:length(date)) {
    ifelse(year_day[i] >= yday(day_start), 
           water_year[i] <- years[i], # if true
           water_year[i] <- years[i] - 1 # else
    )
  }
  water_year
}

# water year day function
water_year_day <- function(date, day_start) {
  year_day <- yday(date)
  water_year_day <- vector(length=length(date))
  for(i in 1:length(date)) {
    ifelse(year_day[i] >= yday(day_start), 
           water_year_day[i] <- year_day[i] - yday(day_start) + 1 , # if ttrue
           water_year_day[i] <- year_day[i] - yday(day_start) + 366 # else
    )
  }
  water_year_day
}

# Weather at Blind channel
stations_search("Blind")
stns <- stations_search(coords= c(50.883868, -125.608020), dist=100, interval="day")

ggplot(stns, y=station_name, x=seq(min(start), max(end),1)) + 
  geom_linerange(aes(y=station_name, xmin=start, xmax=end))  

leaflet() %>%
  addTiles() %>%
  addMarkers(lng=stns$lon, lat = stns$lat, popup=stns$station_name)


wt <- weather_dl(station_ids = 140, interval = "day")
wt$yday<- yday(wt$date)
wt$water_year <- water_year(wt$date, day_start="2020-09-01")
wt$water_yday <- water_year_day(wt$date, day_start="2020-09-01")

test <- wt[,c("date","yday", "water_year", "water_yday")]
test

wt$roll_rain <- rollsum(wt$total_rain, k=2 , fill=NA, align="right")
head(wt)
names(wt)
goodyrs <- c(1957)
badyrs <- c(1958)
ggplot(wt, aes(y=roll_rain, x=water_yday, group=water_year)) + 
  geom_line(alpha=0.1) +
  geom_line(data=wt[wt$water_year %in% goodyrs,], aes(y=roll_rain, x=water_yday), colour="dodgerblue", size=1.1) +
  geom_line(data=wt[wt$water_year %in% badyrs,], aes(y=roll_rain, x=water_yday), colour="orange", size=1.1) +
  #facet_wrap(~water_year) +
  theme_classic()



# Hydrometric data- Klinaklini River (head of KNight inlet)

kli <- hy_daily_flows(station_number = "08GE002")
kli$yday <- yday(kli$Date)
kli$year <- year(kli$Date)
kli$water_year <- water_year(kli$Date, day_start="2020-09-01")
kli$water_yday <- water_year_day(kli$Date, day_start="2020-09-01")

head(kli)
range(kli$Date)

goodyrs <- c(1994)
badyrs <- c(1991)
ggplot(kli, aes(y=Value, x=water_yday, group=water_year)) + 
  geom_line(alpha=0.5) +
  geom_line(data=kli[kli$water_year %in% goodyrs,], aes(y=Value, x=water_yday), colour="dodgerblue", size=1.1) +
  geom_line(data=kli[kli$water_year %in% badyrs,], aes(y=Value, x=water_yday), colour="orange", size=1.1) +
  theme_classic()

# merge flow and weather data
cd <- merge(wt, kli, by.x=c("date", "yday", "water_yday", "water_year", "yday"), by.y=c("Date", "yday", "water_yday", "water_year", "yday"), all=TRUE)

ggplot(cd, aes(y=Value, x=roll_rain, colour=water_year)) + 
  geom_point() +
  facet_wrap(~ month) +
  stat_smooth(method="lm") +
  theme_classic()

ggplot(cd, aes(y=Value, x=, colour=water_year)) + 
  geom_point() +
  facet_wrap(~ month) +
  stat_smooth(method="lm") +
  theme_classic()

#fs <- kli %>% filter(yday<120) %>% group_by(year) %>% summarise(max_jan_apr_flow = max(Value, na.rm=TRUE))
fs <- kli %>% filter(yday>300) %>% group_by(year) %>% summarise(max_nov_dec_flow = max(Value, na.rm=TRUE))

dat$RS <- dat$recruits/dat$spawners
datf <- merge(dat[dat$CU=="3 - Upper Knight",], fs, by="year", all.x=TRUE)

plot(log(datf$RS) ~ datf$max_nov_dec_flow)
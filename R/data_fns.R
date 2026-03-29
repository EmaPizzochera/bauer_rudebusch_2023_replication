loadShortRate <- function() {
    require(lubridate)
    ## real short rate: 1-year CM Treasury yield - one-year survey inflation
    data <- setNames(read.csv("data/GS1.csv", na.strings="."), c("date", "y1")) %>% # last month of the year
        mutate(date = as.Date(date),
               year = year(date)) %>%
        select(-date)
    ## Livingston median CPI growth
    inf <- read.csv("data/MedianGrowthRate_CPI.csv", sep=";", dec=",", na.string="#N/A") %>%
        rename(date = Date, pi = G_BP_To_12M) %>%
        select(date, pi) %>%
        mutate(date = as.Date(date),
               year = year(date)) %>%
        filter(month(date) == 12) %>%
        select(-date)
    inner_join(data, inf) %>%
        mutate(rr = y1 - pi) %>%
        filter(year <= 2019)
}

loadLongRate <- function() {
    require(lubridate)
    ## real long rate: 10-year CM Treasury yield - PTR
    data <- setNames(read.csv("data/GS10.csv", na.strings="."), c("date", "y"))
    data <- data %>%
        mutate(date = as.Date(date),
               year = year(date)) %>%
        select(year, y)
    ## PTR from PRB/US
    ptr <- read.csv("data/HISTDATA_FRBUS.TXT") %>% as_tibble %>%
        select(OBS, PTR) %>%
        rename(pi = PTR) %>%
        mutate(year = as.numeric(substr(OBS,1,4))) %>%
        group_by(year) %>%
        slice_tail(n=1) %>%
        select(year, pi) %>% ungroup
    pi <- ptr
    data <- inner_join(data, pi) %>%
        mutate(rr = y - pi) %>%
        filter(year <= 2019)
}

## other r-star estimates (all are smoothed estimates from state space models)

loadLW <- function() {
    ## Laubach-Williams
    lw <- read.csv("data/rstar_lw.csv" , na.strings = "#N/A")
    lw$Date <- as.Date(lw$Date, format="%m/%d/%Y")
    lw$yyyymm <- as.numeric(format(lw$Date, "%Y%m"))
    names(lw)[2:3] <- c("rstar.lw", "rstar.lw.sm")
    lw[c("yyyymm", "rstar.lw", "rstar.lw.sm")]
}

Matlab2yyyymm <- function(date) {
    ## transform Matlab's dates (Kiley/JM), 2012.00 = 2012:Q1, 2012.75 = 2012:Q4
    year <- floor(date)
    x <- date - year
    quarter <- match(x, x[1:4])
    year*100 + quarter*3
}

loadKiley <- function() {
    ## Kiley
    kiley <- read.csv("data/rstar_kiley.csv", header=FALSE)
    names(kiley) <- c("date", "rstar.kiley.sm", "rstar.kiley")
    kiley$yyyymm <- Matlab2yyyymm(kiley$date)
    kiley <- kiley[c("yyyymm", "rstar.kiley", "rstar.kiley.sm")]
    kiley
}

loadJM <- function() {
    ## Johannsen and Mertens
    jm <- read.csv("data/rstar_jm.csv")
    jm$yyyymm <- Matlab2yyyymm(jm$Date)
    names(jm)[1] <- "rstar.jm.sm"
    jm <- jm[c("yyyymm", "rstar.jm.sm")]
}

loadDN <- function() {
    ## Del Negro et al
    dn <- read.csv("data/rstar_delnegro.csv")
    names(dn) <- c("yyyymm", "rstar.dn.sm", "rstar.dn.lb", "rstar.dn.ub")
    dn[c("yyyymm", "rstar.dn.sm")]
}

library(tidyverse)
# code to read in example graph sizes
# from a bunch of studies

study1 <- bind_rows(
  qin2020 = tibble::tribble(~u, ~v, ~e,
                            120,500,11500,
                            120,1000,23000,
                            120,2000,46000,
                            120,4000,92000,
                            120,8000,184000,
                            120,8124,186852),
  
  lu2020 = tibble::tribble(~u, ~v, ~e,
                           829, 551, 1476,
                           16726, 22015, 58595),
  
  abidi2020 = tibble::tribble(~u, ~v, ~e,
                              20,44,99,
                              254,868,1255,
                              899,1421,33720,
                              10106,16730,50632,
                              4009,20537,95580,
                              4009,11610,95580,
                              16528,24129,95580,
                              1408,26546,193618,
                              2884,30997,201727),
  
  makino = tibble::tribble(~u, ~v, ~e,
                           22677,18484,247003,
                           33347,32757,233450,
                           20433,4297,127713),
  
  shaham2016 = tibble::tribble(~u, ~v, ~e,
                               943,1682,100000,
                               94238,30087,293360,
                               172091,53407,293697,
                               56519,120867,440237,
                               21607,94756,549210,
                               10764,163008,901416,
                               6040,3706,1000209,
                               32583,134942,1164576,
                               17122,82035,2298816,
                               69878,10677,10000054),
  lyu2020 = tibble::tribble(~u, ~v, ~e,
                            89355, 46213,144340,
                            124325,94238,293360,
                            56519,120867,440237,
                            105278,340523,1149739,
                            545195,96678,1301942,
                            901130,34461,1366466,
                            127823,383640,1470404,
                            64415,87678,3232134,
                            2036440,1853493,3795796,
                            499610,395979,8545307,
                            1425813,4000150,8649016,
                            3201203,7489073,112307385,
                            27665730,12756244,140613762,
                            78582023,23827661,184265522,
                            141839807,65589796,1307950593,
                            272227190,75350951,1319706942),
  
  sozdinler2018 = tribble(~u, ~v, ~e,
                          961,2096,3634,
                          4026,96,40197,
                          500,19,1465,
                          1693,10,3114,
                          1024,7,1674,
                          1747,12,3746,
                          21167,10,11115,
                          4582,82,22279,
                          3891,55,15747,
                          34,568,5165,
                          734,69,3810,
                          943,1682,100000),
  
  .id="study")

study2 <- bind_rows(
  mukherjee2014 = tibble(v=c(8114,125551,334863,262111,251226,175944,3928,49142,150615,5021,50000,60000,70000,80000,90000,100000,250000,500000,150000),
                         e=c(26013,168087,925872,1234877,419573,252075,35397,171421,300398,17409,
                             275659,330015,393410,448289,526943,600038,1562707,3751823,1999002)),
  
  das2019 = tibble(v=c(225486,19428,124325, 1199919, 641873, 445801),
                   e=c(293697, 96662, 293360, 3782463, 1301942, 1149739)),
  .id='study')

ratio <- study1 |> summarise(rat = median(pmin(u,v)/(u+v))) |> pull(rat)

study2 |> mutate(u = round(v*ratio), v = v - u)

study1 |> bind_rows(study2) |>
  write_csv(here::here("data/meb_example_graphs.csv"))

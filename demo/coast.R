require(PrecipStat)

data(mu.qp)
mu.qp <- dist2coast(mu.qp)
save(file="~/R/qqplotter/data/mu.qp.rda",mu.qp)

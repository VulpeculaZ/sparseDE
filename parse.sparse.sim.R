load("pen2d01.RData")
res.pen <- sim.res
rm(list = c("sim.res"))
load("lars2d01.RData")
res.lars <- sim.res
rm(list = c("sim.res"))
load("nnls2d01.RData")
res.nnls <- sim.res
rm(list = c("sim.res"))
load("adlars2d02.RData")
res.adlars <- sim.res
rm(list = c("sim.res"))

## load("res.all.RData")

pars.true <- 0.5
beta.true <- rep(0,16)
beta.true[c(4,7)] <- 2
div.nls <- fdp.nls <- fnp.nls <- div.sel <- fdp.sel <- fnp.sel <- 0
div.al <- fdp.al <- fnp.al <- 0
div.lars <- fdp.lars <- fnp.lars <- 0
mpars.nnls <- mpars.lars <- mpars.al <- mpars.sel <- 0
mb.nnls <- mb.lars <- mb.al <- mb.sel <- rep(0, 16)
f.nnls <- f.lars <- f.pen <- f.al <- f.oracle <- 0

for(i in 1:length(res.pen)){
    beta.nnls <- res.nnls[[i]]$res$beta
    pars.nnls <- res.nnls[[i]]$res$pars
    beta.lars <- res.lars[[i]]$beta
    pars.lars <- res.lars[[i]]$pars
    beta.al <- res.adlars[[i]]$beta
    pars.al <- res.adlars[[i]]$pars
    beta.sel <- res.pen[[i]]$beta.pen
    pars.sel <- res.pen[[i]]$pars.pen
    f.nnls <- res.nnls[[i]]$res$f + f.nnls
    f.lars <- res.lars[[i]]$f + f.lars
    f.pen <- res.pen[[i]]$f + f.pen
    f.al <- res.adlars[[i]]$f + f.al
    ## f.oracle <- sim.res[[i]]$res$f + f.oracle
    mpars.nnls <- pars.nnls + mpars.nnls
    mpars.lars <- pars.lars + mpars.lars
    mpars.al <- pars.al+mpars.al
    mpars.sel <- pars.sel + mpars.sel
    mb.nnls <- mb.nnls + beta.nnls
    mb.lars <- mb.lars + beta.lars
    mb.al <- mb.al + beta.al
    mb.sel <- mb.sel + beta.sel
    div.nls <- div.nls + sum(beta.nnls)
    div.sel <- div.sel + sum(beta.sel)
    div.al <- div.al + sum(beta.al)
    div.lars <- div.sel + sum(beta.lars)
    fdp.nls <- fdp.nls + sum(beta.nnls[-c(4,7)] != 0) / sum(beta.true == 0)
    fnp.nls <- fnp.nls + sum(beta.nnls[c(4,7)] == 0) / 2
    fdp.sel <- fdp.sel + sum(beta.sel[-c(4,7)] != 0) / sum(beta.true == 0)
    fnp.sel <- fnp.sel + sum(beta.sel[c(4,7)] == 0) / 2
    fdp.al <- fdp.al + sum(beta.al[-c(4,7)] != 0) / sum(beta.true == 0)
    fnp.al <- fnp.al + sum(beta.al[c(4,7)] == 0) / 2
    fdp.lars <- fdp.sel + sum(beta.lars[-c(4,7)] != 0) / sum(beta.true == 0)
    fnp.lars <- fnp.sel + sum(beta.lars[c(4,7)] == 0) / 2
}

"nnls"
div.nls / length(res.pen)
fdp.nls / length(res.pen)
fnp.nls / length(res.pen)
mpars.nnls /100
## > [1] 3.891313
## > [1] 0.2314286
## > [1] 0.585
## >     gamma
## 0.4793769

"pen"
div.sel / length(res.pen)
fdp.sel / length(res.pen)
fnp.sel  / length(res.pen)
mpars.sel / 100
## > [1] 3.958938
## > [1] 0.1878571
## > [1] 0.62
## >     gamma
## 0.4930961

"a l"
div.al / length(res.pen)
fdp.al / length(res.pen)
fnp.al  / length(res.pen)
mpars.al / 100
## > [1] 3.926803
## > [1] 0.1135714
## > [1] 0.785
## >     gamma
## 0.4793769

"lars"
div.lars / length(res.pen)
fdp.lars / length(res.pen)
fnp.lars  / length(res.pen)
mpars.lars / 100
## > [1] 3.998047
## > [1] 0.19
## > [1] 0.63
## >     gamma
## 0.4793769

barplot(mb.sel/100)

for(i in 1:length(res.al)){
    print(res.al[[i]]$select$beta)
}

hist()

o2 <- load("oracle2d.RData")


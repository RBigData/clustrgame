suppressMessages(library(float))
suppressMessages(library(kazaam))
suppressMessages(library(clustrgame))
comm.set.seed(1234, diff=TRUE)

# rnorm_fun = stats::rnorm
rnorm_fun = curand::rnorm

.pbd_env$SPMD.CT$print.quiet = TRUE

m = 3e5
m.local = kazaam:::get_local_dim(m)
n = 250
k = 3
maxiter = 2

centers = c(0, 2, 10)

generator = function(m, n, centers)
{
  centers_rows = sample(3, size=m, replace=TRUE)
  data = matrix(0.0, m, n)
  for (i in 1:length(centers))
  {
    rows = which(centers_rows == i)
    data[rows, ] = rnorm_fun(length(rows)*n, mean=centers[i])
  }
  
  data
}

generator = compiler::cmpfun(generator)
t_gen = system.time(data <- generator(m.local, n, centers))
comm.cat("data gen (dbl):   ", t_gen[3], "\n")

x = shaq(data, m, n, checks=FALSE)
s = shaq(fl(data), m, n, checks=FALSE)


t_kazaam = system.time(out1 <- kazaam::km(x, k=k, maxiter=maxiter))
comm.cat("kazaam (dbl):     ", t_kazaam[3], "\n")
t_clustrgame = system.time(out2 <- clustrgame::km_game(x, k=k, maxiter=maxiter))
comm.cat("clustrgame (dbl): ", t_clustrgame[3], "\n")
t_clustrgamef = system.time(out3 <- clustrgame::km_game(s, k=k, maxiter=maxiter))
comm.cat("clustrgame (flt): ", t_clustrgamef[3], "\n")

comm.print(c(out1$iterations, out2$iterations, out3$iterations))


finalize()

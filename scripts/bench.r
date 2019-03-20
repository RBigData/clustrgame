suppressMessages(library(kazaam))
suppressMessages(library(clustrgame))

m = 3e5
m.local = kazaam:::get_local_dim(m)
n = 250
k = 3
maxiter = 2

centers = c(0, 2, 10)

generator = function(m, n, centers)
{
  data = matrix(0.0, m, n)
  for (i in 1:m)
    data[i, ] = stats::rnorm(n, mean=sample(centers, size=1))
  
  data
}

generator = compiler::cmpfun(generator)
data = generator(m.local, n, centers)

x = shaq(data, m, n, checks=FALSE)

t_old = system.time(out1 <- kazaam::km(x, k=k, maxiter=maxiter))
comm.print(t_old)
t_new = system.time(out2 <- clustrgame::km_game(x, k=k, maxiter=maxiter))
comm.print(t_new)
# t_new = system.time(out2 <- clustrgame::km_game(x, k=k, maxiter=maxiter))
# comm.print(t_new)

out2$iterations

# comm.print(out1)
# comm.print(out2)

finalize()

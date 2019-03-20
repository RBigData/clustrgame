suppressMessages(library(clustrgame))

if (comm.rank() == 0){
  x = as.matrix(iris[, 1:4])
  y = as.numeric(iris[, 5])
} else {
  x = NULL
  y = NULL
}

dx = expand(x)
dy = expand(y)

dl = km_game(dx, k=3)$labels
l = collapse(dl)

if (comm.rank() == 0){
  rm = rand_measure(y, l[, 1, drop=TRUE])
  print(rm)
}

finalize()

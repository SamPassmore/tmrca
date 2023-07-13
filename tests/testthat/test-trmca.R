
test_that("Single TMRCA call", {
  set.seed(1234)
  tree = ape::rcoal(40, br = runif)
  states = ape::rTraitDisc(tree)

  out = .tmrca_onetree(tree, states)

  expect_equal(out$Node, c("56", "59", "62", "67", "76", "77"))

  expect_equal(out, structure(
    list(
      Node = c("56", "59", "62", "67", "76", "77"),
      Date = c(10.257, 9.755, 9.02, 6.414, 1.969, 1.345),
      n_tips = c(2L,
                 6L, 8L, 4L, 2L, 2L),
      tips = c(
        "t19; t22",
        "t12; t8; t27; t33; t3; t20",
        "t18; t25; t11; t5; t26; t28; t10; t38",
        "t32; t31; t35; t7",
        "t21; t39",
        "t6; t34"
      ),
      state = structure(
        c(2L, 1L, 1L, 2L,
          1L, 1L),
        levels = c("A", "B"),
        class = "factor"
      )
    ),
    row.names = c(NA,-6L),
    class = "data.frame"
  )
  ),
    row.names = c(NA,-6L),
    class = "data.frame"
  ))

})

# plot(tree, label.offset = 0.25)
# ape::tiplabels(pch = 19, col = states)
# ape::nodelabels()


test_that("Looped TMRCA call", {
  set.seed(1234)
  tree = ape::rcoal(40, br = runif)
  states = ape::rTraitDisc(tree)

  trees = rep(tree, 10)
  class(trees) = "multiPhylo"

  out = tmrca(trees, states)

})

trees = rep(tree, 10)
tic()
test = tmrca(trees, states)
toc()

tic()
#loop over each tree
I <- 0
mono.dates <- foreach(tree = trees, .combine=rbind) %do% {
  I <- I + 1
  print(names(trees)[I])

  #get monophyletic clades
  all.sub.trees <- clade.members.list(phy=tree, tip.labels=T)
  tip.states <- lapply(all.sub.trees, function(x) sapply(strsplit(x, "_"), "[", 3))
  samples <- lapply(all.sub.trees, function(x) sapply(strsplit(x, "_"), "[", 1))
  clades <- sapply(tip.states, function(x) length(unique(x))==1)
  clade.tips <- all.sub.trees[which(clades==T)]

  #node dates
  node.height <- nodeHeights(tree)
  tree.height <- max(node.height)
  nh <- data.table(Node=c(tree$edge[, 1], tree$edge[, 2]),
                   Date=c(round(tree.height - node.height[, 1], 3),
                          round(tree.height - node.height[, 2], 3)))
  nh <- nh[!duplicated(Node)]

  #find largest clades
  l <- length(clade.tips)
  L <- names(clade.tips)
  j <- L[1]
  node.mat <- NULL
  while(length(L)>0){
    c.now <- clade.tips[[which(names(clade.tips)==j)]]
    ol <- sapply(clade.tips, function(x){length(intersect(c.now, x))})
    clade.remove <- as.numeric(names(ol[ol>0]))
    big.clade <- as.numeric(names(ol[ol==max(ol)]))
    out <- data.table(Tree=names(trees)[I], nh[Node==big.clade],
                      #State=unique(ts.now), N.tips=length(ts.now),
                      State=unique(c.now), N.tips=length(c.now),
                      tips=paste(sort(samples[[which(names(samples)==big.clade)]]),
                                 collapse=","))
    node.mat <- rbind(node.mat, out)
    L <- setdiff(L, clade.remove)
    j <- L[1]
  }
  node.mat
}
toc()

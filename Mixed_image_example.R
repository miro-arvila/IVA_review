library(imager)
library(ivaBSS)

# Load the images
img1 <- load.image("duckboard.jpg")
img2 <- load.image("fruits.jpg")
img3 <- load.image("car.jpg")
img4 <- load.image("hut.jpg")
img5 <- load.image("bird.jpg")
imgs <- list(img1, img2, img3, img4, img5)

par(mfrow = c(4, 5),    
    oma = c(0, 0, 0, 0), 
    mar = c(1, 1, 0, 0), 
    mgp = c(1, 0.5, 0), 
    xpd = NA)

set.seed(3)

# Plots of the images
plot(img1, axes=FALSE)
plot(img2, axes=FALSE)
plot(img3, axes=FALSE)
plot(img4, axes=FALSE)
plot(img5, axes=FALSE)

# P is the amount of components, t samples and K datasets
P <- 5; t <- 10000; K <- 3;
S <- array(NA, c(P,t,K))
for (i in 1:P) {
  S[i,,1] <- as.data.frame(R(imgs[[i]]))$value
  S[i,,2] <- as.data.frame(G(imgs[[i]]))$value
  S[i,,3] <- as.data.frame(B(imgs[[i]]))$value
}


# Mixing matrices
A <- array(rnorm(P*P*K), c(P,P,K))


# Observed components
X <- array(NaN, c(P,t,K))
for (k in 1:K) {
  X[,,k] <- svd(A[,,k])$u%*%S[,,k]
}

# mimg means mixed image

mimg1 <- imfill(100,100, val=c(1,1,1))
mimg1[,,1] <- X[1,,1]
mimg1[,,2] <- X[1,,2]
mimg1[,,3] <- X[1,,3]

mimg2 <- imfill(100,100, val=c(1,1,1))
mimg2[,,1] <- X[2,,1]
mimg2[,,2] <- X[2,,2]
mimg2[,,3] <- X[2,,3]

mimg3 <- imfill(100,100, val=c(1,1,1))
mimg3[,,1] <- X[3,,1]
mimg3[,,2] <- X[3,,2]
mimg3[,,3] <- X[3,,3]

mimg4 <- imfill(100,100, val=c(1,1,1))
mimg4[,,1] <- X[4,,1]
mimg4[,,2] <- X[4,,2]
mimg4[,,3] <- X[4,,3]

mimg5 <- imfill(100,100, val=c(1,1,1))
mimg5[,,1] <- X[5,,1]
mimg5[,,2] <- X[5,,2]
mimg5[,,3] <- X[5,,3]


# Plots of mixed images
plot(mimg1, axes=FALSE)
plot(mimg2, axes=FALSE)
plot(mimg3, axes=FALSE)
plot(mimg4, axes=FALSE)
plot(mimg5, axes=FALSE)

## IVA-G to the mixed images with Gaussian density

res_G <- NewtonIVA(X, source_density = "gaussian", init = "none", verbose = F, eps = 1e-14, max_iter = 2000, step_size_min = 0.1, alpha = 0.9)

range_G <- max(res_G$S) - min(res_G$S)
sol_G <- res_G$S
sol_G <- sol_G - min(res_G$S)
sol_G <- sol_G / range_G

# rimg means result image

rimg1_G <- imfill(100,100, val=c(1,1,1))
rimg1_G[,,1] <- res_G$S[1,,1]
rimg1_G[,,2] <- res_G$S[1,,2]
rimg1_G[,,3] <- res_G$S[1,,3]

rimg2_G <- imfill(100,100, val=c(1,1,1))
rimg2_G[,,1] <- res_G$S[2,,1]
rimg2_G[,,2] <- res_G$S[2,,2]
rimg2_G[,,3] <- res_G$S[2,,3]

rimg3_G <- imfill(100,100, val=c(1,1,1))
rimg3_G[,,1] <- res_G$S[3,,1]
rimg3_G[,,2] <- res_G$S[3,,2]
rimg3_G[,,3] <- res_G$S[3,,3]

rimg4_G <- imfill(100,100, val=c(1,1,1))
rimg4_G[,,1] <- res_G$S[4,,1]
rimg4_G[,,2] <- res_G$S[4,,2]
rimg4_G[,,3] <- res_G$S[4,,3]

rimg5_G <- imfill(100,100, val=c(1,1,1))
rimg5_G[,,1] <- res_G$S[5,,1]
rimg5_G[,,2] <- res_G$S[5,,2]
rimg5_G[,,3] <- res_G$S[5,,3]

plot(rimg1_G, axes=FALSE)
plot(rimg2_G, axes=FALSE)
plot(rimg3_G, axes=FALSE)
plot(rimg4_G, axes=FALSE)
plot(rimg5_G, axes=FALSE)

############ Fix the colors manually by multiplying the color channels by 1 or -1

fixedrimg1_G <- imfill(100,100, val=c(1,1,1))
fixedrimg1_G[,,1] <- res_G$S[1,,1]
fixedrimg1_G[,,2] <- res_G$S[1,,2]
fixedrimg1_G[,,3] <- res_G$S[1,,3]

fixedrimg2_G <- imfill(100,100, val=c(1,1,1))
fixedrimg2_G[,,1] <- -1*res_G$S[2,,1]
fixedrimg2_G[,,2] <- -1*res_G$S[2,,2]
fixedrimg2_G[,,3] <- -1*res_G$S[2,,3]

fixedrimg3_G <- imfill(100,100, val=c(1,1,1))
fixedrimg3_G[,,1] <- -1*res_G$S[3,,1]
fixedrimg3_G[,,2] <- res_G$S[3,,2]
fixedrimg3_G[,,3] <- -1*res_G$S[3,,3]

fixedrimg4_G <- imfill(100,100, val=c(1,1,1))
fixedrimg4_G[,,1] <- -1*res_G$S[4,,1]
fixedrimg4_G[,,2] <- -1*res_G$S[4,,2]
fixedrimg4_G[,,3] <- -1*res_G$S[4,,3]

fixedrimg5_G <- imfill(100,100, val=c(1,1,1))
fixedrimg5_G[,,1] <- res_G$S[5,,1]
fixedrimg5_G[,,2] <- -1*res_G$S[5,,2]
fixedrimg5_G[,,3] <- res_G$S[5,,3]


# Color fixed images
plot(fixedrimg1_G, axes=FALSE)
plot(fixedrimg2_G, axes=FALSE)
plot(fixedrimg3_G, axes=FALSE)
plot(fixedrimg4_G, axes=FALSE)
plot(fixedrimg5_G, axes=FALSE)

# The joint ISI value
joint_ISI(res_G$W,A)

xy = c(-2,2)
x = c(-1.5, 0, 1.5)
y = c(-1.5*tan(1/6),1.5/cos(1/6),-1.5*tan(1/6))
theta = 2*pi/3
rot = matrix(c(cos(theta), -sin(theta), sin(theta), cos(theta)), 2,2)
x1 = array(0, 3)
y1 = array(0, 3)
for (i in 1:3){
  x1[i] = (rot%*%matrix(c(x[i],y[i]),ncol=1))[1]
  y1[i] = (rot%*%matrix(c(x[i],y[i]),ncol=1))[2]
}


#lines(xy, c(0,0))
#lines(c(0,0), xy)

#1
plot(x,y, cex = 4, xlim = c(-2,5), ylim = c(-2,4), asp = 1, xlab = "", ylab = "", col = "blue")
points(x,y,pch = 65:67, col = "blue")
points(x1 +3 ,y1 + 2, cex = 4, col = "red")
points(x1 +3 ,y1 + 2, pch = 65:67, col = "red")
grid()

#2
plot(x,y, cex = 4, xlim = c(-2,5), ylim = c(-2,4), asp = 1, xlab = "", ylab = "", col = "blue")
points(x,y,pch = 65:67, col = "blue")
points(x1 +3 ,y1 + 2, cex = 4, col = "grey")
points(x1 +3 ,y1 + 2, pch = 65:67, col = "grey")
points(x1 ,y1, cex = 4, col = "red")
points(x1, y1, pch = 65:67, col = "red")
grid()

#3
plot(x,y, cex = 4, xlim = c(-2,5), ylim = c(-2,4), asp = 1, xlab = "", ylab = "", col = "blue")
points(x,y,pch = 65:67, col = "blue")
points(x1 +3 ,y1 + 2, cex = 4, col = "grey")
points(x1 +3 ,y1 + 2, pch = 65:67, col = "grey")
points(x1 ,y1, cex = 4, col = "grey")
points(x1, y1, pch = 65:67, col = "grey")
points(x, y, cex = 4, col = "red")
points(x, y, pch = 65:67, col = "red")
grid()


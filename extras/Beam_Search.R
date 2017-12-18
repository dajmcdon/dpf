library(MASS)
#a = load('~/Desktop/music.Rdata')

multi = function(a,b,n){
	h1 = a[1]
	s1 = a[2]
	m1 = a[3:(n+2)]
	Q1 = matrix(a[(n+3):(n*n+n+2)],n,n)
	h2 = b[1]
	s2 = b[2]
	m2 = b[3:(n+2)]
	Q2 = matrix(b[(n+3):(n*n+n+2)],n,n)
	Q = Q1+Q2
	m = ginv(Q)%*%(Q1%*%m1+Q2%*%m2)
	h = h1+h2+(-.5)*(t(m1)%*%Q1%*%m1+t(m2)%*%Q2%*%m2-t(m)%*%Q%*%m)
	s = s1
	multi = c(h,s,m,Q)
}

ext = function(a,n){
	h1 = a[1]
	s1 = a[2]
	m1 = a[3:(n+2)]
	Q1 = matrix(a[(n+3):(n*n+n+2)],n,n)
	h = h1
	s = s1
	m = c(rep(0,(n/2)),m1[1:(n/2)])
	Q = matrix(0,n,n)
	Q[(n/2+1):n,(n/2+1):n]=Q1[1:(n/2),1:(n/2)]
	ext = c(h,s,m,Q)
}

maxout = function(a,n){
	h1 = a[1]
	s1 = a[2]
	m1 = a[3:(n+2)]
	Q1 = matrix(a[(n+3):(n*n+n+2)],n,n)
	h = h1
	s = s1
	m = c(m1[1:(n/2)],rep(0,(n/2)))
	Q = matrix(0,n,n)
	Q[1:(n/2),1:(n/2)] = Q1[1:(n/2),1:(n/2)]-Q1[1:(n/2),(n/2+1):n]%*%ginv(Q1[(n/2+1):n,(n/2+1):n])%*%Q1[(n/2+1):n,1:(n/2)]
	maxout = c(h,s,m,Q)
}

show = function(a,n){
	h1 = a[1]
	s1 = a[2]
	m1 = a[3:(n+2)]
	Q1 = matrix(a[(n+3):(n*n+n+2)],n,n)
	print(h1)
	print(s1)
	print(m1)
	print(Q1)
}

expand = function(m,s,n){
	m = matrix(m,s,s)
	Q = matrix(0,n,n)
	Q[1:s,1:s] = m
	expand = Q
	expand
}

n = 10
N = n*n+n+2
sq2pid1 = 1/sqrt(2*pi)

##### Setting parameters. These are very important #####

y = tempo.mat[,20]  # choose the performance
x = note.onset
ioi = diff(c(x,61))
#color = rep(c(1,2,3,4,5,6,7),33)
#plot(x,y,type = "b",col = color)

mu_r = mean(y) #tempo prior
sd_r = 30 #standard deviation of tempo
sd_rd1 = 1/sd_r

mu_a = -40 #accel.
sd_a = 20 #standard deviation of accel.
sd_ad1 = 1/sd_a

mu_j = -40 #Stress: a sudden change of tempo
sd_j = 10
sd_jd1 = 1/sd_j

trans = t(matrix(c(.8,.1,0,.1,0,0,.8,.2,0,0,.4,0,.6,0,0,1,0,0,0,0,.5,0,0,0,.5),5,5)) #state 5 is not used now

##### The following parameters work OK #####

mu_l = 0 #observation error
sd_l = 20
sd_ld1 = 1/sd_l

mu_w = 0 #transition error
sd_w = .01
sd_wd1 = 1/sd_w

mu_c = 3 #3rd order accel
sd_c = 1
sd_cd1 = 1/sd_c

#### prepare Gaussian kernal for state transition
trans11 = matrix(c(1,0,0,0),2,2)	
trans11 = expand(trans11,2,n/2)	
		alpha = t(trans11)
		cond_var11 = matrix(0,n,n)
		gammad1 = diag(c(sd_wd1,rep(0,n/2-1)))
		cond_var11[1:(n/2),1:(n/2)] = gammad1
		cond_var11[(n/2+1):n,1:(n/2)] = -alpha%*%gammad1
		cond_var11[1:(n/2),(n/2+1):n] = -gammad1%*%t(alpha)
		cond_var11[(n/2+1):n,(n/2+1):n] = alpha%*%gammad1%*%t(alpha)

trans12 = matrix(c(1,0,1,1),2,2)
trans12 = expand(trans12,2,n/2)
		alpha = t(trans12)
		cond_var12 = matrix(0,n,n)
		gammad1 = diag(c(sd_wd1,sd_wd1,rep(0,n/2-2)))
		cond_var12[1:(n/2),1:(n/2)] = gammad1
		cond_var12[(n/2+1):n,1:(n/2)] = -alpha%*%gammad1
		cond_var12[1:(n/2),(n/2+1):n] = -gammad1%*%t(alpha)
		cond_var12[(n/2+1):n,(n/2+1):n] = alpha%*%gammad1%*%t(alpha)
		
trans14 = matrix(0,n/2,n/2)
trans14[1,1] = 1
trans14[n/2,1] = 1
trans14[1,n/2] = 1
		alpha = t(trans14)
		cond_var14 = matrix(0,n,n)
		gammad1 = diag(c(sd_wd1,rep(0,n/2-2),sd_jd1))
		cond_var14[1:(n/2),1:(n/2)] = gammad1
		cond_var14[(n/2+1):n,1:(n/2)] = -alpha%*%gammad1
		cond_var14[1:(n/2),(n/2+1):n] = -gammad1%*%t(alpha)
		cond_var14[(n/2+1):n,(n/2+1):n] = alpha%*%gammad1%*%t(alpha)
		
trans41 = matrix(0,n/2,n/2)
trans41[1,n/2] = 1
		alpha = t(trans41)
		cond_var41 = matrix(0,n,n)
		gammad1 = diag(c(sd_wd1,rep(0,n/2-1)))
		cond_var41[1:(n/2),1:(n/2)] = gammad1
		cond_var41[(n/2+1):n,1:(n/2)] = -alpha%*%gammad1
		cond_var41[1:(n/2),(n/2+1):n] = -gammad1%*%t(alpha)
		cond_var41[(n/2+1):n,(n/2+1):n] = alpha%*%gammad1%*%t(alpha)

trans15 = matrix(c(1,0,0,0,1,1,0,0,0,1,1,0,0,0,1,1),4,4)
trans15 = expand(trans15,4,n/2)
trans15[n/2,1] = 1
		alpha = t(trans15)
		cond_var15 = matrix(0,n,n)
		gammad1 = diag(c(sd_wd1,sd_wd1,sd_wd1,sd_wd1,rep(sd_wd1,n/2-4)))
		cond_var15[1:(n/2),1:(n/2)] = gammad1
		cond_var15[(n/2+1):n,1:(n/2)] = -alpha%*%gammad1
		cond_var15[1:(n/2),(n/2+1):n] = -gammad1%*%t(alpha)
		cond_var15[(n/2+1):n,(n/2+1):n] = alpha%*%gammad1%*%t(alpha)
		
trans55 = matrix(c(1,0,0,0,1,1,0,0,0,1,1,0,0,0,1,1),4,4)
trans55 = expand(trans55,4,n/2)
trans55[n/2,n/2] = 1
		alpha = t(trans55)
		cond_var55 = matrix(0,n,n)
		gammad1 = diag(c(sd_wd1,sd_wd1,sd_wd1,sd_wd1,rep(sd_wd1,n/2-4)))
		cond_var55[1:(n/2),1:(n/2)] = gammad1
		cond_var55[(n/2+1):n,1:(n/2)] = -alpha%*%gammad1
		cond_var55[1:(n/2),(n/2+1):n] = -gammad1%*%t(alpha)
		cond_var55[(n/2+1):n,(n/2+1):n] = alpha%*%gammad1%*%t(alpha)
		
trans51 = matrix(0,n/2,n/2)
trans51[1,n/2] = 1
		alpha = t(trans51)
		cond_var51 = matrix(0,n,n)
		gammad1 = diag(c(sd_wd1,rep(0,n/2-1)))
		cond_var51[1:(n/2),1:(n/2)] = gammad1
		cond_var51[(n/2+1):n,1:(n/2)] = -alpha%*%gammad1
		cond_var51[1:(n/2),(n/2+1):n] = -gammad1%*%t(alpha)
		cond_var51[(n/2+1):n,(n/2+1):n] = alpha%*%gammad1%*%t(alpha)


prun = 200  #beam width * number of states

branch = matrix(0,N,prun)
branch[1,]=-Inf

curr = temo = prev = matrix(0,length(y),prun)
curr[1,1]=1
curr[1,3]=3

	node_mean = c(mu_r,rep(0,n-1))
	node_var = matrix(c(sd_rd1^2,rep(0,n*n-1)),n,n)
	tempo_dist = c(log(sq2pid1*sd_rd1),1,node_mean,node_var)
	ob_dist = c(log(sq2pid1*sd_ld1),1,y[1],rep(0,n-1),sd_ld1^2,rep(0,n*n-1))
	branch[,1]=multi(tempo_dist,ob_dist,n)

show(branch[,1],n)
show(tempo_dist,n)

	node_mean = c(mu_r+mu_a*1.5,-mu_a,rep(0,n-2))
	node_var = matrix(c(sd_rd1^2,rep(0,n*n-1)),n,n)
	node_var[n+2] = sd_ad1^2
	tempo_dist = c(log(sq2pid1*sd_rd1),3,node_mean,node_var)
	ob_dist = c(log(sq2pid1*sd_ld1),1,y[1],rep(0,n-1),sd_ld1^2,rep(0,n*n-1))
	branch[,3]=multi(tempo_dist,ob_dist,n)


show(branch[,3],n)
show(tempo_dist,n)

#### beam search
len = length(y)
run = len
for(i in 2:run){
	#i = 2
	ob_dist = c(log(sq2pid1*sd_ld1),1,y[i],rep(0,n-1),sd_ld1^2,rep(0,n*n-1))
	new_branch = matrix(0,N,prun)
	new_branch[1,]=-Inf
	
	trans12 = matrix(c(1,0,ioi[i],1),2,2)
	trans12 = expand(trans12,2,n/2)
		alpha = t(trans12)
		cond_var12 = matrix(0,n,n)
		gammad1 = diag(c(sd_wd1,sd_wd1,rep(0,n/2-2)))
		cond_var12[1:(n/2),1:(n/2)] = gammad1
		cond_var12[(n/2+1):n,1:(n/2)] = -alpha%*%gammad1
		cond_var12[1:(n/2),(n/2+1):n] = -gammad1%*%t(alpha)
		cond_var12[(n/2+1):n,(n/2+1):n] = alpha%*%gammad1%*%t(alpha)
	
	for(j in 1:50){if(branch[2,j]!=0){
	#		j = 1
		if(branch[2,j]==1){
			cond_dist = c(log(trans[1,1])+log(sq2pid1*sd_wd1),1,rep(0,n),cond_var11)
			ext_dist = ext(branch[,j],n)
			joint_dist = multi(cond_dist,ext_dist,n)
			max_dist=maxout(joint_dist,n)
			new_branch[,j]=multi(max_dist,ob_dist,n)
			
			cond_dist = c(log(trans[1,2])+log(sq2pid1*sd_wd1),2,rep(0,n),cond_var12)
			temp_branch = branch[,j]
			temp_branch[4] = mu_a
			temp_branch[2*n+4] = sd_ad1
			temp_branch[2*n+3] = 0
			temp_branch[n+4]   = 0
			ext_dist = ext(temp_branch,n)
			joint_dist = multi(cond_dist,ext_dist,n)
			max_dist=maxout(joint_dist,n)
			new_branch[,j+50]=multi(max_dist,ob_dist,n)
			
			cond_dist = c(log(trans[1,4])+log(sq2pid1*sd_wd1),4,rep(0,n),cond_var14)
			temp_branch = branch[,j]
			temp_branch[n/2+2] = mu_j
			temp_branch[n*n/2+n/2+2] = sd_jd1
			ext_dist = ext(temp_branch,n)
			joint_dist = multi(cond_dist,ext_dist,n)
			max_dist=maxout(joint_dist,n)
			new_branch[,j+100]=multi(max_dist,ob_dist,n)
			
			cond_dist = c(log(trans[1,5])+log(sq2pid1*sd_wd1),5,rep(0,n),cond_var15)
			temp_branch = branch[,j]
			temp_branch[4] = 0
			temp_branch[5] = -mu_c
			temp_branch[6] = mu_c/2
			temp_branch[2*n+4] = sd_cd1
			temp_branch[3*n+5] = sd_cd1
			temp_branch[4*n+6] = sd_cd1
			ext_dist = ext(temp_branch,n)
			joint_dist = multi(cond_dist,ext_dist,n)
			max_dist=maxout(joint_dist,n)
			new_branch[,j+150]=multi(max_dist,ob_dist,n)
			new_branch[1,j+150]=-Inf			
		}
		if(branch[2,j]==2){
			cond_dist = c(log(trans[2,2])+log(sq2pid1*sd_wd1),2,rep(0,n),cond_var12)
			ext_dist = ext(branch[,j],n)
			joint_dist = multi(cond_dist,ext_dist,n)
			max_dist=maxout(joint_dist,n)
			new_branch[,j]=multi(max_dist,ob_dist,n)
			
			cond_dist = c(log(trans[2,3])+log(sq2pid1*sd_wd1),3,rep(0,n),cond_var12)
			temp_branch = branch[,j]
			temp_branch[4] = -mu_a
			temp_branch[2*n+4] = sd_ad1
			temp_branch[2*n+3] = 0
			temp_branch[n+4]   = 0
			ext_dist = ext(temp_branch,n)
			joint_dist = multi(cond_dist,ext_dist,n)
			max_dist=maxout(joint_dist,n)
			new_branch[,j+50]=multi(max_dist,ob_dist,n)
		}
		if(branch[2,j]==3){
			cond_dist = c(log(trans[3,3])+log(sq2pid1*sd_wd1),3,rep(0,n),cond_var12)
			ext_dist = ext(branch[,j],n)
			joint_dist = multi(cond_dist,ext_dist,n)
			max_dist=maxout(joint_dist,n)
			new_branch[,j]=multi(max_dist,ob_dist,n)
			#if (new_branch[3,j]<mu_r-20) new_branch[1,j]=-Inf
			
			node_mean = c(mu_r,rep(0,n-1))
			node_var = matrix(c(sd_rd1^2,rep(0,n*n-1)),n,n)
			tempo_dist = c(log(trans[3,1])+log(sq2pid1*sd_rd1),1,node_mean,node_var)
			new_branch[,j+50]=multi(tempo_dist,ob_dist,n)
			new_branch[1,j+50]=new_branch[1,j+50]+branch[1,j]
		}
		if(branch[2,j]==4){
			cond_dist = c(log(trans[4,1])+log(sq2pid1*sd_wd1),1,rep(0,n),cond_var41)
			ext_dist = ext(branch[,j],n)
			joint_dist = multi(cond_dist,ext_dist,n)
			max_dist=maxout(joint_dist,n)
			new_branch[,j]=multi(max_dist,ob_dist,n)
		}
		if(branch[2,j]==5){
			cond_dist = c(log(trans[5,5])+log(sq2pid1*sd_wd1),6,rep(0,n),cond_var55)
			ext_dist = ext(branch[,j],n)
			joint_dist = multi(cond_dist,ext_dist,n)
			max_dist=maxout(joint_dist,n)
			new_branch[,j]=multi(max_dist,ob_dist,n)
		}
		if(branch[2,j]==6){
			cond_dist = c(log(trans[5,1])+log(sq2pid1*sd_wd1),1,rep(0,n),cond_var51)
			ext_dist = ext(branch[,j],n)
			joint_dist = multi(cond_dist,ext_dist,n)
			max_dist=maxout(joint_dist,n)
			new_branch[,j]=multi(max_dist,ob_dist,n)
		}
	}}
	prev[i,] = c(1:prun)
	h_sort = order(new_branch[1,],decreasing = T)
	prev[i,] = prev[i,h_sort]
	for (k in 1:N) new_branch[k,] = new_branch[k,h_sort]
	curr[i,] = new_branch[2,]
	temo[i,] = new_branch[3,]
	branch = new_branch
	print(i)
#	print(branch[c(1:6,n+3),1:5])
}
#### trace back the discrete states
k = 1
while(branch[2,k]>=5) k=k+1
s = branch[2,k]
t = branch[3,k]
prev_index = k
for(i in (run-1):1){
	prev_index = (prev[i+1,prev_index]-1)%%50+1
	s = c(curr[i,prev_index],s)
	print(c(i,prev_index,curr[i,prev_index]))
}


#### compute the continuous variables with discrete states fixed
path = matrix(0,N,len)
if(s[1]==1){
	node_mean = c(mu_r,rep(0,n-1))
	node_var = matrix(c(sd_rd1^2,rep(0,n*n-1)),n,n)
	tempo_dist = c(log(sq2pid1*sd_rd1),1,node_mean,node_var)
	ob_dist = c(log(sq2pid1*sd_ld1),1,y[1],rep(0,n-1),sd_ld1^2,rep(0,n*n-1))
	path[,1]=multi(tempo_dist,ob_dist,n)
}
if(s[1]==3){
	node_mean = c(mu_r+mu_a,-mu_a,rep(0,n-2))
	node_var = matrix(c(sd_rd1^2,rep(0,n*n-1)),n,n)
	node_var[n+2] = sd_ad1^2
	tempo_dist = c(log(sq2pid1*sd_rd1),3,node_mean,node_var)
	ob_dist = c(log(sq2pid1*sd_ld1),1,y[1],rep(0,n-1),sd_ld1^2,rep(0,n*n-1))
	path[,1]=multi(tempo_dist,ob_dist,n)
}
for(i in 2:run){
	#i=5
	ob_dist = c(log(sq2pid1*sd_ld1),1,y[i],rep(0,n-1),sd_ld1^2,rep(0,n*n-1))
	
	trans12 = matrix(c(1,0,ioi[i],1),2,2)
	trans12 = expand(trans12,2,n/2)
		alpha = t(trans12)
		cond_var12 = matrix(0,n,n)
		gammad1 = diag(c(sd_wd1,sd_wd1,rep(0,n/2-2)))
		cond_var12[1:(n/2),1:(n/2)] = gammad1
		cond_var12[(n/2+1):n,1:(n/2)] = -alpha%*%gammad1
		cond_var12[1:(n/2),(n/2+1):n] = -gammad1%*%t(alpha)
		cond_var12[(n/2+1):n,(n/2+1):n] = alpha%*%gammad1%*%t(alpha)
		
	if(s[i-1]==1){
		if(s[i]==1){	
			cond_dist = c(log(trans[1,1])+log(sq2pid1*sd_wd1),1,rep(0,n),cond_var11)
			ext_dist = ext(path[,i-1],n)
			joint_dist = multi(cond_dist,ext_dist,n)
			max_dist=maxout(joint_dist,n)
			path[,i]=multi(max_dist,ob_dist,n)
		}
		if(s[i]==2){	
			cond_dist = c(log(trans[1,2])+log(sq2pid1*sd_wd1),2,rep(0,n),cond_var12)
			temp_branch = path[,i-1]
			temp_branch[4] = mu_a
			temp_branch[2*n+4] = sd_ad1^2
			temp_branch[2*n+3] = 0
			temp_branch[n+4]   = 0
			ext_dist = ext(temp_branch,n)
			joint_dist = multi(cond_dist,ext_dist,n)
			max_dist=maxout(joint_dist,n)
			path[,i]=multi(max_dist,ob_dist,n)
		}
		if(s[i]==4){	
			cond_dist = c(log(trans[1,4])+log(sq2pid1*sd_wd1),4,rep(0,n),cond_var14)
			temp_branch = path[,i-1]
			temp_branch[n/2+2] = mu_j
			temp_branch[n*n/2+n/2+2] = sd_jd1^2
			ext_dist = ext(temp_branch,n)
			joint_dist = multi(cond_dist,ext_dist,n)
			max_dist=maxout(joint_dist,n)
			path[,i]=multi(max_dist,ob_dist,n)
		}
		if(s[i]==5){	
			cond_dist = c(log(trans[1,5])+log(sq2pid1*sd_wd1),5,rep(0,n),cond_var15)
			temp_branch = path[,i-1]
			temp_branch[4] = 0
			temp_branch[5] = -mu_c
			temp_branch[6] = mu_c/2
			temp_branch[2*n+4] = sd_cd1
			temp_branch[3*n+5] = sd_cd1
			temp_branch[4*n+6] = sd_cd1
			ext_dist = ext(temp_branch,n)
			joint_dist = multi(cond_dist,ext_dist,n)
			max_dist=maxout(joint_dist,n)
			path[,i]=multi(max_dist,ob_dist,n)
			#path[1,i]=-Inf	
		}		
	}
	if(s[i-1]==2){
		if(s[i]==2){
			cond_dist = c(log(trans[2,2])+log(sq2pid1*sd_wd1),2,rep(0,n),cond_var12)
			ext_dist = ext(path[,i-1],n)
			joint_dist = multi(cond_dist,ext_dist,n)
			max_dist=maxout(joint_dist,n)
			path[,i]=multi(max_dist,ob_dist,n)
		}
		if(s[i]==3){	
			cond_dist = c(log(trans[2,3])+log(sq2pid1*sd_wd1),3,rep(0,n),cond_var12)
			temp_branch = path[,i-1]
			temp_branch[4] = -mu_a
			temp_branch[2*n+4] = sd_ad1^2
			temp_branch[2*n+3] = 0
			temp_branch[n+4]   = 0
			ext_dist = ext(temp_branch,n)
			joint_dist = multi(cond_dist,ext_dist,n)
			max_dist=maxout(joint_dist,n)
			path[,i]=multi(max_dist,ob_dist,n)
		}
	}
	if(s[i-1]==3){
		if(s[i]==3){
			cond_dist = c(log(trans[3,3])+log(sq2pid1*sd_wd1),3,rep(0,n),cond_var12)
			ext_dist = ext(path[,i-1],n)
			joint_dist = multi(cond_dist,ext_dist,n)
			max_dist=maxout(joint_dist,n)
			path[,i]=multi(max_dist,ob_dist,n)
			#if (new_branch[3,j]<mu_r-20) new_branch[1,j]=-Inf
		}
		if(s[i]==1){	
			node_mean = c(mu_r,rep(0,n-1))
			node_var = matrix(c(sd_rd1^2,rep(0,n*n-1)),n,n)
			tempo_dist = c(log(trans[3,1])+log(sq2pid1*sd_rd1),1,node_mean,node_var)
			path[,i]=multi(tempo_dist,ob_dist,n)
			path[1,i]=path[1,i]+path[1,i-1]
		}
	}
	if(s[i-1]==4){
			cond_dist = c(log(trans[4,1])+log(sq2pid1*sd_wd1),1,rep(0,n),cond_var41)
			ext_dist = ext(path[,i-1],n)
			joint_dist = multi(cond_dist,ext_dist,n)
			max_dist=maxout(joint_dist,n)
			path[,i]=multi(max_dist,ob_dist,n)
	}
	if(s[i-1]==5){
			cond_dist = c(log(trans[5,5])+log(sq2pid1*sd_wd1),6,rep(0,n),cond_var55)
			ext_dist = ext(path[,i-1],n)
			joint_dist = multi(cond_dist,ext_dist,n)
			max_dist=maxout(joint_dist,n)
			path[,i]=multi(max_dist,ob_dist,n)
	}
	if(s[i-1]==6){
			cond_dist = c(log(trans[5,1])+log(sq2pid1*sd_wd1),1,rep(0,n),cond_var51)
			ext_dist = ext(path[,i-1],n)
			joint_dist = multi(cond_dist,ext_dist,n)
			max_dist=maxout(joint_dist,n)
			path[,i]=multi(max_dist,ob_dist,n)
	}
}


t = rep(0,len)
i = 2
while(i<=len){
	if(path[2,i-1]==1 && path[2,i]==2){
		k = i-1
		while(path[2,k]!=3 && k>0) {
			if(path[2,k]<4) t[k]=path[3,i-1]
			k = k - 1
		}
	}
	if(i==len){
		k = i
		while(path[2,k]!=1) k = k-1
		l = k
		while(path[2,l]!=3 && l>0) {
			if(path[2,l]<4) t[l]=path[3,k]
			l = l -1
		}
		if(path[2,i]==2) t[(k+1):i]=t[k-1]+(cumsum(ioi[(k+1):i]))*path[4,i]
	}
	if(path[2,i-1]==2 && path[2,i]==3){
		k = i-1
		while(path[2,k]==2 && k>0) k = k -1
		k = k + 1
		t[k:(i-1)]=t[k-1]+(cumsum(ioi[k:(i-1)]))*path[4,i-1]
	}
	if(path[2,i-1]==3 && path[2,i]==1){
		k = i-1
		while(path[2,k]==3 && k>0) k = k -1
		k = k + 1
		if(k>1) t[k:(i-1)]=t[k-1]+(cumsum(ioi[k:(i-1)]))*path[4,i-1]
		if(k==1) t[k:(i-1)]=path[3,i]-order((cumsum(ioi[k:(i-1)])),decreasing = T)*path[4,i-1]
	}
	if(path[2,i]==4) t[i]=y[i]	
	if(path[2,i-1]==9 && path[2,i]==1){
		k = i-1
		#print(path[2,k])
		while(path[2,k]!=1) k = k-1
		k = k+1
		vec = c(path[3,k-1],0,path[5,i-1]-5*path[6,i-1],path[6,i-1],0)
		for(l in k:(i-1)){
			vec = trans55%*%vec
			t[l] = vec[1]
			print(c(path[3,l],path[4,l],path[5,l],path[6,l],0))
		}
		print(c(path[3,k-1],0,path[5,i-1],path[6,i-1],0))
	}
	i = i+1
}
s[s>5]=5
par(mfrow = c(2,1),mar=c(3,1,0,0))
plot(x*3,y,col= s,type = "b")#,ylim = c(0,600))
plot(x*3,t,col= s,type = "b")#,ylim = c(0,600))

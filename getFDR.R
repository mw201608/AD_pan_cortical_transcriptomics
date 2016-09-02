getFDR=function(x,x0){
#x, a vector of non-negative test statistics from original data.
#x0, a vector of non-negative test statistics from permutation
	if(any(x<0) || any(x0<0)) stop('Input must be non-negative\n')
	n1=length(x)
	n2=length(x0)
	x0=sort(x0)
	order.x=order(x)
	rank.x=rep(0L,n1)
	rank.x[order.x]=1:n1
	p=rep(0,n1)
	p[order.x]=n2-findInterval(x[order.x],x0)
	p=(p/n2)/((n1+1-rank.x)/(n1))
	p[p>1]=1
	if(n1==1) return(p)
	for(i in 1:(n1-1)) if(p[order.x[i+1]] > p[order.x[i]]){p[order.x[i+1]]=p[order.x[i]]} #make sure the larger statistics will have a FDR no bigger than that of the smaller statistics
	p
}

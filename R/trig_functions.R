
### TRIGONOMETRIC FUNCTIONS TO COMPUTE AMPLITUDE AND PHASE

## Function for amplitude
# x is a dataframe / matrix with columns an, bn, cn, dn. It can also be done by x (an, bn) and y (cn, dn)

amplitude <- function(x, coordinate){
	if (coordinate == "x"){
		result <- apply(x,1,function(x) sqrt(x[1]^2+x[2]^2))
	} else if (coordinate == "y"){
		result <- apply(x,1,function(x) sqrt(x[3]^2+x[4]^2))
	} else if (coordinate == "all") {
		result <- apply(x,1,function(x) sqrt(x[1]^2+x[2]^2+x[3]^2+x[4]^2))
	}
	return(result)
}


## Function for phase
## x is a matrix with columns an, bn, cn, dn
## arg coordinate establish with coordinate to compute. If "x", it will compute on an and bn; if "y" it will compute on cn and dn
phase <- function(x, coordinate){
	if (coordinate == "x"){
		obj <- x[,c(1,2)]
	} else if (coordinate == "y"){
		obj <- x[,c(3,4)]
	}
	result <- apply(obj,1,function(x) atan2(x[1],x[2]))
	return(result)
}

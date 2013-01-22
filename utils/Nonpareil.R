
# INTERNAL
Nonpareil.__init_globals <- function(soft=TRUE){
	if(soft && exists('Nonpareil.LastFactor')) return();
	Nonpareil.LastFactor <<- NULL;
	Nonpareil.LastXmax <<- NULL;
	Nonpareil.OpenColors <<- c();
	Nonpareil.OpenNames <<- c();
}

# MODEL FUNCTIONS
#Nonpareil.f <- function(x,a,b){ return( ( 1 - exp(-x/(2^a)) )^b ) }
#Nonpareil.antif <- function(y,a,b){ return( -(2^a)*log(1 - y^(1/b)) ) }
Nonpareil.f <- function(x,a,b){ return(pgamma(log(x+1),a, b)) }
Nonpareil.antif <- function(y,a,b){ return(exp(qgamma(y,a,b))-1) }

# NONPAREIL CURVES
Nonpareil.curve <- function(file,alldata=NULL,factor=NULL,libs=7,modelSD=0,plotError=TRUE,
			xmax=NULL,ymax=1,xlab=NULL,ylab=NULL,r=NULL,g=NULL,b=NULL,
			new=TRUE,plot=TRUE,libname=NULL,model=TRUE,modelOnly=FALSE,
			curve.lwd=2, curve.alpha=0.4, model.lwd=1, model.alpha=1, curve.log='',
			ymin=NULL, xmin=NULL){
	# Create environment
	Nonpareil.__init_globals(!new);
	# Examine consistency
	if(is.null(file)) stop('The file argument is mandatory');
	if(!new && (is.null(Nonpareil.LastFactor) || is.null(Nonpareil.LastXmax)))
		stop('No previous plot found, please use new=TRUE if this is the first plot.');
	if(modelSD!=0 && !is.null(alldata))
		stop('Only one parameter between alldata and modelSD can be passed.');
	if(modelOnly && !model)
		stop('Cannot set modelOnly=T and model=F in the same call');
	# Read input
	a <- read.table(file, sep="\t", h=F);
	Frn <- TRUE;
	if(max(a$V1)>1) Frn <- FALSE;
	if(is.null(factor)){
		if(Frn){	factor <- 1 }
		else if(new){	factor <- 10^(-3*floor(log10(max(a$V1))/3)) }
		else{		factor <- Nonpareil.LastFactor }
	}
	if(is.null(xmax)){
		if(new){	xmax <- max(a$V1)*libs*factor }
		else{		xmax <- Nonpareil.LastXmax }
	}
	if(is.null(libname)) {
	   libname <- basename(file);
	   if(substr(libname, nchar(libname)-3, nchar(libname))==".npo")
	      libname <- substr(libname, 0, nchar(libname)-4);
	}
	
	# For future calls (of this or other Nonpareil.* functions)
	Nonpareil.LastFactor <<- factor;
	
	
	# Draw it
	if(is.null(ymin)) ymin <- min(a$V2[a$V2>0]);
	if(is.null(xmin)) xmin <- min(a$V1[a$V1>0])*factor;
	if(plot){
		if(new){
			if(!is.null(xlab)){	xlab <- as.character(xlab) }
			else if(Frn){		xlab <- 'Number of libraries' }
			else if(factor==1){	xlab <- 'Reads' }
			else if(factor==1e-1){	xlab <- 'Tens of reads' }
			else if(factor==1e-2){	xlab <- 'Hundreds of reads' }
			else if(factor==1e-3){	xlab <- 'Thousands of reads' }
			else if(factor==1e-6){	xlab <- 'Millions of reads' }
			else if(factor==1e-9){	xlab <- 'Billions of reads' }
			else if(factor==1e-12){	xlab <- 'Trillions of reads' }
			else if(factor==1e-15){	xlab <- 'Zillions of reads' }
			else if(factor==1e-18){	xlab <- 'Gazillions of reads' }
			else{			xlab <- paste('Reads (by ', factor, ')', sep='') }
			if(is.null(ylab)) ylab <- 'Redundant fraction';

			plot(1, t='n', xlim=c(xmin, xmax), ylim=c(ymin, ymax), xlab=xlab, ylab=ylab, log=curve.log);
			abline(h=c(1,0.05), lty=2, col='red');
		}

		if(is.null(r)) r <- runif(1);
		if(is.null(g)) g <- runif(1);
		if(is.null(b)) b <- runif(1);
		if(r>1) r <- r/256;
		if(g>1) g <- g/256;
		if(b>1) b <- b/256;
		
		if(!modelOnly){
			if(plotError){
				err.x <- c(a$V1, rev(a$V1))*factor;
				err.y <- c(a$V2+a$V3, rev(a$V2-a$V3));
				polygon(err.x, ifelse(err.y<=ymin*0.1, ymin*0.1, err.y), col=rgb(r,g,b,.2), border=NA);
			}
			lines(a$V1*factor, a$V2, col=rgb(r,g,b,curve.alpha), lwd=curve.lwd);
		}
		
		# Save some info
		Nonpareil.LastXmax   <<- xmax;
		Nonpareil.OpenColors <<- c(Nonpareil.OpenColors, rgb(r,g,b));
		Nonpareil.OpenNames  <<- c(Nonpareil.OpenNames, libname);
	}
	
	# Model it
	if(model){
	   if(is.null(alldata)){
	      data <- data.frame(x=a$V1*factor, y=a$V2+modelSD*a$V3);
	   }else{
	      all  <- read.table(alldata, sep="\t", h=F);
	      data <- data.frame(x=a$V1*factor, y=a$V3);
	   }
	   model <- nls(y ~ Nonpareil.f(x, a, b), data=data, start=list(a=1, b=0.1), lower=c(a=0, b=0), algorithm='port', 
	   		control=nls.control(minFactor=1e-25000, tol=1e-5, maxiter=102400, warnOnly=T));
	   #model <- nls(y ~ Nonpareil.f(x, a, b), data=data, start=list(a=1, b=1), 
	   #		control=nls.control(minFactor=1e-250000000, tol=1e-150000000, maxiter=1024000));
	   model.lty=2;
	   if(modelOnly) model.lty=1;
	   if(plot)
	      lines(seq(0, xmax, length.out=1e4), predict(model, list(x=seq(0, xmax, length.out=1e4))),
	      		col=rgb(r,g,b,model.alpha), lty=model.lty, lwd=model.lwd);
	   return(model);
	}
}

Nonpareil.legend <- function(x=NULL, y=.3, ...){
	if(is.null(Nonpareil.LastFactor) || is.null(Nonpareil.LastXmax))
		stop('There must be at least one active plot to draw a legend.  Use Nonpareil.curve().');
	if(is.null(x)) x <- 0.75*Nonpareil.LastXmax;
	legend(x=x, y=y, legend=Nonpareil.OpenNames, fill=Nonpareil.OpenColors, ...);
}

## PREDICTIONS
Nonpareil.predict <- function(model, y=0.05, x=1e6*Nonpareil.LastFactor, stderr=0, absolute=FALSE, conv.only=TRUE, antif=TRUE){
	if(is.null(model)) stop('The model argument is mandatory.');
	if(y>=1)stop('The y argument must be strictly lesser than 1.');
	if(y<0)	stop('The y argument must be positive.');
	if(x<0)stop('The x argument must be positive.');
	if(conv.only && !summary(model)$convInfo$isConv) stop('The model did not converge.');
	if(absolute && is.null(Nonpareil.LastFactor))
		stop('Cannot set absolute=TRUE: Cannot find the Nonpareil.LastFactor global variable');

	a <- as.numeric(summary(model)$parameters['a', c(1, 2)]);
	b <- as.numeric(summary(model)$parameters['b', c(1, 2)]);
	if(stderr==0){
		if(antif){
			out <- Nonpareil.antif(y, a[1], b[1]);
		}else{
			out <- Nonpareil.f(x, a[1], b[1]);
		}
	}else{
		if(antif){
			out <- c(
				Nonpareil.antif(y, a[1]-a[2]*stderr, b[1]-b[2]*stderr),
				Nonpareil.antif(y, a[1]+a[2]*stderr, b[1]+a[2]*stderr)
			);
		}else{
			out <- c(
				Nonpareil.f(x, a[1]-a[2]*stderr, b[1]-b[2]*stderr),
				Nonpareil.f(x, a[1]+a[2]*stderr, b[1]+a[2]*stderr)
			);
		}
	}
	if(absolute) return(out/Nonpareil.LastFactor);
	return(out);
}

Nonpareil.estimate <- function(file, overlap=75, plot=FALSE, factor=1e-3, ...){
	if(is.null(file)) stop('The file argument is mandatory.');
	# Fit parameters
	if(overlap==100){
		rho.a <- 5.564128;
		rho.b <- 1.142516;
		Rstar.b <- 1.196262;
	}else if(overlap==75){
		rho.a <- 1.921118;
		rho.b <- 1.000554;
		Rstar.b <- 1.524259;
	}else if(overlap==50){
		rho.a <- 1.420594;
		rho.b <- 1.028292;
		Rstar.b <- 1.615895;
	}else if(overlap==25){
		rho.a <- 1.238295;
		rho.b <- 1.143099;
		Rstar.b <- 1.671156;
	}else{
		stop('Unsupported overlap.  Supported overlaps are: 100, 75, 50, and 25.');
	}
	
	# Sequencing depth
	data <- read.table(file, sep="\t", h=F);
	kappa <- tail(data$V2, n=1)
	rho <- exp(rho.a + rho.b*log(kappa));
	
	# Sequencing effort for nearly complete coverage
	model <- Nonpareil.curve(file, factor=factor, plot=plot, ...);
	Sstar <- Nonpareil.predict(model)/factor;
	Rstar <- exp(Rstar.b*log(Sstar));

	# Return
	out <- data.frame(rho, Rstar);
	return(out);
}



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
#Nonpareil.f <- function(x,a,b){ lx<-log(max(1, x)); return(exp((lx^a)/(b + lx^a))) }
#Nonpareil.antif <- function(y,a,b){ ly<-log(max(1e-10, y)); return(exp((b*ly/(1-ly))^(1/a))) }

# NONPAREIL CURVES
Nonpareil.curve.batch <- function(files, overlap, r=NA, g=NA, b=NA, libnames=NA, ...){
   if(!is.vector(files)) files = as.vector(files);
   new=TRUE;
   for(i in 1:length(files)){
      o = Nonpareil.curve(files[i], overlap, r=r[i], g=g[i], b=b[i], libname=ifelse(is.null(libnames), NA, libnames[i]), new=new, ...);
      if(new){
	 out.m = matrix(NA, ncol=length(o), nrow=length(files));
         colnames(out.m) <- rownames(as.matrix(o));
	 if(!is.na(libnames)) rownames(out.m) <- libnames;
      }
      out.m[i, ] <- as.numeric(o);
      new = FALSE;
   }
   return(out.m);
}

Nonpareil.curve <- function(file,overlap,
			factor=1,plotDispersion=NA,
			xmax=10e12,ymax=1,xmin=1e3,ymin=1e-6,xlab=NULL,ylab=NULL,
			r=NA,g=NA,b=NA,
			new=TRUE,plot=TRUE,libname=NA,modelOnly=FALSE, plotModel=TRUE,
			curve.lwd=2, curve.alpha=0.4, model.lwd=1, model.alpha=1, log='x',
			data.consistency=TRUE, useValue='mean', star=95){
	# Create environment
	Nonpareil.__init_globals(!new);
	# Examine consistency
	if(is.null(file)) stop('The file argument is mandatory');
	if(is.null(overlap)) stop('The overlap argument is mandatory');
	if(!new && (is.null(Nonpareil.LastFactor) || is.null(Nonpareil.LastXmax)))
		stop('No previous plot found, please use new=TRUE if this is the first plot.');
	
	# Read input
	out <- list(kappa=0, C=0, Sstar=0, Rstar=0, modelR2=0)
	a <- read.table(file, sep="\t", h=F);
	for(i in 2:6){
	   a[, i] <- a[, i]^Nonpareil.coverageFactor(overlap);
	}
	if(useValue=='median'){
	   values = a[, 5]
	}else if(useValue=='ub'){
	   values = pmin(1, a[, 2] + a[, 3]*1.9)
	}else if(useValue=='lb'){
	   values = pmax(0, a[, 2] - a[, 3]*1.9)
	}else if(useValue=='q1'){
	   values = a[, 4]
	}else if(useValue=='q3'){
	   values = a[, 6]
	}else{
	   values = a[, 2]
	}
	out$kappa <- values[nrow(a)];
	a$V1 = exp(max(log(a$V1)) + max(values^0.27)*(log(a$V1) - max(log(a$V1))));
	a$V1 = exp(log(a$V1)*0.61 + 10);
	if(is.na(libname)) {
	   libname <- basename(file);
	   if(substr(libname, nchar(libname)-3, nchar(libname))==".npo")
	      libname <- substr(libname, 0, nchar(libname)-4);
	}

	# Check data
	out$C <- max(values);
	if(data.consistency){
	   if(a[0.2*nrow(a), 5]==0){
	      warning('Median of the curve is zero at 20% of the reads, check parameters and re-run (e.g., decrease value of -L in nonpareil).');
	      return(out);
	   }
	   if(a[nrow(a), 2]<=1e-5){
	      warning('Curve too low, check parameters and re-run (e.g., decrease value of -L in nonpareil).');
	      return(out);
	   }
	   if(a[2, 2]>=1-1e-5){
	      warning('Curve too steep, check parameters and re-run (e.g., decrease value of -i in nonpareil).');
	      return(out);
	   }
	   #if(min(a[a[,2]>0, 2])<0.01){
	   #   warning('Curve undefined around S* (kappa=05%), check parameters and re-run (e.g., decrease the value of -i in nonpareil).');
	   #   return(0);
	   #}
	}

	
	# For future calls (of this or other Nonpareil.* functions)
	Nonpareil.LastFactor <<- factor;
	
	# Draw it
	if(plot){
		if(new){
			if(!is.null(xlab)){	xlab <- as.character(xlab) }
			else if(factor==1){	xlab <- 'Sequencing effort (bp)' }
			else if(factor==1e-1){	xlab <- 'Sequencing effort (Dbp)' }
			else if(factor==1e-2){	xlab <- 'Sequencing effort (Hbp)' }
			else if(factor==1e-3){	xlab <- 'Sequencing effort (Kbp)' }
			else if(factor==1e-6){	xlab <- 'Sequencing effort (Mbp)' }
			else if(factor==1e-9){	xlab <- 'Sequencing effort (Gbp)' }
			else if(factor==1e-12){	xlab <- 'Sequencing effort (Tbp)' }
			else if(factor==1e-15){	xlab <- 'Sequencing effort (Pbp)' }
			else if(factor==1e-18){	xlab <- 'Sequencing effort (Ebp)' }
			else if(factor==1e-21){	xlab <- 'Sequencing effort (Zbp)' }
			else{			xlab <- paste('Sequencing effort (by ', factor, 'bp)', sep='') }
			if(is.null(ylab)) ylab <- 'Estimated average coverage';

			plot(1, t='n', xlim=c(xmin, xmax), ylim=c(ymin, ymax), xlab=xlab, ylab=ylab, log=log);
			abline(h=c(1,star/100), lty=2, col='red');
			abline(v=10^seq(0,15,by=3), lty=2, col='gray80')
		}

		if(is.na(r)) r <- sample(200,1);
		if(is.na(g)) g <- sample(200,1);
		if(is.na(b)) b <- sample(200,1);
		if(r>1) r <- r/256;
		if(g>1) g <- g/256;
		if(b>1) b <- b/256;
		
		if(!modelOnly){
			if(!is.na(plotDispersion)){
				if(plotDispersion == 'sd'){
				   err.y <- c(values+a$V3, rev(values-a$V3));
				}else if(plotDispersion == 'ci95'){
				   err.y <- c(values+a$V3*1.9, rev(values-a$V3*1.9));
				}else if(plotDispersion == 'ci90'){
				   err.y <- c(values+a$V3*1.64, rev(values-a$V3*1.64));
				}else if(plotDispersion == 'ci50'){
				   err.y <- c(values+a$V3*.67, rev(values-a$V3*.67));
				}else if(plotDispersion == 'iq'){
				   err.y <- c(a$V4, rev(a$V6));
				}
				polygon(c(a$V1, rev(a$V1))*factor, ifelse(err.y<=ymin*0.1, ymin*0.1, err.y), col=rgb(r,g,b,.2), border=NA);
			}
			lines(a$V1*factor, values, col=rgb(r,g,b,curve.alpha), lwd=curve.lwd);
		}
		
		# Save some info
		Nonpareil.LastXmax   <<- xmax;
		Nonpareil.OpenColors <<- c(Nonpareil.OpenColors, rgb(r,g,b));
		Nonpareil.OpenNames  <<- c(Nonpareil.OpenNames, libname);
	}
	
	# Model it
	sel  <- values>0 & values<0.9;
	x <- a$V1[sel];
	if(length(x)>10){
	   y <- values[sel];
	   data <- list(x=x, y=y)
	   model <- nls(y ~ Nonpareil.f(x, a, b), data=data, weights=(a$V3[sel]^-1.1), start=list(a=1, b=0.1), lower=c(a=0, b=0), algorithm='port', 
			control=nls.control(minFactor=1e-25000, tol=1e-15, maxiter=1024, warnOnly=T));
	   #model <- nls(y ~ Nonpareil.f(x, a, b), data=data, start=list(a=1, b=1), 
	   #		control=nls.control(minFactor=1e-250000000, tol=1e-150000000, maxiter=1024000));
	   if(summary(model)$convInfo$isConv){
	      model.lty=2;
	      if(modelOnly) model.lty=1;
	      if(plot & plotModel){
		 model.x <- exp(seq(log(xmin), log(xmax), length.out=1e3));
		 lines(model.x*factor, predict(model, list(x=model.x)),
			   col=rgb(r,g,b,model.alpha), lty=model.lty, lwd=model.lwd);
		 if(modelOnly) points(max(a$V1), predict(model, list(x=max(a$V1))), col=rgb(r,g,b), pch=21, bg='white');
		 #if(modelOnly) points(min(a$V1[a$V1>0]), predict(model, list(x=min(a$V1[a$V1>0]))), col=rgb(r,g,b), pch=8, bg='white');
	      }
	      pa <- as.numeric(summary(model)$parameters['a', 1])
	      pb <- as.numeric(summary(model)$parameters['b', 1])
	      out$Sstar <- Nonpareil.antif(star/100, pa, pb);
	      out$Rstar <- out$Sstar; #^1.6
	      out$modelR2 <- cor(y, predict(model, list(x=x)));
	   }else{
	      warning('Model didn\'t converge.');
	   }
	}else{
	   warning('Insufficient resolution below 90% coverage, skipping model');
	}
	return(out)
}

Nonpareil.legend <- function(x=NULL, y=.3, ...){
	if(is.null(Nonpareil.LastFactor) || is.null(Nonpareil.LastXmax))
		stop('There must be at least one active plot to draw a legend.  Use Nonpareil.curve().');
	if(is.null(x)) x <- 0.75*Nonpareil.LastXmax;
	legend(x=x, y=y, legend=Nonpareil.OpenNames, fill=Nonpareil.OpenColors, ...);
}

## PREDICTIONS
Nonpareil.coverageFactor <- function(overlap){
   if(overlap==25){
      return(.845);
   }else if(overlap==50){
      return(.757);
   }else if(overlap==75){
      return(.633);
   }else if(overlap==100){
      return(.137);
   }else{
      stop('Unsupported overlap.  Supported values are: 100, 75, 50, and 25.');
   }
}


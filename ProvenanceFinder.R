### FUNCTIONS
small_log_number = -100
library(scales)

get_std_CI <- function(m,M){
	"Std that yields 95% CI = {m, M}"
	std = (M-m)/(2*1.96) 
	std = max(std,0.0001)
	return(std)
}

plot_LSP <- function(source_i,target_samples,unif_samples,source_name,base_prob,low_Y_axis=-1500,high_Y_axis=10,log_prob=T,plot_target=T,return_lik=F){
	# plot Likelihood surface
	sample_Pr = c(base_prob)
	for (i in 1:dim(source_i)[1]){
		a = source_i$Age[i]
		m = source_i$lower_confidence[i]
		M = source_i$upper_confidence[i]
		std=get_std_CI(m,M)
		pdf_sample = dnorm(unif_samples,a,std)
		sample_Pr = sample_Pr + pdf_sample
	}
	sample_Pr=log(sample_Pr/(dim(source_i)[1]+1))
	sample_Pr[sample_Pr<small_log_number]=small_log_number
	if (return_lik==T){
		return( sample_Pr )
	}else{
		if (log_prob==F){
			plot(exp(sample_Pr)~ unif_samples,type="l",main=paste("Likelihood surface plot -", source_name),ylab="Probability",xlab="Time",ylim=c(low_Y_axis,high_Y_axis))
		}else{
			plot(sample_Pr~ unif_samples,type="l",main=paste("Likelihood surface plot -", source_name), ylim=c(low_Y_axis,high_Y_axis),ylab="Log Likelihood",xlab="Time")	
		}
		if (plot_target==T){
			points(x=target_samples,y=rep(0,length(target_samples)),col="darkred",pch=20)
		}
		for (i in 1:(dim(source_i)[1])){
			m = source_i$lower_confidence[i]
			M = source_i$upper_confidence[i]
			segments(x0=m,x1=M,y0=low_Y_axis,y1=low_Y_axis,lwd=4,col="darkblue")
		}
	}
	
}

get_prob_analytical_uncertainty <- function(source_i,ms_m,ms_s,unif_samples,source_name,base_prob,return_prob=F){
	sample_Pr = c(base_prob) # background uniform probability
	for (i in 1:dim(source_i)[1]){
		a = source_i$Age[i]
		m = source_i$lower_confidence[i]
		M = source_i$upper_confidence[i]
		std=get_std_CI(m,M)
		pdf_sample = dnorm(unif_samples,a,std)
		sample_Pr = sample_Pr + pdf_sample
	}
	target_sample_Pr = dnorm(unif_samples,ms_m,ms_s)	
	# integration 
	bin_size= unif_samples[2]-unif_samples[1]
	# integral of combined MN source and N target functions
	# rectangular itegration
	# f_integral = sum(target_sample_Pr*sample_Pr/dim(source_i)[1]+1) * bin_size 
	# trapezoidal integration
	y = target_sample_Pr*sample_Pr/(dim(source_i)[1]+1)
	f_integral = sum((diff(y)*bin_size)/2 + y[1:length(y)-1]*bin_size)
	# PROB UNDER UNIFORM ONLY (no contribution fromsamples in source)
	y_null = target_sample_Pr*base_prob/(dim(source_i)[1]+1)
	f_null_integral = sum((diff(y_null)*bin_size)/2 + y_null[1:length(y_null)-1]*bin_size)
	#
	#
	#temp= sort(sample_Pr/(dim(source_i)[1]+1))* bin_size # probability (sums up to 1)
	#cumsum(temp)
	#print(c("JHGFD",sum(temp) ))
	
	# avoid underflow
	if (f_integral>exp(small_log_number)) {
		h1_prob = log(f_integral) 
	}else{	
		h1_prob = small_log_number
	}
	if (return_prob==T){
		return(list(h1_prob,y))
	}else{
		return(list(h1_prob,log(f_null_integral),y))
	}
}

calcLikMixedModel <- function(w_vec,dim1,dim2,exp_Pr_tbl){
	w_vec=abs(w_vec)
	w= w_vec/sum(w_vec)
	W=rep(w,dim1)
	W <- data.frame(t(matrix(W, ncol = dim1, nrow = dim2)))
	pr=sum(log(apply((exp_Pr_tbl*W),1,FUN=sum)))
	return (-pr)
}

get_likelihood_mixed_model <- function(exp_Pr_tbl,output_summary_file){
	cat("\n\nOptimizing mixed model...")
	dim1= dim(exp_Pr_tbl)[1]
	dim2= dim(exp_Pr_tbl)[2]
	resMixed = optim(fn=calcLikMixedModel,par=c(rgamma(n=dim2,shape=rep(1,dim2))),dim1=dim1,dim2=dim2,exp_Pr_tbl=exp_Pr_tbl)
	likMix = -resMixed$value
	parMix = abs(resMixed$par)/sum(abs(resMixed$par))
	cat("\nlogLikelihood Mixed model",likMix,"\t\n",file= output_summary_file,append=T)
	cat("\nSource\tcontribution\n",file= output_summary_file,append=T)
	for (i in 1:dim2){
		cat(colnames(exp_Pr_tbl)[i],parMix[i],"\n",file= output_summary_file,append=T,sep="\t")
	}
}

run_provenance_estimation <- function(tbl,name_target_sample,output_name="provenance",n_points=10000,plot_PDF=T,plot_CORR=T,calcProb=T,runJK=1000){
	# name output files
	LSP_file            = paste(output_name,"_prob_density_plots.pdf",sep="")
	CORR_file           = paste(output_name,"_source_correlation.pdf",sep="")
	PROB_file           = paste(output_name,"_prob_data.pdf",sep="")
	log_lik_file        = paste(output_name,"_logLikelihood_table.txt",sep="")
	rel_prob_file       = paste(output_name,"_relProbability_table.txt",sep="")
	output_summary_file = paste(output_name,"_Summary.txt",sep="")
	
	obs_time_range = c(min(tbl$lower_confidence), max(tbl$upper_confidence))
	sources = unique(tbl$Source)
	# points delimiting the bins
	unif_samples = seq(obs_time_range[1],obs_time_range[2],length.out=n_points)
	# set baseline uniform probability
	baseline_uniform_pdf = rep(1/diff(obs_time_range),length(unif_samples))

	### ACCOUNT FOR CI IN TARGET SAMPLE
	target_sample = tbl[tbl$Source==name_target_sample,]
	#target_sample = target_sample[order(target_sample$Age),]
	target_sample_m = c()
	target_sample_s = c()

	for (i in 1:dim(target_sample)[1]){
		target_sample_m = c(target_sample_m,target_sample$Age[i])
		m = target_sample$lower_confidence[i]
		M = target_sample$upper_confidence[i]
		target_sample_s =c(target_sample_s, get_std_CI(m,M))
	}	

	# PLOT LIK SURFACE
	if (plot_PDF==T){
		cat("\nGenerating likelihood surface plots...")
		pdf(file=LSP_file,height=6,width=9)
		for (s in sources){	
			plot_LSP(tbl[tbl$Source==s,],target_sample_m,unif_samples,s,low_Y_axis=-15,high_Y_axis=2,log=T,base_prob = baseline_uniform_pdf)
		}
		n <- dev.off()
		cat("done.")
	}
	
	# PLOT CORRELATION BETWEEN SOURCES
	if (plot_CORR==T){
		cat("\nGenerating correlation plots...")
		j=1
		sources_temp = sources[sources != name_target_sample]
		for (s in sources_temp){	
			temp = plot_LSP(tbl[tbl$Source==s,],target_sample_m,unif_samples,s,base_prob = baseline_uniform_pdf,return_lik=T)
			temp = as.data.frame(temp)
			colnames(temp)=s
			if (j==1) {prob_tbl=temp}
			if (j>1) {prob_tbl = cbind(prob_tbl,temp)}
			j= j+1
		}
		pdf(file=CORR_file,height=12,width=12)
		pairs(prob_tbl,pch=19,col="darkblue",lower.panel=NULL,main="Correlation between likelihood surfaces")
		dev.off()
		cat("done.")

		
	}
	
	# CALC LIK of DATA GIVEN EACH SOURCE
	if (calcProb==T){
		cat("\nCalculating likelihood of the data for each source...")
		sources = sources[sources != name_target_sample]
		Pr_tbl = NULL
		Pr_tbl_null = NULL	
		cat("Source\tlogLikelihood\n",file= output_summary_file,append=F)
		pdf(file=PROB_file,height=6,width=9)
		for (s in sources){
			source_i = tbl[tbl$Source==s,]
			Pr_source_i = c()
			Pr_source_i_null = c()
			Pr_for_plot =rep(0,n_points)
			for (j in 1:dim(target_sample)[1]){
				ms_m = target_sample_m[j]
				ms_s = target_sample_s[j]
				Pr_source_i_list = get_prob_analytical_uncertainty(source_i,ms_m,ms_s,unif_samples,s,baseline_uniform_pdf)
				Pr_source_i = c(Pr_source_i, Pr_source_i_list[[1]])
				# Probability of a sample based only on the background uniform PDF, i.e. no contribution from the source's samples
				Pr_source_i_null = c(Pr_source_i_null,Pr_source_i_list[[2]])
				Pr_for_plot = Pr_for_plot + Pr_source_i_list[[3]]
			}
			minProb_plotted = 12
			Pr_for_plot[Pr_for_plot< exp(-minProb_plotted)]=exp(-minProb_plotted)
			logPr_for_plot = log(Pr_for_plot)
			plot(log(Pr_for_plot)~unif_samples,type="n",main=s, ylim=c(-(minProb_plotted+2),-4),xlab="Time",ylab="Log Probability")
			polygon( x =  c(unif_samples,rev(unif_samples)), y = c(logPr_for_plot,rep(-minProb_plotted,length(unif_samples))),col=alpha("darkblue",0.3),border=NA)
			points(x=target_sample_m,y=rep(-(minProb_plotted+1),length(target_sample_m)),col="darkred",pch=20)
			points(x=source_i$Age,y=rep(-(minProb_plotted+1.5),length(source_i$Age)),col="darkblue",pch=20)
			Pr_tbl = cbind(Pr_tbl,Pr_source_i)
			Pr_tbl_null =  cbind(Pr_tbl_null,Pr_source_i_null)
			cat(s, sum(Pr_source_i), "\n",file= output_summary_file,append=T,sep="\t")
		}
		colnames(Pr_tbl)=c(sources)
		dev.off()
		# count number of samples with (likely) unknown source
		#exp_Pr_tbl_temp = exp_Pr_tbl                               # background probability for each source (uniform pdf)
		#for (i in 1:dim(exp_Pr_tbl)[2]){                           # samples with probablity ~ bkg prob for all sources 
		#	exp_Pr_tbl_temp[exp_Pr_tbl_temp<=exp(bkgProb[i])]=0  # are likely to come from an unknown source
		#}
		tbl_null=Pr_tbl*0
		Pr_tbl_null = log(exp(Pr_tbl_null)+exp(Pr_tbl_null)*0.05)
		tbl_null[Pr_tbl>(Pr_tbl_null)]=1
	
		zero_prob = which(apply(tbl_null,FUN=sum,1)==0)
		cat("\n\nNote:",length(zero_prob), "samples are likely to come from unknown sources (5% threshold):\n", zero_prob,file=output_summary_file,append=T)	
		cat("done.\nOutput saved in:",output_summary_file,sep=" ")	
	
	
	
		# save table with log likelihoods per sample
		write.table(Pr_tbl,file=log_lik_file,row.names = F, sep="\t",quote=F)
	
		# save table with relative probabilities
		exp_Pr_tbl = exp(Pr_tbl)
		den = apply(exp_Pr_tbl,1,FUN=sum)
		rel_Pr_tbl = exp_Pr_tbl/den
		write.table(round(rel_Pr_tbl,5),file=rel_prob_file,row.names = F, sep="\t",quote=F)
		cat("done.\nOutput saved as:",output_summary_file,log_lik_file, rel_prob_file,sep=" ")
		# run mixed model 
		get_likelihood_mixed_model(exp(Pr_tbl),output_summary_file)
		cat("done.\nOutput saved in:",output_summary_file,sep=" ")	
	}


	# RUN JACK KNIFE
	if (runJK > 0){
		cat("\nJackknifing data...\n")
		sources = sources[sources != name_target_sample]
		min_samples_per_source = min(table(tbl$Source))
		Pr_tbl = NULL
		if (calcProb==T){
			cat("\nSource\tJackknife support\n",file= output_summary_file,append=T)
		}else{
			cat("\nSource\tJackknife support\n",file= output_summary_file,append=F)
		}
		for (bs_rep in 1:runJK){
			if (bs_rep %% 10==0){print(bs_rep)}
			# sub-sample also target sample
			indx_target_samples = sort(sample(1:dim(target_sample)[1],min_samples_per_source))
			prob_per_replicate = c()
			for (s in sources){
				source_i_complete = tbl[tbl$Source==s,]
				source_i = source_i_complete[sort(sample(1:dim(source_i_complete)[1],min_samples_per_source)),]
				Pr_source_i = c()
				for (j in indx_target_samples){
					ms_m = target_sample_m[j]
					ms_s = target_sample_s[j]
					Pr_source_i_list = get_prob_analytical_uncertainty(source_i,ms_m,ms_s,unif_samples,s,baseline_uniform_pdf)
					Pr_source_i = c(Pr_source_i, Pr_source_i_list[[1]])
				}
				prob_per_replicate = c(prob_per_replicate,sum(Pr_source_i))
			}
		Pr_tbl = rbind(Pr_tbl, prob_per_replicate)
		
		}
	Pr_tbl = as.data.frame(Pr_tbl)
	JK_res =apply(Pr_tbl,1,which.max)
	for (j in 1:length(sources)){
		cat(sources[j],length(JK_res[JK_res==j])/runJK,"\n",sep="\t",file =output_summary_file,append=T)
	}
	cat("done.\nOutput saved in:",output_summary_file,sep=" ")	
	}



}




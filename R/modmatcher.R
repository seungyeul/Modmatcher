gender_sig<-function(mat=mat,gender=ann_gender){
	### Gender check based on molecular profile
	### mat: data matrix with row of features and with columns of the samples
	### gender: a vector of number (1 female, 2=male) indicating annotated genders for the samples in the data matrix,"mat". Order should be matched with the column name.
	female<-mat[,which(gender==1)]
	f_na<-apply(female,1,function(i) length(which(is.na(i))))
	male<-mat[,which(gender==2)]
	m_na<-apply(female,1,function(i) length(which(is.na(i))))
	na_ind<-which(f_na>0.7*dim(female)[[2]]|m_na>0.7*dim(female)[[2]])
	if(length(na_ind)>0){
		female<-female[-(na_ind),]
		male<-male[-(na_ind),]
	}
	f_mean<-apply(female,1,function(x) mean(x,na.rm=T))
	m_mean<-apply(male,1, function(x) mean(x,na.rm=T))
	diff<-f_mean-m_mean
	t.pva<-sapply(seq.int(length(female[,1])), function(i) t.test(as.numeric(female[i,]),as.numeric(male[i,]))$p.value)
	t.qva<-p.adjust(t.pva,method="BH",n=length(t.pva))
	gen_tab<-cbind(rownames(female),f_mean,m_mean,diff,t.pva)
	colnames(gen_tab)<-c("Name","Female_mean","Male_mean","Female-Male","p-value")
	len_rep<-ifelse(length(which(t.pva<0.001))<20,length(which(t.pva<0.001)),20)
	gen_tab<-gen_tab[sort.list(t.pva),]
	gen_tab<-gen_tab[1:len_rep,]
	return(gen_tab)
}


gender_check<-function(mat=mat,feature1=gene1,feature2=gene2,gender=ann_gender){
	### Plot for the two features separated between male and female
	feat1<-as.numeric(mat[feature1,])
	feat2<-as.numeric(mat[feature2,])
	sam_col<-ifelse(gender==1,"red",ifelse(gender==2,"blue","grey"))
	plot(feat1,feat2,xlim=range(feat1,na.rm=T),ylim=range(feat2,na.rm=T),xlab=feature1,ylab=feature2,col=sam_col,cex.lab=1.5,cex.axis=1.2)
	if(length(table(gender))==3){
	  legend("top",c("Female","Male","Unknown"),pch=1,col=c("red","blue","grey"))
	} else {
	  legend("top",c("Female","Male"),pch=1,col=c("red","blue"))
	}
}


cis_pair<-function(data1=mat1,data2=mat2){
	# Cis gene identification based on correlation test
	# The two matrix should have the same number of features (i.e genes) and samples
	com_gene<-rownames(data1)
	cis_pair<-sapply(seq.int(length(com_gene)),function(i) cor.test(as.numeric(data1[i,]),as.numeric(data2[i,])))
	cis_rho<-sapply(seq.int(length(com_gene)), function(i) cis_pair[,i]$estimate)
	cis_pva<-sapply(seq.int(length(com_gene)), function(i) cis_pair[,i]$p.value)
	cis_qva<-p.adjust(cis_pva,method="BH",n=length(cis_pva))
	cis_tab<-cbind(com_gene,cis_rho,cis_pva,cis_qva)
	colnames(cis_tab)<-c("GeneName","rho","p-value","q-value")
	cis_tab<-cis_tab[sort.list(cis_qva),]
	return(cis_tab)
}


sam_score<-function(data1=mat1,data2=mat2,type1="Type1",type2="Type2"){
	# Generation of sample similarity score based on ranked-transformed values of selected cis-genes.
	# The direction of the correlation should be considered before calling the function
	rdata1<-t(apply(data1,1,function(x) rank(x,na.last="keep")))
	rdata2<-t(apply(data2,1,function(x) rank(x,na.last="keep")))
	score<-t(sapply(seq.int(dim(rdata1)[[2]]), function(i) sapply(seq.int(dim(rdata2)[[2]]), function(j) cor.test(as.numeric(rdata1[,i]),as.numeric(rdata2[,j]),use="pairwise.complete.obs")$estimate)))
	rownames(score)<-colnames(data1)
	colnames(score)<-colnames(data2)
	return(score)
}


self_align<-function(score=mat,type1="Type1",type2="Type2",top=5,rnd=1){
	# Determination of cutoff based on z-transformed distribution of scores
	# Top 5% or top 2.5% (ci; 5 or 2.5)
	cut1<-rep(0,length(rownames(score)))
	cut2<-cut1
	ci<-ifelse(top==5,1.645,1.960)
	for(i in 1:length(cut1)){
		sim_row<-as.numeric(score[i,])
		fis_mat<-0.5*(log(1+sim_row)-log(1-sim_row))
		cut<-mean(fis_mat)+ci*sd(fis_mat)
		cut1[i]<-length(which(fis_mat>cut))
		sim_col<-as.numeric(score[,i])
		fis_mat<-0.5*(log(1+sim_col)-log(1-sim_col))
		cut<-mean(fis_mat)+ci*sd(fis_mat)
		cut2[i]<-length(which(fis_mat>cut))
	}
	topN1<-round(mean(cut1))
	topN2<-round(mean(cut2))

        # Initial alignment for the pre-aligned pairs
	sam1<-rownames(score)
	sam2<-colnames(score)
	diff1<-rep(0,length(sam1))
	diff2<-diff1
	secd1<-rep(0,length(sam2))
	secd2<-secd1
	self_tab<-matrix(NA,nr=length(sam1),nc=8)
	colnames(self_tab)<-c("Sample1","Sample2","selfscore","cut1","cut2","rank1","rank2","Note")
	for(i in 1:length(sam1)){
		name1<-sam1[i]
		name2<-sam2[i]
		self<-score[i,i]
		row_scr<-as.numeric(score[i,])
		scr1<-sort(row_scr,decreasing=T)
		mean1<-mean(scr1)
		sec1<-scr1[2]
		diff1[i]<-self-mean1
		secd1[i]<-self-sec1
		cut1<-scr1[topN1]
		rank1<-which(scr1==self)
		col_scr<-as.numeric(score[,i])
		scr2<-sort(col_scr,decreasing=T)
		mean2<-mean(scr2)
		sec2<-scr2[2]
 		cut2<-scr2[topN2]
		rank2<-which(scr2==self)
		diff2[i]<-self-mean2
		secd2[i]<-self-sec2
		sam_col<-rep("grey",length(sam1))
		sam_col[i]<-"red"
		if(self>=cut1&self>=cut2){
			note<-"Ok"
		} else if(self>=cut1|self>=cut2){
			note<-"Poor"
		} else {
			note<-"No"
			row_ind<-which(row_scr==max(scr1))
			col_ind<-which(col_scr==max(scr2))
			sam_ind<-c(row_ind,col_ind)
			no_sam<-c(sam1[row_ind],sam2[col_ind])
			sam_col[sam_ind]<-"blue"
		}
		pdf(paste(type1,"_",type2,"_Round_",rnd,"/Pairwise/",type1,"_",name1,"_",type2,"_",name2,".pdf",sep=""))
      		plot(row_scr,col_scr,col=sam_col,xlab=paste(type1,": ",name1," -> ",type2,": All",sep=""),ylab=paste(type2,": ",name2," -> ",type2,": All",sep=""),main=paste(type1,": ",name1," - ",type2,": ",name2,sep=""))
		if(note=="No"){
			text(row_scr[sam_ind]-0.1,col_scr[sam_ind]-0.02,no_sam)
		}
		dev.off()
		self_tab[i,]<-c(name1,name2,self,cut1,cut2,rank1,rank2,note)
	}
	pdf(paste(type1,"_",type2,"_Round_",rnd,"/Compare_self_score_mean.pdf",sep=""))
	plot(diff1,diff2,xlab=paste(type1,"->",type2),ylab=paste(type2,"->",type1),main="Self_score - mean(other)",cex.axis=1.2,cex.lab=1.5,cex.main=2)
	abline(0,1,lty=2,col="red")
	dev.off()
	pdf(paste(type1,"_",type2,"_Round_",rnd,"/Compare_self_score_second_best.pdf",sep=""))
	plot(secd1,secd2,xlab=paste(type1,"->",type2),ylab=paste(type2,"->",type1),main="Self_score - second_best_score",cex.axis=1.2,cex.lab=1.5,cex.main=2)
	abline(0,1,lty=2,col="red")
	dev.off()
	return(self_tab)
}


pairwise_alignment<-function(data1=mat1,data2=mat2,type1="Type1",type2="Type2",rnd=round){
	# Gene the common genes and initial samples for the identification of cis-gene
	dir.create(paste(type1,"_",type2,"_Round_",rnd,sep=""))
	dir.create(paste(type1,"_",type2,"_Round_",rnd,"/Pairwise/",sep=""))
	if(is.null(rnd) | rnd==1){
		sam1<-intersect(colnames(data1),colnames(data2))
		sam2<-sam1
	} else {
		pre_align<-read.table(paste(type1,"_",type2,"_Round_",(rnd-1),"/Alignment_pairs.txt",sep=""),header=T,sep="\t")
		sam1<-as.character(pre_align[,1])
		sam2<-as.character(pre_align[,2])
	}
	sam_pair<-paste(sam1,sam2)
	# Prepare the data matrix for cis gene identification
	gene1<-rownames(data1)
	gene2<-rownames(data2)
	com_gene<-intersect(gene1,gene2)
	ndata1<-data1[match(com_gene,gene1),match(sam1,colnames(data1))]
	ndata2<-data2[match(com_gene,gene2),match(sam2,colnames(data2))]
	odata1<-data1[match(com_gene,gene1),which(!(colnames(data1)%in%sam1))]
	odata2<-data2[match(com_gene,gene2),which(!(colnames(data2)%in%sam2))]

	# Cis gene identification based on correlation test
	cis_tab<-cis_pair(ndata1,ndata2)
	write.table(cis_tab,paste(type1,"_",type2,"_Round_",rnd,"/",type1,"_",type2,"_cis_pairs.txt",sep=""),row.names=F,col.names=T,quote=F,sep="\t")

	# Choose cis-gene up to 1000 or feautures more than 4 times of samples or
	sel_num<-max(2.5*min(dim(data1)[[2]],dim(data2)[[2]]),min(1000,4*min(dim(data1)[[2]],dim(data2)[[2]]),length(which(as.numeric(cis_tab[,4])<0.01)),length(which(as.numeric(cis_tab[,2])>0.4))))
	print(paste("Round",rnd,":",sel_num,"selected_cis_genes"))
	cis_tab<-cis_tab[1:sel_num,]
	cis_dir<-as.numeric(cis_tab[,2])
	cis_gene<-as.character(cis_tab[,1])

	sam_1<-c(sam1,colnames(odata1))
	sam_2<-c(sam2,colnames(odata2))

	# Generating scoreing matrix
	sdata1<-data1[match(cis_gene,gene1),match(sam_1,colnames(data1))]
	sdata2<-data2[match(cis_gene,gene2),match(sam_2,colnames(data2))]
	if(length(which(cis_dir<0))>0){
	  sdata2[which(cis_dir<0),]<- -(sdata2[which(cis_dir<0),])
	}
	score<-sam_score(sdata1,sdata2,type1,type2)
	write.table(score,paste(type1,"_",type2,"_Round_",rnd,"/",type1,"_",type2,"_score.txt",sep=""),row.names=T,col.names=T,quote=F,sep="\t")

	# Best matched pairs
	match1<-as.numeric(apply(score,1,function(x) which(x==max(x))))
	match2<-as.numeric(apply(score,2,function(x) which(x==max(x))))
	best1<-rep(0,length(match1))
	best2<-rep(0,length(match2))

	# self_alignment
	self_tab<-self_align(score[1:length(sam1),1:length(sam2)],type1,type2,5,rnd)
	write.table(self_tab,paste(type1,"_",type2,"_Round_",rnd,"/",type1,"_",type2,"_self_align.txt",sep=""),row.names=F,col.names=T,quote=F,sep="\t")
	ok_sam<-self_tab[which(self_tab[,8]=="Ok"),]
	sam1_ok<-ok_sam[,1]
	sam2_ok<-ok_sam[,2]
	best1[match(sam1_ok,rownames(score))]<-match(sam2_ok,colnames(score))
	best2[match(sam2_ok,colnames(score))]<-match(sam1_ok,rownames(score))

	# reciprocal alignment
	no_sam1<-sam_1[!(sam_1%in%sam1_ok)]
	no_sam2<-sam_2[!(sam_2%in%sam2_ok)]
	if(length(no_sam1)>0){
		for(i in 1:length(no_sam1)){
			name1<-no_sam1[i]
			ind1<-which(sam_1==name1)
			part1<-match1[ind1]
			name2<-sam_2[part1]
			ind2<-which(sam_2==name2)
			part2<-match2[ind2]
			if(ind1==part2&name2%in%no_sam2){
				best1[ind1]<-ind2
				best2[ind2]<-ind1
				print(paste(type1,":",name1,"->",type2,":",name2))
			}
		}
	}

	# Summary of pairwise alignment
	final1<-sam_1[best1>0]
	final2<-sam_2[best1[best1>0]]
	final_pair<-paste(final1,final2)
	flag<-ifelse(final1==final2,"Self-aligned","Cross-aligned")
	final_tab<-cbind(final1,final2,flag)
	colnames(final_tab)<-c(type1,type2,"Note")
	write.table(final_tab,paste(type1,"_",type2,"_Round_",rnd,"/Alignment_pairs.txt",sep=""),row.names=F,col.names=T,quote=F,sep="\t")
	token<-ifelse(length(sam_pair)==length(intersect(final_pair,sam_pair)),1,0)
	return(token)
}


quality_evaluation<-function(data1=mat1,data2=mat2,type1="Type1",type2="Type2",rnd=1){
	# Draw plot for alignment quality
	if(rnd==1){
		stop("Round should be greater than 1")
	}
	num_sam<-rep(0,rnd)
	num_cis<-rep(0,rnd)
	for(i in 1:(rnd-1)){
		cis_tab<-read.table(paste(type1,"_",type2,"_Round_",i,"/",type1,"_",type2,"_cis_pairs.txt",sep=""),header=T,sep="\t")
		num_cis[i]<-length(which(cis_tab$q.value<0.01))
		if(i == 1){
			num_sam[i]<-length(intersect(colnames(data1),colnames(data2)))
			cis_first<-cis_tab
		} else {
			pre_align<-read.table(paste(type1,"_",type2,"_Round_",(i-1),"/Alignment_pairs.txt",sep=""),header=T,sep="\t")
			num_sam[i]<-dim(pre_align)[[1]]
		}
	}
	num_sam[rnd]<-num_sam[rnd-1]
	num_cis[rnd]<-num_cis[rnd-1]

	# Relationship between the number of sample and the number of cis genes along the sample alignment rounds
	pdf(paste(type1,"_",type2,"_sample_alignment_round_change.pdf",sep=""))
	par(mar=c(5, 4, 4, 5) + 0.1)
	plot(seq(1,rnd),num_cis,pch=16,xlab="",ylab="",axes=F,type="b",col="red",main="Sample alignment quality")
	axis(4,col="red",col.axis="red",las=0)
	mtext("Number of cis associated genes",col="red",side=4,line=2.5)
	box()

	par(new=T)

	plot(seq(1,rnd),num_sam,pch=15,xlab="",ylab="",axes=F,type="b",col="blue")
	mtext("Numbe of aligned sample",col="blue",side=2,line=4)
	axis(2,col="blue",col.axis="blue",las=0)

	axis(1,1:rnd,labels=as.character(c(0:(rnd-1))))
	mtext("Alignment round",side=1,col="black",line=2.5)
	dev.off()

	# Comparison of cis strength from the first round to the last
	cis_tab<-cis_tab[match(cis_first[,1],cis_tab[,1]),]
	first_qva<- -log10(cis_first[,4])
	last_qva<- -log10(cis_tab[,4])

	pdf(paste(type1,"_",type2,"_sample_alignment_cis_compare.pdf",sep=""))
	plot(first_qva,last_qva,xlab="-log10(FDR) - First round",ylab="-log10(FDR) - Last round",main="Cis assoication strength")
	abline(0,1,lty=2,col="red")
	dev.off()
}

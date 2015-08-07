#!/usr/bin/env Rscript

topics <- read.table('0100_topics.txt',sep='\t',header=FALSE)
vocab <- read.table('ap/vocab.txt',header=FALSE,stringsAsFactors=FALSE)

for(k in 1:nrow(topics)) {
	N <- 20
	top_words <- order(topics[k,],decreasing=TRUE)[1:N]
	cat('\n=====================\n')
	cat(sprintf('Topic %3d\n',k))
	cat(paste(vocab[top_words,1],topics[k,top_words],sep='\t',collapse='\n'))
	cat('\n')
}


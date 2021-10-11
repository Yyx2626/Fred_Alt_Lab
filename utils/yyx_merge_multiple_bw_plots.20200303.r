
options(stringsAsFactors=FALSE)


require("rtracklayer")


`%.%` = function(x,y) paste0(x,y)

join = function(sep, vec, ...) paste(collapse=sep, c(vec, ...))
cat0 = function(...) cat(sep="", ...)

echo_str <- function(x, sep=" =\t", collapse=", ") deparse(substitute(x)) %.% sep %.% join(collapse, x)
echo <- function(x, sep=" =\t", collapse=", ") cat0(deparse(substitute(x)), sep, join(collapse, x), "\n");


merge_bw_plot = function(bw_filenames, chr, start, end, num_windows, y_max=NA, mean_col="gray", se_col="black", xaxis=FALSE, yaxis=TRUE, should_annotate_V=FALSE, name=""){
#	bw_input_list = list()
	coordinate_range = chr %.% ":" %.% start %.% "-" %.% end
	cat("  coordinate_range = " %.% coordinate_range %.% "\n")
 
	cat("Request overlapping genomic coordinate range =\t" %.% chr %.% " : " %.% start %.% " - " %.% end %.% "\n")
	cat("\tcount in " %.% num_windows %.% " continuous even-length non-overlapping windows\n")
	coordinate_GRange = GRanges(seqnames=chr, ranges=IRanges(start=as.integer(start), end=as.integer(end)))

	tiles = tile(range(coordinate_GRange), n=num_windows)
	tiles.gr = unlist(tiles)

	NF = length(bw_filenames)
	for(repidx in 1:NF){
		cat("Now read and parse " %.% bw_filenames[repidx] %.% " ...\n");
		now_input = import(bw_filenames[repidx], selection=BigWigSelection(coordinate_GRange))
	
		### ref: https://support.bioconductor.org/p/100222/
		hits = findOverlaps(tiles.gr, now_input)
	#	agg = aggregate(now_input, hits, score=max(score))   # Error: may be less than 5000 !?
	#	mcols(tiles.gr)[["rep" %.% repidx]] = agg$score
		mcols(tiles.gr)[["rep" %.% repidx]] = 0
		for(qi in 1:num_windows){
			idx = (queryHits(hits) == qi)
			if(sum(idx) > 0){
				now_score_vec = now_input$score[subjectHits(hits)[idx]]
				if((!is.na(y_max) && y_max < 0) || all(now_score_vec <= 0)){
					mcols(tiles.gr)[["rep" %.% repidx]][qi] = min(now_score_vec, na.rm=TRUE)
				}else{
					mcols(tiles.gr)[["rep" %.% repidx]][qi] = max(now_score_vec, na.rm=TRUE)
				}
			}
		}
	}

	mcols(tiles.gr)$mid = (start(tiles.gr)+end(tiles.gr))/2
	if(NF > 1){
		mcols(tiles.gr)$mean = apply(mcols(tiles.gr)[, "rep" %.% 1:NF], 1, mean)
		mcols(tiles.gr)$sd = apply(mcols(tiles.gr)[, "rep" %.% 1:NF], 1, sd)
		mcols(tiles.gr)$se = mcols(tiles.gr)$sd / sqrt(NF)
	}else{
		mcols(tiles.gr)$mean = mcols(tiles.gr)$rep1
		mcols(tiles.gr)$sd = 0
		mcols(tiles.gr)$se = 0
	}
	mcols(tiles.gr)$meanPlusSe = mcols(tiles.gr)$mean + mcols(tiles.gr)$se
	if((!is.na(y_max) && y_max < 0) || all(mcols(tiles.gr)$mean <= 0)){   #  <-- 2020-02-28 debugged !!
		mcols(tiles.gr)$meanPlusSe = mcols(tiles.gr)$mean - mcols(tiles.gr)$se
	}
	head(tiles.gr)

	cat("Range of signal mean =\t" %.% min(mcols(tiles.gr)$mean, na.rm=TRUE) %.% " - " %.% max(mcols(tiles.gr)$mean, na.rm=TRUE) %.% "\n")
	cat("Range of signal mean + se =\t" %.% min(mcols(tiles.gr)$meanPlusSe, na.rm=TRUE) %.% " - " %.% max(mcols(tiles.gr)$meanPlusSe, na.rm=TRUE) %.% "\n\n")

	if(is.na(y_max)){
		y_max = max(mcols(tiles.gr)$meanPlusSe, na.rm=TRUE)
		if(all(mcols(tiles.gr)$mean <= 0)){
			y_max = min(mcols(tiles.gr)$meanPlusSe, na.rm=TRUE)
		}
	}
	x_lim = c(start(ranges(coordinate_GRange)), end(ranges(coordinate_GRange)))
	if(should_reverse_coord){
		x_lim = range(-x_lim)
	}
	adjust_x_lim = (x_lim - mean(x_lim)) * 1.01 + mean(x_lim)
	y_lim = range(c(0, y_max))
	adjust_y_lim = (y_lim - mean(y_lim)) * 1.01 + mean(y_lim)
	plot(0, type="n", xlim=adjust_x_lim, ylim=adjust_y_lim, axes=FALSE, xlab="", ylab="Signal", xaxs="i", yaxs="i")
	if(xaxis)  axis(1, xpd=NA)
	if(yaxis)  axis(2, xpd=NA)
	if(should_reverse_coord){
		if(should_annotate_V){
			with(V_family_range[-1,], rect(-end, y_max*1.05, -start, y_max*1.15, col=hsv((0:3)/4, 0.2, 1), border=NA, xpd=NA))
			with(all_V_locations_DF, rect(-end, y_max*1.05, -start, y_max*1.15, col="black", border=NA, xpd=NA))
#			text(mean(adjust_x_lim), y_max*1.04, adj=c(0.5,1), name)
		}
		rect(-end(tiles.gr), mcols(tiles.gr)$mean, -start(tiles.gr), mcols(tiles.gr)$meanPlusSe, col=se_col, border=NA)
		rect(-end(tiles.gr), 0, -start(tiles.gr), mcols(tiles.gr)$mean, col=mean_col, border=NA)
	}else{
		if(should_annotate_V){
			with(V_family_range[-1,], rect(start, y_max*1.05, end, y_max*1.15, col=hsv((0:3)/4, 0.2, 1), border=NA, xpd=NA))
			with(all_V_locations_DF, rect(start, y_max*1.05, end, y_max*1.15, col="black", border=NA, xpd=NA))
#			text(mean(adjust_x_lim), y_max*1.04, adj=c(0.5,1), name)
		}
		rect(start(tiles.gr), mcols(tiles.gr)$mean, end(tiles.gr), mcols(tiles.gr)$meanPlusSe, col=se_col, border=NA)
		rect(start(tiles.gr), 0, end(tiles.gr), mcols(tiles.gr)$mean, col=mean_col, border=NA)
	}
}





if(!interactive()){
	args = commandArgs(TRUE)
}
if(length(args) < 12){
	stop("no enough command-line arguments

Usage: Rscript this.r
	<output_prefix> <pdf_width> <pdf_height>
	<cex> <xxx> <xxx> <se_col>
	<should_reverse_coord>
	<each_V.bed> <V_family.bed> <target.bed> <num_windows>
	<mean_col_1:ymax_1:1.1.bw,1.2.bw,...>
	[mean_col_2:ymax_2:2.1.bw,2.2.bw,...] ...

Version: 0.2.0 (2020-03-03)
Author: Adam Yongxin Ye @ BCH

")
}

output_prefix = args[1]
pdf_width = as.numeric(args[2])
pdf_height = as.numeric(args[3])
cex = as.numeric(args[4])
y_max = as.numeric(args[5])
mean_col = args[6]
se_col = args[7]
should_reverse_coord = as.logical(args[8])
each_V_bed_filename = args[9]
V_family_bed_filename = args[10]
target_bed_filename = args[11]
num_windows = as.integer(args[12])
bw_file_groups = args[13:length(args)]

echo(output_prefix)
echo(y_max)
echo(each_V_bed_filename)
echo(V_family_bed_filename)
echo(target_bed_filename)
echo(num_windows)


all_V_locations_DF = read.delim(each_V_bed_filename, header=FALSE)
names(all_V_locations_DF) = c("chr", "start", "end", "name")


V_family_range = read.delim(V_family_bed_filename, header=FALSE)
V_family_range = V_family_range[, 1:4]
names(V_family_range) = c("chr", "start", "end", "family")
tmp_list = strsplit(V_family_range$family, "_")
V_family_range$family = sapply(tmp_list, function(one) one[1])

all_V_DF = data.frame(chr=V_family_range$chr[1], start=min(V_family_range$start), end=max(V_family_range$end), family="all_V")
V_family_range = rbind(all_V_DF, V_family_range)

V_family_range$family2 = sub("/", "_", V_family_range$family)


target_range = read.delim(target_bed_filename, header=FALSE)
names(target_range) = c("chr", "start", "end", "family")
tmp_list = strsplit(target_range$family, "_")
target_range$family = sapply(tmp_list, function(one) one[1])

all_V_DF = data.frame(chr=target_range$chr[1], start=min(target_range$start), end=max(target_range$end), family="all_V")
target_range = rbind(all_V_DF, target_range)

target_range$family2 = gsub("[ /()]", "_", target_range$family)


for(i in 1:nrow(target_range)){
	chr = target_range$chr[i]
	start = target_range$start[i]
	end = target_range$end[i]
	coordinate_range = chr %.% ":" %.% start %.% "-" %.% end
	output_filename = output_prefix %.% "." %.% target_range$family2[i] %.% "." %.% sub("[:-]", "_", coordinate_range) %.% ".pdf"
	cat("V family = " %.% target_range$family[i] %.% "\n")
	cat("  coordinate_range = " %.% coordinate_range %.% "\n")
 
	cat("Now output pdf " %.% output_filename %.% " ...\n");
	cat("\tset cex =\t" %.% cex %.% "\n") 

	NG = length(bw_file_groups)
	pdf(output_filename, width=pdf_width, height=pdf_height)
	par(mar=c(0,1.8,0,0), mgp=c(2,0.8,0), cex=cex, mfrow=c(NG, 1), oma=c(3,0,3,0))

	for(k in 1:NG){
		now_bw_file_group = bw_file_groups[k]
		tmp = strsplit(now_bw_file_group, ":", fixed=TRUE)[[1]]
		mean_col = tmp[1]
		y_max = as.numeric(tmp[2])
		now_bw_filenames = strsplit(tmp[3], ",", fixed=TRUE)[[1]]
		cat("merge_bw_plot for " %.% k %.% "-th group of " %.% length(now_bw_filenames) %.% " input files: " %.% join(", ", now_bw_filenames) %.% "\n")
		cat("\tmean color =\t" %.% mean_col %.% "\n\tSE color =\t" %.% se_col %.% "\n") 
		cat("\tset y_max =\t" %.% y_max %.% "\n")

		merge_bw_plot(now_bw_filenames, chr, start, end, num_windows, y_max=y_max, mean_col=mean_col, se_col=se_col, xaxis=(k==NG), yaxis=TRUE, should_annotate_V=(k==1))
	}

	dev.off()

}   # end for(i in 1:nrow(target_range))

cat("All done. Congratulations!\n");

library(ape)
library(scales)

### Output directory structure
### /scratch/aubkbk001_01_slim/
### |__ slim_output_files_m${m}_r${r}_p${p}
### |__ sample_vcf_files_m${m}_r${r}_p${p}
### |__ sample_fasta_files_m${m}_r${r}_p${p} **
### |__ sample_id_list_m${m}_r${r}_p${p}.txt
### |__ vcf_file_list_m${m}_r${r}_p${p}.txt

### Loop over different settings for mutation and recombination rates and population sizes
## to summarize output
pdf('/scratch/aubkbk001_01_slim/all_roh_output.pdf', width=15, height=6)
## final settings
for(m in c('5e-07')){
  for(r in c('1e-8')){
    for(p in c('500')) {
      ## select parameter set directory
      setwd(paste0('/scratch/aubkbk001_01_slim/sample_fasta_files_m',m,'_r',r,'_p',p,'/')) ##  
      print(paste0('m',m,'_r',r,'_p',p,'...'))

      ## get list of files in the directory to test
      fns <- list.files()
      nums <- do.call(rbind, strsplit(fns, split='_'))[,1]
      nums <- gsub('i', '', nums)
      nums <- unique(nums) ## unique list of individual IDs
      nums <- as.numeric(nums)
      
      for(n in nums){
        print(n)
        f.1 <- list.files('./', pattern=paste0('i',n,'_1'))
        f.2 <- list.files('./', pattern=paste0('i',n,'_2'))
        
        seq.1 <- read.FASTA(f.1)
        seq.2 <- read.FASTA(f.2)
        names(seq.2) <- '2'
        
        temp <- cbind(seq.1$`1`, seq.2$`2`)
        rownames(temp) <- c(1:length(seq.1$`1`))
        print('Sequences read in...')
        ## heterozygous positions
        locs <- as.numeric(rownames(temp[which(temp[,1] != temp[,2]),]))
        print(paste0(length(locs),' het sites...'))
        
        OUT <- NULL
        for(i in 1:length(locs)){
          if(i %% 25000 == 0){
            print(paste0(i,' - ',Sys.time()))
            write.table(OUT, paste0('/scratch/aubkbk001_01_slim/roh_results_m',m,'_r',r,'_p',p,'.txt'),
                        sep='\t', row.names=FALSE, quote=FALSE, append=TRUE, col.names=FALSE)
            OUT <- NULL
          }
          if(i == 1 & locs[i] != 1){ ## if it's the first heterozygous site on the chromosome and it's not the first position,
            s <- 1 ## save the first position as the ROH start,
            e <- locs[i]-1 ## and save the position just before the heterozygous site as the ROH end.
            save <- c(n, s, e, e-s+1) ## and save the information.
            OUT <- rbind(OUT, save)
            if(locs[i]+1 != locs[i+1]){ ## if the next site isn't heterozygous
              s <- locs[i]+1 ## save the next site as the next ROH start
            }
          }
          else if(i == 1 & locs[i] == 1){ ## if it's the first heterozygous site AND the first position,
            next ## just go to the next heterozygous site and record nothing
          }
          else if(i != 1 & i != length(locs)){ ## if it's not the first or last heterozygous site on the chrom,
            if(locs[i-1] != (locs[i] - 1) & locs[i+1] != (locs[i] + 1)){ ## and if it's not in the middle of a run of heterozygosity
              e <- locs[i]-1 ## save the position just before the heterozygous site as a ROH end,
              save <- c(n, s, e, e-s+1) ## and save the information.
              OUT <- rbind(OUT, save)
              s <- locs[i]+1 ## and start a new ROH at the next position.
            }
            else if(locs[i-1] == (locs[i] - 1) & locs[i+1] == (locs[i] + 1)){ ## and if it's in the middle of a run of heterozygosity,
              next ## just go to the next site
            }
            else if(locs[i-1] != (locs[i] - 1) & locs[i+1] == (locs[i] + 1)){ ## and if it's at the beginning of a run of heterozygosity,
              e <- locs[i]-1 ## save the position just before the heterozygous site as a ROH end,
              save <- c(n, s, e, e-s+1) ## and save the information.
              OUT <- rbind(OUT, save)
            }
            else if(locs[i-1] == (locs[i] - 1) & locs[i+1] != (locs[i] + 1)){ ## and if it's at the end of a run of heterozygosity,
              s <- locs[i]+1 ## start a new ROH.
            }
            else if(i == length(locs)){ ## if it's the last site and heterozygous,
              e <- locs[i]-1 ## save the position just before the heterozygous site as a ROH end,
              save <- c(n, s, e, e-s+1) ## and save the information.
              OUT <- rbind(OUT, save)
            }
          }
          if(i == length(locs)){ ## if it's the last heterozygous site on the chrom,
            if(locs[i] != 30e6){ ## and it's not the last chrom position,
              e <- locs[i]-1 ## save the position just before the heterozygous site as a ROH end,
              save <- c(n, s, e, e-s+1) ## and save the information.
              OUT <- rbind(OUT, save)
              s <- locs[i]+1 ## and save the information for the last ROH on the chrom
              e <- 30e6
              save <- c(n, s, e, e-s+1) ## and save the information.
              OUT <- rbind(OUT, save)
            }
            if(locs[i] == 30e6){ ## and it's the last chrom position,
              e <- locs[i]-1 ## save the position just before the heterozygous site as a ROH end,
              save <- c(n, s, e, e-s+1) ## and save the information.
              OUT <- rbind(OUT, save)
            }
          }
        }
        
        # lens <- OUT[,3] - OUT[,2]
        print('ROH lengths calculated...')
        print(paste0('... ',Sys.time()))
        write.table(OUT, paste0('/scratch/aubkbk001_01_slim/roh_results_m',m,'_r',r,'_p',p,'.txt'),
                    sep='\t', row.names=FALSE, quote=FALSE, append=TRUE, col.names=FALSE)
      }
      
      print('Results written...')
      ## visualize ROH distributions
      rohs <- read.table(paste0('/scratch/aubkbk001_01_slim/roh_results_m',m,'_r',r,'_p',p,'.txt'), 
                         sep='\t', header=FALSE)
      colnames(rohs) <- c('samp','start','end')
      rohs$length <- rohs$end - rohs$start
      frohs.all <- NULL
      for(i in unique(rohs$samp)){
        frohs.all <- c(frohs.all, sum(rohs[rohs$samp == i, 'length'])/30e6)
      }
      print('Beginning writing plots...')
      par(mfrow=c(1,3))
      ## plot 1
      hist(frohs.all, main=paste0('m ',m,' - r',r,' - p',p), xlab='f(ROH) (All ROHs)')
      
      long.rohs <- rohs[rohs$length >= 100e3,]
      frohs.long <- NULL
      for(i in unique(long.rohs$samp)){
        frohs.long <- c(frohs.long, sum(long.rohs[long.rohs$samp == i, 'length'])/30e6)
      }
      ## plot 2
      if(length(frohs.long) > 0){
        hist(frohs.long, xlab='f(ROH) (ROHs >= 100kb)', main=paste0('m ',m,' - r',r,' - p',p))
        ## plot 3
        hist(long.rohs$length, breaks=30, main=paste0('m ',m,' - r',r,' - p',p), xlab='Lenghs ROHs >= 100kb')
      }
        
      ## write output to summarize some stats across parameter combos
      # if(length(frohs.long > 0)){
        save <- c(m, r, p,
                  min(30e6 - (frohs.all*30e6)), mean(30e6 - (frohs.all*30e6)), max(30e6 - (frohs.all*30e6)),
                  mean(rohs$length), sd(rohs$length), mean(frohs.all), mean(long.rohs$length), sd(long.rohs$length), mean(frohs.long))
      # }
      # if(length(frohs.long == 0)){
      #   save <- c(m, r, p,
      #             min(30e6 - (frohs.all*30e6)), mean(30e6 - (frohs.all*30e6)), max(30e6 - (frohs.all*30e6)),
      #             mean(rohs$length), sd(rohs$length), mean(frohs.all), NA, NA, NA)
      # }
      print('Writing summary stats...')
      write.table(save, '/scratch/aubkbk001_01_slim/param_combo_sum_stats.txt',
                  sep='\t', row.names=FALSE, quote=FALSE, append=TRUE, col.names=FALSE)
      
      ## plot cumulative ROH lengths for all samples
      print('Writing final plot...')
      par(mfrow=c(1,1))
      ## plot 1
      if(length(frohs.long) > 0){
        plot(x=c(min(long.rohs$length), max(long.rohs$length)), y=c(0, max(frohs.long)), col='transparent', log='x',
             xlab='ROH length', ylab='Cumulative f(ROH)', main=paste0('m ',m,' - r',r,' - p',p))
          for(i in unique(long.rohs$samp)){
            temp <- long.rohs[long.rohs$samp == i,]
            temp <- temp[order(temp$length),]
            temp$csum <- cumsum(temp$length)
            temp$cfroh <- temp$csum/30e6
            lines(x=temp$length, y=temp$cfroh)
            points(max(long.rohs$length), max(temp$cfroh), pch=19, cex=0.5, col=alpha('deepskyblue3', 0.4))
          }
    	  legend('topleft', legend=c(paste0('Mean # het sites = ',mean(30e6 - (frohs.all*30e6)))),
          	bty = 'n', fill = 'transparent')
        }
    }
  }
}
dev.off()

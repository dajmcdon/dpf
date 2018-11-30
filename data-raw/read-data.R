file.names = list.files(tmp, pattern='svaed')
n.recordings = length(file.names)
splits = strsplit(file.names,'-')
front = sapply(splits, function(x) paste(x[1],x[2],sep='-'))
ids = substr(front,4,nchar(front))

meta = read.csv(paste0(tmp,'/','mazurka_ids.csv')) # should work
meta$performer = as.character(meta$performer)
meta$pid = as.character(meta$pid)
loc = match(ids, meta$pid)
performers = as.character(meta$perf[loc])
years = meta$year[loc]
performance_title = paste(performers, years, sep = '_')
performance_title[ids == '9082-11'] = paste(performance_title[ids == '9082-11'], 'a', sep = '')
performance_title[ids == '9091-17'] = paste(performance_title[ids == '9091-17'], 'b', sep = '')
#performers[is.na(performers)] = ids[is.na(performers)]

recordings = list()
for(i in 1:n.recordings){
  d1 = read.delim(paste0(tmp,'/',file.names[i]), comment='#', head=FALSE)
  n = nrow(d1)
  d1[[3]] = as.character(d1[[3]])
  meas.beat = strsplit(d1[[3]], ':')
  d1$measure = as.numeric(sapply(meas.beat, function(x) x[1]))
  d1$beat = as.numeric(sapply(meas.beat, function(x) x[2]))
  d1[[1]] = d1[[1]] - d1[[1]][1]
  d1$inter = c(0, d1[[1]][-1] - d1[[1]][-n])
  names(d1)[1:3] = c('time','dynamic','tick')
  recordings[[i]] = d1
}
rm(i, meas.beat)
names(recordings) = ids


inter.mat = matrix(NA, n, n.recordings)
for(i in 1:n.recordings) inter.mat[,i] = recordings[[i]]$inter
dynamic.mat = matrix(NA, n, n.recordings)
for(i in 1:n.recordings) dynamic.mat[,i] = recordings[[i]]$dynamic
time.mat = matrix(NA, n, n.recordings)
for(i in 1:n.recordings) time.mat[,i] = recordings[[i]]$time


#Save data
dynamics = data.frame(cbind(d1$measure, d1$beat, d1$measure + (d1$beat - 1)/3, dynamic.mat))
colnames(dynamics) = c('meas_num', 'beat', 'note_onset', performance_title)
devtools::use_data(dynamics) # puts it in the package in the right format

tempos = array(dim = c(nrow(dynamics), ncol(dynamics)))
tempos[,1:3] = as.matrix(dynamics[,1:3])
for(i in 4:ncol(tempos)){
    tempos[,i] = diff(c(dynamics$note_onset,61))*3*60/diff(c(recordings[[i-3]]$time,113))
    tempos[nrow(tempos),i] = mean(c(tempos[nrow(tempos) - 1,i], tempos[nrow(tempos) - 2,i]))
}
colnames(tempos) = colnames(dynamics)
tempos = data.frame(tempos)
devtools::use_data(tempos)

remove(recordings)
recordings = data.frame(performers, years)
colnames(recordings) = c('performer', 'year')
devtools::use_data(recordings)
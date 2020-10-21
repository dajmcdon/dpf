pars = unlist(pvec_ml["Richter 1976",])
yt = matrix(tempos$Richter_1976,nrow=1)
lt = diff(c(tempos$note_onset, 61))
pmats = musicModel(
  lt, pars[1], pars[2:4], c(pars[5],1,1),
  pars[6:12], c(132,0), c(400,10))
beam = with(
  pmats, beamSearch(
    a0, P0, c(1,0,0,0,0,0,0,0,0,0), 
    dt, ct, Tt, Zt, HHt, GGt, yt, transMat, 400))
bestpath = beam$paths[which.max(beam$weights),]
kal = kalman(pmats, bestpath, yt)

richt = read_csv("aoas/richter.csv", col_names = FALSE)
blah = as.numeric(richt$X4[richt$X3=="Tempo"])
csv_onset = as.numeric(richt$X2[richt$X3=="Tempo"])
m_onset = (csv_onset %/% 120) /3 + 1
m_onset = as.integer(m_onset * 12)
ref_temps = as.integer(tempos$note_onset * 12)
m_temps = 60000000/blah 
replace_temps = kal$ests[match(m_onset, ref_temps)]
replace_temps[is.na(replace_temps)] = m_temps[is.na(replace_temps)]

richt$X4[richt$X3=="Tempo"] = round(60000000/replace_temps)
write_csv(richt, "aoas/richter_smooth.csv", col_names = FALSE)

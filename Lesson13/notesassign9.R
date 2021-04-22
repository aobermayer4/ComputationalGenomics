#----Question 1 code----#
#I wrote this code to list the columns I was ranking to find the top observed

data.frame(count=sort(table(demo_asv_taxa[,2]), decreasing=TRUE))
data.frame(count=sort(table(demo_asv_taxa[,6]), decreasing=TRUE))


#----Question 2 code----#

percloss <- expression(abs(((n-o)/o)*100))
n <- demo_track_steps[1:4,1]
o <- demo_track_steps[1:4,2]
loss2filterd <- eval(percloss)
n <- demo_track_steps[1:4,2]
o <- demo_track_steps[1:4,3]
loss2denoisedF <- eval(percloss)
n <- demo_track_steps[1:4,3]
o <- demo_track_steps[1:4,4]
loss2denoisedR <- eval(percloss)
n <- demo_track_steps[1:4,4]
o <- demo_track_steps[1:4,5]
loss2merged <- eval(percloss)
n <- demo_track_steps[1:4,5]
o <- demo_track_steps[1:4,6]
loss2nochim <- eval(percloss)
prlossall <- cbind(loss2filterd,loss2denoisedF,loss2denoisedR,loss2merged,loss2nochim)
coln <- c('%loss2filterd', '%loss2denoisedF','%loss2denoisedR','%loss2merged','%loss2nochim')
colnames(prlossall) <- coln
prlossall




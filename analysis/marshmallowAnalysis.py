#Basic analyses of Marshmallow data
import numpy as np 
import matplotlib.pyplot as plt 
import scipy.stats as stats
import scipy.io as sio
import pandas as pd

##### column name meanings #####
#testData
# % run
# % trial
# % cb: 1 = beaches nback, 2 = gazebos nback
# % category_cue: category of the cued scene and competitor (1 = beach, 2 = gazebo)
# % category_noncomp: category of the non-competitor item (1 = beach, 2 = gazebo)
# % pairnum_cue: which pair within a category the cue and competitor comes from
# % pairnum_noncomp: which pair within a category the non-competitor comes from
# % pairitem_cue: which item is being shown on this trial (pairmate 1 or pairmate 2) for the cue
# % pairitem_noncomp: which item is being shown on this trial (pairmate 1 or pairmate 2) for the non-competitor
# % pairid_cue: ID # of this pair (based on all possible combinations) for the cue and competitor
# % pairid_noncomp: ID # of this pair (based on all possible combinations) for the non-competitor
# % scenenum: which scene is being shown for the cue
# % objnum_obj1: which object is being shown for the left item
# % objnum_obj2 which object is being shown for the middle item
# % objnum_obj3: which object is being shown for the right item
# % condition_cue: similarity condition (0-4) of the cue and competitor
# % condition_noncomp: similarity condition (0-4) of the non-competitor
# % objpos1: which item appeared on the left (1 = target, 2 = competitor, 3 = non-competitor)
# % objpos2: which item appeared on the middle (1 = target, 2 = competitor, 3 = non-competitor)
# % objpos3:which item appeared on the right (1 = target, 2 = competitor, 3 = non-competitor)
# % resp: actual response
# % score: 1 = target, 2 = competitor, 3 = non-competitor
# % rt: response time (s)

#nback
# % run
# % trial
# % category: 1 = beach, 2 = gazebo
# % stim: stimulus number shown
# % cresp: correct response
# % resp: actual response
# % acc: accuracy
# % rt: response time (s)

##### Helpers #####
#helper to load data
def loadData(sn,task="both",robertData=0):
	datadir = "../data/"
	nbackData = np.empty([1000,8])
	testData = np.empty([240,23])
	subj_str = "marshmallow_"+str(sn)
	if robertData==1:
		if sn<10:
			subj_str = "marshmallow_00"+str(sn)
		else:
			subj_str = "marshmallow_0"+str(sn)
	fullpath_str = datadir+subj_str+"/"+subj_str
	if np.logical_and(robertData==1,sn==4):
		fullpath_str = datadir+subj_str+"/marshmallow_003" #dunno how this happened but I suspect it was just an error while running the code
	header = sio.loadmat(fullpath_str+"_header.mat")
	for irun in range(20):
		data = sio.loadmat(fullpath_str+"_nback_1_"+str(irun+1)+".mat")
		nbackData[(irun*50):(irun+1)*50] = data["data"]["rundata"][0][0]
	nbackData = pd.DataFrame(nbackData,columns=["run","trial","category","stim","cresp","resp","acc","rt"])
	for irun in range(3):
		data = sio.loadmat(fullpath_str+"_studytest_"+str(irun+1)+".mat")
		testData[(irun*80):(irun+1)*80] = data["data"]["testdata"][0][0]
	testData = pd.DataFrame(testData,columns=["run","trial","cat_nb","cat_cue","cat_noncomp",
												"pairnum_cue","pairnum_noncomp","pairitem_cue","pairitem_noncomp",
												"pairid_cue","pairid_noncomp",
												"scenenum","objnum_obj1","objnum_obj2","objnum_obj3","cond_cue",
												"cond_noncomp","objpos1","objpos2","objpos3","resp","score","rt"])
	return header, nbackData, testData

#do some descriptive stats on the nback data
def nbackSummary(nbackdata,sn,printSubj=0):
	mean_nback_acc = np.nanmean(nbackData.acc)
	mean_nback_rt = np.nanmean(nbackData.rt)
	if printSubj:
		print("subj"+str(sn)+" mean nBack acc = "+str(mean_nback_acc))
		print("subj"+str(sn)+" mean nBack rt = "+str(mean_nback_rt))
	meanRT = np.zeros(20)
	meanAcc = np.zeros(20)
	for i in range(20):
		meanRT[i] = np.nanmean(nbackData.rt[nbackData.run==(i+1)])
		meanAcc[i] = np.nanmean(nbackData.acc[nbackData.run==(i+1)])
	return meanRT,meanAcc

def studytestSummary(testData,sn,printSubj=0):
	#first, we fix this bullshit
	testData.cat_nb = testData.cat_nb - 1 #who tf uses 1 and 2 as their binary codes????
	testData.cat_cue = testData.cat_cue - 1
	testData.cat_noncomp = testData.cat_noncomp - 1
	#assign catagory labels
	learned_cat = testData.cat_nb[0]
	novel_cat = (learned_cat*-1)+1

	meanAcc = np.nansum(testData.score==1)/240
	#break it up by stimulus condition
	learned_meanAcc = np.nansum(testData.score[testData.cat_cue==learned_cat]==1)/np.nansum(testData.cat_cue==learned_cat)
	novel_meanAcc = np.nansum(testData.score[testData.cat_cue==novel_cat]==1)/np.nansum(testData.cat_cue==novel_cat)
	if printSubj:
		print("subj"+str(sn)+" mean ST acc = "+str(meanAcc))
		print("subj"+str(sn)+" learned stimulus mean acc = "+str(learned_meanAcc))
		print("subj"+str(sn)+" novel stimulus mean acc = "+str(novel_meanAcc))
		#document nan responses to make sure people are actually doing the task
		print("subj"+str(sn)+" nan responses = "+str(np.sum(np.isnan(testData.rt))))
		print("\n")
	#break it up by block
	#basic stats
	learned_acc_block = np.zeros(3)
	novel_acc_block = np.zeros(3)
	#histogram stats
	learned_num_resps = np.zeros((3,3))
	novel_num_resps = np.zeros((3,3))
	learned_rt_resps = np.zeros((3,3))
	novel_rt_resps = np.zeros((3,3))
	#rt stats
	learned_rt_block = np.zeros(3)
	novel_rt_block = np.zeros(3)
	for i in range(3):
		learned_rt_block[i] = np.nanmean(testData.rt[np.logical_and(testData.cat_cue==learned_cat,testData.run==(i+1))])
		novel_rt_block[i] = np.nanmean(testData.rt[np.logical_and(testData.cat_cue==novel_cat,testData.run==(i+1))])
		learned_acc_block[i] = np.nansum(testData.score[np.logical_and(testData.cat_cue==learned_cat,testData.run==(i+1))]==1)/np.nansum(np.logical_and(testData.cat_cue==learned_cat,testData.run==(i+1)))
		novel_acc_block[i] = np.nansum(testData.score[np.logical_and(testData.cat_cue==novel_cat,testData.run==(i+1))]==1)/np.nansum(np.logical_and(testData.cat_cue==novel_cat,testData.run==(i+1)))
		learned_num_resps[0,i] = np.nansum(testData.score[np.logical_and(testData.cat_cue==learned_cat,testData.run==(i+1))]==1)
		learned_num_resps[1,i] = np.nansum(testData.score[np.logical_and(testData.cat_cue==learned_cat,testData.run==(i+1))]==2)
		learned_num_resps[2,i] = np.nansum(testData.score[np.logical_and(testData.cat_cue==learned_cat,testData.run==(i+1))]==3)
		novel_num_resps[0,i] = np.nansum(testData.score[np.logical_and(testData.cat_cue==novel_cat,testData.run==(i+1))]==1)
		novel_num_resps[1,i] = np.nansum(testData.score[np.logical_and(testData.cat_cue==novel_cat,testData.run==(i+1))]==2)
		novel_num_resps[2,i] = np.nansum(testData.score[np.logical_and(testData.cat_cue==novel_cat,testData.run==(i+1))]==3)
		#these sometimes throw a warning because some subjects did not chose the competitor/noncompetior for whole blocks, so theres nothing to average
		learned_rt_resps[0,i] = np.nanmean(testData.rt[np.logical_and(np.logical_and(testData.cat_cue==learned_cat,testData.run==(i+1)),testData.score==1)])
		learned_rt_resps[1,i] = np.nanmean(testData.rt[np.logical_and(np.logical_and(testData.cat_cue==learned_cat,testData.run==(i+1)),testData.score==2)])
		learned_rt_resps[2,i] = np.nanmean(testData.rt[np.logical_and(np.logical_and(testData.cat_cue==learned_cat,testData.run==(i+1)),testData.score==3)])
		novel_rt_resps[0,i] = np.nanmean(testData.rt[np.logical_and(np.logical_and(testData.cat_cue==novel_cat,testData.run==(i+1)),testData.score==1)])
		novel_rt_resps[1,i] = np.nanmean(testData.rt[np.logical_and(np.logical_and(testData.cat_cue==novel_cat,testData.run==(i+1)),testData.score==2)])
		novel_rt_resps[2,i] = np.nanmean(testData.rt[np.logical_and(np.logical_and(testData.cat_cue==novel_cat,testData.run==(i+1)),testData.score==3)])
	return learned_acc_block,novel_acc_block,learned_num_resps,novel_num_resps,learned_rt_block,novel_rt_block,learned_rt_resps,novel_rt_resps,learned_cat

#assumes levels of responses are balanced (wont work with nback)
#Adapted from this source: https://psychology.stackexchange.com/questions/13386/in-a-forced-choice-task-what-proportion-of-responses-is-above-chance-level
def binomial_cutoff(ntrials,npermutations,p_cutoff,num_levels):
	#randomly generate a 'correct' set of responses
	expected = np.random.choice(num_levels,ntrials)
	#randomly generate a responses to many runs of the experiment
	random_resps = np.random.choice(num_levels,[npermutations,ntrials])
	#get mean accuracy for each random run
	percent_correct = np.nansum(random_resps==expected,1)/ntrials
	#find the point on the distribution of accuracies where anything less would be statistically significant
	cutoff = np.quantile(percent_correct,1-p_cutoff)
	return cutoff

############# DRIVER ##############
#set constants
exclusion = 1 #1 = exclude subjects below chance perfomance on the association task 0 = no exclusion
printSubj = 0 #toggle whether to print individual subject's behavioral performance
show_plot = 0 #toggle whether to show plots or just save as pdfs
save_plot = 1 #toggle whether to save plots (useful if you're fucking around)
nsubj = 52
nRuns_nback = 20
nRuns_ST = 3
figure_dir = "../figures/"
#preallocate data frames
meanRT_nback = np.zeros((nsubj,nRuns_nback))
meanAcc_nback = np.zeros((nsubj,nRuns_nback))
meanAcc_learned = np.zeros((nsubj,nRuns_ST))
meanAcc_novel = np.zeros((nsubj,nRuns_ST))
meanRT_learned = np.zeros((nsubj,nRuns_ST))
meanRT_novel = np.zeros((nsubj,nRuns_ST))
histData_learned = np.zeros((nsubj,3,nRuns_ST))
histData_novel = np.zeros((nsubj,3,nRuns_ST))
rtData_learned = np.zeros((nsubj,3,nRuns_ST))
rtData_novel = np.zeros((nsubj,3,nRuns_ST))
category_learned = np.zeros(nsubj)
#determine what "chance" performance cutoff is
cutoff = binomial_cutoff(240,10000,0.05,3)
print("acc cutoff = ",cutoff)
#load data by subject
for sn in range(nsubj):
	subjID = sn+1
	if subjID>40:
		subjID = subjID-40
		[header,nbackData,testData] = loadData(subjID,robertData=1)
	else:
		[header,nbackData,testData] = loadData(subjID)

	#nback analysis
	[meanRT_nback[sn], meanAcc_nback[sn]] = nbackSummary(nbackData,subjID,printSubj)

	#testData analysis
	[meanAcc_learned[sn], meanAcc_novel[sn],histData_learned[sn],histData_novel[sn],
	meanRT_learned[sn],meanRT_novel[sn],rtData_learned[sn],rtData_novel[sn],category_learned[sn]] = studytestSummary(testData,subjID,printSubj)

	#apriori below chance exclusion criteria
	if exclusion:
		if np.nansum(testData.score==1)/240 < cutoff: 
			 meanAcc_learned[sn,:] = np.nan
			 meanAcc_novel[sn,:] = np.nan
			 histData_learned[sn,:,:] = np.nan
			 histData_novel[sn,:,:] = np.nan
			 meanAcc_nback[sn,:] = np.nan
			 meanRT_nback[sn,:] = np.nan
			 meanRT_learned[sn,:] = np.nan
			 meanRT_novel[sn,:] = np.nan
			 rtData_learned[sn,:] = np.nan
			 rtData_novel[sn,:] = np.nan

#recalculate sample size after exclusion
nsubj = nsubj-np.sum(np.isnan(np.mean(meanAcc_learned,1)))
print('nsubj after exclusion = ',nsubj)
#because of the error when originally running this the stimuli arent balanced so I'm just excluding all of roberts data for this single analysis since thats the most principled way I can balance this
category_learned[40:52] = np.nan 
print(np.nansum(category_learned==1))

#print the basics
print("Familiar accuracy = "+str(np.nanmean(meanAcc_learned))," Unfamiliar accuracy = "+str(np.nanmean(meanAcc_novel)))

#ttest of mean accuracy distributions between familiar and unfamiliar stimuli
print("Familiar vs unfamiliar acc ttest: ",np.round(stats.ttest_rel(np.nanmean(meanAcc_learned,1),np.nanmean(meanAcc_novel,1),nan_policy='omit'),6))
#plot difference between familiar and unfamiliar accuracy means
fig = plt.figure()
plt.bar([0,1],[np.nanmean(meanAcc_learned),np.nanmean(meanAcc_novel)],color=['blue','m'])
#error term is SEM
plt.errorbar([0,1],[np.nanmean(meanAcc_learned),np.nanmean(meanAcc_novel)],yerr=[np.nanstd(meanAcc_learned)/np.sqrt(nsubj),np.nanstd(meanAcc_novel)/np.sqrt(nsubj)],color='k')
plt.xlabel("stimulus condition")
plt.xticks([0,1],["familiar","unfamiliar"])
plt.ylabel("mean accuracy (%)")
plt.yticks([0,.2,.4,.6,.8],[0,20,40,60,80])
plt.title("Mean accuracy by condition")
#plt.scatter([0.5],[.8],marker='*',color='k')
if save_plot:
	plt.savefig(figure_dir+'acc_diff_familiar_vs_nonfamiliar.pdf',format='pdf')
if show_plot:
	plt.show()

####### Plot distribution of responses by condition #######
#plot histograms of responses for each condition
#ttest between familiar and unfamiliar target #dunno why I did this but it looks pretty I guess
print("Familiar vs unfamiliar target resps ttest: ",np.round(stats.ttest_rel(np.sum(histData_learned,2)[:,0],np.sum(histData_novel,2)[:,0],nan_policy='omit'),6))
#ttest between familiar and unfamiliar competitor
print("Familiar vs unfamiliar competitor resps ttest: ",np.round(stats.ttest_rel(np.sum(histData_learned,2)[:,1],np.sum(histData_novel,2)[:,1],nan_policy='omit'),6))
#ttest between familiar and unfamiliar noncompetitor
print("Familiar vs unfamiliar noncompetitor resps ttest: ",np.round(stats.ttest_rel(np.sum(histData_learned,2)[:,2],np.sum(histData_novel,2)[:,2],nan_policy = 'omit'),6))
#anova between familiar responses 
print("Familiar anova: ",stats.f_oneway(np.nansum(histData_learned,2)[:,0],np.nansum(histData_learned,2)[:,1],np.nansum(histData_learned,2)[:,2],axis=0))
#anova between unfamiliar responses 
print("Unfamiliar anova: ",stats.f_oneway(np.nansum(histData_novel,2)[:,0],np.nansum(histData_novel,2)[:,1],np.nansum(histData_novel,2)[:,2],axis=0))
#plot
fig = plt.figure()
#collapse across subjects
histAllsubj_learned = np.nansum(histData_learned,0)
histAllsubj_novel = np.nansum(histData_novel,0)
c1 = plt.bar([0,1.5,3],[np.nansum(histAllsubj_learned[0]),np.nansum(histAllsubj_learned[1]),np.nansum(histAllsubj_learned[2])],width=.5,color='b')
c2 = plt.bar([.5,2,3.5],[np.nansum(histAllsubj_novel[0]),np.nansum(histAllsubj_novel[1]),np.nansum(histAllsubj_novel[2])],width=.5,color='m')
#compute error term (will be len(3) array)
learned_hist_err = np.std(np.nansum(histData_learned,2),0)
novel_hist_err = np.std(np.nansum(histData_novel,2),0)
#plot
plt.errorbar([0,.5],[np.nansum(histAllsubj_learned[0]),np.nansum(histAllsubj_novel[0])],yerr=[learned_hist_err[0],novel_hist_err[0]],color='k')
plt.errorbar([1.5,2],[np.nansum(histAllsubj_learned[1]),np.nansum(histAllsubj_novel[1])],yerr=[learned_hist_err[1],novel_hist_err[1]],color='k')
plt.errorbar([3,3.5],[np.nansum(histAllsubj_learned[2]),np.nansum(histAllsubj_novel[2])],yerr=[learned_hist_err[2],novel_hist_err[2]],color='k')
plt.legend([c1,c2],['familiar','unfamiliar'])
plt.xticks([.25,1.75,3.25],['target','competitor','noncompetitor'])
plt.ylabel('# of responses')
plt.title("Histogram of responses by condition")
if save_plot:
	plt.savefig(figure_dir+'response_histograms.pdf',format='pdf')
if show_plot:
	plt.show()

#plot histogrm but as a mean
totTrials = (240*nsubj)/2
fig = plt.figure()
#collapse across subjects
histAllsubj_learned = np.nansum(histData_learned,0)
histAllsubj_novel = np.nansum(histData_novel,0)
c1 = plt.bar([0,1.5,3],[np.nansum(histAllsubj_learned[0])/totTrials,np.nansum(histAllsubj_learned[1])/totTrials,np.nansum(histAllsubj_learned[2])/totTrials],width=.5,color='b')
c2 = plt.bar([.5,2,3.5],[np.nansum(histAllsubj_novel[0])/totTrials,np.nansum(histAllsubj_novel[1])/totTrials,np.nansum(histAllsubj_novel[2])/totTrials],width=.5,color='m')
#compute error term (will be len(3) array)
learned_hist_err = np.std(np.nansum(histData_learned,2)/totTrials,0)
novel_hist_err = np.std(np.nansum(histData_novel,2)/totTrials,0)
#plot
plt.errorbar([0,.5],[np.nansum(histAllsubj_learned[0])/totTrials,np.nansum(histAllsubj_novel[0])/totTrials],yerr=[learned_hist_err[0],novel_hist_err[0]],color='k')
plt.errorbar([1.5,2],[np.nansum(histAllsubj_learned[1])/totTrials,np.nansum(histAllsubj_novel[1])/totTrials],yerr=[learned_hist_err[1],novel_hist_err[1]],color='k')
plt.errorbar([3,3.5],[np.nansum(histAllsubj_learned[2])/totTrials,np.nansum(histAllsubj_novel[2])/totTrials],yerr=[learned_hist_err[2],novel_hist_err[2]],color='k')
plt.legend([c1,c2],['familiar','unfamiliar'])
plt.xticks([.25,1.75,3.25],['target','competitor','noncompetitor'])
plt.ylabel('# of responses')
plt.title("Mean responses by condition")
if save_plot:
	plt.savefig(figure_dir+'response_means.pdf',format='pdf')
if show_plot:
	plt.show()

###### nback performance over blocks #######
#plot nback acc as a function of blocks
#ttest between blocks 1 and 20
print("nback block acc ttest: ",stats.ttest_rel(meanAcc_nback[:,0],meanAcc_nback[:,19],nan_policy='omit'))
#plot
fig = plt.figure()
plt.plot(np.arange(0,20),np.nanmean(meanAcc_nback,0),color='m',linewidth=3)
#error is sem
plt.errorbar(np.arange(0,20),np.nanmean(meanAcc_nback,0),yerr=np.nanstd(meanAcc_nback,0)/np.sqrt(nsubj),color='k',capsize=3)
plt.yticks([0.5,0.6,0.7,0.8,0.9,1],[50,60,70,80,90,100])
#plt.scatter([0.5,2.5,9.5],[1,1,1],marker='*',color='k')
plt.ylabel('accuracy (%)')
plt.xlabel('block number')
plt.xticks([0,5,10,15,20])
plt.title('nback mean accuracy by block')
if save_plot:
	plt.savefig(figure_dir+'nback_acc_blockwise.pdf',format='pdf')
if show_plot:
	plt.show()

#RTs, cause why not
#ttest between blocks 1 and 20
print("nback block rt ttest: ",stats.ttest_rel(meanRT_nback[:,0],meanRT_nback[:,19],nan_policy='omit'))
fig = plt.figure()
plt.plot(np.arange(0,20),np.nanmean(meanRT_nback,0),color='m',linewidth=3)
#error is std
plt.errorbar(np.arange(0,20),np.nanmean(meanRT_nback,0),yerr=np.nanstd(meanRT_nback,0)/np.sqrt(nsubj),color='k',capsize=3)
plt.yticks([0.5,0.6,0.7,0.8],[500,600,700,800])
#plt.scatter([0.5,2.5,9.5],[1,1,1],marker='*',color='k')
plt.ylabel('RT (ms)')
plt.xlabel('block number')
plt.xticks([0,5,10,15,20])
plt.title('nback mean response time by block')
if save_plot:
	plt.savefig(figure_dir+'nback_rt_blockwise.pdf',format='pdf')
if show_plot:
	plt.show()

###### correlations between nback performance and association acc ########
#across subject correlation of nback acc and learned acc
#remove nans from data
clean_meanAcc_nback = meanAcc_nback[~np.isnan(meanAcc_nback[:,1])]
clean_meanAcc_novel = meanAcc_novel[~np.isnan(meanAcc_novel[:,1])]
clean_meanAcc_learned = meanAcc_learned[~np.isnan(meanAcc_learned[:,1])]
[r,p] = stats.pearsonr(np.nanmean(clean_meanAcc_nback,1),np.nanmean(clean_meanAcc_learned,1))
print("Pearson r of across subject correlation of nback and familiar acc: ",np.round(r,5),np.round(p,5))
[m,b] = np.polyfit(np.nanmean(clean_meanAcc_nback,1),np.nanmean(clean_meanAcc_learned,1),1)
fig = plt.figure()
plt.scatter(np.nanmean(clean_meanAcc_nback,1),np.nanmean(clean_meanAcc_learned,1),color='gray')
plt.plot(np.nanmean(clean_meanAcc_nback,1),m*np.nanmean(clean_meanAcc_nback,1)+b,color='m')
plt.xlabel('mean nback accuracy (%)')
plt.ylabel('mean learned stimulus accuracy (%)')
plt.yticks([.2,.3,.4,.5,.6,.7,.8,.9],[20,30,40,50,60,70,80,90])
plt.xticks([.65,.7,.75,.8,.85,.9,.95],[65,70,75,80,85,90,95])
plt.title('nback and learned stimulus accuracy correlation\nr = '+str(np.round(r,2)))
if save_plot:
	plt.savefig(figure_dir+'famAcc_correlation.pdf',format='pdf')
if show_plot:
	plt.show()

#across subject correlation of nback and novel acc
[r,p] = stats.pearsonr(np.nanmean(clean_meanAcc_nback,1),np.nanmean(clean_meanAcc_novel,1))
print("Pearson r of across subject correlation of nback and unfamiliar acc: ",np.round(r,5),np.round(p,5))
[m,b] = np.polyfit(np.nanmean(clean_meanAcc_nback,1),np.nanmean(clean_meanAcc_novel,1),1)
fig = plt.figure()
plt.scatter(np.nanmean(clean_meanAcc_nback,1),np.nanmean(clean_meanAcc_novel,1),color='gray')
plt.plot(np.nanmean(clean_meanAcc_nback,1),m*np.nanmean(clean_meanAcc_nback,1)+b,color='m')
plt.xlabel('mean nback accuracy (%)')
plt.ylabel('mean novel stimulus accuracy (%)')
plt.yticks([.2,.3,.4,.5,.6,.7,.8,.9],[20,30,40,50,60,70,80,90])
plt.xticks([.65,.7,.75,.8,.85,.9,.95],[65,70,75,80,85,90,95])
plt.title('nback and novel stimulus accuracy correlation\nr = '+str(np.round(r,2)))
if save_plot:
	plt.savefig(figure_dir+'unfamAcc_correlation.pdf',format='pdf')
if show_plot:
	plt.show()

###### Look at response times for the association task #######
#plot rts for study test
print("Familiar rts = "+str(np.nanmean(meanRT_learned))," Unfamiliar rts = "+str(np.nanmean(meanRT_novel)))
#ttest of mean accuracy distributions between familiar and unfamiliar stimuli
print("Familiar vs unfamiliar rt ttest: ",np.round(stats.ttest_rel(np.nanmean(meanRT_learned,1),np.nanmean(meanRT_novel,1),nan_policy='omit'),6))
#plot difference between familiar and unfamiliar accuracy means
fig = plt.figure()
plt.bar([0,1],[np.nanmean(meanRT_learned),np.nanmean(meanRT_novel)],color=['blue','m'])
#error term is SEM
plt.errorbar([0,1],[np.nanmean(meanRT_learned),np.nanmean(meanRT_novel)],yerr=[np.nanstd(meanRT_learned)/np.sqrt(nsubj),np.nanstd(meanRT_novel)/np.sqrt(nsubj)],color='k')
plt.xlabel("stimulus condition")
plt.xticks([0,1],["familiar","unfamiliar"])
plt.ylabel("mean RTs (ms)")
plt.yticks([0,.3,.6,.9,1.2,1.5],[0,300,600,900,1200,1500])
plt.title("Mean RT by condition")
plt.scatter([0.5],[1.5],marker='*',color='k')
if save_plot:
	plt.savefig(figure_dir+'rt_diff_familiar_vs_nonfamiliar.pdf',format='pdf')
if show_plot:
	plt.show()

#plot mean respone times split by responses
print("Familiar vs unfamiliar target rt ttest: ",np.round(stats.ttest_rel(np.mean(rtData_learned,2)[:,0],np.mean(rtData_novel,2)[:,0],nan_policy='omit'),6))
#ttest between familiar and unfamiliar competitor
print("Familiar vs unfamiliar competitor rt ttest: ",np.round(stats.ttest_rel(np.mean(rtData_learned,2)[:,1],np.mean(rtData_novel,2)[:,1],nan_policy='omit'),6))
#ttest between familiar and unfamiliar noncompetitor
print("Familiar vs unfamiliar noncompetitor rt ttest: ",np.round(stats.ttest_rel(np.mean(rtData_learned,2)[:,2],np.mean(rtData_novel,2)[:,2],nan_policy = 'omit'),6))
#anova between familiar responses 
print("Familiar anova: ",stats.f_oneway(np.nansum(histData_learned,2)[:,0],np.nanmean(rtData_learned,2)[:,1],np.nanmean(rtData_learned,2)[:,2],axis=0))
#anova between unfamiliar responses 
print("Unfamiliar anova: ",stats.f_oneway(np.nansum(histData_novel,2)[:,0],np.nanmean(rtData_novel,2)[:,1],np.nanmean(rtData_novel,2)[:,2],axis=0))
fig = plt.figure()
#collapse across subjects
rtData_learned_blockmean = np.nanmean(rtData_learned,2)
rtData_novel_blockmean = np.nanmean(rtData_novel,2)
c1 = plt.bar([0,1.5,3],[np.nanmean(rtData_learned_blockmean[:,0]),np.nanmean(rtData_learned_blockmean[:,1]),np.nanmean(rtData_learned_blockmean[:,2])],width=.5,color='b')
c2 = plt.bar([.5,2,3.5],[np.nanmean(rtData_novel_blockmean[:,0]),np.nanmean(rtData_novel_blockmean[:,1]),np.nanmean(rtData_novel_blockmean[:,2])],width=.5,color='m')
#compute error term (will be len(3) array)
learned_rt_err = np.std(rtData_learned_blockmean,0)/np.sqrt(nsubj)
novel_rt_err = np.std(rtData_novel_blockmean,0)/np.sqrt(nsubj)
#plot
plt.errorbar([0,.5],[np.nanmean(rtData_learned_blockmean[:,0]),np.nanmean(rtData_novel_blockmean[:,0])],yerr=[learned_rt_err[0],novel_rt_err[0]],color='k')
plt.errorbar([1.5,2],[np.nanmean(rtData_learned_blockmean[:,1]),np.nanmean(rtData_novel_blockmean[:,1])],yerr=[learned_rt_err[1],novel_rt_err[1]],color='k')
plt.errorbar([3,3.5],[np.nanmean(rtData_learned_blockmean[:,2]),np.nanmean(rtData_novel_blockmean[:,2])],yerr=[learned_rt_err[2],novel_rt_err[2]],color='k')
plt.legend([c1,c2],['familiar','unfamiliar'])
plt.xticks([.25,1.75,3.25],['target','competitor','noncompetitor'])
plt.ylabel('mean rt (s)')
plt.title("Mean response times by condition")
if save_plot:
	plt.savefig(figure_dir+'responsetime_means.pdf',format='pdf')
if show_plot:
	plt.show()

######## Look at blockwise performance on the association task ########
#ttests between condition means per block
print("Familiar vs unfamiliar block 1 acc ttest: ",np.round(stats.ttest_rel(meanAcc_learned[:,0],meanAcc_novel[:,0],nan_policy='omit'),6))
print("Familiar vs unfamiliar block 2 acc ttest: ",np.round(stats.ttest_rel(meanAcc_learned[:,1],meanAcc_novel[:,1],nan_policy='omit'),6))
print("Familiar vs unfamiliar block 3 acc ttest: ",np.round(stats.ttest_rel(meanAcc_learned[:,2],meanAcc_novel[:,2],nan_policy='omit'),6))
#plot the mean accuracy by block
fig = plt.figure()
c1 = plt.bar([0,1.5,3],np.nanmean(meanAcc_learned,0),width=.5,color='b')
c2 = plt.bar([.5,2,3.5],np.nanmean(meanAcc_novel,0),width=.5,color='m')
plt.errorbar([0,.5],[np.nanmean(meanAcc_learned,0)[0],np.nanmean(meanAcc_novel,0)[0]],yerr=[np.nanstd(meanAcc_learned,0)[0]/np.sqrt(nsubj),np.nanstd(meanAcc_novel,0)[0]/np.sqrt(nsubj)],color='k')
plt.errorbar([1.5,2],[np.nanmean(meanAcc_learned,0)[1],np.nanmean(meanAcc_novel,0)[1]],yerr=[np.nanstd(meanAcc_learned,0)[1]/np.sqrt(nsubj),np.nanstd(meanAcc_novel,0)[1]/np.sqrt(nsubj)],color='k')
plt.errorbar([3,3.5],[np.nanmean(meanAcc_learned,0)[2],np.nanmean(meanAcc_novel,0)[2]],yerr=[np.nanstd(meanAcc_learned,0)[2]/np.sqrt(nsubj),np.nanstd(meanAcc_novel,0)[2]/np.sqrt(nsubj)],color='k')
plt.xlabel('block #')
plt.xticks([.25,1.75,3.25],[1,2,3])
plt.ylabel('accuracy (%)')
plt.yticks([.1,.2,.3,.4,.5,.6,.7,.8],[10,20,30,40,50,60,70,80])
plt.legend([c1,c2],['familiar','unfamiliar'])
if save_plot:
	plt.savefig(figure_dir+'assoc_acc_by_block.pdf',format='pdf')
if show_plot:
	plt.show()

#look at competitor vs noncompetitor
blockTrials = (80*nsubj)/2
fig = plt.figure()
c1 = plt.bar([0,2.5,5],histAllsubj_learned[1]/blockTrials,color='b',width=.5)
c2 = plt.bar([1,3.5,6],histAllsubj_learned[2]/blockTrials,color='b',width=.5)
c3 = plt.bar([.5,3,5.5],histAllsubj_novel[1]/blockTrials,color='m',width=.5)
c4 = plt.bar([1.5,4,6.5],histAllsubj_novel[2]/blockTrials,color='m',width=.5)
learned_hist_err = np.nanstd(histData_learned/40,0)/np.sqrt(nsubj)
novel_hist_err = np.nanstd(histData_novel/40,0)/np.sqrt(nsubj)
plt.errorbar([0,.5],[histAllsubj_learned[1,0]/blockTrials,histAllsubj_novel[1,0]/blockTrials],yerr=[learned_hist_err[1,0],novel_hist_err[1,0]],color='k')
plt.errorbar([1,1.5],[histAllsubj_learned[2,0]/blockTrials,histAllsubj_novel[2,0]/blockTrials],yerr=[learned_hist_err[2,0],novel_hist_err[2,0]],color='k')
plt.errorbar([2.5,3],[histAllsubj_learned[1,1]/blockTrials,histAllsubj_novel[1,1]/blockTrials],yerr=[learned_hist_err[1,1],novel_hist_err[1,1]],color='k')
plt.errorbar([3.5,4],[histAllsubj_learned[2,1]/blockTrials,histAllsubj_novel[2,1]/blockTrials],yerr=[learned_hist_err[2,1],novel_hist_err[2,1]],color='k')
plt.errorbar([5,5.5],[histAllsubj_learned[1,2]/blockTrials,histAllsubj_novel[1,2]/blockTrials],yerr=[learned_hist_err[1,2],novel_hist_err[1,2]],color='k')
plt.errorbar([6,6.5],[histAllsubj_learned[2,2]/blockTrials,histAllsubj_novel[2,2]/blockTrials],yerr=[learned_hist_err[2,2],novel_hist_err[2,2]],color='k')
plt.xticks([.25,1.25,2.75,3.75,5.25,6.25],['comp','noncomp','comp','noncomp','comp','noncomp'])
plt.yticks([0,.05,.1,.15,.2,.25,.3],[0,5,10,15,20,25,30])
plt.ylabel("accuracy (%)")
plt.legend([c1,c3],['familiar','unfamiliar'])
if save_plot:
	plt.savefig(figure_dir+'comp_vs_nocomp_blockwise.pdf',format='pdf')
if show_plot:
	plt.show()


####### check that gazebos arent just easier to learn than beaches ######
print('beachs vs gaxebos ttest: ',stats.ttest_rel(np.nanmean(meanAcc_learned[category_learned==0],1),np.nanmean(meanAcc_learned[category_learned==1],1),nan_policy='omit'))
fig = plt.figure()
plt.bar([0,1],[np.nanmean(meanAcc_learned[category_learned==0]),np.nanmean(meanAcc_learned[category_learned==1])],color=['b','g'])
plt.errorbar([0,1],[np.nanmean(meanAcc_learned[category_learned==0]),np.nanmean(meanAcc_learned[category_learned==1])],yerr=[np.nanmean(meanAcc_learned[category_learned==0]),np.nanmean(meanAcc_learned[category_learned==1]/20)])
plt.xticks([0,1],['beaches','gazebos'])
if show_plot:
	plt.show()

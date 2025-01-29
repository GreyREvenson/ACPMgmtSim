#Create delta TN plots
for name, acp_trajectory_set in dt_acp_trajectory_set.items():
    xs = list(numpy.arange(0,len(acp_trajectory_set[0][0])+1,1))
    fig, ax1 = matplotlib.pyplot.subplots(1, 1)
    fig.set_size_inches(10, 10)
    for acp_trajectory in acp_trajectory_set:
        tmp = acp_trajectory[1].copy()
        for i in range(1,len(tmp)):
            tmp[i] = tmp[i]+tmp[i-1]
        ax1.plot(xs,[0]+tmp)
    ax1.set_title('ACP Impacts on TN ('+name+')')
    ax1.set_xlabel('# of ACPs Implemented')
    ax1.set_ylabel('Cumulative TN Load Reduction (kg/yr)')
    fig.tight_layout()
    #fig.savefig(os.path.join(dirOutput,'delta_tn_'+name+'.png'),dpi=600)

#Create delta TP plots
for name, acp_trajectory_set in dt_acp_trajectory_set.items():
    xs = list(numpy.arange(0,len(acp_trajectory_set[0][0])+1,1))
    fig, ax1 = matplotlib.pyplot.subplots(1, 1)
    fig.set_size_inches(10, 10)
    for acp_trajectory in acp_trajectory_set:
        tmp = acp_trajectory[2].copy()
        for i in range(1,len(tmp)):
            tmp[i] = tmp[i]+tmp[i-1]
        ax1.plot(xs,[0]+tmp)
    ax1.set_title('ACP Impacts on TP ('+name+')')
    ax1.set_xlabel('# of ACPs Implemented')
    ax1.set_ylabel('Cumulative TP Load Reduction (kg/yr)')
    fig.tight_layout()
    #fig.savefig(os.path.join(dirOutput,'delta_tp_'+name.replace('-','_')+'.png'),dpi=600)

#Create baseline TN yield boxplots
for name, acp_trajectory_set in dt_acp_trajectory_set.items():
    xs = list(numpy.arange(0,len(acp_trajectory_set[0][0])+1,1))
    bvs = list()
    for acp_trajectory in acp_trajectory_set:
        tmp = [acp_trajectory[3][i]/acp_trajectory[0][i] for i in range(len(acp_trajectory[3])) if acp_trajectory[0][i] > 0]
        bvs.extend(tmp)
    fig, ax1 = matplotlib.pyplot.subplots(1, 1)
    fig.set_size_inches(7.5, 3)
    ax1.boxplot(bvs,vert=False,whis=3.0)
    ax1.set_title('Baseline TN Yield (kg/ha/yr) ('+name+') (n='+str(len(bvs))+')')
    ax1.set_yticklabels('')
    fig.tight_layout()
    #fig.savefig(os.path.join(dirOutput,'baseline_yields_tn_'+name.replace('-','_')+'.png'),dpi=600)

###### Create baseline TP yield boxplots
for name, acp_trajectory_set in dt_acp_trajectory_set.items():
    xs = list(numpy.arange(0,len(acp_trajectory_set[0][0])+1,1))
    bvs = list()
    for acp_trajectory in acp_trajectory_set:
        tmp = [acp_trajectory[4][i]/acp_trajectory[0][i] for i in range(len(acp_trajectory[4])) if acp_trajectory[0][i] > 0]
        bvs.extend(tmp)
    fig, ax1 = matplotlib.pyplot.subplots(1, 1)
    fig.set_size_inches(7.5, 3)
    ax1.boxplot(bvs,vert=False,whis=3.0)
    ax1.set_title('Baseline TP Yield (kg/ha/yr) ('+name+') (n='+str(len(bvs))+')')
    ax1.set_yticklabels('')
    fig.tight_layout()
    #fig.savefig(os.path.join(dirOutput,'baseline_yields_tp_'+name.replace('-','_')+'.png'),dpi=600)

###### Create TN effectiveness boxplots
for name, acp_trajectory_set in dt_acp_trajectory_set.items():
    xs = list(numpy.arange(0,len(acp_trajectory_set[0][0])+1,1))
    evs = list()
    for acp_trajectory in acp_trajectory_set:
        tmp = [(-1*acp_trajectory[1][i])/acp_trajectory[3][i] for i in range(len(acp_trajectory[1])) if acp_trajectory[3][i] > 0]
        evs.extend(tmp)
    evs = numpy.array(evs)*100
    fig, ax1 = matplotlib.pyplot.subplots(1, 1)
    fig.set_size_inches(7.5, 3)
    ax1.boxplot(evs,vert=False,whis=3.0)
    ax1.set_title('Percent of TN Load Removed ('+name+') (n='+str(len(evs))+')')
    ax1.set_yticklabels('')
    fig.tight_layout()
    #fig.savefig(os.path.join(dirOutput,'effectiveness_coefficients_tn_'+name.replace('-','_')+'.png'),dpi=600)

###### Create TP effectiveness boxplots
for name, acp_trajectory_set in dt_acp_trajectory_set.items():
    xs = list(numpy.arange(0,len(acp_trajectory_set[0][0])+1,1))
    evs = list()
    for acp_trajectory in acp_trajectory_set:
        tmp = [(-1*acp_trajectory[2][i])/acp_trajectory[4][i] for i in range(len(acp_trajectory[2])) if acp_trajectory[4][i] > 0]
        evs.extend(tmp)
    evs = numpy.array(evs)*100
    fig, ax1 = matplotlib.pyplot.subplots(1, 1)
    fig.set_size_inches(7.5, 3)
    ax1.boxplot(evs,vert=False,whis=3.0)
    ax1.set_title('Percent of TP Load Removed ('+name+') (n='+str(len(evs))+')')
    ax1.set_yticklabels('')
    fig.tight_layout()
    #fig.savefig(os.path.join(dirOutput,'effectiveness_coefficients_tp_'+name.replace('-','_')+'.png'),dpi=600)

###### Create ACP pie chart
for name, acp_trajectory_set in dt_acp_trajectory_set.items():
    dtacpcount = dict()
    for acp_trajectory in acp_trajectory_set:
        acp_counts = {acp:[acp_trajectory[5].count(acp)] for acp in set(acp_trajectory[5])}
        for acp in acp_counts:
            if acp not in dtacpcount: dtacpcount[acp] = list()
            dtacpcount[acp].extend(acp_counts[acp])
    dtacpcount = {acp:numpy.average(dtacpcount[acp]) for acp in dtacpcount}
    fig, ax1 = matplotlib.pyplot.subplots(1, 1)
    fig.set_size_inches(10, 10)
    ax1.pie(dtacpcount.values(),labels=[key+' ('+"{:.2f}".format(dtacpcount[key])+')' for key in dtacpcount])
    ax1.set_title('ACP Trajectory Makeup (Average #)')
    fig.tight_layout()
    #fig.savefig(os.path.join(dirOutput,'acp_pie'+name.replace('-','_')+'.png'),dpi=600)
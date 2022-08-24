import numpy as np
import pandas as pd
import re
from statistics import mean
import glob
import os
import sys


def getcontigsinregion(xmap):
    xmap_headercols = ['XmapEntryID', 'QryContigID', 'RefContigID', 'QryStartPos', 'QryEndPos', 'RefStartPos',
                       'RefEndPos', 'Orientation', 'Confidence', 'HitEnum', 'QryLen', 'RefLen', 'LabelChannel',
                       'Alignment', 'MapWt']
    contigx = pd.read_csv(xmap, sep='\t', comment='#', header=None, names=xmap_headercols, index_col=None)
    contigsneeded = list(set(contigx['QryContigID']))
    return contigsneeded


def getrelevantcontigcoordinates(start: int, stop, cmapID, alignment, contigqmap, contigrmap):
    # open contig files
    qmap_headercols = ['CMapId', 'ContigLength', 'NumSites', 'SiteID', 'LabelChannel', 'Position', 'StdDev', 'Coverage',
                       'Occurrence', 'ChimQuality', 'SegDupL', 'SegDupR', 'FragileL', 'FragileR', 'OutlierFrac',
                       'ChimNorm', 'Mask']
    contigq = pd.read_csv(contigqmap, sep='\t', comment='#', header=None, names=qmap_headercols)

    rmap_headercols = ['CMapId', 'ContigLength', 'NumSites', 'SiteID', 'LabelChannel', 'Position', 'StdDev', 'Coverage',
                       'Occurrence', 'ChimQuality', 'SegDupL', 'SegDupR', 'FragileL', 'FragileR', 'OutlierFrac',
                       'ChimNorm']
    contigr = pd.read_csv(contigrmap, sep='\t', comment='#', header=None, names=rmap_headercols)

    # get relevant positions on contig
    sub_qmap = contigq.query("CMapId == @cmapID")
    querydict = pd.Series(sub_qmap.Position.values, index=sub_qmap.SiteID).to_dict()
    chrom_ref = contigr.query("CMapId == @chrom")
    refdict = pd.Series(chrom_ref.SiteID.values, index=chrom_ref.Position).to_dict()
    keylist = list(refdict.keys())
    refstart_nicksite = closest(keylist, start)
    refend_nicksite = closest(keylist, stop)
    alignmentlist = re.split('\(|\)', alignment)
    alignmentlist = list(filter(None, alignmentlist))
    alignment_df = pd.DataFrame(alignmentlist)
    alignment_df[['refsite', 'querysite']] = alignment_df[0].str.split(',', expand=True)
    alignment_df.drop(columns=0, inplace=True)
    alignment_df = alignment_df.set_index('refsite', drop=False)
    # alignment_df.to_csv('temp_alignmentdf.tsv', sep='\t')
    reflist = alignment_df.index.values.tolist()
    refdictstart = refdict[refstart_nicksite]
    refdictend = refdict[refend_nicksite]
    # print(refdictstart, refdictend)
    while str(refdictstart) not in reflist:
        refdictstart = refdictstart - 1
    while str(refdictend) not in reflist:
        refdictend = refdictend + 1
    # print(refdictstart, refdictend)
    contigqstart = alignment_df.loc[str(refdictstart), 'querysite']
    contigqend = alignment_df.loc[str(refdictend), 'querysite']
    # print(contigqstart, contigqend)
    contigstartpos = min([querydict[int(contigqstart)], querydict[int(contigqend)]])
    contigendpos = max([querydict[int(contigqstart)], querydict[int(contigqend)]])

    return contigstartpos, contigendpos


def closest(lst, K):
    return lst[min(range(len(lst)), key=lambda i: abs(lst[i] - K))]


def getredlabelindices(labellist):
    indices = []
    for i in range(len(labellist)):
        if labellist[i] == 1:
            indices.append(i)
    return indices


def createalignmentdf(alignment, col1, col2):
    contigalignmentlist = re.split('\(|\)', alignment)
    contigalignmentlist = list(filter(None, contigalignmentlist))
    temp_df = pd.DataFrame(contigalignmentlist)
    temp_df[[col1, col2]] = temp_df[0].str.split(',', expand=True)
    temp_df.drop(columns=0, inplace=True)
    return temp_df


def getunderthresholdlocs(distlist):
    indices = []
    for i in range(len(distlist)):
        if labellist[i] <= 1250:
            indices.append(i)
    return indices


def main(sampleID):
    print("Analyzing " + sampleID)
    # pd.set_option('display.max_columns', None)
    # sampleID = 'NA12878_Inversion-1_assembly_pipeline_results'

    # # making output directory
    outdir = sampleID + '_redlabeloutput'
    os.makedirs(outdir)

    mapdir = sampleID + '/contigs/exp_refineFinal1_sv/merged_smaps/'
    contigqmap = mapdir + 'exp_refineFinal1_merged_q.cmap'
    contigxmap = mapdir + 'exp_refineFinal1_merged.xmap'
    contigrmap = mapdir + 'exp_refineFinal1_merged_r.cmap'
    cmapheadercols = ['CMapId', 'ContigLength', 'NumSites', 'SiteID',
                      'LabelChannel', 'Position', 'StdDev', 'Coverage',
                      'Occurrence', 'ChimQuality', 'SegDupL', 'SegDupR',
                      'FragileL', 'FragileR', 'OutlierFrac', 'ChimNorm']
    contigrmapdf = pd.read_csv(contigrmap, sep='\t', header=None, names=cmapheadercols, comment='#')
    headercols = ['XmapEntryID', 'QryContigID', 'RefContigID', 'QryStartPos',
                  'QryEndPos', 'RefStartPos', 'RefEndPos', 'Orientation',
                  'Confidence', 'HitEnum', 'QryLen', 'RefLen', 'LabelChannel', 'Alignment']
    contigxmapdf = pd.read_csv(contigxmap, sep='\t', header=None, names=headercols, comment='#')

    molmapdir = sampleID + '/contigs/exp_refineFinal1/alignmol/merge/'
    headercols.append('MapWt')
    cmapheadercols = cmapheadercols + ['Mask', 'GmeanSNR', 'lnSNRsd']

    # donefiles = []
    # for filename in glob.glob('NA12878_Inversion-1_assembly_pipeline_results_redlabeloutput/*.tsv'):
    #     contigID = int(filename.split('/')[-1].split('.')[0].split('_')[0].split('g')[-1])
    #     donefiles.append(contigID)
    # print(list(set(relevantcontigs) - set(donefiles)))

    relevantcontigs = []
    for filename in glob.glob(molmapdir + '*.xmap'):
        contigID = int(filename.split('/')[-1].split('.')[0].split('contig')[-1])
        relevantcontigs.append(contigID)
    # # print(len(set(relevantcontigs)))
    # relevantcontigs = [322]
    for contig in relevantcontigs:
        print('contig', contig)
        # if contig == 1590:
        molxmap = molmapdir + 'exp_refineFinal1_contig' + str(contig) + '.xmap'
        molqmap = molmapdir + 'exp_refineFinal1_contig' + str(contig) + '_q.cmap'
        molrmap = molmapdir + 'exp_refineFinal1_contig' + str(contig) + '_r.cmap'

        molxmapdf = pd.read_csv(molxmap, sep='\t', header=None, names=headercols, comment='#')
        molqmapdf = pd.read_csv(molqmap, sep='\t', header=None, names=cmapheadercols, comment='#')
        molrmapdf = pd.read_csv(molrmap, sep='\t', header=None, names=cmapheadercols, comment='#')
        relevantdf = contigxmapdf.query("QryContigID == @contig")
        if relevantdf.shape[0] > 0:
            relevantdf['alignedlen'] = relevantdf['RefEndPos'] - relevantdf['RefStartPos']
            longestalignment = max(relevantdf['alignedlen'])
            relevantdf = relevantdf.query("alignedlen == @longestalignment")
            contigorientation = list(relevantdf['Orientation'])[0]
            contigalignment = relevantdf['Alignment']
            chromosome = list(relevantdf['RefContigID'])[0]

            contigsubxmap_df = contigxmapdf.query("QryContigID == @contig")
            # # print(contigsubxmap_df)
            alignmentstrings = contigsubxmap_df['Alignment']
            contigalignment_df = pd.DataFrame()
            for alignstr in alignmentstrings:
                temp_df = createalignmentdf(alignstr, 'refsite', 'querysite')
                contigalignment_df = contigalignment_df.append(temp_df, ignore_index=True)
                alignedcontigsites = list(contigalignment_df['querysite'])
            contigalignment_df.to_csv('contigalignmentdf.tsv', sep='\t')

            molslist = list(set(molxmapdf['QryContigID']))
            # print(molslist)
            dflist = []
            for mol in molslist:
                # print(mol)
                # if mol == 1091114:
                ## get indices of red labels in list
                relevantqmaprows = molqmapdf.query("CMapId == @mol")
                relevantqmaprows = relevantqmaprows.reset_index(drop=True)
                labels = list(relevantqmaprows['LabelChannel'])
                redindices = getredlabelindices(labels)
                if len(redindices) >= 1:
                    # print(mol)
                    relevantmoldf = molxmapdf.query("QryContigID == @mol")
                    molorientation = list(relevantmoldf['Orientation'])[0]
                    molalignment = list(relevantmoldf['Alignment'])[0]
                    molalignment_df = createalignmentdf(molalignment, 'contigsite', 'molsite')
                    molalignment_df = molalignment_df.set_index('molsite', drop=False)
                    molalignment_df.to_csv('tempmolalignment.tsv', sep='\t')
                    alignedgreensites = list(molalignment_df['molsite'])
                    # print(alignedgreensites)
                    c = 0
                    # print(redindices)
                    for index in redindices:
                        c += 1
                        if index > 0 and c < index:
                            dfdict = {}
                            lastgreenindex = index - 1
                            distancefromgreen = relevantqmaprows.loc[index, 'Position'] - relevantqmaprows.loc[
                                lastgreenindex, 'Position']
                            lastgreensiteID = relevantqmaprows.loc[lastgreenindex - (c-1), 'SiteID']
                            nextgreenindex = index + 1
                            distancetogreen = relevantqmaprows.loc[nextgreenindex, 'Position'] - relevantqmaprows.loc[
                                index, 'Position']
                            nextgreensiteID = relevantqmaprows.loc[nextgreenindex - c, 'SiteID']
                            # print(c, 'index = ', index, 'lastgreenindex = ', lastgreenindex,
                            #       'lastgreensiteID = ', lastgreensiteID, 'nextgreenindex = ', nextgreenindex,
                            #       'nextgreensiteID = ', nextgreensiteID)
                            if distancefromgreen < distancetogreen or str(nextgreensiteID) not in alignedgreensites:
                                # print('distancefromgreen < distancetogreen')
                                if str(lastgreensiteID) in alignedgreensites:
                                    # print('lastgreen in aligned green')
                                    anchorsite = lastgreensiteID
                                    contigsiteID = molalignment_df.loc[str(lastgreensiteID), 'contigsite']
                                    # print(contigsiteID)
                                    try:
                                        contigpos = list(molrmapdf.query("CMapId == @contig and SiteID == @contigsiteID")[
                                                'Position'])[0]
                                        refsiteID = list(contigalignment_df.query('querysite == @contigsiteID')['refsite'])[
                                            0]
                                        refpos = list(contigrmapdf.query("CMapId == @chromosome and SiteID == @refsiteID")[
                                                'Position'])[0]
                                    except:
                                        # print('contigsite not aligned to ref')
                                        # unalignedgreencontigsite.append(mol)
                                        continue
                                    print('chromosome: ', chromosome, 'index: ', index,
                                          'c: ', c, 'distancefromgreen: ', distancefromgreen,
                                          'distancetogreen: ', distancetogreen,
                                          'lastgreenindex: ', lastgreenindex,
                                          'lastgreensiteID: ', lastgreensiteID, 'contigsiteID: ',
                                          contigsiteID, 'refsiteID: ', refsiteID,'refpos: ', refpos,
                                          'orientation: ', molorientation, 'contigpos: ', contigpos)
                                    if contigorientation == '+':
                                        if molorientation == '+':
                                            redlabelonref = refpos + distancefromgreen
                                        if molorientation == '-':
                                            redlabelonref = refpos - distancefromgreen
                                    elif contigorientation == '-':
                                        if molorientation == '+':
                                            redlabelonref = refpos - distancefromgreen
                                        elif molorientation == '-':
                                            redlabelonref = refpos + distancefromgreen
                                    dfdict['MolID'] = mol
                                    dfdict['ContigID'] = contig
                                    dfdict['LabelChannel'] = 1
                                    dfdict['ClosestGreenOnContig'] = contigpos
                                    dfdict['Position'] = redlabelonref
                                    dflist.append(dfdict)
                                    print('redlabelonref: ', redlabelonref)

                            if distancetogreen < distancefromgreen or str(lastgreensiteID) not in alignedgreensites:
                                # print('distancetogreen < distancefromgreen')
                                if str(nextgreensiteID) in alignedgreensites:
                                    # print('nextgreeninalignedgreen')
                                    anchorsite = nextgreensiteID
                                    contigsiteID = molalignment_df.loc[str(nextgreensiteID), 'contigsite']
                                    try:
                                        contigpos = list(molrmapdf.query("CMapId == @contig and SiteID == @contigsiteID")[
                                                         'Position'])[0]
                                        refsiteID = list(contigalignment_df.query('querysite == @contigsiteID')['refsite'])[
                                            0]
                                        refpos = list(contigrmapdf.query("CMapId == @chromosome and SiteID == @refsiteID")[
                                                'Position'])[0]
                                    except:
                                        # print('contigsite not aligned to ref')
                                        # unalignedgreencontigsite.append(mol)
                                        continue
                                    print('chromosome: ', chromosome, 'index: ', index,
                                          'c: ', c, 'distancefromgreen: ', distancefromgreen,
                                          'distancetogreen: ', distancetogreen,
                                          'nextgreenindex: ', nextgreenindex,
                                          'nextgreensiteID: ', nextgreensiteID, 'contigsiteID: ',
                                          contigsiteID, 'refsiteID: ', refsiteID, 'refpos: ', refpos,
                                          'orientation: ', molorientation, 'contigpos: ', contigpos)

                                    # print(molorientation, refpos)
                                    if contigorientation == '+':
                                        if molorientation == '+':
                                            redlabelonref = refpos - distancetogreen
                                        if molorientation == '-':
                                            redlabelonref = refpos + distancetogreen
                                    elif contigorientation == '-':
                                        if molorientation == '+':
                                            redlabelonref = refpos + distancetogreen
                                        elif molorientation == '-':
                                            redlabelonref = refpos - distancetogreen
                                    dfdict['MolID'] = mol
                                    dfdict['ContigID'] = contig
                                    dfdict['LabelChannel'] = 1
                                    dfdict['ClosestGreenOnContig'] = contigpos
                                    dfdict['Position'] = redlabelonref
                                    dflist.append(dfdict)
                                    print('redlabelonref: ', redlabelonref)

            if len(dflist) >= 1:
                finaldf = pd.DataFrame(dflist)
                # print(finaldf.head(5))
                try:
                    finaldf = finaldf.sort_values(by=['Position'])
                except:
                    print("couldn't sort dataframe for contig", contig)
                    continue
                finaldf.to_csv(outdir + '/Contig' + str(contig) + '_redlabelpositions.tsv', sep='\t', index=None)

    print("Running Part 2 for " + sampleID)

    finaldir = sampleID + '_groupedredlabelfiles'
    os.makedirs(finaldir)

    for outfile in glob.glob(outdir + '/*_redlabelpositions.tsv'):
        contig = outfile.split('/')[-1].split('_')[0].split('g')[-1]
        print(contig)
        contigxmap = mapdir + 'exp_refineFinal1_merged.xmap'
        headercols = ['XmapEntryID', 'QryContigID', 'RefContigID', 'QryStartPos',
                      'QryEndPos', 'RefStartPos', 'RefEndPos', 'Orientation',
                      'Confidence', 'HitEnum', 'QryLen', 'RefLen', 'LabelChannel', 'Alignment']
        contigxmapdf = pd.read_csv(contigxmap, sep='\t', header=None, names=headercols, comment='#')


        headercols.append('MapWt')
        molxmap = molmapdir + 'exp_refineFinal1_contig' + str(contig) + '.xmap'
        molxmapdf = pd.read_csv(molxmap, sep='\t', header=None, names=headercols, comment='#')

        relevantdf = contigxmapdf.query("QryContigID == @contig")

        relevantdf['alignedlen'] = relevantdf['RefEndPos'] - relevantdf['RefStartPos']
        longestalignment = max(relevantdf['alignedlen'])
        relevantdf = relevantdf.query("alignedlen == @longestalignment")
        chromosome = list(relevantdf['RefContigID'])[0]

        finaldf = pd.read_csv(outfile, sep='\t')
        # finaldf = finaldf.query("18000000 > Position > 17000000")
        finaldf = finaldf.sort_values(by=['Position']).reset_index(drop=True)
        # print(finaldf)
        poslist = finaldf['Position']
        distlist = list(np.diff(poslist))
        distlist.append(0)
        finaldf['consecutivedist'] = distlist
        # # print(finaldf.head(5))
        positionlist = []
        contigposlist = []
        distancelist = []
        inagroup = False
        grouplist = []
        molids = []
        for index, row in finaldf.iterrows():
            if row['consecutivedist'] <= 1250 and inagroup == False:
                firstredlabel = row['Position']
                firstredcontigpos = row['ClosestGreenOnContig']
                positionlist = []
                contigposlist = []
                distancelist = []
                molids = []
                groupdict = {}
                positionlist.append(row['Position'])
                contigposlist.append(row['ClosestGreenOnContig'])
                distancelist.append(row['consecutivedist'])
                molids.append(row['MolID'])
                inagroup = True
            elif row['consecutivedist'] <= 1250 and inagroup == True:
                positionlist.append(row['Position'])
                contigposlist.append(row['ClosestGreenOnContig'])
                distancelist.append(row['consecutivedist'])
                molids.append(row['MolID'])
                # currentredlabel = index
            elif row['consecutivedist'] > 1250 and inagroup == True:
                molids.append(row['MolID'])
                lastredlabel = row['Position']
                lastredcontigpos = row['ClosestGreenOnContig']
                positionlist.append(row['Position'])
                contigposlist.append(row['ClosestGreenOnContig'])
                molids.append(row['MolID'])
                groupdict['ContigID'] = contig
                groupdict['firstpos'] = firstredlabel
                groupdict['lastpos'] = lastredlabel
                groupdict['firstposoncontig'] = firstredcontigpos
                groupdict['lastposoncontig'] = lastredcontigpos
                groupdict['positions_in_group'] = positionlist
                groupdict['distances_in_group'] = distancelist
                groupdict['mols_in_group'] = molids
                grouplist.append(groupdict)
                inagroup = False
                positionlist = []
                distancelist = []
                molids = []
            else:
                inagroup = False
        # print(grouplist)
        if len(grouplist) >= 1:
            writetodf = False
            for group in grouplist:
                if len(group['positions_in_group']) > 3:
                    writetodf = True
                    distlist = group['distances_in_group']
                    totaldist = sum(distlist)
                    group['total_distance'] = totaldist
                    # print(group)
                    contigstartpos = min(group['firstposoncontig'], group['lastposoncontig'])
                    contigendpos = max(group['firstposoncontig'], group['lastposoncontig'])
                    # print(group['firstpos'], group['lastpos'], contigstartpos, contigendpos)
                    molsdf = molxmapdf.query("@contigstartpos <= RefEndPos and @contigendpos >= RefStartPos")
                    molsinregion = list(molsdf['QryContigID'])
                    group['all_mols_in_region'] = molsinregion
                    totalmolnum = len(molsinregion)
                    group['total_molnum'] = totalmolnum
                    redsitenum = len(group['positions_in_group'])
                    group['redsitenum'] = redsitenum
                    percredmols = redsitenum / totalmolnum
                    group['perc_redmols'] = percredmols
                    group['Chromosome'] = chromosome

            groupdf = pd.DataFrame(grouplist)
            # print(groupdf[['firstpos', 'lastpos']])
            groupdf = groupdf.dropna()
            # print(groupdf)
            if writetodf:
                groupdf = groupdf[['Chromosome', 'firstpos', 'lastpos', 'ContigID',
                                   'redsitenum', 'total_molnum', 'perc_redmols',
                                   'mols_in_group', 'all_mols_in_region',
                                   'positions_in_group', 'distances_in_group']]
                groupdf.to_csv(finaldir + '/Contig' + str(contig) + '_redlabelgroupstats.tsv', sep='\t', index=None)


if __name__ == '__main__':
    sampleID = sys.argv[1]
    main(sampleID)

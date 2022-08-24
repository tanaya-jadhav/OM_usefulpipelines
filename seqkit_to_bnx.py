import pandas as pd
import numpy as np
import random
import matplotlib.pyplot as plt


def create_bnx_file(inputdf, chromosome, centromerestart, centromereend, chromregion):
    outfile = 'chr' + str(chromosome) + '_' + chromregion + 'centromere_fakemolecules.bnx'
    # outfile = 'chr1_aftercentromere_fakemolecules.bnx'
    output_file = open(outfile, "w")
    headerfile = 'newbnxheader.txt'
    with open(headerfile, 'r') as hfile:
        output_file.write(hfile.read())

    df = pd.read_csv(inputdf, sep='\t')
    if chromregion == 'after':
        subdf = df.query("patternName == 'BstNBI' and start > @centromereend")
        molcounter = int(str(chromosome) + '002' + '0000001')
    elif chromregion == 'before':
        subdf = df.query("patternName == 'BstNBI' and start < @centromerestart")
        molcounter = int(str(chromosome) + '001' + '0000001')
    elif chromregion == 'centromere':
        subdf = df.query("patternName == 'BstNBI' and start > @centromerestart and start < @centromereend")
        molcounter = int(str(chromosome) + '003' + '0000001')
    else:
        print("chromregion needs to be specified as either 'before' or 'after'")
        return

    subdf = subdf.sort_values(by=['start'])
    subdf = subdf.reset_index()
    # df['dist'] = np.diff(df['start'])
    # subdf.to_csv('temp.tsv', sep='\t')
    # breakpoints = []
    counter = 0
    fragment = []
    molecules = []
    for index, row in subdf.iterrows():
        counter += 1
        # print(index, counter)
        if counter == 1:
            currentstrand = row['strand']
            currentpos = row['start']
            fragment.append(currentpos)
        # breakpoints.append(currentpos)
        elif counter == len(subdf):
            currentpos = row['start']
            # print('last pos', currentpos)
            # breakpoints.append(currentpos)
            fragment.append(currentpos)
            # print(fragment)
            molecules.append(np.diff(fragment))
            fragment = []
        else:
            currentstrand = row['strand']
            previousstrand = subdf.loc[index - 1, 'strand']
            currentpos = row['start']
            previouspos = subdf.loc[index - 1, 'start']
            if currentstrand == previousstrand:
                fragment.append(currentpos)
            elif currentstrand != previousstrand and currentpos - previouspos > 500:
                fragment.append(currentpos)
            else:
                fragment.append(currentpos)
                # print(fragment)
                molecules.append(np.diff(fragment))
                fragment = []
                # breakpoints.append(previouspos)
    # print(breakpoints)

    distlist = [list(mol) for mol in molecules if len(mol) > 1]
    # print(len(distlist))

    with open(outfile, 'a') as o:
        # molcounter = 10020000001
        for mol in distlist:
            mollength = sum(mol) + 2000
            mol = [x for x in mol if x > 0]
            labelpos = [0]
            for x in mol:
                labelpos.sort()
                pos = labelpos[-1] + x
                labelpos.append(pos)
            labelpos = labelpos[1:]
            hrow = [0, molcounter, mollength, round(random.uniform(1500, 4000), 2), round(random.uniform(90, 600), 3),
                    len(labelpos), molcounter, 1, -1,
                    'chips,SN_NT5FCKWLPQYERNWU,Run_00f10215-19fa-40f1-ae2d-9c05a53fd7d3,0',
                    1, 37, random.randrange(1, 200), random.randrange(1, 5), random.randrange(40, 900),
                    random.randrange(1, 1500),
                    random.randrange(1, 5), random.randrange(1, 1000), random.randrange(250, 2000)]
            intensities = random.sample(range(10000, 40000), len(labelpos))
            intensities = [x / 100 for x in intensities]
            snr = random.sample(range(1000000, 2500000), len(labelpos))
            snr = [x / 100000 for x in snr]
            labelpos.append(mollength)
            hrow = [str(x) for x in hrow]
            hrow = '\t'.join(hrow)
            o.write(hrow + '\n')

            labelpos.insert(0, 1)
            labelpos = [str(x) for x in labelpos]
            labelpos = '\t'.join(labelpos)
            o.write(labelpos + '\n')

            snr.insert(0, 'QX11')
            snr = [str(x) for x in snr]
            snr = '\t'.join(snr)
            o.write(snr + '\n')

            intensities.insert(0, 'QX12')
            intensities = [str(x) for x in intensities]
            intensities = '\t'.join(intensities)
            o.write(intensities + '\n')
            molcounter += 1


def main():
    # # creates a fake bnx file by fragmenting if 2 sites are on opposite strands and <=500bp apart.
    create_bnx_file('wholegenomefragmentsizes/chr1_vs_chosenenzymes.txt', 1, 122026459, 124932724, 'centromere')


    # molfile = 'allcentromeremols_over1Mb.bnx'
    # molfile = 'refaligner_chr1aftercentromereNEW_repeat.bnx'
    # with open(molfile, 'r') as mf:
    #     lines = mf.readlines()
    #
    # counter = 0
    # outlist = []
    # for line in lines:
    #     linesplit = line.strip('\n').split('\t')
    #     if linesplit[0] == '0':
    #         molID = linesplit[1]
    #     elif linesplit[0] == '1':
    #         greenpositions = [float(pos.strip('\n')) for pos in linesplit[1:-1]]
    #     elif linesplit[0] == 'QX12':
    #         # print(molID)
    #         # if molID == '1041969':
    #             # print(molID, 'FOUND')
    #         greendistances = np.diff(greenpositions)
    #         # reddistances = np.diff(redpositions)
    #         # x = reddistances
    #         y = greendistances
    #         bins = list(range(0, 10000, 250))
    #         # n, bins, patches = plt.hist(x, bins=bins)
    #         n2, bins2, patches2 = plt.hist(y, bins=bins)
    #         plt.clf()
    #         # plt.plot(bins[1:], n, marker='o', color='r')
    #
    #         plt.plot(bins2[1:], n2, marker='o', color='g')
    #         # figname = 'moleculeprobehistograms/' + molID + '_histogram.png'
    #         # plt.savefig(figname)
    #         # plt.clf()
    #         plt.show()
    #         break


if __name__ == '__main__':
    main()
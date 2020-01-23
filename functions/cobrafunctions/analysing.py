from __future__ import absolute_import
import cobra, os, matplotlib, collections, itertools, errno, math

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from tqdm import tqdm_notebook as tqdm

from cobrafunctions.buildingediting import *
from cobrafunctions.adjustedcobrafunctions import *


def nameexperiment(experimentname):
    '''
    This function is used to define the name of the experiment in the current jupyter notebook and create a results folder into which results will be deposited.
    It takes a single argument, experimentname, a string.
    '''
    savelocation = createfolder("Results/" + experimentname)
    inputlocation = ("Inputs/")
    try:
        saveinputs(inputlocation, savelocation)
    except:
        pass
    return savelocation, inputlocation

def createfolder(path):
    '''
    This function is used by nameexperiment to create the folder for a new experiment, and ask the user
    if they want to overwrite if there is already an experiment with that name.
    It takes a single argument, path, a string, which defines where the new folder is created.
    When used with nameexperiment path is the name of a folder inside the "Results" folder, defined by the argument to the nameexperiment function.
    '''
    if os.path.isdir(path):
        while True:
            query = raw_input('This experiment already exists, continue and overwrite?')
            Fl = query[0].lower()
            if query == '' or not Fl in ['y','n']:
                print('Please answer with yes or no!')
            else:
                break
        if Fl == 'y':
            print("OK, Continuing and overwriting experiment")
            try:
                os.makedirs(path)
            except OSError as exception:
                if exception.errno != errno.EEXIST:
                    raise
            return path
        if Fl == 'n':
            experimentname = raw_input("Ok, please enter the name of the new experiment")
            path = "Results/" + experimentname
            try:
                os.makedirs(path)
            except OSError as exception:
                if exception.errno != errno.EEXIST:
                    raise
    else:
        try:
                os.makedirs(path)
        except OSError as exception:
            if exception.errno != errno.EEXIST:
                raise
        return path

def saveinputs(inputlocation, savelocation):
    '''
    This function is used by nameexperiment in order to copy the inputs folder into the results folder of that experiment.
    This enables the experimentor to track the inputs they used to get those results.
    It takes two arguments, the path to the input location, and the path to the savelocation where to copy this folder.
    '''
    import shutil
    shutil.copytree(inputlocation, savelocation + "/Inputs")
    return os.path.join(savelocation, "/Inputs")

def getreactioncompartment(reaction):
    '''
    This is a very simple function which takes a reaction id and returns the last three tags to give a compartment string
    '''
    compartment = "_".join(reaction.id.split("_")[-3:])
    return compartment

def plotlinkers(model, reactionsdf, coloursdict = "", savelocation=None, name=None):
    '''
    This function returns a matplotlib figure containing a plot of how linker concentrations change over the phases of the model.
    This function takes two arguments and 3 optional arguments:
        model: The model that is to be used for plotting
        reactionsdf: a dataframe containing the results of a solution of the model, with index reaction ids and columns for reaction
            flux then maximum flux and minimum flux according to fva.
        coloursdict: an optional argument to specify the colours for the lines in the plot
        savelocation: path defining where the plot should be saved
        name: a name for the plot, otherwise the default plot name will be linkersplot. This will be numbered if there already
            exists a linkersplot in the specified path
    '''

    import itertools

    def flip(items, ncol):
        return itertools.chain(*[items[i::ncol] for i in range(ncol)])

    reactionsdict = {}

    if coloursdict == "":
        coloursdict = {
            'Cl': '#5C6BC0',
            'STARCH':'#8D6E63' ,
            'K': '#AB47BC',
            'MAL': '#2E7D32',
            'SUCROSE': '#EF5350',
            'CIT': '#FF7043',
            'FRU': '#26C6DA',
            'GLC': '#9CCC65',
            'Other': '#BDBDBD',
         }
    for reaction, row in reactionsdf.iterrows():
        if reaction.__contains__("Linker"):
            if "_a_" not in reaction:
                try:
                    reactionsdict[reaction[-13:-9]][reaction[:-14]][int(reaction[-1:])] = (row["pfba"], row["minimum"], row["maximum"])
                except:
                    try:
                        reactionsdict[reaction[-13:-9]][reaction[:-14]] = {}
                        reactionsdict[reaction[-13:-9]][reaction[:-14]][int(reaction[-1:])] = (row["pfba"], row["minimum"], row["maximum"])
                    except:
                        reactionsdict[reaction[-13:-9]] = {}
                        reactionsdict[reaction[-13:-9]][reaction[:-14]] = {}
                        reactionsdict[reaction[-13:-9]][reaction[:-14]][int(reaction[-1:])] = (row["pfba"], row["minimum"], row["maximum"])

    coloursdict = collections.OrderedDict(sorted(coloursdict.items()))

    matplotlib.rcParams.update({'font.size': 30})

    matplotlib.rcParams['lines.linewidth'] = 5

    matplotlib.rcParams['axes.linewidth'] = 3

    plt.rcParams["font.family"] = "Arial"

    titlelist = ["Guard Cell Plastid", "Guard Cell Vacuole", "Guard Cell Cytosol", "Mesophyll Cell Plastid",  "Mesophyll Cell Vacuole"]
    unitslist = ["(mmoles)", "(M)", "(M)", "(mmoles)", "(mM)"]
    complist = ["p_gc", "v_gc", "c_gc", "p_me", "v_me",]

    for compartment in ["v_me", "v_gc"]:
        amal = reactionsdict[compartment].pop("aMAL")
        for i in range(1, checknumberofmodels(model)+1):
            reactionsdict[compartment]["MAL"][i] = reactionsdict[compartment]["MAL"][i] + amal[i]
        acit = reactionsdict[compartment].pop("aCIT")
        for i in range(1, checknumberofmodels(model)+1):
            reactionsdict[compartment]["CIT"][i] = reactionsdict[compartment]["CIT"][i] + acit[i]
    apertures = [6,17,6,6]
    gcvolume = []
    for i in range(checknumberofmodels(model)):
        gcvolume.append((apertures[i]*200+4000)*17.2*(10.0**(-8.0)))

    mevolume = 0.284
    v_mevolume = mevolume*0.9
    linkerconcs = reactionsdict.copy()
    for compartment in linkerconcs:
        for reaction in linkerconcs[compartment]:
            for phase in linkerconcs[compartment][reaction]:
                if "v_me" in compartment:
                    linkerconcs[compartment][reaction][phase] = ((linkerconcs[compartment][reaction][phase][0]*0.001)/v_mevolume*1000, (linkerconcs[compartment][reaction][phase][1]*0.001)/v_mevolume*1000, (linkerconcs[compartment][reaction][phase][2]*0.001)/v_mevolume*1000)
                elif "v_gc" in compartment:
                    linkerconcs[compartment][reaction][phase] = ((linkerconcs[compartment][reaction][phase][0]*0.001)/(gcvolume[phase-1]*0.9), (linkerconcs[compartment][reaction][phase][1]*0.001)/(gcvolume[phase-1]*0.9), (linkerconcs[compartment][reaction][phase][2]*0.001)/(gcvolume[phase-1]*0.9))
                elif "c_gc" in compartment:
                    linkerconcs[compartment][reaction][phase] = ((linkerconcs[compartment][reaction][phase][0]*0.001)/(gcvolume[phase-1]*0.1), (linkerconcs[compartment][reaction][phase][1]*0.001)/(gcvolume[phase-1]*0.1), (linkerconcs[compartment][reaction][phase][2]*0.001)/(gcvolume[phase-1]*0.1))
    reactionsdict = linkerconcs


    legendset = set()

    fig, axs = plt.subplots(6,1,figsize=(10, 25), sharex = True)
    i=0
    c=0

    phaseconversion = [0,6,6.5,18,24]
    for compartment in complist:
        nonzero = False
        m = 0
        for reaction in reactionsdict[compartment]:
            try:
                colour = coloursdict[reaction]
            except:
                colour = coloursdict['Other']
            x = phaseconversion
            y0 = []
            y1 = []
            y2 = []
            for phase in reactionsdict[compartment][reaction].values():
                y0.append(phase[0])
                y1.append(phase[1])
                y2.append(phase[2])

            y0.insert(0, reactionsdict[compartment][reaction][4][0])
            y1.insert(0, reactionsdict[compartment][reaction][4][1])
            y2.insert(0, reactionsdict[compartment][reaction][4][2])
            if max(y2) > 0.001 and ("c_gc" in compartment or "v_gc" in compartment) or max(y2) > 0.0001 and "c_gc" not in compartment and "v_gc" not in compartment:
                if "sum" not in reaction:
                    axs[i].plot(x,y0, color = colour,)
                    axs[i].fill_between(x, y1, y2, color = colour, alpha = 0.2)
                    legendset.add(reaction)
                    nonzero = True
            else:
                axs[i].plot(x,[0,0,0,0,0], color = colour,)
            m += 1
        axs[i].set_ylabel(unitslist[c], labelpad = 20)
        axs[i].set_title(titlelist[c], pad = 10)
        axs[i].tick_params(axis='y', which='major',length = 15, width = 3)
        axs[i].tick_params(axis='x', which='major',length = 10, width = 3)
        plt.xlabel("Time of Day")

        if nonzero == False:
            axs[i].set_ylim(-0.05, 1)


        c += 1
        i += 1
    osmolarity = [0.234091748752047, 1.9583538166251, 0.234091748752047, 0.23409174875204702]
    aperture = [6,6,17,6,6]
    axs[i].plot(phaseconversion, aperture, color = 'black')
    axs[i].set_ylabel(u"Aperture Diameter \n (\u03bcm)", labelpad = 20)
    axs[i].set_title("Stomatal Aperture", pad = 10)
    axs[i].set_yticks([6,17])
    axs[i].set_xticks([0,6,18,24])
    axs[i].set_xticklabels(["Midnight", "Sunrise", "Sunset", "Midnight",], rotation = 90)
    axs[i].set_ylim(0,20)
    axs[i].tick_params(axis='y', which='major',length = 15, width = 3)
    axs[i].tick_params(axis='x', which='major',length = 10, width = 3)


    fig.subplots_adjust(hspace=0.3, wspace = 0.3, left=0.2, )

    namesdict = {
        "CIT" : "Citrate",
        "Cl" : r"Cl$^{-}$",
        "STARCH" : "Starch",
        "K" : r"K$^{+}$",
        "MAL" : "Malate",
        "GLC" : "Glucose",
        "FRU" : "Fructose",
        "SUCROSE" : "Sucrose"
    }

    from matplotlib.lines import Line2D
    legendlist = sorted(legendset)
    legendlines = []
    legendnames = []
    for reaction in legendlist:
        includeother = False
        try:
            legendlines.append(Line2D([0], [0], color=coloursdict[reaction]))
            legendnames.append(namesdict[reaction])
        except:
            includeother = True
    if includeother == True:
        legendnames.append("Amino Acids")
        legendlines.append(Line2D([0], [0], color='#BDBDBD'))

    fig.text(0.0, 0.8, "Concentration or amount of ion/metabolite in titular compartment", rotation = 90)

    plt.legend(flip(legendlines, 3), flip(legendnames,3), loc = 'upper center', bbox_to_anchor=(0.5, -1), ncol = 3, prop={'size': 25})

    if name != None and savelocation != None:
        fig.savefig(os.path.join(savelocation, name +  ".svg"), format='svg',  bbox_inches='tight')
    return fig

def linkersplotfrommodel(model, savelocation = None, name = None):
    '''
    This function will return an optimised model and a matplotlib figure containing a plot of osmolytes over the diel cycle.
    It takes as an argument an unoptimised model, an optional location to save the figure, and an optional name for the figure. The defaults for
    the latter two are as described in plotlinkers.
    '''
    gclinkers = []
    for reaction in model.reactions:
        if "gc_Linker" in reaction.id:
            gclinkers.append(reaction)
    model, solution = FBA_FVA_run(model, obj = "Phloem_tx_overall", rxnlist = gclinkers)
    reactionsdf = pd.DataFrame(columns = ["pfba", "minimum", "maximum"])
    for reaction in tqdm(model.reactions):
        try:
            reactionsdf.loc[reaction.id] = [reaction.flux, model.fva["minimum"][reaction.id], model.fva["maximum"][reaction.id]]
        except:
            reactionsdf.loc[reaction.id] = [reaction.flux, reaction.flux, reaction.flux]
    fig = plotlinkers(model, reactionsdf, coloursdict = "", savelocation=savelocation, name=name)
    return model, fig

def comparelinkers(modellist, savelocation, coloursdict = "", name=None):

    '''
    This function will return a matplotlib figure containing plots of osmolytes over the diel cycle for a series of models.
    It takes as an argument a list of unoptimised models, a location to save the figure, an optional dictionary of colours
    for the lines and an optional name for the figure. The defaults for the latter two are as described in plotlinkers.
    '''

    reactionsdflist = []
    for model in modellist:
        gclinkers = []
        for reaction in model.reactions:
            if "gc_Linker" in reaction.id:
                gclinkers.append(reaction)
        model, solution = FBA_FVA_run(model, obj = "Phloem_tx_overall", rxnlist = gclinkers)
        reactionsdf = pd.DataFrame(columns = ["pfba", "minimum", "maximum"])
        for reaction in tqdm(model.reactions):
            try:
                reactionsdf.loc[reaction.id] = [reaction.flux, model.fva["minimum"][reaction.id], model.fva["maximum"][reaction.id]]
            except:
                reactionsdf.loc[reaction.id] = [reaction.flux, reaction.flux, reaction.flux]
        reactionsdflist.append(reactionsdf)

    import itertools

    def flip(items, ncol):
        return itertools.chain(*[items[i::ncol] for i in range(ncol)])

    reactionsdictlist = []
    for reactionsdf in reactionsdflist:
        reactionsdict = {}
        for reaction, row in reactionsdf.iterrows():
            if reaction.__contains__("Linker"):
                if "_a_" not in reaction:
                    try:
                        reactionsdict[reaction[-13:-9]][reaction[:-14]][int(reaction[-1:])] = (row["pfba"], row["minimum"], row["maximum"])
                    except:
                        try:
                            reactionsdict[reaction[-13:-9]][reaction[:-14]] = {}
                            reactionsdict[reaction[-13:-9]][reaction[:-14]][int(reaction[-1:])] = (row["pfba"], row["minimum"], row["maximum"])
                        except:
                            reactionsdict[reaction[-13:-9]] = {}
                            reactionsdict[reaction[-13:-9]][reaction[:-14]] = {}
                            reactionsdict[reaction[-13:-9]][reaction[:-14]][int(reaction[-1:])] = (row["pfba"], row["minimum"], row["maximum"])
        reactionsdictlist.append(reactionsdict)

    if coloursdict == "":
        coloursdict = {
            'Cl': '#5C6BC0',
            'STARCH':'#8D6E63' ,
            'K': '#AB47BC',
            'MAL': '#2E7D32',
            'SUCROSE': '#EF5350',
            'CIT': '#FF7043',
            'FRU': '#26C6DA',
            'GLC': '#9CCC65',
            'Other': '#BDBDBD',
         }
    coloursdict = collections.OrderedDict(sorted(coloursdict.items()))

    matplotlib.rcParams.update({'font.size': 30})
    matplotlib.rcParams['lines.linewidth'] = 5
    matplotlib.rcParams['axes.linewidth'] = 3
    plt.rcParams["font.family"] = "Arial"

    titlelist = ["Guard Cell Plastid", "Guard Cell Vacuole", "Guard Cell Cytoplasm", "Mesophyll Cell Plastid",  "Mesophyll Cell Vacuole"]
    unitslist = ["(mmoles)", "(M)", "(M)", "(mmoles)", "(mM)"]
    complist = ["p_gc", "v_gc", "c_gc", "p_me", "v_me",]

    apertures = [6,17,6,6]
    gcvolume = []
    for i in range(checknumberofmodels(modellist[0])):
        gcvolume.append((apertures[i]*200+4000)*17.2*(10.0**(-8.0)))

    mevolume = 0.284
    v_mevolume = mevolume*0.9

    for index, reactionsdict in enumerate(reactionsdictlist):
        for compartment in ["v_me", "v_gc"]:
            amal = reactionsdict[compartment].pop("aMAL")
            for i in range(1, checknumberofmodels(modellist[index])+1):
                reactionsdict[compartment]["MAL"][i] = reactionsdict[compartment]["MAL"][i] + amal[i]
            acit = reactionsdict[compartment].pop("aCIT")
            for i in range(1, checknumberofmodels(modellist[index])+1):
                reactionsdict[compartment]["CIT"][i] = reactionsdict[compartment]["CIT"][i] + acit[i]
        linkerconcs = reactionsdict.copy()
        for compartment in linkerconcs:
            for reaction in linkerconcs[compartment]:
                for phase in linkerconcs[compartment][reaction]:
                    if "v_me" in compartment:
                        linkerconcs[compartment][reaction][phase] = ((linkerconcs[compartment][reaction][phase][0]*0.001)/v_mevolume*1000, (linkerconcs[compartment][reaction][phase][1]*0.001)/v_mevolume*1000, (linkerconcs[compartment][reaction][phase][2]*0.001)/v_mevolume*1000)
                    elif "v_gc" in compartment:
                        linkerconcs[compartment][reaction][phase] = ((linkerconcs[compartment][reaction][phase][0]*0.001)/(gcvolume[phase-1]*0.9), (linkerconcs[compartment][reaction][phase][1]*0.001)/(gcvolume[phase-1]*0.9), (linkerconcs[compartment][reaction][phase][2]*0.001)/(gcvolume[phase-1]*0.9))
                    elif "c_gc" in compartment:
                        linkerconcs[compartment][reaction][phase] = ((linkerconcs[compartment][reaction][phase][0]*0.001)/(gcvolume[phase-1]*0.1), (linkerconcs[compartment][reaction][phase][1]*0.001)/(gcvolume[phase-1]*0.1), (linkerconcs[compartment][reaction][phase][2]*0.001)/(gcvolume[phase-1]*0.1))
        reactionsdict = linkerconcs


    legendset = set()

    fig, axs = plt.subplots(6,len(modellist),figsize=(20, 25), sharex = True, sharey = 'row')
    i=0
    j=0
    c=0

    phaseconversion = [0,6,6.5,18,24]

    nonzerolist = []
    for reactionsdict in reactionsdictlist:
        for compartment in complist:
            nonzero = False
            m = 0
            for reaction in reactionsdict[compartment]:
                try:
                    colour = coloursdict[reaction]
                except:
                    colour = coloursdict['Other']
                x = phaseconversion
                y0 = []
                y1 = []
                y2 = []
                for phase in reactionsdict[compartment][reaction].values():
                    y0.append(phase[0])
                    y1.append(phase[1])
                    y2.append(phase[2])

                y0.insert(0, reactionsdict[compartment][reaction][4][0])
                y1.insert(0, reactionsdict[compartment][reaction][4][1])
                y2.insert(0, reactionsdict[compartment][reaction][4][2])
                if max(y2) > 0.001 and ("c_gc" in compartment or "v_gc" in compartment) or max(y2) > 0.0001 and "c_gc" not in compartment and "v_gc" not in compartment:
                    if "sum" not in reaction:
                        axs[i][j].plot(x,y0, color = colour,)
                        axs[i][j].fill_between(x, y1, y2, color = colour, alpha = 0.2)
                        axs[i][j].tick_params(axis='y', which='major',length = 15, width = 3)
                        axs[i][j].tick_params(axis='x', which='major',length = 10, width = 3)
                        legendset.add(reaction)
                        nonzero = True
                else:
                    axs[i][j].plot(x,[0,0,0,0,0], color = colour,)
                    axs[i][j].tick_params(axis='y', which='major',length = 15, width = 3)
                    axs[i][j].tick_params(axis='x', which='major',length = 10, width = 3)
                m += 1
            if j == 0:
                axs[i][j].set_ylabel(unitslist[c], labelpad = 20)
                axs[i][j].set_title(titlelist[c], pad = 15, x = 1.1, ha = "center",)
                axs[i][j].tick_params(axis='y', which='major')
                nonzerolist.append(nonzero)
            elif nonzerolist[i] == False:
                nonzerolist[i] = nonzero

            c += 1
            i += 1
            if i == 5:
                osmolarity = [0.234091748752047, 1.9583538166251, 0.234091748752047, 0.23409174875204702]
                aperture = [6,6,17,6,6]
                axs[i][j].plot(phaseconversion, aperture, color = 'black')
                if j == 0:
                    axs[i][j].set_ylabel(u"Aperture Diameter \n (\u03bcm)", labelpad = 30, size = 30)
                    axs[i][j].set_title("Stomatal Aperture", pad = 15, x = 1.1, ha = "center",)
                axs[i][j].set_yticks([6,17])
                axs[i][j].set_xticks([0,6,18,24])
                axs[i][j].set_xticklabels(["Midnight", "Sunrise", "Sunset", "Midnight",], rotation = 90)
                axs[i][j].set_ylim(0,20)
                axs[i][j].set_xlabel("Time of Day")

                j = 1
                c = 0
                i = 0

    for index, lim in enumerate(nonzerolist):
        if lim == False:
            axs[index][j].set_ylim(-0.05, 1)


    fig.subplots_adjust(hspace=0.4, wspace = 0.1, left=0.2, )

    namesdict = {
        "CIT" : "Citrate",
        "Cl" : r"Cl$^{-}$",
        "STARCH" : "Starch",
        "K" : r"K$^{+}$",
        "MAL" : "Malate",
        "GLC" : "Glucose",
        "FRU" : "Fructose",
        "SUCROSE" : "Sucrose"
    }

    from matplotlib.lines import Line2D
    legendlist = sorted(legendset)
    legendlines = []
    legendnames = []
    for reaction in legendlist:
        includeother = False
        try:
            legendlines.append(Line2D([0], [0], color=coloursdict[reaction]))
            legendnames.append(namesdict[reaction])
        except:
            includeother = True
    if includeother == True:
        legendnames.append("Amino Acids")
        legendlines.append(Line2D([0], [0], color='#BDBDBD'))

    axs[0][0].text(0.5, 1.4, "+Starch", horizontalalignment='center', verticalalignment='center', transform=axs[0][0].transAxes, size = 40, weight = "bold")
    axs[0][1].text(0.5, 1.4, "-Starch", horizontalalignment='center', verticalalignment='center', transform=axs[0][1].transAxes, size = 40, weight = "bold")
    fig.text(0.1, 0.8, "Concentration or amount of ion/metabolite in titular compartment", rotation = 90)

    plt.legend(flip(legendlines, 4), flip(legendnames,4), loc = 'upper center', bbox_to_anchor=(-0.1, -1), ncol = 4, prop={'size': 25})



    if name != None and savelocation != None:
        fig.savefig(os.path.join(savelocation, name + ".svg"), format='svg',  bbox_inches='tight')
    return fig

def getphloemoutput(model):
    '''
    This function returns the overall phloem output of a model over the diel cycle.
    It takes one argument, which is an optimised model.
    '''
    phloemtotal = 0
    for i, j in zip([1,2,3,4,], [6,0.5,11.5,6]):
        phloemtotal = phloemtotal + model.reactions.get_by_id("Phloem_output_tx_me_" + str(i)).flux*j
    return phloemtotal

def modelsolutionstocsv(modellist, modelnames, savelocation, csvname = "modelsolutions.csv"):
    '''
    This function takes a list of optimised models and saves these solutions to a csv. It returns the dataframe that is written to csv.
    It takes 4 arguments:
    modellist: a list of optimised models which might have had fva run on them, not necessarily for all reactions of the models
    modelnames: a list of names matching the list of models
    savelocation: where to save the csv
    csvname: an optional name for the csv, default is modelsolutions
    '''

    resultsdf = pd.DataFrame(columns = ["Reaction Name", "Equation", "Compartment"])
    resultsdf.index.name = "Reaction ID"
    newdf = True
    for model, modelname in tqdm(zip(modellist, modelnames)):
        rxndict = {}
        maxdict = {}
        mindict = {}
        for reaction in model.reactions:
            if newdf == True:
                resultsdf.loc[reaction.id] = [reaction.name, reaction.reaction, getreactioncompartment(reaction)]
            rxndict[reaction.id] = reaction.flux
            try:
                maxdict[reaction.id] = model.fva["maximum"][reaction.id]
                mindict[reaction.id] = model.fva["minimum"][reaction.id]
            except:
                maxdict[reaction.id] = "NA"
                mindict[reaction.id] = "NA"
        resultsdf[modelname + "_pfba"] = resultsdf.index.to_series().map(rxndict)
        resultsdf[modelname + "_max"] = resultsdf.index.to_series().map(maxdict)
        resultsdf[modelname + "_min"] = resultsdf.index.to_series().map(mindict)
        newdf = False
    resultsdf.to_csv(os.path.join(savelocation, csvname))
    return resultsdf

def multifactorialconstraintscan(model, scanarray, savelocation, csvname = "mastersolutions.csv", aperture = True):
    '''
    This function runs a multifactorial constraint scan of the model based on a scanarray generated by pyDOE
    It adds the results to a csv, and returns a list of the results.
    It takes five arguments:
    model: the unoptimised model to be used for the scan
    scanarray: the array with the constraints to be used in the scan
    savelocation: the location to save the csv to
    csvname: the name of the csv, if the csv already exists addsolutiontomastercsv will keep adding the solutions to those that are already there
    aperture: boolean that allows the user to define if the aperture is constrained as part of the scan or not
    '''

    doescanresults = []
    osmolytestoconstrain = ["Cl", "K", "SUCROSE", "FRU", "GLC", "MAL"]
    for constraints in tqdm(scanarray):
        tempmodel = model.copy()
        osmolarity = ((((2.5*math.exp(0.16*(15+constraints[6]*5)))*((((15+constraints[6]*5)*200+4000)*17.2*10**-8)-((33*(10**-6))+(565*(10**3)*(0)))))/(0.082*293))*10**3)
        totalosmo = 0
        for i, osmolyte in enumerate(osmolytestoconstrain):
            closedosmo = (((constraints[i])*osmolarity)+20)/26
            if osmolyte == "MAL":
                setflux(tempmodel, osmolyte + "_c_gc_Linker", 0, (constraints[i])*closedosmo*0.1, multi = True)
                setflux(tempmodel, osmolyte + "_v_gc_Linker", 0, (constraints[i])*closedosmo*0.9*0.7, multi = True)
                setflux(tempmodel, "a" + osmolyte + "_v_gc_Linker", 0, (constraints[i])*closedosmo*0.9*0.3, multi = True)
                setflux(tempmodel, osmolyte + "_c_gc_Linker_2", 0, (constraints[i])*osmolarity*0.1, multi = False)
                setflux(tempmodel, osmolyte + "_v_gc_Linker_2", 0, (constraints[i])*osmolarity*0.9*0.7, multi = False)
                setflux(tempmodel, "a" + osmolyte + "_v_gc_Linker_2", 0, (constraints[i])*osmolarity*0.9*0.3, multi = False)
            else:
                setflux(tempmodel, osmolyte + "_c_gc_Linker", 0, (constraints[i])*closedosmo*0.1, multi = True)
                setflux(tempmodel, osmolyte + "_v_gc_Linker", 0, (constraints[i])*closedosmo*0.9, multi = True)
                setflux(tempmodel, osmolyte + "_c_gc_Linker_2", 0, (constraints[i])*osmolarity*0.1, multi = False)
                setflux(tempmodel, osmolyte + "_v_gc_Linker_2", 0, (constraints[i])*osmolarity*0.9, multi = False)
            print "Max " + osmolyte + " constrained to " + str(round(constraints[i]*(17*200+4000)*17.2*(10.0**(-8.0))*10**3, 3)) + "M (" + str((constraints[i])*osmolarity)+ ") " + str(round((constraints[i])*osmolarity/osmolarity*100)) + "% of total osmolarity required"
            totalosmo = totalosmo + constraints[i]*osmolarity
        print "Constraints can account for %s percent of total osmolarity" % (totalosmo/osmolarity*100,)
        if aperture == True:
            apertureconstraint = [6,15+constraints[6]*5,6,6]
            tempmodel = constrainaperture(tempmodel, apertureconstraint) #constrain open aperture
            print "Aperture constrained to:" + str(apertureconstraint)
        else:
            pass
        setflux(tempmodel, "Photon_tx_gc_2", 0, constraints[7]*8.64, multi = False) #constrain photon import opening
        setflux(tempmodel, "Photon_tx_gc_3", 0, constraints[7]*8.64, multi = False) #constrain photon import closing
        print "Photon influx constrained to:" + str(constraints[7]*8.64)
        try:
            print "Running pFBA"
            tempmodel, tempsolution = pfba_Weighted(tempmodel, objective="Phloem_tx_overall")
        except:
            print "Constraint " + str(constraints) + "gives infeasible solution."
            continue
        print "pFBA completed, creating results dictionary"
        resultsdict = pd.DataFrame(columns = ["Reaction Name", "Equation", "Compartment", "pfba", "upb", "lob"])
        resultsdict.index.name = "Reaction ID"
        for reaction in tempmodel.reactions:
            resultsdict.loc[reaction.id] = [reaction.name, reaction.reaction, getreactioncompartment(reaction), reaction.flux, reaction.upper_bound, reaction.lower_bound]
        print "Constraints completed: " +str(constraints)
        resultslist = [tempmodel, resultsdict, constraints]
        doescanresults.append(resultslist)
        try:
            addsolutiontomastercsv(resultsdict, filename = os.path.join(savelocation, csvname))
        except:
            raw_input('Error writing to file, check file is not open and press enter')
            addsolutiontomastercsv(resultsdict, filename = os.path.join(savelocation, csvname))
    return doescanresults

def addsolutiontomastercsv(solutiondf, filename = "dfmodelsolutions.csv"):
    '''
    This function is used by multifactorialconstraintscan to add model solutions to a solution csv.
    If the csv does not yet exist it creates it and adds the index columns (reaction id, name, equation, compartment), otherwise
    the model solutions are added as new columns, with the fba solution, lower bound, then upper bound of the reaction.
    It takes as an argument a dataframe with the aforementioned index columns and fba, lower bound, upper bound, and a filename
    '''
    import pandas as pd
    try:
        newdf = solutiondf.rename(columns={'pfba': 'fba'})
        masterdf = pd.read_csv(filename, index_col = 0)
        modelnumber = int(masterdf.columns[-1][3:])+1
        for i in ['fba', "upb", "lob"]:
            masterdf[i + str(modelnumber)] = pd.Series(newdf[i].to_dict())
    except IOError:
        print "Starting new file"
        newdf = solutiondf.rename(columns={'pfba': 'fba1', "lob" : "lob1", "upb": "upb1"})
        masterdf = newdf
    masterdf.to_csv(filename)
    return masterdf

def parsesolutionsforpcplot(csvname, savelocation, limitlist):
    '''
    This function takes the solutions from the multifactorial constraint scan and parses them to be used for
    the parallel coordinate plot. It performs transformations and changes names ready for plotting, and also
    calculates whether certain reactions have hit their upper or lower bounds.
    It returns scaleddf, the dataframe to be used for plotting, and the columns of the dataframe.
    It takes three arguments:
    csvname: the name of the csv file to be parsed
    savelocation: the location of the csv file to be parsed
    limitlist: a list of the reactions to check to see if they are limiting
    '''
    masterdf = pd.read_csv(os.path.join(savelocation, csvname + ".csv"), index_col = 0)
    newdf = pd.DataFrame()
    for i in range(1, int(masterdf.columns[-1][3:])+1):
        rxndict = {}
        for rxnid in masterdf.index:
            rxndict[rxnid] = masterdf["fba" + str(i)][rxnid]
        for reactionid in limitlist:
            if (abs(masterdf["upb" + str(i)][reactionid] - masterdf["fba" + str(i)][reactionid])<0.00001 and masterdf["upb" + str(i)][reactionid] != 0) or (abs(masterdf["lob" + str(i)][reactionid] - masterdf["fba" + str(i)][reactionid])<0.00001 and masterdf["lob" + str(i)][reactionid] != 0):
                rxndict[reactionid + "_limit"] = 0
            else:
                rxndict[reactionid + "_limit"] = 1
        newdf[i] = pd.Series(rxndict)
    newdf =newdf.transpose()
    newdf["SUCROSE_importfvaconstraint_log"] = newdf["SUCROSE_importfvaconstraint"].apply(lambda x: np.log10(x) if x >= 0.0001 else -4)
    newdf["Cl_importfvaconstraint_log"] = newdf["Cl_importfvaconstraint"].apply(lambda x: np.log10(x) if x >= 0.0001 else -4)
    newdf["FRU_importfvaconstraint_log"] = newdf["FRU_importfvaconstraint"].apply(lambda x: np.log10(x) if x >= 0.0001 else -4)
    newdf["GLC_importfvaconstraint_log"] = newdf["GLC_importfvaconstraint"].apply(lambda x: np.log10(x) if x >= 0.0001 else -4)
    newdf["STARCH_linkerfvaconstraint_abs"] = newdf["STARCH_linkerfvaconstraint"].abs()

    scaleddf = pd.DataFrame()

    newdf["Volume"] = ((200*17+4000)*17.2*10**(-8))-(33*10**(-6))

    scaleddf["Osmo/2 (M)"] = (newdf["pseudoOs_constraint_gc_2"]-newdf["pseudoOs_constraint_gc_1"])/newdf["Volume"]*10**(-3)/2
    scaleddf['Photons/5 ($mmolm^{-2}leafh^{-1}$)'] = newdf["Photon_tx_gc_3"]/5
    scaleddf[r'Cl$^{-}$ Acc (M)'] = newdf["Cl_linkerfvaconstraint"]/newdf["Volume"]*10**(-3)
    scaleddf[r'K$^{+}$ Acc (M)'] = newdf["K_linkerfvaconstraint"]/newdf["Volume"]*10**(-3)
    scaleddf['Mal Acc (M)']  = newdf["MAL_linkerfvaconstraint"]/newdf["Volume"]*10**(-3)
    scaleddf['Suc Acc (M)'] = newdf["SUCROSE_linkerfvaconstraint"]/newdf["Volume"]*10**(-3)
    scaleddf['Fru Acc (M)']  = newdf["FRU_linkerfvaconstraint"]/newdf["Volume"]*10**(-3)
    scaleddf['Glc Acc (M)'] = newdf["GLC_linkerfvaconstraint"]/newdf["Volume"]*10**(-3)
    scaleddf['Starch Deg $(mmolm^{-2})$'] = newdf["STARCH_linkerfvaconstraint_abs"]

    columns = []
    for xtick in scaleddf:
        columns.append(xtick)

    for reaction in limitlist:
        scaleddf[reaction + "_limit"] = newdf[reaction + "_limit"]
    return scaleddf, columns

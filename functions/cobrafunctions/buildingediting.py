import cobra

def splitmodel (model, labels):
    '''
    This function takes two arguments: the model to be split and the labels for the models to split it into.
    len(labels) number of models will be created, with the suffix "_label" to define them.
    '''

    numberofmodels = len(labels)
    oldmodel = model.copy()

    # Create new metabolites for models
    oldmetabolites = []
    for i in range(len(model.metabolites)):
        oldmetabolite = model.metabolites[i]
        for j in labels:
            newmetabolite = oldmetabolite.copy()
            newmetabolite.name = newmetabolite.name + "_" + str(j)
            newmetabolite.id = newmetabolite.id + "_" + str(j)
            model.add_metabolites(newmetabolite)
        oldmetabolites.append(oldmetabolite)

    # Create new reactions for models
    oldreactions = []
    for i in range(len(model.reactions)):
        oldreaction = model.reactions[i]
        for j in labels:
            newreaction = oldreaction.copy()
            newreaction.name = newreaction.name + "_" + str(j)
            newreaction.id = newreaction.id + "_" + str(j)
            newmetabolites = {}
            for metabolite, value in newreaction.metabolites.items():
                    newmetabolites[str(metabolite) + "_" + str(j)] = value
            model.add_reaction(newreaction)
            newreaction = model.reactions.get_by_id(newreaction.id)
            newreaction.subtract_metabolites(newreaction.metabolites)
            newreaction.add_metabolites(newmetabolites)
        oldreactions.append(oldreaction)

    # Remove old reactions and metabolites
    model.remove_reactions(oldreactions)
    model.remove_metabolites(oldmetabolites)

#tests

    # Check no. of metabolites and reactions is equal
    if (len(model.reactions) == numberofmodels*len(oldmodel.reactions) and len(model.metabolites) == numberofmodels*len(oldmodel.metabolites)) == False:
        raise UserWarning("Number of metabolites and reactants does not match")

    emptyreactions = []
    for reaction in model.reactions:
        if len(reaction.metabolites) == 0:
            emptyreactions.append(reaction.name)
    if len(emptyreactions) != 0:
        raise UserWarning("There are some reactions with no metabolites, they are: " + str(emptyreactions))

def addlinkers(model, linkersfile, compartments, cells, phasetimes):
    '''
    This function adds linker reactions for osmolytes to a phased model of guard cell, and adds osmolarity and charge linker pseudoreactions.
    It takes five arguments:
    model: the model to which linkers are to be added
    linkersfile: a .csv file with a list of osmolytes and their osmotic coefficients
    compartments: a dictionary with the tags of compartments as keys and their relative volumes as values
    cells: a list of cell tags to be used
    phasetimes: a list of the lengths of each phase, for scaling of the linker reactions
    '''

    numberofmodels = checknumberofmodels(model)
    linkers = deflinkers(linkersfile)
    cellcompartments = []
    for compartment in compartments:
        if compartment == "a":
            cellcompartments.append(compartment)
        else:
            for cell in cells:
                cellcompartments.append(compartment + "_" + cell)

    #create the linker reactions from the 'osmolytes' list
    newreactions = []
    for linker in linkers:
        addreaction(model, linker + "_Linker")


    #add pseudometabolites to represent accumulation reactions
    for compartment in cellcompartments:
        addmetabolite(model, "pseudoOs_"+ compartment, compartment)
        addmetabolite(model, "pseudoCharge_" + compartment, compartment)

    #add metabolites to linker reactions
    for compartment in cellcompartments:
        for linker in linkers:
            if "_"+ compartment in linker:
                for i in range(1, numberofmodels):
                    model.reactions.get_by_id(linker + "_Linker_" + str(i)).add_metabolites({
                        linker + "_" + str(i) : (1.0/phasetimes[i-1])*-1,
                        linker + "_" + str(i+1) : (1.0/phasetimes[i])*1,
                        "pseudoOs_" + compartment + "_" + str(i): int(linkers[linker]),
                        "pseudoCharge_" +  compartment + "_" + str(i): float(model.metabolites.get_by_id(linker + "_" + str(i)).charge) * (1.0/phasetimes[i-1])*-1

                    })
                model.reactions.get_by_id(linker + "_Linker_" + str(numberofmodels)).add_metabolites({
                        linker + "_" + str(numberofmodels) : (1.0/phasetimes[numberofmodels-1])*-1,
                        linker + "_" + str(1) : (1.0/phasetimes[0])*1,
                        "pseudoOs_" + compartment + "_" + str(numberofmodels): int(linkers[linker]),
                        "pseudoCharge_" +  compartment + "_" + str(numberofmodels): float(model.metabolites.get_by_id(linker + "_" + str(i)).charge) * (1.0/phasetimes[numberofmodels-1])* -1
                    })

    #add osmotic constraint reactions

    for cell in cells:
        for i in range(1, numberofmodels+1):
            addreaction(model, "pseudoOs_constraint" + "_" + cell + "_" + str(i), multi = "")
            for compartment in compartments:
                if compartment == "a":
                    addreaction(model, "pseudoOs_constraint_a_" + str(i), multi = "")
                    model.reactions.get_by_id("pseudoOs_constraint_a_" + str(i)).add_metabolites({
                        "pseudoOs_" + compartment + "_" + str(i): -1 * compartments[compartment],
                    })
                else:
                    model.reactions.get_by_id("pseudoOs_constraint" + "_" + cell + "_" + str(i)).add_metabolites({
                        "pseudoOs_" + compartment + "_" + cell + "_" + str(i): -1 * compartments[compartment],
                    })
#tests

    #run a basic test to ensure new reactions have metabolites attached to them
    emptyreactions = []
    for reaction in model.reactions:
        if len(reaction.metabolites) == 0:
            emptyreactions.append(reaction.name)
    if len(emptyreactions) != 0:
        raise UserWarning ("There are some reactions with no metabolites, they are: " + str(emptyreactions))

def deflinkers(file):
    '''
    This function creates a dict of osmolytes and their osmotic coefficients from the .csv file, and is used by addlinkers
    '''
    import csv
    with open(file, mode='r') as csv_file:
        csvlinkers = csv.DictReader(csv_file)
        linkers = {}
        line_count = 0
        for linker in csvlinkers:
            if line_count == 0:
                line_count += 1
            linkers[linker["Linker"]] = linker["Osmotic Coefficient"]
            line_count += 1
    return linkers

def constrainaperture(model, apertures):
    '''
    This function will constrain the osmolarity of the model relative to aperture constraints.
    It takes a model to constrain, and a list of aperture sizes in um to use, and returns the constrained model.
    '''
    import math
    osmoconstraints = []
    co2constraints = []
    for i in range(checknumberofmodels(model)):
        osmoconstraints.append(model.problem.Constraint(model.reactions.get_by_id("pseudoOs_constraint_gc_"+str(i+1)).flux_expression - ((((2.5*math.exp(0.16*apertures[i]))*(((apertures[i]*200+4000)*17.2*10**-8)-((33*(10**-6))+(565*(10**3)*(0)))))/(0.082*293))*10**3),
            lb=0,
            ub=0))
        model.add_cons_vars(osmoconstraints[i])
    return model

def constrainmaintenance(model):
    '''
    This function constrains the maintenance reactions in the model relative to the input of photons into the model.
    '''
    for i in range(1, checknumberofmodels(model) + 1):
        memaintenance = model.problem.Constraint((model.reactions.get_by_id("ATPase_tx_me_" + str(i)).flux_expression-(model.reactions.Photon_tx_me_3.flux_expression*0.0049+2.7851)), lb=0, ub=1000)
        model.add_cons_vars(memaintenance)
        gcmaintenance = model.problem.Constraint((model.reactions.get_by_id("ATPase_tx_gc_" + str(i)).flux_expression-(model.reactions.Photon_tx_gc_3.flux_expression*0.0049+(2.7851/500))), lb=0, ub=1000)
        model.add_cons_vars(gcmaintenance)
    return model

def addmetabolite(model, name, compartment = "default", multi=True):
    '''
    This function adds a metabolite to the given model, and returns this metabolite.
    If multi is set to true it will be added to every phase of the model.
    '''
    if multi == True:
        for i in range(1,checknumberofmodels(model)+1):
            newname = name + "_" + str(i)
            metabolite = cobra.core.Metabolite(newname)
            metabolite.id = newname
            metabolite.name = newname
            metabolite.compartment = compartment
            model.add_metabolites(metabolite)
    else:
        metabolite = cobra.core.Metabolite(name)
        metabolite.id = name
        metabolite.name = name
        metabolite.compartment = compartment
        model.add_metabolites(metabolite)

    return metabolite

def addreaction(model, name, multi=True):
    '''
    This function adds a reaction to the given model, and returns this reaction.
    If multi is set to true it will be added to every phase of the model.
    '''
    if multi == True:
        for i in range(1,checknumberofmodels(model)+1):
            newname = name + "_" + str(i)
            reaction = cobra.core.Reaction(newname)
            reaction.id = newname
            reaction.name = newname
            model.add_reaction(reaction)
    else:
        reaction = cobra.core.Reaction(name)
        reaction.id = name
        reaction.name = name
        model.add_reaction(reaction)
    return reaction

def addmetabolitestoreactionmulti(model, reaction, metabolitesdict):
    '''
    This function adds a dictionary of metabolites to a reaction in every phase of the model
    '''
    for i in range(1,checknumberofmodels(model)+1):
        for metabolite in metabolitesdict:
            model.reactions.get_by_id(reaction + "_" + str(i)).add_metabolites({
                metabolite + "_" + str(i): metabolitesdict[metabolite]
            })

def checknumberofmodels(model):
    '''
    This function returns the number of phases of the model based on the maximal number of the last tag
    '''
    submodelnumbers = set()
    for metabolite in model.metabolites:
        try:
            submodelnumbers.add(int(metabolite.name[-1:]))
        except:
            pass
    try:
        modelnumber = max(submodelnumbers)
    except:
        modelnumber = 1
    return modelnumber

def setflux(model, reactionid, lowerflux, upperflux, multi = True):

    '''
    Use setflux in order to constrain a reaction's flux.
    If multi == True, the function will constrain that reaction in every phase.
    '''

    if multi == True:
        for i in range(1, checknumberofmodels(model)+1):
            model.reactions.get_by_id(reactionid + "_" + str(i)).lower_bound = lowerflux
            model.reactions.get_by_id(reactionid + "_" + str(i)).upper_bound = upperflux
    else:
        model.reactions.get_by_id(reactionid).lower_bound = lowerflux
        model.reactions.get_by_id(reactionid).upper_bound = upperflux

    return model

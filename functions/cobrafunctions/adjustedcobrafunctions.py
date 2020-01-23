from __future__ import absolute_import
from cobrafunctions.buildingediting import *
from cobrafunctions.analysing import *

import logging
from warnings import warn
from itertools import chain
from optlang.symbolics import add, Zero
from numpy import zeros
from cobra.util import solver as sutil
from cobra.core.solution import get_solution

def pfba_Weighted(model, weightings = None, fraction_of_optimum=1.0, objective=None, reactions=None):
    """Perform basic pFBA (parsimonious Enzyme Usage Flux Balance Analysis)
    to minimize total flux.
    pFBA [1] adds the minimization of all fluxes the the objective of the
    model. This approach is motivated by the idea that high fluxes have a
    higher enzyme turn-over and that since producing enzymes is costly,
    the cell will try to minimize overall flux while still maximizing the
    original objective function, e.g. the growth rate.
    Parameters
    ----------
    model : cobra.Model
        The model
    fraction_of_optimum : float, optional
        Fraction of optimum which must be maintained. The original objective
        reaction is constrained to be greater than maximal_value *
        fraction_of_optimum.
    objective : dict or model.problem.Objective
        A desired objective to use during optimization in addition to the
        pFBA objective. Dictionaries (reaction as key, coefficient as value)
        can be used for linear objectives.
    reactions : iterable
        List of reactions or reaction identifiers. Implies `return_frame` to
        be true. Only return fluxes for the given reactions. Faster than
        fetching all fluxes if only a few are needed.
    Returns
    -------
    cobra.Solution
        The solution object to the optimized model with pFBA constraints added.
    References
    ----------
    .. [1] Lewis, N. E., Hixson, K. K., Conrad, T. M., Lerman, J. A.,
       Charusanti, P., Polpitiya, A. D., Palsson, B. O. (2010). Omic data
       from evolved E. coli are consistent with computed optimal growth from
       genome-scale models. Molecular Systems Biology, 6,
       390. doi:10.1038/msb.2010.47
    """
    reactions = model.reactions if reactions is None \
        else model.reactions.get_by_any(reactions)
    tempmodel = model.copy()
    with tempmodel as m:
        add_pfba_Weighted(m, weightings, objective=objective,
                 fraction_of_optimum=fraction_of_optimum)
        m.slim_optimize(error_value=None)
        solution = get_solution(m, reactions=reactions)
    return m, solution


def add_pfba_Weighted(model, weightings = None, objective=None, fraction_of_optimum=1.0):
    """
    #################################################################################
    # This function is a modified version of cobrapy add_pfba function			#
    #										#
    #################################################################################
    Add pFBA objective
    Add objective to minimize the summed flux of all reactions to the
    current objective.
    See Also
    -------
    pfba
    Parameters
    ----------
    model : cobra.Model
        The model to add the objective to
    objective :
        An objective to set in combination with the pFBA objective.
    fraction_of_optimum : float
        Fraction of optimum which must be maintained. The original objective
        reaction is constrained to be greater than maximal_value *
        fraction_of_optimum.
    """
    if weightings == None:
        weightings = getweightings(model)
    if objective is not None:
        model.objective = objective
    if model.solver.objective.name == '_pfba_objective':
        raise ValueError('The model already has a pFBA objective.')
    sutil.fix_objective_as_constraint(model, fraction=fraction_of_optimum)
    reaction_variables = ((rxn.forward_variable, rxn.reverse_variable)
                          for rxn in model.reactions)
    variables = chain(*reaction_variables)
    model.objective = model.problem.Objective(
        Zero, direction='min', sloppy=True, name="_pfba_objective")
    #print([v for v in variables])
    tempDict = dict()
    for v in variables:
        w = str(v).split("=")[1].replace(" ","").replace("<","")
        found=False
        for rxn in weightings.keys():
            if w.__contains__(rxn):
                tempDict[v]=weightings[rxn]
                found=True
                break
        if not found:
            print("Weightings for reaction "+w+" not found, so assuming weighting = 1")
            tempDict[v] = 1
    model.objective.set_linear_coefficients(tempDict)


def getweightings(model):
    '''
    This function is used by pfba_weighted to generate default weightings for the guard cell model
    It takes the model as an argument and returns the weightings based on the phase lengths of the model.
    '''
    weightings = {}
    numberofmodels = checknumberofmodels(model)
    for i in range(1,numberofmodels+1):
        lengthofphase = 1/(-model.reactions.get_by_id("SUCROSE_v_gc_Linker_" +  str(i)).get_coefficient("SUCROSE_v_gc_" + str(i)))
        for reaction in model.reactions:
            if "constraint" in reaction.id or "overall" in reaction.id or "sum" in reaction.id:
                weightings[reaction.id] = 0
            elif "_" + str(i) in reaction.id:
                weightings[reaction.id] = lengthofphase
    return weightings

def pfba_rxns2avoid(model,rxns2avoid, fraction_of_optimum=1.0, objective=None, reactions=None):

    temp = model.copy()
    reactions = temp.reactions if reactions is None \
        else temp.reactions.get_by_any(reactions)
    with temp as m:
        add_pfba_rxns2avoid(m, rxns2avoid, objective=objective,
                 fraction_of_optimum=fraction_of_optimum)
        m.slim_optimize(error_value=None)
        solution = cobra.core.solution.get_solution(m, reactions=reactions)
    return temp, solution

def add_pfba_rxns2avoid(model, rxns2avoid, objective=None, fraction_of_optimum=1.0):

    if objective is not None:
        model.objective = objective
    if model.solver.objective.name == '_pfba_objective':
        raise ValueError('The model already has a pFBA objective.')

    cobra.util.solver.fix_objective_as_constraint(model, fraction=fraction_of_optimum)
    rxnlist = [r for r in model.reactions if r.id  not in rxns2avoid]
    reaction_variables = ((rxn.forward_variable, rxn.reverse_variable)
                         for rxn in rxnlist)

    variables = itertools.chain(*reaction_variables)
    model.objective = model.problem.Objective(
        Zero, direction='min', sloppy=True, name="_pfba_objective")
    model.objective.set_linear_coefficients({v: 1.0 for v in variables})

def FBA_FVA_run(cobra_model,obj,rxnlist=[]):
  from cobra import flux_analysis


  if len(rxnlist)==0:
    rxnlist = cobra_model.reactions
  #print("Rxn list ="+str(rxnlist))
  print("Running pFBA")
  cobra_model, solution = pfba_Weighted(cobra_model, objective=obj)
  objvalue = solution.get_primal_by_id(obj)

  sumOfFluxes = getsumoffluxes(cobra_model)

  weightings = getweightings(cobra_model)
  cobra_model2 = cobra_model.copy()
  irr_model = rev2irrev(cobra_model2, weightings)
  print("Setting SOF model")

  for reaction in irr_model.reactions:
      if reaction.id.__contains__("_reverse"):
          id = reaction.id
          originalreaction = id.replace("_reverse","")
          weightings[reaction.id] = weightings[originalreaction]

  coefficients = {}
  for reaction in irr_model.reactions:
      coefficients[reaction.forward_variable] = weightings[reaction.id]
      coefficients[reaction.reverse_variable] = weightings[reaction.id]
  sofconstraint = irr_model.problem.Constraint(0, lb=sumOfFluxes, ub=sumOfFluxes)
  irr_model.add_cons_vars(sofconstraint)
  irr_model.solver.update()
  sofconstraint.set_linear_coefficients(coefficients=coefficients)

  phloemconstraint = irr_model.problem.Constraint(irr_model.reactions.Phloem_tx_overall.flux_expression
      , lb = objvalue, ub = objvalue)
  irr_model.add_cons_vars(phloemconstraint)


  sfmodel = irr_model.copy()

  rxnlist2 = list()
  if rxnlist == cobra_model.reactions:
    rxnlist2 = sfmodel.reactions
  else:
    for rxn in rxnlist:
      if rxn.lower_bound<0 and rxn.upper_bound>0 and weightings[rxn.id] != 0:
        rxnlist2.append(sfmodel.reactions.get_by_id(rxn.id+"_reverse"))
      rxnlist2.append(sfmodel.reactions.get_by_id(rxn.id))
  #print("Rxn list ="+str(rxnlist2))
  print("Running FVA")


  fva = flux_variability_analysis(sfmodel,reaction_list = rxnlist2)
  print("Processing results")

  fva2=dict()
  for mode in fva.keys():
    if mode == "maximum":
      tempdict = dict()
      FVArxnSet = set()
      for rxn in fva[mode].keys():
        if rxn.__contains__("_reverse"):
          rxn = rxn.replace("_reverse","")
        if FVArxnSet.__contains__(rxn):
          continue
        FVArxnSet.add(rxn)
        if not fva[mode].keys().__contains__(rxn+"_reverse"):
          maxi = fva[mode][rxn]
        else:
          maxi = fva[mode][rxn]+fva[mode][rxn+"_reverse"]
        tempdict[rxn]=maxi
    else:
      tempdict=dict()
      FVArxnSet = set()
      for rxn in fva[mode].keys():
        if rxn.__contains__("_reverse"):
          rxn = rxn.replace("_reverse","")
        if FVArxnSet.__contains__(rxn):
          continue
        FVArxnSet.add(rxn)
        if not fva[mode].keys().__contains__(rxn+"_reverse"):
          mini = fva[mode][rxn]
        else:
          mini = fva[mode][rxn]+fva[mode][rxn+"_reverse"]
        tempdict[rxn]=mini
    fva2[mode]=tempdict

  sfmodel.fva = fva
  cobra_model.fva = fva2
  return cobra_model, sfmodel

def getsumoffluxes(model):
     weightings = getweightings(model)
     sumoffluxes = 0
     for reactionid in weightings.keys():
         sumoffluxes = sumoffluxes + (abs(model.reactions.get_by_id(reactionid).flux)*weightings[reactionid])
     return sumoffluxes

def rev2irrev(cobra_model, weightings):
  '''
  #Function to convert any model with reversible reactions to a copy of the same m-
  #-odel with only irreversible reactions. ID of reverse reactions are generated by
  #suffixing "_reverse" to the ID of the orignal reaction.
  #args: 1) a cobra model
  #output: a cobra model with only irreversible reactions
  '''
  exp_model=cobra_model.copy()
  for RXN in cobra_model.reactions:
    if weightings[RXN.id] == 0:
        continue
    else:
        rxn=exp_model.reactions.get_by_id(RXN.id)
        if (rxn.lower_bound < 0):
          rxn_reverse = rxn.copy()
          rxn_reverse.id = "%s_reverse" %(rxn.id)
          rxn.lower_bound = 0
          rxn_reverse.upper_bound = 0
          exp_model.add_reaction(rxn_reverse)

  return exp_model

from cobra.flux_analysis.loopless import loopless_fva_iter
from cobra.flux_analysis.parsimonious import add_pfba
from cobra.flux_analysis.deletion import (
    single_gene_deletion, single_reaction_deletion)

def flux_variability_analysis(model, reaction_list=None, loopless=False, fraction_of_optimum=1.0, pfba_factor=None, showprogress = True):
    """
    This is identical to the default fva function, but with an added boolean parameter for if one would like to show a progress bar.

    Determine the minimum and maximum possible flux value for each reaction.

    Parameters
    ----------
    model : cobra.Model
        The model for which to run the analysis. It will *not* be modified.
    reaction_list : list of cobra.Reaction or str, optional
        The reactions for which to obtain min/max fluxes. If None will use
        all reactions in the model (default).
    loopless : boolean, optional
        Whether to return only loopless solutions. This is significantly
        slower. Please also refer to the notes.
    fraction_of_optimum : float, optional
        Must be <= 1.0. Requires that the objective value is at least the
        fraction times maximum objective value. A value of 0.85 for instance
        means that the objective has to be at least at 85% percent of its
        maximum.
    pfba_factor : float, optional
        Add an additional constraint to the model that requires the total sum
        of absolute fluxes must not be larger than this value times the
        smallest possible sum of absolute fluxes, i.e., by setting the value
        to 1.1 the total sum of absolute fluxes must not be more than
        10% larger than the pFBA solution. Since the pFBA solution is the
        one that optimally minimizes the total flux sum, the ``pfba_factor``
        should, if set, be larger than one. Setting this value may lead to
        more realistic predictions of the effective flux bounds.

    Returns
    -------
    pandas.DataFrame
        A data frame with reaction identifiers as the index and two columns:
        - maximum: indicating the highest possible flux
        - minimum: indicating the lowest possible flux

    Notes
    -----
    This implements the fast version as described in [1]_. Please note that
    the flux distribution containing all minimal/maximal fluxes does not have
    to be a feasible solution for the model. Fluxes are minimized/maximized
    individually and a single minimal flux might require all others to be
    suboptimal.

    Using the loopless option will lead to a significant increase in
    computation time (about a factor of 100 for large models). However, the
    algorithm used here (see [2]_) is still more than 1000x faster than the
    "naive" version using ``add_loopless(model)``. Also note that if you have
    included constraints that force a loop (for instance by setting all fluxes
    in a loop to be non-zero) this loop will be included in the solution.

    References
    ----------
    .. [1] Computationally efficient flux variability analysis.
       Gudmundsson S, Thiele I.
       BMC Bioinformatics. 2010 Sep 29;11:489.
       doi: 10.1186/1471-2105-11-489, PMID: 20920235

    .. [2] CycleFreeFlux: efficient removal of thermodynamically infeasible
       loops from flux distributions.
       Desouki AA, Jarre F, Gelius-Dietrich G, Lercher MJ.
       Bioinformatics. 2015 Jul 1;31(13):2159-65.
       doi: 10.1093/bioinformatics/btv096.
    """
    if reaction_list is None:
        reaction_list = model.reactions
    else:
        reaction_list = model.reactions.get_by_any(reaction_list)

    prob = model.problem
    fva_results = pd.DataFrame({
        "minimum": zeros(len(reaction_list), dtype=float),
        "maximum": zeros(len(reaction_list), dtype=float)
    }, index=[r.id for r in reaction_list])
    with model:
        # Safety check before setting up FVA.
        model.slim_optimize(error_value=None,
                            message="There is no optimal solution for the "
                                    "chosen objective!")
        # Add the previous objective as a variable to the model then set it to
        # zero. This also uses the fraction to create the lower/upper bound for
        # the old objective.
        if model.solver.objective.direction == "max":
            fva_old_objective = prob.Variable(
                "fva_old_objective",
                lb=fraction_of_optimum * model.solver.objective.value)
        else:
            fva_old_objective = prob.Variable(
                "fva_old_objective",
                ub=fraction_of_optimum * model.solver.objective.value)
        fva_old_obj_constraint = prob.Constraint(
            model.solver.objective.expression - fva_old_objective, lb=0, ub=0,
            name="fva_old_objective_constraint")
        model.add_cons_vars([fva_old_objective, fva_old_obj_constraint])

        if pfba_factor is not None:
            if pfba_factor < 1.:
                warn("The 'pfba_factor' should be larger or equal to 1.",
                     UserWarning)
            with model:
                add_pfba(model, fraction_of_optimum=0)
                ub = model.slim_optimize(error_value=None)
                flux_sum = prob.Variable("flux_sum", ub=pfba_factor * ub)
                flux_sum_constraint = prob.Constraint(
                    model.solver.objective.expression - flux_sum, lb=0, ub=0,
                    name="flux_sum_constraint")
            model.add_cons_vars([flux_sum, flux_sum_constraint])

        model.objective = Zero  # This will trigger the reset as well
        if showprogress == True:
            with tqdm(total=len(reaction_list*2)) as pbar:
                pbar.set_description("FVA Progress:")
                for what in "minimum", "maximum":
                    sense = "min" if what == "minimum" else "max"
                    model.solver.objective.direction = sense
                    for rxn in reaction_list:
                        # The previous objective assignment already triggers a reset
                        # so directly update coefs here to not trigger redundant resets
                        # in the history manager which can take longer than the actual
                        # FVA for small models
                        model.solver.objective.set_linear_coefficients(
                            {rxn.forward_variable: 1, rxn.reverse_variable: -1})
                        model.slim_optimize()
                        sutil.check_solver_status(model.solver.status)
                        if loopless:
                            value = loopless_fva_iter(model, rxn)
                        else:
                            value = model.solver.objective.value
                        fva_results.at[rxn.id, what] = value
                        model.solver.objective.set_linear_coefficients(
                            {rxn.forward_variable: 0, rxn.reverse_variable: 0})
                        pbar.update(1)
        else:
            for what in "minimum", "maximum":
                sense = "min" if what == "minimum" else "max"
                model.solver.objective.direction = sense
                for rxn in reaction_list:
                    # The previous objective assignment already triggers a reset
                    # so directly update coefs here to not trigger redundant resets
                    # in the history manager which can take longer than the actual
                    # FVA for small models
                    model.solver.objective.set_linear_coefficients(
                        {rxn.forward_variable: 1, rxn.reverse_variable: -1})
                    model.slim_optimize()
                    sutil.check_solver_status(model.solver.status)
                    if loopless:
                        value = loopless_fva_iter(model, rxn)
                    else:
                        value = model.solver.objective.value
                    fva_results.at[rxn.id, what] = value
                    model.solver.objective.set_linear_coefficients(
                        {rxn.forward_variable: 0, rxn.reverse_variable: 0})


    return fva_results[["minimum", "maximum"]]

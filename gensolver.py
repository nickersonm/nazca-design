"""
Submodule for integration with with GenSol package for liner photonic circuit simulation based on Scattering Matrix
"""
import numpy as np
from copy import copy
from copy import deepcopy
from itertools import tee

import nazca as nd


solver_flag = 1
try:
    import solver as sv
    solver_flag = 0
except ImportError:
    solver_flag = 1



def get_solver(Cell,
    fullreturn=False,
    drc=True,
    prune=True,
    infolevel=0,
    atol=1e-4,
    ampl_model='ampl',
    optlength_model='optlen',
    loss_model='optloss',
    allowed=None,
):
    """Build a Solver Object from a Nazca Cell

    Args:
        Cell (Cell): Nazca cell
        fullreturn (bool): If True, the returned object is a dictionary containing all the solver of all sub-cells in the hierarchy. If False (default) only the top solver is returned.
        drc (bool): if True (default) pin connections violating drc will not be connected
        prune (bool): Remove empty branches from the solver structure. Default is True
        infolevel (int): Regulates the amount of info written to stdout. Default is 0 (no info added). The higher the number the more the information
        atol (flost): absolute tolerance of pin connection
        ampl_model (string): Name of the tracker to get the scattering amplidtude commpact model. It takes precedence over the loss and optical length model.
        optlength_model (string): Name of the tracker to get the optical length compact model. It is only used if ampl_model is not found in the cell.
        loss_model (string): Name of the tracker to get the loss compact model. It is only used if ampl_model is not found in the cell.

    Returns:
        Solver or dict: see fullreturn
    """
    global solver_flag
    if solver_flag == 1:  # Only show this message once.
        nd.main_logger('No Gensol package installed. Circuit simulation will not be possible', 'warning')
        solver_flag += 1
    return None

    allowed = {'': dict(pol=0,mode=0)} if allowed is None else allowed

    models={}
    structures={}
    cells_with_model=[]

    it1,it2=tee(nd.cell_iter(Cell, hierarchy='full', topdown=False))

    for params in it1:
        if params.cell_create:
            if params.cell.auxiliary: continue
            if 'bbox' in params.cell.cell_name: continue
            if 'icon' in params.cell.cell_name: continue
            if 'stub' in params.cell.cell_name: continue
            #if set([t.cell.cnode for t in params.branch]).intersection(cells_with_model)!=set(): continue

            if infolevel>0: sv.logger.debug(f'{params.cell.cnode}')
            #sv.logger.debug(params.cell.model_info)
            if params.cell.model_info['model'] is None:
                models[params.cell.cnode] = sv.Model_from_NazcaCM.nazca_init(params.cell,
                                                                             ampl_model=ampl_model,
                                                                             loss_model=loss_model,
                                                                             optlength_model=optlength_model,
                                                                             allowed=allowed
                                                                             )
            else:
                models[params.cell.cnode]=params.cell.model_info['model']
            if not models[params.cell.cnode].is_empty():
                cells_with_model.append(params.cell.cnode)
                if infolevel>1: sv.logger.debug(f'  {models[params.cell.cnode]}')
                if infolevel>0: sv.logger.debug('')
                continue


            params.iters['instance'],itcopy=tee(params.iters['instance'])
            for inode, xya, flip in itcopy:
                if inode.instance.cell.auxiliary: continue
                if 'bbox' in inode.instance.cell.cell_name: continue
                if 'stub' in inode.instance.cnode.up.cell.cell_name: continue
                if 'icon' in inode.instance.cnode.up.cell.cell_name: continue
                param_mapping = inode.instance.model_info['param_mapping'] if 'param_mapping' in inode.instance.model_info else {}
                if infolevel>1: sv.logger.debug(f'  {inode}')
                if isinstance(models[inode.instance.cell.cnode], sv.Model):
                    ST=sv.Structure(model=models[inode.instance.cell.cnode], param_mapping=param_mapping)
                if isinstance(models[inode.instance.cell.cnode],sv.Solver):
                    ST=sv.Structure(solver=models[inode.instance.cell.cnode], param_mapping=param_mapping)
                structures[inode] = ST
                models[params.cell.cnode].add_structure(ST)


            for inode, xya, flip in params.iters['instance']:
                if inode.instance.cell.auxiliary: continue
                if 'bbox' in inode.instance.cell.cell_name: continue
                if 'stub' in inode.instance.cnode.up.cell.cell_name: continue
                if 'icon' in inode.instance.cnode.up.cell.cell_name: continue
                if infolevel>0: sv.logger.debug(f'  {inode}')
                for name,pin in inode.instance.pin.items():
                    #if infolevel>1: sv.logger.debug(f'    {name} : {pin}')
                    if name not in structures[inode].get_pin_basenames(): continue
                    if name in ['org'] : continue
                    if pin.type == 'bbox' : continue
                    if 'stub' in pin.cnode.up.cell.cell_name: continue
                    if 'icon' in pin.cnode.up.cell.cell_name: continue
                    if infolevel>0: sv.logger.debug(f'    {name} : {pin}')
                    new_nb_geo=[]
                    for (tnode,pointer) in pin.nb_geo:
                        if infolevel>1: sv.logger.debug(f'      {tnode}')
                        if tnode.type == 'bbox' : continue
                        if tnode.name.endswith('org'): continue
                        if 'bbox' in tnode.name: continue
                        if tnode.io is None: continue
                        new_nb_geo.append((tnode,pointer))
                    connected = False
                    for (tnode,pointer) in new_nb_geo:
                        if tnode.name.isdigit() and len(new_nb_geo)>1: continue
                        if drc==True and pin.xs!=tnode.xs: continue
                        xya=list(nd.diff(pin, tnode))
                        xya[2] = xya[2] % 360.0
                        if infolevel>0: sv.logger.debug(f'      {tnode} : {xya}')
                        mode_out = structures[inode].get_pin_modenames(name)
                        if infolevel>2: sv.logger.debug(f'      {mode_out}, {xya}')
                        if np.all(np.isclose(xya, [0.0,0.0,180.0], atol=atol)) and tnode.up.name in structures[tnode.cnode].get_pin_basenames():
                            models[params.cell.cnode].connect_all(structures[inode], name , structures[tnode.cnode], tnode.up.name)
                            connected = True
                    if connected: continue
                    for (tnode,pointer) in new_nb_geo:
                        if tnode.name.isdigit() and len(new_nb_geo)>1: continue
                        if drc==True and pin.xs!=tnode.xs: continue
                        xya=list(nd.diff(pin, tnode))
                        xya[2] = xya[2] % 360.0
                        if infolevel>0: sv.logger.debug(f'      {tnode} : {xya}')
                        mode_out = structures[inode].get_pin_modenames(name)
                        if infolevel>2: sv.logger.debug(f'      {mode_out}, {xya}')
                        if np.all(np.isclose(xya, [0.0,0.0,0.0], atol=atol)) or np.all(np.isclose(xya, [0.0,0.0,360.0], atol=atol)):
                            for mi in mode_out:
                                name_out = name if mi=='' else '_'.join([name,mi])
                                name_in  = tnode.name if mi=='' else '_'.join([tnode.name,mi])
                                models[params.cell.cnode].map_pins( {name_in : structures[inode].pin[name_out]})
                                if infolevel>1: sv.logger.debug(f'        {name_in} : {name_out}')
            if infolevel>0: sv.logger.debug('\n')

    if prune:
        for cnode, solver in copy(models).items():
            if solver.prune():
                models.pop(cnode)

    if len(models) ==0:
        nd.main_logger(f'Impossible to buld solver for cell {Cell.name}', 'error')
        return None

    if fullreturn:
        return models
    else:
        return models[Cell.cnode]



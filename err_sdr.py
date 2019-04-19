import io
import os
import yaml
import shutil
import argparse
import numpy as np
from sets import Set
from pyne import nucname
from pyne.mesh import Mesh, NativeMeshTag
from pyne.bins import pointwise_collapse
from pyne.material import Material, MaterialLibrary
from pyne.partisn import write_partisn_input, isotropic_vol_source, mesh_to_isotropic_source
from pyne.dagmc import discretize_geom, load, cell_material_assignments
from pyne.alara import calc_eta, calc_T
from pyne.cccc import Atflux
from pyne.variancereduction import cadis
from pyne.mcnp import Wwinp

config_filename = 'config.yml'


def step5(cfg, cfg2, cfg5):
    """
    This function creates the adjoint neutron source and writes the
    PARTISN input file for adjoint neutron transport.
 
    Parameters
    ----------
    cfg : dictionary
        User input from 'general' section of config.yml file 
    cfg2 : dictionary
        User input for step 2 from the config.yml file
    cfg5 : dictionary 
        User input for step 3 from the config.yml file 
    """
    # Get user-input from config file
    num_n_groups = cfg['n_groups']
    num_p_groups = cfg['p_groups']
    n_geom = cfg2['n_geom_file']
    decay_times = str(cfg2['decay_times']).split(' ')
    meshflux = cfg5['meshflux']
    
    # Base of geometry file name
    basename = n_geom.split("/")[-1].split(".")[0]

    # The total number of photon and neutron energy groups 
    total_num_groups = num_n_groups + num_p_groups

    # Read given flux file and tag to a mesh
    if meshflux:
      # Adjoint flux file is an hdf5 mesh file 
      fw_n_err = meshflux
      m = Mesh(structured=True, mesh=fw_n_err, mats=None)
    else:
      raise RuntimeError("No neutron flux file given")

    # Size of flux tag is equal to the total number of energy groups
    m.ERROR_TAG = NativeMeshTag(num_n_groups, name="ERROR_TAG")
    fw_n_err = m.ERROR_TAG[:]

    # Load geometry and get material assignments
    load(n_geom)
    ml = MaterialLibrary(n_geom)
    mat_names = list(ml.keys())
    cell_mats = cell_material_assignments(n_geom)
    
    # Load T matrix
    if not os.path.exists('step2_T.npy'):
        raise RuntimeError("T matrix from step 2 (step2_T.npy) not found")
    T = np.load('step2_T.npy')
    
    # Loop over mesh voxels and calculate the error in the photon source by multiplying neutron flux and T matrix 
    dg = discretize_geom(m)
    for t, dt in enumerate(decay_times): 
        temp = np.zeros(shape=(len(m), num_p_groups))
        for i in range(len(m)):
            for row in np.atleast_1d(dg[dg["idx"] == i]):
                cell = row[1]
                if not cell_mats[cell] in mat_names:
                    continue
                vol_frac = row[2]
                mat = mat_names.index(cell_mats[cell])
                for h in range(num_p_groups):
                    for g in range(num_n_groups):
                        temp[i, h] += (fw_n_err[i, g]**2)*T[mat, t, g, h]*vol_frac
                    print("err ", temp[i,h])
        # Tag the mesh with the squared error in the photon source values
        tag_name = "sq_err_q_src_{0}".format(dt)
        m.err_q_src = NativeMeshTag(num_p_groups, name=tag_name)
        m.err_q_src[:] = temp
    
    # Save adjoint neutron source mesh file tagged with values for all decay times   
    m.save("sq_err_q_src.h5m")

#def get_mesh_elements(filename, tag_name, element_dict, num_e_groups):
#def step6(cfg, cfg6):
#    """
#    """
#    num_n_groups = cfg['n_groups']
#    num_p_groups = cfg['p_groups']
#    total_num_groups = num_n_groups + num_p_groups
#
#    adj_flux_mesh = Mesh(structured=True, mesh='adj_n_mesh.h5m', mats=None)
#    adj_flux_mesh.flux = NativeMeshTag(total_num_groups, name="flux")
#    adj_p_flux = adj_flux_mesh[:]
#    adj_flux_mesh.flux_42 = NativeMeshTag(total_num_groups, name="flux_42")
#    adj_flux_mesh.flux_42[:] = adj_p_flux[174:]
##    adj_flux_tag = "flux_42"
#
#    p_src_mesh = Mesh(structured=True, mesh='p_src_mesh.h5m', mats=None)
#    p_src_tag_name = cfg6['p_src_tag']
#    p_src_mesh.p_src_tag_name = NativeMeshTag(num_p_groups)
#    p_src = p_src_mesh[:]


def main():
    """ 
    This function manages the setup and steps 1-5 for the GT-CADIS workflow.
    """
    
    parser = argparse.ArgumentParser()
    #subparsers = parser.add_subparsers(help=gtcadis_help, dest='command')
    subparsers = parser.add_subparsers(dest='command')
    step5_parser = subparsers.add_parser('step5')

    args, other = parser.parse_known_args()
    if args.command == 'setup':
        setup()
    else:
        with open(config_filename, 'r') as f:
            cfg = yaml.load(f)
            
    if args.command == 'step5':
        step5(cfg['general'], cfg['step2'], cfg['step5'])    
    elif args.command == 'get_mesh_elements':
        step5(cfg['general'])    

if __name__ == '__main__':
    main()

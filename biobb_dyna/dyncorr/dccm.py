from biobb_common.generic.biobb_object import BiobbObject
from biobb_common.configuration import settings
from biobb_common.tools import file_utils as fu
from biobb_common.tools.file_utils import launchlogger
from mdigest.core.correlation import DynCorr
from mdigest.core.dcorrelation import DihDynCorr
import MDAnalysis as mda
import numpy as np
import argparse
try:
    from mdigest.core.parsetrajectory import MDS
except ImportError:
    import importlib
    MDS = importlib.import_module('mdigest.core.parsetrajectory').MDS

class Dccm(BiobbObject):
    """
    | biobb_dccm Dccm
    | Computes per-atom covariance, correlation or lmi correlation matrices from an MD trajectory.
    | Computes per-atom covariance, correlation or lmi correlation matrices from an MD trajectory.

    Args:
        input_traj_path (str): Path to trajectory file (e.g., XTC, DCD). File type: input. Accepted formats: xtc, dcd (edam:format_1476).
        input_top_path (str): Path to topology file (e.g., PDB). File type: input. Accepted formats: pdb (edam:format_1476).
        output_matrix_path (str): Path to output matrix file. File type: output. Accepted formats: npz (edam:format_1476).
        properties (dict - Python dictionary object containing the tool parameters, not input/output files):
            * **matrix_type** (*str*) - ("correlation") Options: 'covariance', 'correlation', 'mi', 'lmi'.
            * **selection** (*str*) - ("name CA") Atom selection (MDAnalysis syntax, e.g., "name CA" for alpha carbons). Used for displacements; for dihedrals, defaults to 'backbone'.
            * **dynamics_type** (*str*) - ("displacements") Options: 'displacements', 'dihedrals'.
            * **num_replicas** (*int*) - (1) Number of replicas in the trajectory.
            * **system_selstr** (*str*) - ("protein") MDAnalysis selection for the system.
            * **align_selection** (*str*) - (None) MDAnalysis selection for alignment. If None, defaults to "protein and {selection}" for displacements or "protein and backbone" for dihedrals.
            * **in_mem** (*bool*) - (True) Load and align trajectory in memory.
            * **reference** (*any*) - (None) Reference for alignment.
            * **stride_initial** (*int*) - (0) Starting frame for striding.
            * **stride_final** (*int*) - (-1) Ending frame for striding.
            * **stride_step** (*int*) - (1) Step size for striding.
            * **scale** (*bool*) - (True) For parse_dynamics: Remove mean from coordinates.
            * **normalize** (*bool*) - (True) For parse_dynamics and parse_dih_dynamics (as **kwargs): Normalize cross-correlation matrices.
            * **lmi** (*str*) - (None) Linearized mutual information method ('gaussian' or None).
            * **mi** (*str*) - ("None") Mutual information method. None to skip computation of MI based correlation. 'knn_arg1_arg2' to compute MI, with k = arg1, and estimator= arg2, default is 'knn_5_1'
            * **dcc** (*bool*) - (False) Compute dynamical cross-correlation.
            * **pcc** (*bool*) - (False) Compute Pearson correlation.
            * **cov_disp** (*bool*) - (False) Compute covariance of displacements/dihedrals.
            * **verbose** (*bool*) - (True) Enable verbose output.
            * **mean_center** (*bool*) - (True) For parse_dih_dynamics: Subtract mean from dihedral data.
            * **center** (*str*) - (None) For parse_dih_dynamics **kwargs: Covariance computation method ('mean' or 'square_disp').
            * **subset** (*list*) - (None) For parse_dih_dynamics **kwargs: Indices for nodes in MI computation.
            * **remove_tmp** (*bool*) - (True) Remove temporary files.
            * **restart** (*bool*) - (False) Do not execute if output files exist.

    Examples:

        from biobb_dccm.biobb_dccm.dccm import dccm
        prop = {'matrix_type': 'lmi', 'selection': 'name CA', 'dynamics_type': 'dihedrals'}
        dccm(input_traj_path='/path/to/traj.xtc',
             input_top_path='/path/to/top.pdb',
             output_matrix_path='/path/to/matrix.npz',
             properties=prop)

    Info:
        * wrapped_software:
            * name: MDiGest
            * version: >=0.1.0
            * license: Apache-2.0
        * ontology:
            * name: EDAM
            * schema: http://edamontology.org/EDAM.owl
    """

    def __init__(self, input_traj_path, input_top_path, output_matrix_path, properties=None, **kwargs) -> None:
        properties = properties or {}
        super().__init__(properties, **kwargs)
        self.locals_var_dict = locals().copy()

        # Input/Output files
        self.io_dict = {
            "in": {"input_traj_path": input_traj_path, 
                   "input_top_path": input_top_path},
            "out": {"output_matrix_path": output_matrix_path}
        }

        # Properties
        self.matrix_type = properties.get('matrix_type', 'correlation')
        self.selection = properties.get('selection', 'name CA')
        self.dynamics_type = properties.get('dynamics_type', 'displacements')
        self.num_replicas = properties.get('num_replicas', 1)
        self.system_selstr = properties.get('system_selstr', 'protein')
        self.align_selection = properties.get('align_selection', None)
        self.in_mem = properties.get('in_mem', True)
        self.reference = properties.get('reference', None)
        self.stride_initial = properties.get('stride_initial', 0)
        self.stride_final = properties.get('stride_final', -1)
        self.stride_step = properties.get('stride_step', 1)
        self.scale = properties.get('scale', True)
        self.normalize = properties.get('normalize', True)
        self.lmi = properties.get('lmi', None)
        self.mi = properties.get('mi', 'None')
        self.dcc = properties.get('dcc', False)
        self.pcc = properties.get('pcc', False)
        self.cov_disp = properties.get('cov_disp', False)
        self.verbose = properties.get('verbose', True)
        self.mean_center = properties.get('mean_center', True)
        self.center = properties.get('center', None)
        self.subset = properties.get('subset', None)
        self.properties = properties
        
        # Check the properties
        self.check_properties(properties)
        self.check_arguments()

    @launchlogger
    def launch(self) -> int:
        """Execute the computation using MDiGest."""
        fu.log(f'Executing Dccm with matrix_type: {self.matrix_type}', self.out_log, self.global_log)

        # Load universe with MDAnalysis
        u = mda.Universe(self.io_dict['in']['input_top_path'], self.io_dict['in']['input_traj_path'])

        # Setup MDS
        mds = MDS()
        mds.set_num_replicas(self.num_replicas)
        mds.load_system(self.io_dict['in']['input_top_path'], self.io_dict['in']['input_traj_path'], inmem=self.in_mem)
        
        if self.align_selection is None:
            if self.dynamics_type == 'dihedrals':
                align_sel = "protein and backbone"
            else:
                align_sel = f"protein and {self.selection}"
        else:
            align_sel = self.align_selection
        mds.align_traj(inmem=self.in_mem, reference=self.reference, selection=align_sel)
        
        if self.dynamics_type == 'dihedrals':
            atom_group = "protein and backbone"
        else:
            atom_group = f"protein and {self.selection}"
        mds.set_selection(atom_group, self.system_selstr)
        mds.stride_trajectory(initial=self.stride_initial, final=self.stride_final, step=self.stride_step)

        # Configure parameters based on matrix_type if not overridden
        lmi = self.lmi
        mi  = self.mi
        dcc = self.dcc
        pcc = self.pcc
        cov_disp = self.cov_disp

        if self.matrix_type == 'covariance':
            cov_disp = True
        elif self.matrix_type == 'correlation':
            dcc = True
        elif self.matrix_type == 'lmi':
            if lmi is None:
                lmi = 'gaussian'
        else:
            raise ValueError(f"Unsupported matrix_type: {self.matrix_type}")

        # Setup Corr object and parse
        if self.dynamics_type == 'displacements':
            dyncorr = DynCorr(mds)
            dyncorr.parse_dynamics(scale=self.scale, normalize=self.normalize, LMI=lmi, MI=mi, DCC=dcc, PCC=pcc, VERBOSE=self.verbose, COV_DISP=cov_disp)
            d = 3.0
        elif self.dynamics_type == 'dihedrals':
            dyncorr = DihDynCorr(mds)
            dyncorr.parse_dih_dynamics(mean_center=self.mean_center, LMI=lmi, MI=mi, DCC=dcc, PCC=pcc, COV_DISP=cov_disp, normalized=self.normalize, subset=self.subset, center=self.center)
            d = 4.0
        else:
            raise ValueError(f"Unsupported dynamics_type: {self.dynamics_type}")

        # Extract the matrix
        if self.matrix_type == 'covariance':
            matrix = dyncorr.cov_disp_allreplicas['rep_0']
        elif self.matrix_type == 'correlation':
            matrix = dyncorr.dcc_allreplicas['rep_0']
        elif self.matrix_type == 'lmi':
            matrix = dyncorr.gcc_allreplicas['rep_0']['gcc_lmi']

        # Set diagonal for correlation/lmi if necessary
        if self.matrix_type in ['correlation', 'lmi']:
            np.fill_diagonal(matrix, 1.0)

        # Save as NPZ
        np.savez(self.io_dict['out']['output_matrix_path'], dccm=matrix)

        fu.log('Computation complete.', self.out_log, self.global_log)
        return 0
    
def compute_dccm(input_traj_path: str, input_top_path: str, output_matrix_path: str, properties: dict = None, **kwargs) -> int:
    """Create :class:`Dccm <biobb_dyna.dyncorr.dccm.Dccm>` class and
    execute the :meth:`launch() <biobb_dyna.dyncorr.dccm.Dccm.launch>` method."""

    return Dccm(input_traj_path=input_traj_path,
                 input_top_path=input_top_path,
                 output_matrix_path=output_matrix_path,
                 properties=properties, **kwargs).launch()

compute_dccm.__doc__ = Dccm.__doc__


def main():
    """Command line execution of this building block. Please check the command line documentation."""

    parser = argparse.ArgumentParser(
        description="",
        formatter_class=lambda prog: argparse.RawTextHelpFormatter(prog, width=99999),
    )

    parser.add_argument(
    "-c",
    "--config",
    required=False,
    help="This file can be a YAML file, JSON file or JSON string",
    )

    required_args = parser.add_argument_group("required arguments")
    required_args.add_argument(
        "-i",
        "--input_topology_path",
        required=True,
        help="Input topology file path"
    )
    
    required_args.add_argument(
        "-f",
        "--input_trajectory_path",
        required=True,
        help="Input trajectory file path"
    )

    required_args.add_argument(
        "-o",
        "--output_matrix_path",
        required=True,
        help="Output matrix file path"
    )

    args = parser.parse_args()
    config = args.config if args.config else None
    properties = settings.ConfReader(config=config).get_prop_dic()

    compute_dccm(input_trajectory_path=args.input_trajectory_path,
                  input_topology_path=args.input_topology_path,
                  output_matrix_path=args.output_matrix_path,
                  properties=properties,
                  )

if __name__ == "__main__":
    main()

# # Usage example

# prop = {
#     'matrix_type': 'lmi',  # Options: 'covariance', 'correlation', 'lmi'
#     'selection': 'name CA',  # MDTraj selection, e.g., alpha carbons for protein backbone
#     'dynamics_type': 'displacements',  # Options: 'displacements', 'dihedrals'
#     'num_replicas': 1,  # Number of replicas in the trajectory
#     'system_selstr': 'protein',  # MDAnalysis selection for the system
#     'align_selection': None,  # MDAnalysis selection for alignment
#     'in_mem': True,  # Load and align trajectory in memory
#     'reference': None,  # Reference for alignment
#     'stride_initial': 0,  # Starting frame for striding
#     'stride_final': -1,  # Ending frame for striding
#     'stride_step': 1,  # Step size for striding
#     'scale': True,  # For parse_dynamics: Remove mean from coordinates
#     'normalize': True,  # For parse_dynamics and parse_dih_dynamics (as **kwargs): Normalize cross-correlation matrices
#     'lmi': None,  # Linearized mutual information method ('gaussian' or None)
#     'mi': None,  # Mutual information method. None to skip computation of MI based correlation
#     'dcc': False,  # Compute dynamical cross-correlation
#     'pcc': False,  # Compute Pearson correlation
#     'cov_disp': False,  # Compute covariance of displacements/dihedrals
#     'verbose': True,  # Enable verbose output
#     'mean_center': True,  # For parse_dih
#     'center': None,  # For parse_dih **kwargs: Covariance computation method ('mean' or 'square_disp')
#     'subset': None,  # For parse_dih **kwargs: Indices for nodes in MI computation
# }

# # Run the class
# dccm = Dccm(
#     input_top_path='/content/drive/MyDrive/AI/Data/PDB/WT_apo_md_1_CA_ChainsA_DRY_imagedFit_first.pdb',  # e.g., protein structure
#     input_traj_path='/content/drive/MyDrive/AI/Data/XTC/WT_apo_CA_ChainsA_2500.xtc',  # e.g., MD simulation output
#     output_matrix_path='output_matrix.npz',  # Saved as NumPy archive
#     properties=prop
# )
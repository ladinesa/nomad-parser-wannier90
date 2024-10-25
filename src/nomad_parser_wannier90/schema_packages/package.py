from nomad.config import config
from nomad.metainfo import SchemaPackage

from nomad_simulations.schema_packages import general, model_system, model_method, numerical_settings, outputs, properties
from nomad.parsing.file_parser.mapping_parser import MappingAnnotationModel


configuration = config.get_plugin_entry_point(
    'nomad_parser_wannier90.schema_packages:schema_package_entry_point'
)

m_package = SchemaPackage()


OUT_ANNOTATION_KEY = 'out'
IN_ANNOTATION_KEY = 'in'
BAND_ANNOTATION_KEY = 'band'


class Program(general.Program):

    general.Program.version.m_annotations[OUT_ANNOTATION_KEY] = MappingAnnotationModel(mapper='.version')

class AtomsState(model_system.AtomsState):

    model_system.AtomsState.chemical_symbol.m_annotations[OUT_ANNOTATION_KEY] = MappingAnnotationModel(mapper='.@')


class AtomicCell(model_system.AtomicCell):

    model_system.AtomicCell.atoms_state.m_annotations[OUT_ANNOTATION_KEY] = MappingAnnotationModel(mapper='.labels')

    model_system.AtomicCell.positions.m_annotations[OUT_ANNOTATION_KEY] = MappingAnnotationModel(mapper='.positions', unit='angstrom')

    model_system.AtomicCell.lattice_vectors.m_annotations[OUT_ANNOTATION_KEY] = MappingAnnotationModel(mapper=('get_lattice_vectors', ['lattice_vectors']), unit='angstrom')


class ModelSystem(model_system.ModelSystem):

    model_system.AtomicCell.m_def.m_annotations[OUT_ANNOTATION_KEY] = MappingAnnotationModel(mapper='.@')

    model_system.ModelSystem.model_system.m_annotations[IN_ANNOTATION_KEY] = MappingAnnotationModel(mapper=('get_projections', ['.projections']))

    model_system.ModelSystem.branch_label.m_annotations[IN_ANNOTATION_KEY] = MappingAnnotationModel(mapper=('get_branch_label_indices', ['.atom', 'structure.positions', 'structure.labels','lattice_vectors']), search='label', cache=True)

    model_system.ModelSystem.atom_indices.m_annotations[IN_ANNOTATION_KEY] = MappingAnnotationModel(mapper=('get_branch_label_indices', ['.atom', 'structure.positions', 'structure.labels','lattice_vectors']), search='indices')


class KMesh(numerical_settings.KMesh):

    numerical_settings.KMesh.n_points.m_annotations[OUT_ANNOTATION_KEY] = MappingAnnotationModel(mapper='.n_points')


class KSpace(numerical_settings.KSpace):

    numerical_settings.KSpace.k_mesh.m_annotations[OUT_ANNOTATION_KEY] = MappingAnnotationModel(mapper='.k_mesh')


class ModelMethod(model_method.ModelMethod):

    numerical_settings.KSpace.m_def.m_annotations[OUT_ANNOTATION_KEY] = MappingAnnotationModel(mapper='.@')


class Wannier(model_method.Wannier):

    model_method.Wannier.is_maximally_localized.m_annotations[OUT_ANNOTATION_KEY] = MappingAnnotationModel(mapper=('is_maximally_localized', ['.Niter'], dict(default=0)))

    model_method.Wannier.energy_window_outer.m_annotations[OUT_ANNOTATION_KEY] = MappingAnnotationModel(mapper='.energy_windows.outer')

    model_method.Wannier.n_orbitals.m_annotations[OUT_ANNOTATION_KEY] = MappingAnnotationModel(mapper='.Wannier')


class ElectronicBandStructure(properties.ElectronicBandStructure):

    properties.ElectronicBandStructure.n_bands.m_annotations[OUT_ANNOTATION_KEY] = MappingAnnotationModel(mapper='.Nwannier')

    properties.ElectronicBandStructure.value.m_annotations[BAND_ANNOTATION_KEY] = MappingAnnotationModel(mapper=('get_data', ['.data']))


class Outputs(outputs.Outputs):

    outputs.Outputs.electronic_band_structures.m_annotations[OUT_ANNOTATION_KEY] = MappingAnnotationModel(mapper='.@')

    outputs.Outputs.electronic_band_structures.m_annotations[BAND_ANNOTATION_KEY] = MappingAnnotationModel(mapper='.@')


class Simulation(general.Simulation):

    general.Simulation.program.m_annotations[OUT_ANNOTATION_KEY] = MappingAnnotationModel(mapper='.@')

    general.Simulation.model_system.m_annotations[OUT_ANNOTATION_KEY] = MappingAnnotationModel(mapper='.structure')

    general.Simulation.model_system.m_annotations[IN_ANNOTATION_KEY] = MappingAnnotationModel(mapper='.@')

    model_method.Wannier.m_def.m_annotations[OUT_ANNOTATION_KEY] = MappingAnnotationModel(mapper='.@')

    general.Simulation.outputs.m_annotations[OUT_ANNOTATION_KEY] = MappingAnnotationModel(mapper='.@')

    general.Simulation.outputs.m_annotations[BAND_ANNOTATION_KEY] = MappingAnnotationModel(mapper='.@')


Simulation.m_def.m_annotations[OUT_ANNOTATION_KEY] = MappingAnnotationModel(mapper='@')

Simulation.m_def.m_annotations[IN_ANNOTATION_KEY] = MappingAnnotationModel(mapper='@')

Simulation.m_def.m_annotations[BAND_ANNOTATION_KEY] = MappingAnnotationModel(mapper='@')


m_package.__init_metainfo__()
